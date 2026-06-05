"""
Sparse solver backend abstraction.

Supports three backends:
- 'scipy'    : default, uses scipy.sparse.linalg (no extra dependencies)
- 'pardiso'  : uses pypardiso (Intel MKL PARDISO), ``pip install pypardiso``
- 'cupy'     : uses CuPy (NVIDIA GPU via cuSPARSE/cuSOLVER), ``pip install cupy-cuda12x``

Optional eigensolver:
- LOBPCG with PyAMG preconditioning, ``pip install pyamg``
  Supports warm-starting from previous eigenvectors (X0).

Dependencies
------------
Only scipy is required. The 'pardiso', 'cupy', and 'pyamg' packages are
optional — they are imported lazily when selected, and raise ImportError with
install instructions if the package is missing.

Usage
-----
    from pyNastran.utils.solver_utils import get_solver, set_solver

    set_solver('pardiso')  # or 'cupy' or 'scipy'
    backend = get_solver()
    x = backend.spsolve(K, F)
    eigenvalues, eigenvectors = backend.eigsh(K, k=10, M=M)

    # LOBPCG with warm-start (optional, all backends)
    eigenvalues, eigenvectors = backend.lobpcg(K, k=40, M=M, X0=prev_modes)
"""

from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from scipy.sparse import spmatrix

_ACTIVE_BACKEND: str = "scipy"


def set_solver(name: str) -> None:
    """Set the active sparse solver backend.

    Parameters
    ----------
    name : str
        One of 'scipy', 'pardiso', 'cupy'.
    """
    global _ACTIVE_BACKEND
    name = name.lower().strip()
    valid = ("scipy", "pardiso", "cupy")
    if name not in valid:
        raise ValueError(f"Unknown solver backend {name!r}. Must be one of {valid}.")
    _ACTIVE_BACKEND = name


def get_solver() -> "SolverBackend":
    """Return the active solver backend instance."""
    if _ACTIVE_BACKEND == "scipy":
        return _ScipyBackend()
    elif _ACTIVE_BACKEND == "pardiso":
        return _PardisoBackend()
    elif _ACTIVE_BACKEND == "cupy":
        return _CuPyBackend()
    raise RuntimeError(f"Invalid backend state: {_ACTIVE_BACKEND!r}")


class SolverBackend:
    """Base class defining the solver interface."""

    name: str = ""

    def spsolve(self, A: spmatrix, b: np.ndarray) -> np.ndarray:
        """Solve A x = b for sparse A."""
        raise NotImplementedError

    def eigsh(
        self,
        A: spmatrix,
        k: int,
        M: np.ndarray | spmatrix | None = None,
        sigma: float | None = None,
        which: str = "SM",
        **kwargs,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Sparse symmetric eigenvalue solve (shift-invert Lanczos)."""
        raise NotImplementedError

    def lobpcg(
        self,
        A: spmatrix,
        k: int,
        M: np.ndarray | spmatrix | None = None,
        X0: np.ndarray | None = None,
        tol: float = 1e-6,
        maxiter: int = 400,
        largest: bool = False,
    ) -> tuple[np.ndarray, np.ndarray]:
        """LOBPCG eigensolver with optional AMG preconditioning and warm-start.

        Parameters
        ----------
        A : spmatrix
            Stiffness matrix (sparse, symmetric).
        k : int
            Number of eigenvalues/eigenvectors to compute.
        M : array_like, optional
            Mass matrix (sparse or dense). If None, standard eigenproblem.
        X0 : ndarray, optional
            Initial guess for eigenvectors, shape (n, k). Pass previous
            converged modes for warm-start (3-5x faster on parameter sweeps).
        tol : float
            Convergence tolerance.
        maxiter : int
            Maximum iterations.
        largest : bool
            If False (default), find smallest eigenvalues.

        Returns
        -------
        eigenvalues : ndarray of shape (k,)
        eigenvectors : ndarray of shape (n, k)

        Notes
        -----
        If pyamg is installed (``pip install pyamg``), a smoothed-aggregation
        AMG preconditioner is built automatically. This makes convergence
        nearly mesh-independent. Without pyamg, LOBPCG still works but
        converges slower on cold starts.
        """
        raise NotImplementedError

    def factorize(self, A: spmatrix) -> "callable[[np.ndarray], np.ndarray]":
        """Factor A once, return a callable that solves A x = b via back-substitution.

        Parameters
        ----------
        A : spmatrix
            Sparse matrix to factor (LU or LDL^T depending on backend).

        Returns
        -------
        solve : callable
            Function with signature solve(b) -> x that reuses the stored factors.
        """
        raise NotImplementedError

    def inv(self, A: np.ndarray) -> np.ndarray:
        """Dense matrix inverse (for frequency response)."""
        raise NotImplementedError


# ---------------------------------------------------------------------------
# Helper: build AMG preconditioner for LOBPCG
# ---------------------------------------------------------------------------
def _build_amg_preconditioner(A: spmatrix) -> object | None:
    """Build a smoothed-aggregation AMG preconditioner if pyamg is available."""
    try:
        import pyamg
        from scipy.sparse import csr_matrix

        A_csr = csr_matrix(A)
        ml = pyamg.smoothed_aggregation_solver(A_csr)
        return ml.aspreconditioner()
    except ImportError:
        return None


def _run_lobpcg_scipy(
    A: spmatrix,
    k: int,
    M=None,
    X0: np.ndarray | None = None,
    tol: float = 1e-6,
    maxiter: int = 400,
    largest: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Run scipy LOBPCG with optional AMG preconditioner."""
    from scipy.sparse.linalg import lobpcg
    from scipy.sparse import issparse, csr_matrix

    n = A.shape[0]

    if X0 is None:
        rng = np.random.default_rng(42)
        X0 = rng.standard_normal((n, k))

    if X0.shape[1] < k:
        rng = np.random.default_rng(42)
        extra = rng.standard_normal((n, k - X0.shape[1]))
        X0 = np.column_stack([X0, extra])
    elif X0.shape[1] > k:
        X0 = X0[:, :k]

    precond = _build_amg_preconditioner(A)

    M_arg = None
    if M is not None:
        M_arg = csr_matrix(M) if not issparse(M) else M

    eigenvalues, eigenvectors = lobpcg(
        A, X0, B=M_arg, M=precond, tol=tol, maxiter=maxiter, largest=largest
    )

    # Sort by eigenvalue magnitude (smallest first)
    idx = np.argsort(eigenvalues)
    return eigenvalues[idx], eigenvectors[:, idx]


# ---------------------------------------------------------------------------
# SciPy backend (default)
# ---------------------------------------------------------------------------
class _ScipyBackend(SolverBackend):
    name = "scipy"

    def spsolve(self, A: spmatrix, b: np.ndarray) -> np.ndarray:
        from scipy.sparse.linalg import spsolve

        return spsolve(A, b)

    def factorize(self, A: spmatrix):
        from scipy.sparse.linalg import splu
        from scipy.sparse import csc_matrix

        lu = splu(csc_matrix(A))
        return lu.solve

    def eigsh(
        self, A: spmatrix, k: int, M=None, sigma=None, which="SM", **kwargs
    ) -> tuple[np.ndarray, np.ndarray]:
        from scipy.sparse.linalg import eigsh

        return eigsh(A, k=k, M=M, sigma=sigma, which=which, **kwargs)

    def lobpcg(self, A, k, M=None, X0=None, tol=1e-8, maxiter=200, largest=False):
        return _run_lobpcg_scipy(A, k, M=M, X0=X0, tol=tol, maxiter=maxiter, largest=largest)

    def inv(self, A: np.ndarray) -> np.ndarray:
        from scipy.linalg import inv

        return inv(A)


# ---------------------------------------------------------------------------
# pypardiso backend (Intel MKL PARDISO)
# ---------------------------------------------------------------------------
class _PardisoBackend(SolverBackend):
    name = "pardiso"

    def __init__(self) -> None:
        try:
            from pypardiso import spsolve as _pardiso_spsolve  # noqa: F401
        except ImportError as e:
            raise ImportError(
                "pypardiso is not installed. Install with: pip install pypardiso"
            ) from e

    def spsolve(self, A: spmatrix, b: np.ndarray) -> np.ndarray:
        from pypardiso import spsolve
        from scipy.sparse import csr_matrix

        A_csr = csr_matrix(A)
        return spsolve(A_csr, b)

    def factorize(self, A: spmatrix):
        from pypardiso import factorized
        from scipy.sparse import csc_matrix

        return factorized(csc_matrix(A))

    def eigsh(
        self, A: spmatrix, k: int, M=None, sigma=None, which="SM", **kwargs
    ) -> tuple[np.ndarray, np.ndarray]:
        from scipy.sparse import csc_matrix, eye
        from scipy.sparse.linalg import eigsh, LinearOperator

        # pypardiso doesn't provide eigsh; use PARDISO as shift-invert operator
        if sigma is not None:
            from pypardiso import factorized

            A_shifted = A - sigma * (M if M is not None else eye(A.shape[0]))
            A_shifted_csc = csc_matrix(A_shifted)
            solve_shifted = factorized(A_shifted_csc)
            OPinv = LinearOperator(A.shape, matvec=solve_shifted)
            return eigsh(A, k=k, M=M, sigma=sigma, which=which, OPinv=OPinv, **kwargs)
        return eigsh(A, k=k, M=M, sigma=sigma, which=which, **kwargs)

    def lobpcg(self, A, k, M=None, X0=None, tol=1e-8, maxiter=200, largest=False):
        return _run_lobpcg_scipy(A, k, M=M, X0=X0, tol=tol, maxiter=maxiter, largest=largest)

    def inv(self, A: np.ndarray) -> np.ndarray:
        from scipy.linalg import inv

        return inv(A)


# ---------------------------------------------------------------------------
# CuPy backend (NVIDIA GPU)
# ---------------------------------------------------------------------------
class _CuPyBackend(SolverBackend):
    name = "cupy"

    def __init__(self) -> None:
        try:
            import cupy  # noqa: F401
            import cupyx.scipy.sparse  # noqa: F401
            import cupyx.scipy.sparse.linalg  # noqa: F401
        except ImportError as e:
            raise ImportError(
                "CuPy is not installed. Install with: pip install cupy-cuda12x "
                "(adjust cuda version as needed)"
            ) from e

    def spsolve(self, A: spmatrix, b: np.ndarray) -> np.ndarray:
        import cupy as cp
        import cupyx.scipy.sparse as cp_sparse
        import cupyx.scipy.sparse.linalg as cp_linalg

        A_gpu = cp_sparse.csc_matrix(A.tocsc())
        b_gpu = cp.asarray(b)
        x_gpu = cp_linalg.spsolve(A_gpu, b_gpu)
        return cp.asnumpy(x_gpu)

    def factorize(self, A: spmatrix):
        from scipy.sparse.linalg import splu
        from scipy.sparse import csc_matrix

        lu = splu(csc_matrix(A))
        return lu.solve

    def eigsh(
        self, A: spmatrix, k: int, M=None, sigma=None, which="SM", **kwargs
    ) -> tuple[np.ndarray, np.ndarray]:
        import cupy as cp
        import cupyx.scipy.sparse as cp_sparse
        import cupyx.scipy.sparse.linalg as cp_linalg
        from scipy.sparse import issparse, csc_matrix

        A_gpu = cp_sparse.csc_matrix(A.tocsc())
        M_gpu = None
        if M is not None:
            if issparse(M):
                M_gpu = cp_sparse.csc_matrix(M.tocsc())
            else:
                M_gpu = cp_sparse.csc_matrix(csc_matrix(M))

        # CuPy eigsh does not support sigma; fall back to scipy for shift-invert
        if sigma is not None:
            from scipy.sparse.linalg import eigsh

            return eigsh(A, k=k, M=M, sigma=sigma, which=which, **kwargs)

        eigenvalues_gpu, eigvecs_gpu = cp_linalg.eigsh(A_gpu, k=k, M=M_gpu, which=which)
        return cp.asnumpy(eigenvalues_gpu), cp.asnumpy(eigvecs_gpu)

    def lobpcg(self, A, k, M=None, X0=None, tol=1e-8, maxiter=200, largest=False):
        import cupy as cp
        import cupyx.scipy.sparse as cp_sparse
        from cupyx.scipy.sparse.linalg import lobpcg as cp_lobpcg
        from scipy.sparse import issparse, csr_matrix

        n = A.shape[0]
        A_gpu = cp_sparse.csr_matrix(csr_matrix(A))

        M_gpu = None
        if M is not None:
            M_sparse = csr_matrix(M) if not issparse(M) else csr_matrix(M)
            M_gpu = cp_sparse.csr_matrix(M_sparse)

        if X0 is None:
            X0_gpu = cp.random.randn(n, k, dtype="float64")
        else:
            if X0.shape[1] < k:
                extra = np.random.default_rng(42).standard_normal((n, k - X0.shape[1]))
                X0 = np.column_stack([X0, extra])
            elif X0.shape[1] > k:
                X0 = X0[:, :k]
            X0_gpu = cp.asarray(X0)

        eigenvalues_gpu, eigvecs_gpu = cp_lobpcg(
            A_gpu, X0_gpu, B=M_gpu, tol=tol, maxiter=maxiter, largest=largest
        )

        eigenvalues = cp.asnumpy(eigenvalues_gpu)
        eigenvectors = cp.asnumpy(eigvecs_gpu)
        idx = np.argsort(eigenvalues)
        return eigenvalues[idx], eigenvectors[:, idx]

    def inv(self, A: np.ndarray) -> np.ndarray:
        import cupy as cp
        from cupy.linalg import inv as cp_inv

        A_gpu = cp.asarray(A)
        result_gpu = cp_inv(A_gpu)
        return cp.asnumpy(result_gpu)
