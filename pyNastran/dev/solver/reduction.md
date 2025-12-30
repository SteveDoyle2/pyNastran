M-Set Reduction
---------------

$$ \begin{bmatrix}
    K_{mm}       & K_{mn} \\
    K_{nm}       & K_{nn} \\
\end{bmatrix} \begin{Bmatrix}
    u_{m} \\
    u_{n} \\
\end{Bmatrix} = \begin{Bmatrix}
    F_{m} \\
    F_{n} \\
\end{Bmatrix} $$

$$ K_{mm} u_m + K_{mn} u_n = F_m $$

$$ K_{mn} u_m + K_{nn} u_n = F_n $$

Let:

$$ u_m = G_{mn} u_n $$

$$ (K_{mm} G_{mn} + K_{mn}) u_n = F_m $$

$$ (K_{mn} G_{mn} + K_{nn}) u_n = F_n $$

MPC Definitinon
---------------
$$ \Sigma A_j u_j = 0 $$

RBE2 Definition
---------------
Independent (n) = Dependent (m)

$$ u_{m,2} = u_{n,1} $$

$$ G_mn = [1] $$

RBE3 Definition
---------------
```
       /2
#-k---1---3
       \4
```
$$ F_n1 = F_{d2} + F_{d3} + F_{d4} $$

$ K11 = 2
$ F1 = 10
K = [2 0 0 0]; F=[0]
    [0 0 0 0]    [Fd2]
    [0 0 0 0]    [Fd3]
    [0 0 0 0]    [Fd4]

$$ (K_{mn} G_{mn} + K_{nn}) u_1 = F_1 $$

$$ (K_{mn} G_{mn} + K_{nn}) u_n = F_n $$

(K_{mn} G_{mn} + K_{nn}) u_n

$$ K_{mm} u_m + K_{mn} u_n = F_m $$
$$ K_{mn} u_m + K_{nn} u_n = F_n $$

Per Elements
------------
$$ [RG_m, RG_n] [nm, xn]^T = 0 $$

$$ [RG_m] x_m + [RG_n] x_n = 0 $$

$$ x_m = -[RG_m]^{-1} [RG_n] x_n = [GM] x_n $$

So:

$$ x_m = [GM] x_n $$

$$ [GM] = -[RG_m]^{-1} [RG_n] $$

Prior to reduction:

$$ \begin{bmatrix}
    M_{mm}  & 0      \\
    0       & M_{nn} \\
\end{bmatrix} \begin{Bmatrix}
    \ddot u_{m} \\
    \ddot u_{n} \\
\end{Bmatrix} + \begin{bmatrix}
    K_{mm}  & 0      \\
    0       & K_{nn} \\
\end{bmatrix} \begin{Bmatrix}
    u_{m} \\
    u_{n} \\
\end{Bmatrix} = \begin{Bmatrix}
    F_{m} \\
    F_{n} \\
\end{Bmatrix} $$

$$ M_{mm} \ddot x_m + K_{mm} x_m = F_m $$

$$ M_{nn} \ddot x_n + K_{nn} x_n = F_n $$

$$ M_{mm} [GM] \ddot x_n + K_{mm} [GM] x_n = F_m $$

Pre-multiplying:

$$ [GM]^T M_{mm} [GM] \ddot x_n + [GM]^T K_{mm} [GM] x_n = [GM]^T F_m $$

Adding:

$$ ( M_{nn} + [GM]^T M_{mm} [GM]) \ddot x_n + (K_{nn} + [GM]^T K_{mm} [GM]) x_n = F_n + [GM]^T F_m $$

S-Set Reduction
---------------

$$ \begin{bmatrix}
    K_{ss}  &  K_{sf} \\
    K_{fs}  &  K_{ff} \\
\end{bmatrix} \begin{Bmatrix}
    u_{s} \\
    u_{f} \\
\end{Bmatrix} = \begin{Bmatrix}
    F_{s} \\
    F_{f} \\
\end{Bmatrix} $$

$$ K_{ss} u_s + K_{sf} u_f = F_s $$

$$ K_{ss} u_s = F_s - K_{sf} u_f $$

$$ u_s = K_{ss}^{-1} (F_s - K_{sf} u_f) $$


$$ K_{sf} u_s + K_{ff} u_f = F_f $$

$$ K_{ff} u_f = F_f - K_{sf} u_s $$

$$ u_f = K_{ff}^{-1} (F_f - K_{sf} u_s) $$

O-Set Reduction
---------------
