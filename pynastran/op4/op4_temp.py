
def __backup(self, name, matrix, form=2, precision='default'):
	"""
	Put this documentation somewhere else...

	Parameters
	----------
	name : str
		the name of the matrix
	matrix : ndarray
		a two-dimensional NUMPY.NDARRAY
	form : int (default=2)
		Form is defined as one of the following:
	precision : str; default=True
		{'default', 'single', 'double'}

	+======+================+
	| Form |   Definition   |
	+======+================+
	|  1   | Square         |
	+------+----------------+
	|  2   | Rectangular    |
	+------+----------------+
	|  3   | Diagonal       |
	+------+----------------+
	|  6   | Symmetric      |
	+------+----------------+
	|  8   | Id entity      |
	+------+----------------+
	|  9   | Pseudoidentity |
	+------+----------------+

	Not Supported by all OP4s (this is not a restriction of the OP4
	reader/writer)

	+======+================================+
	| Form |         Definition             |
	+======+================================+
	|  4   | Lower triangular factor        |
	+------+--------------------------------+
	|  5   | Upper triangular factor        |
	+------+--------------------------------+
	|  10  | Cholesky factor                |
	+------+--------------------------------+
	|  11  | Trapezoidal factor             |
	+------+--------------------------------+
	|  13  | Sparse lower triangular factor |
	+------+--------------------------------+
	|  15  | Sparse upper triangular factor |
	+------+--------------------------------+

	.. note:: form defaults to 2, but 1 can be easily determined.
			  Any others must be specified.
	"""
	assert isinstance(name, string_types), name
	assert isinstance(form, int), form
