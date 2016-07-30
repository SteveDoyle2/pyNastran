
import os
from op4 import OP4

def op4_matrices_get():
	"""creates dummy matrices"""
	strings = array([
		[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[1, 0, 3, 0, 5, 0, 7, 0, 9, 0, 11, 0, 13, 0, 15, 0, 17, 0, 19, 0],
		[1, 0, 3, 0, 5, 0, 7, 0, 9, 0, 11, 0, 13, 0, 15, 0, 17, 0, 19, 0],
		[1, 0, 3, 0, 5, 0, 7, 0, 9, 0, 11, 0, 13, 0, 15, 0, 17, 0, 19, 0],
		[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[1, 0, 3, 0, 5, 0, 7, 0, 9, 0, 11, 0, 13, 0, 15, 0, 17, 0, 19, 0],
		[1, 0, 3, 0, 5, 0, 7, 0, 9, 0, 11, 0, 13, 0, 15, 0, 17, 0, 19, 0],
		[1, 0, 3, 0, 5, 0, 7, 0, 9, 0, 11, 0, 13, 0, 15, 0, 17, 0, 19, 0],
		[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
		[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
		[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
		[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
		[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
		[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
		[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
		[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
		[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
		[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]], dtype='float32') # f?
	return strings

def op4_read(op4_filename=None, matrix_names=None, precision='default', debug=False, log=None):
	"""
	Reads a NASTRAN OUTPUT4 file, and stores the
	matrices as the output arguments.  The number of
	matrices read is defined by the list matrix_names.  By default, all
	matrices will be read.  The resulting output is a dictionary of
	matrices that are accessed by their name.

	.. code-block:: python

	  >>> from pyNastran.op4.op4 import OP4
	  >>> op4 = OP4()

	  # get all the matrices
	  >>> matrices = op4.read_op4(op4_filename)
	  >>> (formA, A) = matrices['A']
	  >>> (formB, B) = matrices['B']
	  >>> (formC, C) = matrices['C']

	  # or to reduce memory usage
	  >>> matrices = op4.read_op4(op4_filename, matrix_names=['A', 'B'])
	  >>> (formA, A) = matrices['A']
	  >>> (formB, B) = matrices['B']

	  # or because you only want A
	  >>> matrices = op4.read_op4(op4_filename, matrix_names='A')
	  >>> (formA, A) = matrices['A']

	  # get all the matrices, but select the file using a file dialog
	  >>> matrices = op4.read_op4()
	  >>>

	Parameters
	----------
	op4_filename : str / None
		an OP4 filename.  Type=STRING.
	matrix_names : List[str], str / None
		matrix name(s) (None -> all)
	precision : str; {'default', 'single', 'double'}
		specifies if the matrices are in single or double precsion
		which means the format will be whatever the file is in

	Returns
	-------
	matricies : dict[str] = (int, matrix)
		dictionary of matrices where the key is the name and the value is [form, matrix]

		+------+----------------+
		| Form |   Definition   |
		+======+================+
		|  1   | Square		 |
		+------+----------------+
		|  2   | Rectangular	|
		+------+----------------+
		|  3   | Diagonal	   |
		+------+----------------+
		|  6   | Symmetric	  |
		+------+----------------+
		|  8   | Id entity	  |
		+------+----------------+
		|  9   | Pseudoidentity |
		+------+----------------+

		+--------+-------------------------+
		|  Type  | Object				  |
		+========+=========================+
		| Dense  | NUMPY.NDARRAY		   |
		+--------+-------------------------+
		| Sparse | SCIPY.SPARSE.COO_MATRIX |
		+--------+-------------------------+

	.. note:: based off the MATLAB code SAVEOP4 developed by ATA-E and
			  later UCSD.
	.. note:: it's strongly recommended that you convert sparse matrices to
			  another format before doing math on them.  This is standard
			  with sparse matrices.
	"""
	op4 = OP4(log=log, debug=debug)
	return op4.read_op4(op4_filename, matrix_names, precision)


def main():
        op4_filename = os.path.abspath(os.path.join('../..', 'models', 'mat_t_dn.op4'))
        print(op4_filename)
        #assert os.path.exists(op4_filename), print_bad_path(op4_filename)

        op4 = OP4()
        matrices = op4.read(op4_filename)

        """tests various matrices"""
#	from op4_utils import write_DMIG

	#compress_column([14, 15, 16, 20, 21, 22, 26, 27, 28])
#	filenames = [
#		'test/mat_t_dn.op4',
#		'test/mat_t_s1.op4',
#		'test/mat_t_s2.op4',
#		'test/mat_b_dn.op4',
#		'test/mat_b_s1.op4',
#		'test/mat_b_s2.op4',
		#'test/b_sample.op4',
		#'binary.op4',
#	]

	#matrix_names = 'EYE10' # identity
	#matrix_names = 'LOW'
	#matrix_names = 'RND1RS' # real,single
	#matrix_names = 'RND1RD' # real,double
	#matrix_names = 'RND1CS' # complex,single
	#matrix_names = 'RND1CD' # complex,double
	#matrix_names = 'STRINGS'
	#matrix_names = 'EYE5CD' # complex identity
#	matrix_names = None
#	strings = get_matrices()

#	is_big_mat = True
#	if PY2:
#		wb = 'wb'
#	else:
#		wb = 'w'
#	with open('ascii.op4', wb) as op4_file:
#		for fname in filenames:
#			op4 = OP4()
#			op4.endian = '>'
			#if 't' in fname:
			#else:
				#f = open('binary.op4','wb')

#			matrices = op4.read_op4(fname, matrix_names=matrix_names, precision='default')
#			print("keys = %s" % matrices.keys())
#			print("fname=%s" % fname)
#			for name, (form, matrix) in sorted(iteritems(matrices)):
#				op4.write(op4_file, matrices, name_order=name)

#	print("-----------------------------")
#	print("done")
#	print("-----------------------------")

if __name__ == '__main__':  # pragma: no cover
	main()
