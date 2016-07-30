
def matrix_sparse_ascii_write(op4, name, A, form=2, is_big_mat=False, precision='default'):
	"""
	.. todo:: Does this work for complex matrices?
	"""
	msg = ''
	assert isinstance(name, string_types), 'name=%s' % name
	#A = A.tolil() # list-of-lists sparse matrix
	#print dir(A)
	matrix_type, nwords_per_value = matrix_type_nwv_get(A.data[0], precision)
	if matrix_type in [3, 4]:
		complex_factor = 2
	else: # 1, 2
		complex_factor = 1
	(nrows, ncols) = A.shape

	#if nrows == ncols and form == 2:
		#form = 1
	#print("matrix_type=%s" % matrix_type)
	if is_big_mat:
		msg += '%8i%8i%8i%8i%-8s1P,3E23.16\n' % (ncols, -nrows, form, matrix_type, name)
	else:
		msg += '%8i%8i%8i%8i%-8s1P,3E23.16\n' % (ncols, nrows, form, matrix_type, name)

	#print("A.row = ", A.row)
	#print("A.col = ", A.col)

	cols = {}
	for j in A.col:
		cols[j] = []
	for i, jcol in enumerate(A.col):
		cols[jcol].append(i)
	#print("cols = ", cols)

	op4.write(msg)
	msg = ''
	for j, col in iteritems(cols):
		#print("***********")
		#print("j=%s col=%s" % (j, col))
		#col.sort()

		#print('A =', A)
		irows = [A.row[jj] for jj in col]
		#print "irows = ",irows
		(dpacks) = matrix_column_compress(irows)
		#print("dpacks = %s" % (dpacks))

		npacks = len(dpacks)
		nrows = len(irows)
		if is_big_mat:
			#L = complex_factor * (2 * len(irows)) + 1
			L = 2 * npacks * nwords_per_value + nrows
			msg = '%8i%8i%8i\n' % (j + 1, 0, L)
		else:
			L = complex_factor * (2 * len(irows))
			msg = '%8i%8i%8i\n' % (j + 1, 0, L + 1)
		op4.write(msg)

		for (ipack, dpack) in enumerate(dpacks):
			msg = ''
			#print("pack = ",pack)

			irow = A.row[col[dpack[0]]]
			if is_big_mat:
				#L = complex_factor * (2 * len(pack)) + 1
				#L = (nPacks+1) + nRows * complex_factor
				L = (len(dpack) + 1) * nwords_per_value
				#if iPack==0:
					#L+=1

				#L = complex_factor * (2 + npacks) + 1
				#L = len(pack) + complex_factor * 2
				#msg = '%8i%8i%8i\n' % (j+1, 0, L+1)
				msg += '%8i%8i\n' % (L, irow + 1)
			else:
				#L = complex_factor * (2 * len(pack))
				#msg = '%8i%8i%8i\n' % (j+1, 0, L+1)

				IS = irow + 65536 * (L + 1) + 1
				msg += '%8i\n' % IS

			i = 0
			value_str = ''
			#print("ipack=%s rowPack=%s" % (ipack, [A.row[p] for p in dpack]))
			for p in dpack:
				irow = col[p]
				val = A.data[irow]
				irow = A.row[irow]

				if matrix_type in [1, 2]:
					value_str += '%23.16E' % val
					if (i + 1) % 3 == 0:
						msg += value_str + '\n'
						value_str = ''
				else:
					value_str += '%23.16E' % val.real
					if (i + 1) % 3 == 0:
						msg += value_str + '\n'
						value_str = ''
					i += 1
					value_str += '%23.16E' % val.imag
					if (i + 1) % 3 == 0:
						msg += value_str + '\n'
						value_str = ''
				i += 1
			if value_str:
				msg += value_str + '\n'
			op4.write(msg)
	op4.write('%8i%8i%8i\n' % (ncols + 1, 1, 1))
	op4.write(' 1.0000000000000000E+00\n')

def matrix_nrows_get(nrows):
	"""
	Parameters
	----------
	nrows : int
		the number of rows in the matrix

	Returns
	-------
	nrows : int
		the number of rows in the matrix
	BIGMAT : Input-logical-default=FALSE. BIGMAT is applicable only
		when IUNIT < 0. BIGMAT=FALSE selects the format that uses a
		string header as described under Remark 1. But, if the
		matrix has more than 65535 rows, then BIGMAT will
		automatically be set to TRUE regardless of the value
		specified.
	"""
	if nrows < 0:  # if less than 0, big
		is_big_mat = True
		nrows = abs(nrows)
	elif nrows > 0:
		is_big_mat = False
		if nrows > 65535:
			is_big_mat = True
			nrows = abs(nrows)
	else:
		raise RuntimeError('unknown BIGMAT.  nrows=%s' % nrows)
	return is_big_mat, nrows

def dtype_get(matrix_type, precision='default'):
	"""reset the type if 'default' not selected"""
	if precision == 'single':
		if matrix_type in [1, 2]:
			dtype = 'float32'
		else:
			dtype = 'complex64'
	elif precision == 'double':
		if matrix_type in [1, 2]:
			dtype = 'float64'
		else:
			dtype = 'complex128'
	else:  # default
		if matrix_type == 1:
			dtype = 'float32'
		elif matrix_type == 2:
			dtype = 'float64'
		elif matrix_type == 3:
			dtype = 'complex64'
		else:
			dtype = 'complex128'
	return dtype


def matrix_type_nwv_get(A, precision='default'):
	"""
	Determines the Type and number of words per value
	an entry in the matrix takes up.

	Parameters
	----------
	A : matrix
		a matrix or entry in a matrix (to save memory)
	precision : str
		data precision ='default', 'single', 'double'

	Returns
	-------
	matrix_type : int
		the dtype of the matrix as an integer
	nwords_per_value : int
		Number of words per value

	A word is 4 bytes

	+-------------+------------+-----------+--------+--------+
	| Matrix Type |   dtype    | precision | nwords | nbytes |
	+=============+============+===========+========+========+
	|      1      | float32    |   single  |    1   |    4   |
	+-------------+------------+-----------+--------+--------+
	|      2      | complex64  |   single  |    2   |    8   |
	+-------------+------------+-----------+--------+--------+
	|      3      | float64    |   double  |    2   |    8   |
	+-------------+------------+-----------+--------+--------+
	|      4      | complex128 |   double  |    4   |   16   |
	+-------------+------------+-----------+--------+--------+
	"""
	# real
	if isinstance(A.dtype.type(), float32):
		nwords_per_value = 1
		if precision != 'double':
			matrix_type = 1
		else:
			matrix_type = 2
	elif isinstance(A.dtype.type(), float64):
		nwords_per_value = 1
		if precision != 'single':
			matrix_type = 2
		else:
			matrix_type = 1

	# complex
	elif isinstance(A.dtype.type(), complex64):
		nwords_per_value = 2
		if precision != 'double':
			matrix_type = 3
		else:
			matrix_type = 4
	elif isinstance(A.dtype.type(), complex128):
		nwords_per_value = 2
		if precision != 'single':
			matrix_type = 4
		else:
			matrix_type = 3
	else:
		msg = ('invalid matrix_type, only float32, float64, '
			   'complex64, complex128; dtype=%r' % A.dtype)
		raise TypeError(msg)
	return matrix_type, nwords_per_value


def matrix_column_compress(col):
	"""takes a dense matrix column and puts it into OP4 format"""
	packs = []

	n = 0
	i = 0
	packi = []
	while i < len(col):
		#print("i=%s n=%s col[i]=%s" % (i, n, col[i]))
		if col[i] == n + 1:
			#print("i=n=%s" % i)
			packi.append(i)
			n += 1
		else:
			if packi:
				packs.append(packi)
				#print("pack = ", pack)
			packi = [i]
			n = col[i]
		#print("pack = ", pack)
		i += 1

	if packi:
		packs.append(packi)
	#print("packs = ", packs)
	return packs

def DMIG_write(f, name, matrix, form,
			   row_index_to_component_id, col_index_to_component_id,
			   precision='default'):
	"""
	==== ===============
	Form Definition
	==== ===============
	1.   Square
	2.   Rectangular
	3.   Diagonal
	6.   Symmetric
	8.   Id entity
	9.   Pseudoidentity
	==== ===============

	==== =========================
	Type Definition
	==== =========================
	1.   Real, Single Precision
	2.   Real, Double Precision
	3    Complex, Single Precision
	4.   Complex, Double Precision
	==== =========================

	:todo: collapse columns
	:todo: support ndarray
	"""
	assert isinstance(matrix, coo_matrix), 'type(matrix)=%s' % type(matrix)
	matrix = matrix.todense()
	assert precision in ['default', 'single', 'double'], "precison=%r valid=['default','single','double']"
	ifo = form

	#dtype = dtype_get(Type, precision='default').name
	dtype = matrix.dtype.name

	# get the matrix Type
	if dtype == 'float32':
		precision_type = 1
		size = 8
		is_complex = False
		print_card = print_card_8
	elif dtype == 'float64':
		precision_type = 2
		#size = 16
		is_complex = False
		print_card = print_card_double
	elif dtype == 'complex64':
		precision_type = 3
		#size = 8
		is_complex = True
		print_card = print_card_8
	elif dtype == 'complex128':
		precision_type = 4
		#size = 16
		is_complex = True
		print_card = print_card_double
	else:
		raise RuntimeError('dtype = %r' % dtype)

	tin = precision_type
	tout = precision_type

	# real/imaginary vs magnitude/phase
	polar = 0

	ncols = ''
	form = ifo
	if ifo == 9: # Number of columns in a rectangular matrix. Used only for IFO = 9.
		nrows, ncols = A.shape
	f.write('$ type = %s\n' % type(matrix))
	card = ['DMIG', name, 0, ifo, tin, tout, polar, '', ncols]
	f.write(print_card_8(card))

	index_to_node_id = {}
	index_to_component_id = {}

	if isinstance(matrix, coo_matrix):
		if is_complex:
			reals = real(matrix.data)
			imags = imag(matrix.data)
			for i, j, reali, imagi in zip(matrix.row, matrix.col, reals, imags):
				GJ, CJ = col_index_to_node_component_id[i]
				G1, C1 = row_index_to_node_component_id[i]
				card = ['DMIG', name, GJ, CJ, G1, C1, reali, imagi]
				f.write(print_card(card))
		else:
			for i, j, reali in zip(matrix.row, matrix.col, matrix.data):
				GJ, CJ = col_index_to_node_component_id[i]
				G1, C1 = row_index_to_node_component_id[i]
				card = ['DMIG', name, GJ, CJ, G1, C1, reali]
				f.write(print_card(card))
	elif isinstance(matrix, Matrix):
		nrows, ncols = matrix.shape
		#if name == 'EYE5CD':
			#pass
		if is_complex:
			for icol in range(ncols):
				ii = where(matrix[:, icol] != 0.0)[0]
				reals = real(matrix[:, icol])
				imags = imag(matrix[:, icol])
				if ii.shape[1]: #len(ii) > 0:  # Matrix is 2D
					#index_range = [ii.min(), ii.max()]
					GJ, CJ = col_index_to_node_component_id[icol]
					card = ['DMIG', name, GJ, CJ, '']
					ii_max = ii.max()
					ii_min = ii.min()
					if ii_max == ii_min:
						ii_max += 1
					for iirow, irow in enumerate(arange(ii_min, ii_max)):
						G1, C1 = row_index_to_node_component_id[irow]
						card.append(G1)
						card.append(C1)
						card.append(reals[iirow, 0])
						card.append(imags[iirow, 0])
					f.write(print_card(card))
		else:
			for icol in range(ncols):
				ii = where(matrix[:, icol] != 0.0)[0]
				if ii.shape[1]: #len(ii) > 0:  # Matrix is 2D
					#index_range = [ii.min(), ii.max()]
					GJ, CJ = col_index_to_node_component_id[icol]
					card = ['DMIG', name, GJ, CJ, '']
					#for i, j, reali in zip(matrix.row, matrix.col, matrix.data):
					ii_max = ii.max()
					ii_min = ii.min()
					if ii_max == ii_min:
						ii_max += 1
					for iirow, irow in enumerate(arange(ii_min, ii_max)):
						G1, C1 = row_index_to_node_component_id[irow]
						card.append(G1)
						card.append(C1)
						card.append(matrix[irow, icol])  # real number
						card.append('')
					f.write(print_card(card))
	else:
		raise NotImplementedError('type = %s' % type(matrix))
	f.write('$-----------------------------------------------------------------------\n')
		#I = id_to_node_idj]