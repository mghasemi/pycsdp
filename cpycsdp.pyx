from libc.stdlib cimport malloc, free

import numpy

cimport numpy

import gc

cdef extern from 'declarations.h':
	
	# Each block is a diagonal block or a matrix
	cdef enum blockcat:
		DIAG
		MATRIX
		PACKEDMATRIX
	
	# A block data record contains a pointer to the actual 
	# data for the block. Note that matrices are stored in 
	# Fortran form, while vectors are stored using indices 
	# 1, 2, ..., n.
	cdef union blockdatarec:
		double *vec
		double *mat
	
	# A block record describes an individual block within a 
	# matrix. 
	cdef struct blockrec:
		blockdatarec data
		blockcat blockcategory
		int blocksize
	
	# A block matrix contains an entire matrix in block 
	# diagonal form.
	cdef struct blockmatrix:
		int nblocks
		blockrec *blocks
	
	### Definition for constraint matrices. ###
	
	# There's one of these for each nonzero block in each 
	# constraint matrix.
	cdef struct sparseblock:
		sparseblock *next
		sparseblock *nextbyblock
		double *entries
		int *iindices
		int *jindices
		int numentries
		int blocknum
		int blocksize
		int constraintnum
		int issparse
	
	# A constraint matrix is just an array of pointers to 
	# linked lists of the constraint blocks.
	cdef struct constraintmatrix:
		sparseblock *blocks
		
	# Build an initial solution for an SDP problem.
	cdef void initsoln(int n, int k, blockmatrix C, double *a, constraintmatrix *constraints, blockmatrix *pX0, double **py0, blockmatrix *pZ0)
	
	# Free the memory associated with a problem.
	cdef void free_prob(int n, int k, blockmatrix C, double *a, constraintmatrix *constraints, blockmatrix X, double *y, blockmatrix Z)
	
	# This is an easy to call version of the sdp routine.  
	# It takes as input a problem 
	#		(n, k, C, a, constraints, constant_offset), 
	# and an initial solution (X,y,Z), allocates working 
	# storage, and calls sdp() to solve the problem.  
	# The solution is returned in X,y,Z,pobj,dobj, and the 
	# return code from sdp is returned as the return value 
	# from easy_sdp.   
	cdef int easy_sdp(int n, int k, blockmatrix C, double *a, constraintmatrix *constraints, double constant_offset, blockmatrix *pX, double **py, blockmatrix *pZ, double *ppobj, double *pdobj)

CSDP_SOLVER = True

cdef int py_ijtok(i, j, l):
	"""
	Converts indices from column major form to regular format.
	"""
	return ((j-1)*l+i-1)

cdef cy_block_detail(numpy.ndarray B, int blk_size):
	"""
	Extracts information from a block which is required by
	CSDP solver. 
	"""
	cdef int n, num_nonzero = 0
	cdef int i, i_idx, j_idx
	lst = numpy.array([0.0], dtype=numpy.float64)
	i_ind = numpy.array([0])
	j_ind = numpy.array([0])
	for i in numpy.nonzero(B)[0]:
		i_idx = i % blk_size
		j_idx = i / blk_size
		if j_idx >= i_idx:
			num_nonzero += 1
			i_ind = numpy.append(i_ind, i_idx+1)
			j_ind = numpy.append(j_ind, j_idx+1)
			lst = numpy.append(lst, B[i])
	return [num_nonzero, i_ind, j_ind, lst]

cpdef cpycsdp(numpy.ndarray C, numpy.ndarray A, numpy.ndarray b, numpy.ndarray[long int] d):
	"""
	'C' is the objective multiplier in primal form
	'A' contains blocks of constraints
	'a' is the vector consists of the constants at the 
	    right side of constraints in primal form  
	'd' is a list of integers:
	    d[0] is the number of blocks
	    d[1] is the number of constraints
	    d[i] is the size of each block
	"""
	
	## Initiate detailed info
	cdef int sol_size = 0, idx, i, this_block_start = 0, this_block_end
	cdef int total_block_size = 0, block_start = 0, num_non_zeros
	cdef int num_bloacks
	cdef int num_constraints
	num_bloacks = d[0]
	num_constraints = d[1]
	b_size = b.shape[0] -1 
	block_size = d[2:]
	message = ""
	Cblks = <double *> C.data
	for i in block_size:
		sol_size += i
		total_block_size += i*i
	
	## Check data consistency
	if b_size != num_constraints:
		message = "Inconsistency in size of constraints"
		return {'message':message, 'code':10}
	
	## The problem and solution data.
	cdef blockmatrix CC
	cdef double *Cb
	cdef constraintmatrix *constraints
	
	## Storage for the initial and final solutions.
	cdef blockmatrix X,Z
	cdef double *y
	cdef double pobj,dobj
	
	## blockptr will be used to point to blocks in constraint matrices.
	cdef sparseblock *blockptr
	
	## A return code for the call to easy_sdp().
	cdef int ret
	
	## The first major task is to setup the C matrix and right hand side b.
	## First, allocate storage for the C matrix.  We have d[0] blocks, but
	## because C starts arrays with index 0, we have to allocate space for
	## d[0]+1 blocks- we'll waste the 0th block.  Notice that we check to 
	## make sure that the malloc succeeded.
	
	CC.nblocks = num_bloacks
	try:
		CC.blocks = <blockrec *>malloc((num_bloacks+1)*sizeof(blockrec))
		if CC.blocks == NULL:
			raise Exception("Memory allocation failed")
	except:
		message = "Couldn't allocate storage for 'C'"
		return {'message':message, 'code':10}
	
	## Setup blocks.
	for idx in range(num_bloacks):
		CC.blocks[idx+1].blockcategory = MATRIX
		CC.blocks[idx+1].blocksize = block_size[idx]
		
		try:
			CC.blocks[idx+1].data.mat = <double *>malloc(block_size[idx]*block_size[idx]*sizeof(double))
			if CC.blocks[idx+1].data.mat == NULL:
				raise Exception("Memory allocation failed")
		except:
			message = "Couldn't allocate storage for 'C'"
			return {'message':message, 'code':10}
		
		## Put the entries into each block.
		CC.blocks[idx+1].data.mat = &Cblks[block_start]
		block_start += block_size[idx]*block_size[idx]
	
	## Allocate storage for the right hand side, b.
	try:
		Cb = <double *>malloc((num_constraints+1)*sizeof(double))
		if Cb == NULL:
			raise Exception("Memory allocation failed")
	except:
		message = "Failed to allocate storage for 'a'"
		return {'message':message, 'code':10}
	
	## Fill in the entries in b.
	b_csdp = <double *> b.data
	Cb = b_csdp

	## The next major step is to setup the constraint matrices A_i.
	## Again, because C indexing starts with 0, we have to allocate 
	## space for one more constraint. constraints[0] is not used.
	
	try:
		constraints = <constraintmatrix *>malloc((num_constraints+1)*sizeof(constraintmatrix))
		if constraints == NULL:
			raise Exception("Memory allocation failed")
	except:
		message = "Failed to allocate storage for constraints"
		return {'message':message, 'code':10}
	
	for idx in range(num_constraints):
		## Setup the A_i matrices.  Note that we start with last block of A and
		## then proceed backward.  We do this in this order because the blocks 
		## will be inserted into the linked list of each block in reverse order.  
		## Terminate the linked list with a NULL pointer.
		
		constraints[idx+1].blocks = NULL
		
		## Now, we handle each block
		block_no = num_bloacks
		
		## Allocate space for i'th block of A[i].
		this_block_end = (idx + 1)*total_block_size
		
		while block_no > 0:
			try:
				blockptr = <sparseblock *>malloc(sizeof(sparseblock))
				if blockptr == NULL:
					raise Exception("Memory allocation failed")
			except:
				message = "Allocation of constraint block failed"
				return {'message':message, 'code':10}
			
			## Initialize block
			blockptr.blocknum = block_no
			blockptr.blocksize = block_size[block_no - 1]
			this_block_start = this_block_end - block_size[block_no - 1]*block_size[block_no - 1]
			blockptr.constraintnum = idx + 1
			Blk_Det = cy_block_detail(A[this_block_start:this_block_end], block_size[block_no - 1])
			num_non_zeros = Blk_Det[0]
			blockptr.next = NULL
			blockptr.nextbyblock = NULL
			try:
				blockptr.entries = <double *>malloc((num_non_zeros + 1)*sizeof(double))
				if blockptr.entries == NULL:
					raise Exception("Memory allocation failed")
			except:
				message = "Allocation of constraint block failed"
				return {'message':message, 'code':10}
			
			try:
				blockptr.iindices = <int *>malloc((num_non_zeros+1)*sizeof(int))
				if blockptr.iindices == NULL:
					raise Exception("Memory allocation failed")
			except:
				message = "Allocation of constraint block failed"
				return {'message':message, 'code':10}
			
			try:
				blockptr.jindices = <int *>malloc((num_non_zeros+1)*sizeof(int))
				if blockptr.jindices == NULL:
					raise Exception("Memory allocation failed")
			except:
				message = "Allocation of constraint block failed"
				return {'message':message, 'code':10}
			
			## We have num_non_zeros nonzero entry in the upper triangle of current block.
			blockptr.numentries = num_non_zeros
			
			## Filling each entry
			for entry in range(num_non_zeros+1):
				blockptr.iindices[entry] = Blk_Det[1][entry]
				blockptr.jindices[entry] = Blk_Det[2][entry]
				blockptr.entries[entry] = Blk_Det[3][entry]
			
			## Insert current block into the linked list of current constraint.
			blockptr.next = constraints[idx+1].blocks
			constraints[idx+1].blocks = blockptr
			this_block_end = this_block_start
			block_no -= 1
			
	## At this point, we have all of the problem data setup.
	
	## Create an initial solution.  This allocates space for X, y, and Z,
	## and sets initial values.
	initsoln(sol_size, num_constraints, CC, Cb, constraints, &X, &y, &Z)
	
	## Solve the problem.
	ret = easy_sdp(sol_size, num_constraints, CC, Cb, constraints, 0.0, &X, &y, &Z, &pobj, &dobj)
	
	if (ret == 0):
		message = "Success: SDP solved"
	elif (ret == 1):
		message = "Success: SDP is primal infeasible"
	elif (ret == 2):
		message = "Success: SDP is dual infeasible"
	elif (ret == 3):
		message = "Partial Success: SDP solved with reduced accuracy"
	elif (ret >= 4):
		message = "Failure"
	
	# Translate the results back into numpy format
	solX = []
	solZ = []
	num_blocks = X.nblocks
	for idx in range(num_blocks):
		blk_size = X.blocks[idx+1].blocksize
		solX.append([])
		for i in range(blk_size):
			solX[idx].append([])
			for j in range(blk_size):
				solX[idx][i].append(X.blocks[idx+1].data.mat[py_ijtok(i+1,j+1,blk_size)])

	num_blocks = Z.nblocks
	for idx in range(num_blocks):
		blk_size = Z.blocks[idx+1].blocksize
		solZ.append([])
		for i in range(blk_size):
			solZ[idx].append([])
			for j in range(blk_size):
				solZ[idx][i].append(Z.blocks[idx+1].data.mat[py_ijtok(i+1,j+1,blk_size)])
	soly = []
	for idx in range(b_size):
		soly.append(y[idx+1])
	
	## Free storage allocated for the problem and return.
	for i in range(1, X.nblocks):
		free(X.blocks[i].data.mat);
	for i in range(1, Z.nblocks):
		free(Z.blocks[i].data.mat);
	for idx in range(1, num_constraints):
		ptr=constraints[idx].blocks
		while (ptr != NULL):
			free(ptr.entries)
			free(ptr.iindices)
			free(ptr.jindices)
			oldptr = ptr
			ptr=ptr.next
			free(oldptr)
	free(constraints)
	free(y)
	gc.collect()
	
	# Return the generated results as a Python dictionary.
	return {'dual':dobj, 'primal':pobj, 'X':solX, 'Z':solZ, 'y':soly, 'code':ret, 'message':message}
