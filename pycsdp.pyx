from libc.stdlib cimport malloc, free

cdef extern from 'declarations.h':
	
	cdef enum blockcat:
		DIAG
		MATRIX
		PACKEDMATRIX
	
	cdef union blockdatarec:
		double *vec
		double *mat
	
	cdef struct blockrec:
		blockdatarec data
		blockcat blockcategory
		int blocksize
	
	cdef struct blockmatrix:
		int nblocks
		blockrec *blocks
	
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
	
	cdef struct constraintmatrix:
		sparseblock *blocks
	
	cdef int write_prob(char *fname, int n, int k, blockmatrix C, double *a, constraintmatrix *constraints)
	
	cdef void initsoln(int n, int k, blockmatrix C, double *a, constraintmatrix *constraints, blockmatrix *pX0, double **py0, blockmatrix *pZ0)
	
	cdef void free_prob(int n, int k, blockmatrix C, double *a, constraintmatrix *constraints, blockmatrix X, double *y, blockmatrix Z)
	
	cdef int easy_sdp(int n, int k, blockmatrix C, double *a, constraintmatrix *constraints, double constant_offset, blockmatrix *pX, double **py, blockmatrix *pZ, double *ppobj, double *pdobj)

cdef int py_ijtok(i, j, l):
	return ((j-1)*l+i-1)

def cns_block_detail(B):
	n = len(B)
	non_zeros = []
	for i in range(n):
		for j in range(i, n):
			entity = B[i][j]
			if entity != 0:
				non_zeros.append((i,j,entity))
	return non_zeros

CSDP_SOLVER = True

def py_csdp(C, A, b):
	
	## Initiate detailed info
	cdef int gen_int
	cdef double gen_double
	cdef int sol_size
	message = ""
	cdef int num_bloacks
	num_bloacks = len(C)
	cdef int num_constraints
	num_constraints = len(A)
	b_size = len(b)
	block_size = []
	for Cblk in C:
		block_size.append(len(Cblk))
	
	## Check data consistency
	if b_size != num_constraints:
		message = "Inconsistency in size of constraints"
		return {'message':message}
	for Ab in A:
		cns_block_size = []
		for cns_blk in Ab:
			cns_block_size.append(len(cns_blk))
		if cns_block_size != block_size:
			message = "Size of blocks in objective and constraints are not match"
			return {'message':message}
	
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
	## First, allocate storage for the C matrix.  We have three blocks, but
	## because C starts arrays with index 0, we have to allocate space for
	## four blocks- we'll waste the 0th block.  Notice that we check to 
	## make sure that the malloc succeeded.
	
	CC.nblocks = num_bloacks
	CC.blocks = <blockrec *>malloc((num_bloacks+1)*sizeof(blockrec))
	if CC.blocks == NULL:
		message = "Couldn't allocate storage for C"
		return {'message':message}
	
	## Setup blocks.
	for idx in range(num_bloacks):
		CC.blocks[idx+1].blockcategory = MATRIX
		CC.blocks[idx+1].blocksize = block_size[idx]
		CC.blocks[idx+1].data.mat = <double *>malloc(block_size[idx]*block_size[idx]*sizeof(double))
		if CC.blocks[idx+1].data.mat == NULL:
			message = "Couldn't allocate storage for C"
			return {'message':message}
		
		## Put the entries into each block.
		for i in range(block_size[idx]):
			for j in range(block_size[idx]):
				CC.blocks[idx+1].data.mat[py_ijtok(i+1,j+1,block_size[idx])] = C[idx][i][j]
	
	## Allocate storage for the right hand side, b.
	Cb=<double *>malloc((num_constraints+1)*sizeof(double))
	if Cb == NULL:
		message = "Failed to allocate storage for a"
		return {'message':message}
	## Fill in the entries in b.
	for idx in range(b_size):
		Cb[idx+1] = b[idx]
	
	## The next major step is to setup the constraint matrices A_i.
	## Again, because C indexing starts with 0, we have to allocate 
	## space for one more constraint. constraints[0] is not used.
	
	constraints = <constraintmatrix *>malloc((num_constraints+1)*sizeof(constraintmatrix))
	if constraints == NULL:
		message = "Failed to allocate storage for constraints"
		return {'message':message}
	
	for idx in range(num_constraints):
		## Setup the A_i matrices.  Note that we start with last block of A and
		## then proceed backward.  We do this in this order because the blocks 
		## will be inserted into the linked list of each block in reverse order.  
		A[idx].reverse()
		## Terminate the linked list with a NULL pointer.
		constraints[idx+1].blocks = NULL
		## Now, we handle each block
		block_no = num_bloacks
		## Allocate space for block 3 of A1.
		for Atmp in A[idx]:
			blockptr = <sparseblock *>malloc(sizeof(sparseblock))
			if blockptr == NULL:
				message = "Allocation of constraint block failed"
				return {'message':message}
			## Initialize block
			block_detail = cns_block_detail(Atmp)
			num_non_zeros = len(block_detail)
			blockptr.blocknum = block_no
			blockptr.blocksize = block_size[block_no - 1]
			blockptr.constraintnum = idx + 1
			blockptr.next = NULL
			blockptr.nextbyblock = NULL
			blockptr.entries = <double *>malloc((num_non_zeros + 1)*sizeof(double))
			if blockptr.entries == NULL:
				message = "Allocation of constraint block failed"
				return {'message':message}
			blockptr.iindices = <int *>malloc((num_non_zeros+1)*sizeof(int))
			if blockptr.iindices == NULL:
				message = "Allocation of constraint block failed"
				return {'message':message}
			
			blockptr.jindices = <int *>malloc((num_non_zeros+1)*sizeof(int))
			if blockptr.jindices == NULL:
				message = "Allocation of constraint block failed"
				return {'message':message}
			## We have num_non_zeros nonzero entry in the upper triangle of current block.
			blockptr.numentries = num_non_zeros
			## Filling each entry
			for entry in range(num_non_zeros):
				blockptr.iindices[entry+1] = block_detail[entry][0]+1
				blockptr.jindices[entry+1] = block_detail[entry][1]+1
				blockptr.entries[entry+1] = block_detail[entry][2]
			## Insert current block into the linked list of current constraint.
			blockptr.next = constraints[idx+1].blocks
			constraints[idx+1].blocks = blockptr
			block_no -= 1
	
	## At this point, we have all of the problem data setup.
	
	## Write the problem out in SDPA sparse format.
	sol_size = sum(block_size)
	#write_prob("prob.dat-s", sol_size, num_constraints, CC, Cb, constraints)
	
	## Create an initial solution.  This allocates space for X, y, and Z,
	## and sets initial values.
	initsoln(sol_size, num_constraints, CC, Cb, constraints, &X, &y, &Z)
	
	## Solve the problem.
	ret = easy_sdp(sol_size, num_constraints, CC, Cb, constraints, 0.0, &X, &y, &Z, &pobj, &dobj)
	
	if (ret == 0):
		message = "Optimal solution found"
	else:
		message = "SDP failed."
	
	## Write out the problem solution.
	#write_sol("prob.sol",7,2,X,y,Z)
	
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
	free_prob(sol_size, num_constraints, CC, Cb, constraints, X, y, Z)

	return {'dual':dobj, 'primal':pobj, 'X':solX, 'Z':solZ, 'y':soly, 'message':message}
