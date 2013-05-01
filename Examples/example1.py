######################################################################
#	This program shows how to use SDP class to solve the 
#	following semidefinite program, using CvxOpt and CSDP
#	with python arrays:
#
#	maximize tr(C*X)
#	Subject to:
#			tr(A1*X) = 1
#			tr(A2*X) = 2
#			X >= 0 (X is PSD)
#
#	where
#	C = [2 1
#		 1 2
#			3 0 1
#			0 2 0
#			1 0 3
#				 0 0
#				 0 0]
#
#	A1 = [3 1
#		  1 3
#			 0 0 0
#			 0 0 0
#			 0 0 0
#				  1 0
#				  0 0]
#
#	A2 = [0 0
#		  0 0
#			 3 0 1
#			 0 4 0
#			 1 0 5
#				  0 0
#				  0 1]
#
######################################################################

from SDP import *
import numpy

# Initializing C as 3 block matrices
C = [
	 [[2,1],
	  [1,2]], 
	[[3,0,1],
	 [0,2,0],
	 [1,0,3]], 
	[[0,0],
	 [0,0]]
]

# Initializing the first constraint
A1 = [
	  [[3,1],
	   [1,3]], 
	[[0,0,0],
	 [0,0,0],
	 [0,0,0]], 
	[[1,0],
	 [0,0]]
]

# Initializing the second constraint
A2 = [
	  [[0,0],
	   [0,0]], 
	[[3,0,1],
	 [0,4,0],
	 [1,0,5]], 
	[[0,0],
	 [0,1]]
]

# Packing constraints
A = [A1, A2]

# The right hand side of constraints
b = [1,2]

# Creating a SDP object for CvxOpt solver
Prob1= SDP.sdp(solver='cvxopt')

# Solving the problem
Prob1.solve(C, b, A)

# Creating a SDP object for CSDP solver
Prob2= SDP.sdp(solver='csdp')

# Solving the problem
Prob2.solve(C, b, A)

# Printing primal and dual objective values
print "CvxOpt's result:"
print "Primal Objective=",Prob1.Info['PObj'],", Dual Objective=",Prob1.Info['DObj']
print "CSDP's result:"
print "Primal Objective=",Prob2.Info['PObj'],", Dual Objective=",Prob2.Info['DObj']
