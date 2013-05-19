######################################################################
#   This program shows how to use SDP class to solve the 
#   following semidefinite program, using CvxOpt and CSDP
#   with numpy matrices:
#
#   maximize tr(C*X)
#   Subject to:
#           tr(A1*X) = 1
#           tr(A2*X) =-1
#           tr(A3*X) = 1
#           X >= 0 (X is PSD)
#
#   where
#   C = [-33 9
#         9 -26
#               -14 -9 -40
#               -9  -91 -10
#               -40 -10 -15]
#
#   A1 = [7  11
#         11 -3
#               21  11  0
#               11 -10 -8
#                0 -8  -5]
#
#   A2 = [-7 18
#          18 -8
#                0  -10 -16
#               -10  10  10
#               -16  10 -3 ]
#
#   A3 = [2  8
#         8 -1
#                5  -2 17
#               -2   6 -8 
#               17  -8 -6]
#
######################################################################

from SDP import *
import numpy

# Initializing C as 3 block matrices
C = [
    numpy.matrix([[-33,9],
                  [9,-26]]), 
    numpy.matrix([[-14,-9,-40],
                  [-9,-91,-10],
                  [-40,-10,-15]])
]

# Initializing the first constraint
A1 = [
    numpy.matrix([[7, 11],
                  [11, -3]]),
    numpy.matrix([[21, 11, 0],
                  [11, -10, -8],
                  [0, -8, -5]])
]

# Initializing the second constraint
A2 = [
    numpy.matrix([[-7, 18],
                  [18, -8]]),
    numpy.matrix([[0, -10, -16],
                  [-10, 10, 10],
                  [-16, 10, -3]])
]
A3 = [
    numpy.matrix([[2, 8],
                  [8, -1]]),
    numpy.matrix([[5, -2, 17],
                  [-2, 6, -8],
                  [17, -8, -6]])
]

# Packing constraints
A = [A1, A2, A3]

# The right hand side of constraints
b = [1, -1, 1]

# Creating a SDP object for CvxOpt solver
Prob1= SDP.sdp(solver='cvxopt', 
        Settings={'feastol':1e-10, 'detail':True, 'Iterations':150})

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
