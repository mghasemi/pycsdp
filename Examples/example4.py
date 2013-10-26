########################################################################
#   This program shows how to use SDP class to solve the 
#   following dual semidefinite program, using formal 
#   expressions with matrix coefficients:
#
#               [ 2  1
#                 1  2
#                       3  0  1           [ X[1]
#   maximize tr(        0  2  0       ) *   X[2]
#                       1  0  3              X[3] ]
#                               0  0
#                               0  0]
#   subject to:
#               [ 3  1
#                 1  3
#                       0  0  0           [ X[1]
#           tr(         0  0  0       ) *   X[2]    = 1
#                       0  0  0             X[3] ]
#                               1  0
#                               0  0 ]
#
#               [ 0  0
#                 0  0
#                       3  0  1           [ X[1]
#           tr(         0  4  0       ) *   X[2]    = 2
#                       1  0  5             X[3] ]
#                               0  0
#                               0  1 ]
#
########################################################################

# import the package
from SDP import *

# initiate a semidefinite program
q = SemidefiniteProgram(primal = False)

# declare a new formal variable for symbolic expressions
x = q.new_variable()

# define the `objective` as a linear combination of x[1] and x[2]
objective = (matrix([[2,1],[1,2]])*x[1] + matrix([[3,0,1],[0,2,0],[1,0,3]])*x[2] + matrix([[0,0],[0,0]])*x[3]).tr()
# set the objective of the semidefinite program
q.set_objective(objective)

# define a formal inequality with matrix coefficients
constraint1 = (matrix([[3,1],[1,3]])*x[1] + matrix([[0,0,0],[0,0,0],[0,0,0]])*x[2] + matrix([[1,0],[0,0]])*x[3]).tr() == 1
# set the first constraint of the program
q.add_constraint(constraint1)

# define anpther formal inequality with matrix coefficients
constraint2 = (matrix([[0,0],[0,0]])*x[1] + matrix([[3,0,1],[0,4,0],[1,0,5]])*x[2] + matrix([[0,0],[0,1]])*x[3]).tr() == 2
# set the second constraint of the program
q.add_constraint(constraint2)

# solve the program
result = q.solve()

# show the details of solution
print result
