########################################################################
#   This program shows how to use SDP class to solve the 
#   following primal semidefinite program, using formal 
#   expressions with matrix coefficients:
#
#   minimize y[1] + 2*y[2]
#   subject to:
#       [ 3   1             [ 2   1 
#         1   3 ] * y[1] >=   1   2 ]
#
#       [ 3   0   1             [ 3   0   1 
#         0   4   0               0   2   0
#         1   0   5 ] * y[2] >=   1   0   3 ]
#
#       [ 1   0             [ 0   0             [ 0   0
#         0   1 ] * y[1] +    0   1 ] * y[2] >=   0   0 ]
#
########################################################################

# import the package
from SDP import *

# initiate a semidefinite program
p = SemidefiniteProgram(primal = True)

# declare a new formal variable for symbolic expressions
y = p.new_variable()

# define the `objective` as a linear combination of y[1] and y[2]
objective = y[1]+2*y[2]
# set the objective of the semidefinite program
p.set_objective(objective)

# define a formal inequality with matrix coefficients
constraint1 = matrix([[3,1],[1,3]])*y[1] + matrix([[0,0],[0,0]])*y[2]>=matrix([[2,1],[1,2]])
# set the first constraint of the program
p.add_constraint(constraint1)

# define anpther formal inequality with matrix coefficients
constraint2 = matrix([[0,0,0],[0,0,0],[0,0,0]])*y[1] + matrix([[3,0,1],[0,4,0],[1,0,5]])*y[2]>=matrix([[3,0,1],[0,2,0],[1,0,3]])
# set the second constraint of the program
p.add_constraint(constraint2)

# yet another formal inequality with matrix coefficients
constraint3 = matrix([[1,0],[0,0]])*y[1] + matrix([[0,0],[0,1]])*y[2]>=matrix([[0,0],[0,0]])
# set the third constraint of the program
p.add_constraint(constraint3)

# solve the program
result = p.solve()

# show the details of solution
print result
