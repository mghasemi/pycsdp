
====================
cpysdp & SDP
====================

This module provides basic general features for semidefinite programming
in Sage. PyCsdp is a Cython interface to the famous CSDP solver and SDP
aims to provides a unified interface for semidefinite programming using
various solvers that are available for Sage (so far CvxOpt and CSDP). 
One can include these modules like this::

    #!/usr/bin/env python

    from from cpycsdp import cpycsdp
    from SDP import SDP

	
Installation
====================
In Sage -sh mode, simply run 
	
	python setup.py build_ext install

cpysdp
====================

This is an interface for CSDP written in Cython based on the optional 
Sage compatible installation of CSDP. The main purpose of PyCsdp is to
define the 'easy_sdp' function of CSDP which provides a convenient way
of solving semidefinite programs, for Sage. The Sage version of 
'easy_sdp' is the python 'cpycsdp' function.

* cpycsdp
This function takes three arguments as input::

	cpycsdp(C, A, a)
	
and returns a python dictionary consist of output data of 'easy_sdp'.
This function aims to solve a primal or dual semidefinite program in the
following form:

        Primal:
                
            maximize 
                tr(C*X)
            subject to
                tr(A1*X) = a1
                    .
                    .
                    .
                tr(Am*X) = am
                    X is  psd
                    
        Dual:
                
            minimize 
                a^t*y
            subject to
                A^t(y) - C = Z
                    Z is psd
            where
                A^t(y) = y1 A1 + ... + ym Am
            and
                C = [C1			0
					  0	C2		0
					  0		.	0
					  .		.	.
					  0		.	Cm]


Arguments
--------------------
1. 'C' is a python list C = [C1,..., Cm] where 'C1' .... 'Cm' are 2-dim 
numpy matrices defining the block matrices as are shown above.

2. 'A' is a python list of block matrices [A1,..., Ak] where each 'Ai'
is a python list of the similar structure of 'C'. In other words,
Ai = [Ai1,..., Aim], where 'Aij' is a numpy matrix of the size of 'Cj'.

3. 'a' is 1-dim python list. Simply a=[a1,..., am] as is shown above.

Output
--------------------
Is a python dictionary with the following keys:

1. 'primal'
	the final value of primal objective function.
2. 'dual'
	the final value of dual objective function.
3. 'y'
	the vector y in above problem.
4. 'Z'
	the sdp matrix Z in above problem.
5. 'X'
	the sdp matrix X in above problem.
6. 'message'
	a descriptive text explaining the result of CSDP.


SDP
====================
This module aims to provide a unified interface for various sdp solver
that are supported by Sage. Currently CvxOpt is by default available in
Sage and also CSDP can be installed inside Sage. A Sage installable spkg
is available at <https://github.com/dimpase/csdp.git>.
SDP solves a semidefinite program of the form given above. A typical 
usage of SDP looks like the following:

    #!/usr/bin/env python
    from SDP import SDP
    program = SDP.sdp()
    program.solve(C, a, A)
    
The arguments and output of each method are described below.

* Initialization
SDP.sdp() accepts two optional arguments:

1. solver:
	currently there is two possible choices for solver:
	+ 'cvxopt'
	+ 'csdp'
	that are subject to availability of each one under Sage. The value 
	of silver is preset to 'cvxopt'.
2. Settings:
	a python dictionary for fine tuning of the sdp solver. Acceptable 
	keys are the followings:
	+ 'detail':
		is preset to 'True'. Switch to 'False' in order to hide the 
		progress detail of the solver (Currently only works for cvxopt).
	+ 'Iterations':
		the maximum number of iteration. Default is 100 steps (Currently
		only works for cvxopt).
	+ 'refinement':
		A positive integer number. If the solver encounters a singular 
		KKT matrix, it continues to iterate for this given number of 
		iterations (Currently only works for cvxopt).

* solve
The solve method requires three arguments and sets the 'Info' dictionary
of the class based on the output of the solver. Arguments and outputs 
are described bellow:

Arguments
--------------------
1. 'C' is a list of matrices C1,..., Cm. The following types of data for
	Ci's are acceptable:
	+ Python 2-dimensional list,
	+ Sage matrix,
	+ Numpy matrix,
	+ Blas type matrix that is compatible with CvxOpt.
2. 'a' is a list of scalars.
3. 'A' is a list of block matrices. Any type which is admissible for 'C' 
	is also admissible for 'A'.

Output
--------------------
'SDP.solve' method fills a dictionary 'SDP.Info' of following keys:
	+ 'Status'
		shows the final status of the solver,
	+ 'PObj'
		the final value of primal objective function,
	+ 'DObj'
		the final value of dual objective function,
	+ 'y'
		the vector y in above problem as a numpy vector,
	+ 'Z'
		the sdp matrix Z in above problem as a numpy matrix,
	+ 'X'
		the sdp matrix X in above problem as a numpy matrix,
	+ 'Wall'
		the wall time spent by the solver,
	+ 'CPU'
		the CPU time spent by the solver,

The latest version of this code is available at:
<https://github.com/mghasemi/pycsdp.git>
