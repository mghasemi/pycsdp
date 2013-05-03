import numpy
import string

from array import array
from time import time, clock

class sdp:
    """
    sdp class:
        
        Provides a simple unified interface for various SDP solvers.
        Currently supports: 
                'cvxopt'
                'csdp'
    """
    
    def __init__(self, solver = 'CvxOpt', Settings = {'detail':True}):
        """
        Creates a generic semidefinite programming object.
        
        Arguments:
            solver:
                Either 'CvxOpt' or 'CSDP' (case insensitive)
            Settings:
                A python dictionary containing detail for fine configuration
                of SDP solver (requires more detail)
        """
            
        solver = string.lower(solver)
        self.SOLVER = solver
        
        self.Info = {}
        self.Options = Settings
        
        try:
            from sage import matrix
            self.SAGE_MATRIX_TYPE=[matrix.matrix_integer_dense.Matrix_integer_dense,
                            matrix.matrix_rational_dense.Matrix_rational_dense,
                            matrix.matrix_generic_dense.Matrix_generic_dense]#,
                            #sage.matrix.matrix_real_double_dense.Matrix_real_double_dense]
            self.SageAvailable = True
        except:
            self.SAGE_MATRIX_TYPE = []
            self.SageAvailable = False
        
        ## Importing the appropriate library based on input
        if solver == 'cvxopt':
            ## Initializing environment variables as required by CvxOpt
            try:
                from cvxopt import solvers
                from cvxopt.base import matrix as Mtx
                RealNumber = float  # Required for CvxOpt
                Integer = int       # Required for CvxOpt
                self.CvxOpt_Available = True
            except:
                self.CvxOpt_Available = False
                print "CVXOPT is not available."
                return
                
        elif solver == 'csdp':
            CSDP_SOLVER = False
            ## Check for availability of CSDP
            try:
                from pycsdp import py_csdp, CSDP_SOLVER
                if not CSDP_SOLVER:
                    print "CSDP is not available."
                    return
                self.Csdp_Available = True
            except ImportError:
                self.Csdp_Available = False
                print "CSDP is not available"
                return
        else:
            print "'"+solver+"' is not available."
            return
        
    def VEC(self, M):
        """
        Converts the matrix M into a column vector
        
        Argument:
            M:
                Any type of matrix
        """
        
        V = []
        Mt = self.matrix_converter(M, 'numpy')
        n,m = Mt.shape
        for j in range(m):
            for i in range(n):
                V.append(Mt[i,j])
        return V
        
    def matrix_converter(self, M, output_type):
        """
        Converts a matrix of a given type to another type.
        
        Arguments:
            M:
                a sage matrix, double indexed list, numpy matrix or cvxopt base matrix
            output_type:
                'sage'       Output will be a sage matrix
                'numpy'      Output will be a numpy matrix
                'list'       Output will be a double indexed python list
                'cvxopt'     Output will be a cvxopt matrix
        """
        
        from cvxopt.base import matrix as Mtx
        
        l_output_type = string.lower(output_type)
        
        input_type = type(M)
        
        if input_type in self.SAGE_MATRIX_TYPE:
            ## Converting Sage matrix to others
            if l_output_type == 'sage':
                return M
                
            elif l_output_type == 'numpy':
                return M.numpy()
                
            elif l_output_type == 'list':
                n=M.ncols()
                m=M.nrows()
                CM=[]
                for j in range(n):
                    RW=[]
                    for i in range(m):
                        RW.append(M[i,j]);
                    CM.append(RW)
                return CM
                
            elif l_output_type == 'cvxopt':
                n=M.ncols()
                m=M.nrows()
                CM=[]
                for j in range(n):
                    for i in range(m):
                        CM.append(M[i,j])
                CC=Mtx(array('d', CM), (m,n))
                return CC
        
        elif input_type is numpy.matrixlib.defmatrix.matrix:
            ## Converting Numpy matrix to others
            if l_output_type == 'sage':
                if self.SageAvailable:
                    return Matrix(M)
                else:
                    return M
                
            elif l_output_type == 'numpy':
                return M
                
            elif l_output_type == 'list':
                return M.tolist()
                
            elif l_output_type == 'cvxopt':
                n=M.shape[1]
                m=M.shape[0]
                CM=[]
                for j in range(n):
                    for i in range(m):
                        CM.append(M[i,j])
                CC=Mtx(array('d', CM), (m,n))
                return CC
        
        elif input_type is list:
            ## Converting Python list to others
            if l_output_type == 'sage':
                if self.SageAvailable:
                    return Matrix(M)
                else:
                    return M
                
            elif l_output_type == 'numpy':
                return numpy.matrix(M)
                
            elif l_output_type == 'list':
                return M
                
            elif l_output_type == 'cvxopt':
                n=len(M[0])
                m=len(M)
                CM=[]
                for j in range(n):
                    for i in range(m):
                        CM.append(M[i][j])
                CC=Mtx(array('d', CM), (m,n))
                return CC
        elif input_type is type(Mtx()):
            ## Converting CvxOpt matrix to others
            if l_output_type == 'sage':
                if self.SageAvailable:
                    m,n = M.size
                    return matrix(n,list(M)).transpose()
                else:
                    return M
            elif l_output_type == 'numpy':
                m,n = M.size
                return numpy.array(list(M)).reshape(m,n)
            elif l_output_type == 'list':
                m,n = M.size
                return matrix(n,list(M)).transpose().tolist()
            elif l_output_type == 'cvxopt':
                return M
        
    def solve(self, C, a, A):
        """
        Solves a semidefinite programming problem with entered data.
        
        Arguments:
            C:
                a list of square matrices (sage matrix, numpy matrix, 
                double indexed python list or cvxopt matrix)
            a:
                a list of real numbers
            A:
                list of block constraints
                
        Example:
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
                C = [C1,..., Cm]
        
        Returns a dictionary containing following data:
            
            'Status':
                either 'Optimal' or 'Infeasible'
            'PObj':
                value of primal objective function
            'DObj':
                value of dual objective function
            'X':
                the psd matrix X in the primal problem
            'Z':
                the psd matrix Z in the dual problem
            'y':
                the vector y for the dual problem
            'Wall':
                Wall time spent by sdp solver (in seconds)
            'CPU':
                CPU time spent by sdp solver (in seconds)
        
        The output matrix are in numpy format
        """
        
        from time import time, clock
        
        self.num_constraints = len(A)
        self.num_blocks = len(C)
        
        ## setting up the data for cvxopt solver
        if self.SOLVER == 'cvxopt':
            if not self.CvxOpt_Available:
                print "CvxOpt is not available."
                return
            
            from cvxopt import solvers
            Cns = []
            for idx in range(self.num_constraints):
                Cns.append([])
            Acvxopt = []
            Ccvxopt = []
            acvxopt = []
            for M in C:
                Ccvxopt.append(-self.matrix_converter(M,'cvxopt'))
            for blk_no in range(self.num_blocks):
                Ablock = []
                for Cns in A:
                    Ablock.append(self.VEC(Cns[blk_no]))
                Acvxopt.append(-self.matrix_converter(numpy.matrix(Ablock).transpose(), 'cvxopt'))
            aTranspose=[]
            for elmnt in a:
                aTranspose.append([elmnt])
            acvxopt = self.matrix_converter(aTranspose, 'cvxopt')
            
            if 'detail' in self.Options:
                solvers.options['show_progress'] = self.Options['detail']
            else:
                solvers.options['show_progress'] = False
            
            if 'Iterations' in self.Options:
                solvers.options['maxiters'] = max(1, self.Options['Iterations'])
            else:
                solvers.options['maxiters'] = 100
            if 'refinement' in self.Options:
                solvers.options['refinement'] = max(1, self.Options['refinement'])
            else:
                solvers.options['refinement'] = 1
            start1 = time()
            start2= clock()
            
            try:
            #if True:
                sol = solvers.sdp(acvxopt, Gs = Acvxopt, hs = Ccvxopt)
                elapsed1 = (time() - start1)
                elapsed2 = (clock() - start2)
                if sol['status'] != 'optimal':
                    self.Info={'Status':'Infeasible'}
                else:
                    self.Info = {'Status':'Optimal', 'DObj':sol['dual objective'],
                    'PObj':sol['primal objective'], 'Wall':elapsed1, 'CPU':elapsed2}
                    self.Info['y'] = self.matrix_converter(sol['x'], 'numpy')
                    self.Info['Z'] = []
                    for ds in sol['ss']:
                        self.Info['Z'].append(self.matrix_converter(ds, 'numpy'))
                    self.Info['X'] = []
                    for ds in sol['zs']:
                        self.Info['X'].append(self.matrix_converter(ds, 'numpy'))
            except:
                self.Info={'Status':'Infeasible'}
        
        ## setting up the data for csdp solver
        elif self.SOLVER == 'csdp':
            if not self.Csdp_Available:
                print "CSDP is not available."
                return
            from pycsdp import py_csdp
            admissible_types = [type(numpy.matrix([])), type(list())]
            for idx in range(self.num_blocks):
				C[idx] = self.matrix_converter(C[idx], 'list')
            for i in range(self.num_constraints):
                for idx in range(self.num_blocks):
                    A[i][idx] = self.matrix_converter(A[i][idx], 'list')
            b = self.matrix_converter(a,'list')
            start1 = time()
            start2= clock()
            try:
                sol = py_csdp(C, A, b)
                elapsed1 = (time() - start1)
                elapsed2 = (clock() - start2)
                if sol['message'] != 'Optimal solution found':
                    self.Info={'Status':sol['message']}
                else:
                    self.Info = {'Status':'Optimal', 'DObj':sol['dual'],
                    'PObj':sol['primal'], 'Wall':elapsed1, 'CPU':elapsed2}
                    self.Info['y'] = self.matrix_converter(sol['y'], 'numpy')
                    self.Info['Z'] = []
                    for ds in sol['Z']:
                        self.Info['Z'].append(self.matrix_converter(ds, 'numpy'))
                    self.Info['X'] = []
                    for ds in sol['X']:
                        self.Info['X'].append(self.matrix_converter(ds, 'numpy'))
            except:
                self.Info={'Status':'Infeasible'}
        
        self.Info['solver'] = self.SOLVER
