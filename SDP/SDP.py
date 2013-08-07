import numpy
import string
from scipy.linalg import block_diag
from sage.matrix.constructor import Matrix
from sage.rings.real_mpfr import *
from sage.rings.integer import *

from array import array
from time import time, clock

##########################  ##########################
def matrix_converter(M, output_type):
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
    
    try:
        from sage import matrix
        SAGE_MATRIX_TYPE=[matrix.matrix_integer_dense.Matrix_integer_dense,
                        matrix.matrix_rational_dense.Matrix_rational_dense,
                        matrix.matrix_generic_dense.Matrix_generic_dense]#,
                        #sage.matrix.matrix_real_double_dense.Matrix_real_double_dense]
        SageAvailable = True
    except:
        SAGE_MATRIX_TYPE = []
        SageAvailable = False
    
    l_output_type = string.lower(output_type)
    
    input_type = type(M)
    
    if input_type in SAGE_MATRIX_TYPE:
        ## Converting Sage matrix to others
        if l_output_type == 'sage':
            return M
            
        elif l_output_type == 'numpy':
            return M.numpy()
            
        elif l_output_type == 'list':
            return list(M)
            
        elif l_output_type == 'cvxopt':
            return Mtx(M.numpy(), tc='d')
    
    elif (input_type is numpy.matrixlib.defmatrix.matrix) or (input_type is numpy.ndarray):
        ## Converting Numpy matrix to others
        if l_output_type == 'sage':
            if SageAvailable:
                return Matrix(M.tolist())
            else:
                return M
            
        elif l_output_type == 'numpy':
            return M
            
        elif l_output_type == 'list':
            return M.tolist()
            
        elif l_output_type == 'cvxopt':
            return Mtx(M, tc='d')
    
    elif input_type is list:
        ## Converting Python list to others
        if l_output_type == 'sage':
            if SageAvailable:
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
            return Mtx(numpy.array(M).reshape(m*n, order='F'), size=(m, n), tc='d')
    
    elif input_type is type(Mtx()):
        ## Converting CvxOpt matrix to others
        if l_output_type == 'sage':
            if SageAvailable:
                m,n = M.size
                return Matrix(n,list(M)).transpose()
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

def is_constant(c):
    """
    Determines whether the input argument is a constant number or not.
    """
    
    if (type(c) is sage.rings.real_mpfr.RealLiteral) or (type(c) is sage.rings.integer.Integer):
        return True
    return False

def MatlabBlockMat(M):
    """
    Generates a code block to define a MATLAB block matrix.
    
    Argument:
        'M':
            is a python list of square matrices.
            Matrices can be in different formats:
                Sage matrix;
                Numpy matrix;
                Double indexed python list
                CvxOpt matrix
    """
    
    from sage.interfaces.matlab import Matlab
    instance = Matlab()
    n = len(M)
    code = "blkdiag("
    for idx in range(n):
        #print matrix_converter(M[idx], 'sage')
        code += instance.sage2matlab_matrix_string(matrix_converter(M[idx], 'sage'))
        if idx != (n-1):
            code += ", "
    code += ")"
    return code

########################## The generic sdp solver interfacing various solvers ########################## 

class sdp:
    """
    sdp class:
        
        Provides a simple unified interface for various SDP solvers.
        Currently supports: 
                'cvxopt'
                'csdp'
        If MATLAB is installed and it is added to system PATH upon 
        availability:
                'SeDuMi'
                'SDPNAL'
                
    """
    
    def __init__(self, solver = 'CvxOpt', Settings = {'detail':True}):
        """
        Creates a generic semidefinite programming object.
        
        Arguments:
            solver:
                Either 'CvxOpt', 'CSDP', 'SeDuMi' or 'SDPNAL' (case insensitive)
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
                from cpycsdp import cpycsdp, CSDP_SOLVER
                if not CSDP_SOLVER:
                    print "CSDP is not available."
                    return
                self.Csdp_Available = True
            except ImportError:
                self.Csdp_Available = False
                print "CSDP is not available"
                return
        
        elif solver == 'sedumi':
            from sage.interfaces.matlab import Matlab
            try:
                interface = Matlab()
                path_str = interface.get('path')
                self.MATLAB_Available = True
                if path_str.find('SeDuMi') >= 0:
                    self.SeDuMi_Available = True
                else:
                    self.SeDuMi_Available = False
                interface.eval("quit;")
            except:
                self.MATLAB_Available = False
                self.SeDuMi_Available = False
                print "MATLAB SeDuMi is not available"
                return
        
        elif solver == 'sdpnal':
            from sage.interfaces.matlab import Matlab
            try:
                interface = Matlab()
                path_str = interface.get('path')
                self.MATLAB_Available = True
                if path_str.find('SDPNAL') >= 0:
                    self.SDPNAL_Available = True
                else:
                    self.SDPNAL_Available = False
                interface.eval("quit;")
            except:
                self.MATLAB_Available = False
                self.SDPNAL_Available = False
                print "MATLAB SDPNAL is not available"
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
        Mt = matrix_converter(M, 'numpy')
        n,m = Mt.shape
        for j in range(m):
            for i in range(n):
                V.append(Mt[i,j])
        return V
        
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
                
            minimize 
                a^t*y
            subject to
                A^t(y) - C = Z
                    Z is psd
            where
                A^t(y) = y1 A1 + ... + ym Am
            and
                C = [C1,..., Cm]

            Dual:
                
            maximize 
                tr(C*X)
            subject to
                tr(A1*X) = a1
                    .
                    .
                    .
                tr(Am*X) = am
                    X is  psd
        
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
                Ccvxopt.append(-matrix_converter(M,'cvxopt'))
            for blk_no in range(self.num_blocks):
                Ablock = []
                for Cns in A:
                    Ablock.append(self.VEC(Cns[blk_no]))
                Acvxopt.append(-matrix_converter(numpy.matrix(Ablock).transpose(), 'cvxopt'))
            aTranspose=[]
            for elmnt in a:
                aTranspose.append([elmnt])
            acvxopt = matrix_converter(aTranspose, 'cvxopt')
            
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
                    self.Info['y'] = matrix_converter(sol['x'], 'numpy')
                    self.Info['Z'] = []
                    for ds in sol['ss']:
                        self.Info['Z'].append(matrix_converter(ds, 'numpy'))
                    self.Info['X'] = []
                    for ds in sol['zs']:
                        self.Info['X'].append(matrix_converter(ds, 'numpy'))
            except:
                self.Info={'Status':'Infeasible'}
        
        ## setting up the data for csdp solver
        elif self.SOLVER == 'csdp':
            if not self.Csdp_Available:
                print "CSDP is not available."
                return
            for idx in range(self.num_blocks):
				C[idx] = matrix_converter(C[idx], 'numpy')
            for i in range(self.num_constraints):
                for idx in range(self.num_blocks):
                    A[i][idx] = matrix_converter(A[i][idx], 'numpy')
            b = list(a)
            start1 = time()
            start2= clock()
            try:
            #if True:
                sol = self.cpycsdp(C, b, A)
                elapsed1 = (time() - start1)
                elapsed2 = (clock() - start2)
                if sol['code'] >= 4:
                    status = 'Infeasible'
                else:
                    status = 'Optimal'
                self.Info = {'Status':status, 'DObj':sol['dual'],
                    'PObj':sol['primal'], 'Wall':elapsed1, 'CPU':elapsed2}
                self.Info['y'] = matrix_converter(sol['y'], 'numpy')
                self.Info['Z'] = []
                for ds in sol['Z']:
                    self.Info['Z'].append(matrix_converter(ds, 'numpy'))
                self.Info['X'] = []
                for ds in sol['X']:
                    self.Info['X'].append(matrix_converter(ds, 'numpy'))
            except:
                self.Info={'Status':'Infeasible'}
        
        ## setting up the data for SeDuMi solver
        elif self.SOLVER == 'sedumi':
            if not self.SeDuMi_Available:
                print "SeDuMi is not available."
                return
            from sage.interfaces.matlab import Matlab
            from math import sqrt, floor
            instance = Matlab()
            matlab_code = "c = " + instance.sage2matlab_matrix_string(matrix_converter(a, 'sage')) + ";"
            instance.eval(matlab_code)
            matlab_code = "p = length(c);"
            instance.eval(matlab_code)
            matlab_code = "bt = -c;"
            instance.eval(matlab_code)
            matlab_code = "F = -" + MatlabBlockMat(C) + ";"
            instance.eval(matlab_code)
            matlab_code = "ct = vec(F);"
            instance.eval(matlab_code)
            for i in range(self.num_constraints):
                matlab_code = "At(:," + str(i+1) + ") = -vec(" + MatlabBlockMat(A[i]) + ");"
                instance.eval(matlab_code)
            matlab_code = "K.s = size(F,1);"
            instance.eval(matlab_code)
            matlab_code = "[x, y, info] = sedumi(At,bt,ct,K);"
            instance.eval(matlab_code)
            x_size = floor(sqrt(eval(instance.get("length(x)"))))
            y_size = eval(instance.get("length(y)"))
            
            self.Info['X'] = []
            i = 0
            while i < x_size:
                j = 0
                row = []
                while j < x_size:
                    row.append(eval(instance.get("x(" + str(int(i*x_size + j +1)) + ")")))
                    j += 1
                i += 1
                self.Info['X'].append(row)
            self.Info['X'] = matrix_converter(self.Info['X'], 'numpy')
            
            self.Info['y'] = []
            i = 0
            while i < y_size:
                self.Info['y'].append(eval(instance.get("y(" + str(i+1) + ")")))
                i += 1
            
            self.Info['PObj'] = eval(instance.get("-sum(bt'.*y)"))
            self.Info['CPU'] = eval(instance.get("info.cpusec"))
            self.Info['Wall'] = eval(instance.get("info.wallsec"))
            matlab_code = "quit;"
            instance.eval(matlab_code)
            self.Info['DObj'] = None
            self.Info['Z'] = []
            self.Info['Status'] = 'unknown'
        
        ## setting up the data for SDPNAL solver
        elif self.SOLVER == 'sdpnal':
            if not self.SDPNAL_Available:
                print "SDPNAL is not available."
                return
            from sage.interfaces.matlab import Matlab
            instance = Matlab()
            matlab_code = "c = " + instance.sage2matlab_matrix_string(matrix_converter(a, 'sage')) + ";"
            instance.eval(matlab_code)
            matlab_code = "p = length(c);"
            instance.eval(matlab_code)
            matlab_code = "bt = -c;"
            instance.eval(matlab_code)
            matlab_code = "F = -" + MatlabBlockMat(C) + ";"
            instance.eval(matlab_code)
            matlab_code = "ct = vec(F);"
            instance.eval(matlab_code)
            for i in range(self.num_constraints):
                matlab_code = "At(:," + str(i+1) + ") = -vec(" + MatlabBlockMat(A[i]) + ");"
                instance.eval(matlab_code)
            matlab_code = "K.s = size(F,1);"
            instance.eval(matlab_code)
            matlab_code = "opts = []; opts.tol = 1e-7;"
            instance.eval(matlab_code)
            matlab_code = "[blk,At,C,B] = read_sedumi(At,bt,ct,K);"
            instance.eval(matlab_code)
            matlab_code = "[obj,X,y,Z,info,runhist] = sdpnal(blk,At,C,B,opts);"
            instance.eval(matlab_code)

            x_size = eval(instance.get("length(X{1,1})"))
            y_size = eval(instance.get("length(y)"))
            
            self.Info['X'] = []
            i = 0
            while i < x_size:
                j = 0
                row = []
                while j < x_size:
                    row.append(eval(instance.get("X{1,1}(" + str(int(i*x_size + j +1)) + ")")))
                    j += 1
                i += 1
                self.Info['X'].append(row)
            self.Info['X'] = matrix_converter(self.Info['X'], 'numpy')
            
            self.Info['y'] = []
            i = 0
            while i < y_size:
                self.Info['y'].append(eval(instance.get("y(" + str(i+1) + ")")))
                i += 1
            
            self.Info['Z'] = []
            i = 0
            while i < x_size:
                j = 0
                row = []
                while j < x_size:
                    row.append(eval(instance.get("Z{1,1}(" + str(int(i*x_size + j +1)) + ")")))
                    j += 1
                i += 1
                self.Info['Z'].append(row)
            self.Info['Z'] = matrix_converter(self.Info['Z'], 'numpy')
            
            self.Info['PObj'] = eval(instance.get("-obj(1)"))
            self.Info['DObj'] = eval(instance.get("-obj(2)"))
            matlab_code = "quit;"
            instance.eval(matlab_code)
            self.Info['CPU'] = None
            self.Info['Wall'] = None
            self.Info['Status'] = 'unknown'
        
        self.Info['solver'] = self.SOLVER
        
    def cpycsdp(self, C, a, A):
        """
        Translates the input sdp into the format acceptable by the Cython module
        cpycsdp. One should avoid using this method individually, but instead call
        the solver method for 'csdp' solver. 
        """
        
        from cpycsdp import cpycsdp
        d=[]
        # Number of blocks
        d.append(len(C))
        # Number of constraints
        d.append(len(A))
        flag = 0
        a = numpy.insert(a, 0, 0)
        a_csdp = numpy.array(a, dtype = numpy.float64)
        for Cblk in C:
            d.append(len(Cblk))
            temp_block = numpy.matrix(Cblk)
            if flag == 0:
                C_csdp = numpy.array(numpy.reshape(temp_block, temp_block.shape[0]*temp_block.shape[1], order='F'), dtype=numpy.float64)
            else:
                C_csdp = numpy.append(C_csdp, numpy.array(numpy.reshape(temp_block, temp_block.shape[0]*temp_block.shape[1], order='F'), dtype=numpy.float64))
            flag += 1
        d_csdp = numpy.array(d)
        flag = 0
        for Cns in A:
            for cns_blk in Cns:
                temp_block = numpy.matrix(cns_blk)
                if flag == 0:
                    A_csdp = numpy.array(numpy.reshape(temp_block, temp_block.shape[0]*temp_block.shape[1], order='F'))
                else:
                    A_csdp = numpy.append(A_csdp, numpy.array(numpy.reshape(temp_block, temp_block.shape[0]*temp_block.shape[1], order='F')))
                flag+=1
        return cpycsdp(C_csdp, A_csdp, a_csdp, d_csdp)
        
########################## Block Matrix class ##########################

class BlockMat:
    """
    BlockMat class:
        
        A simple interface which implements simple block matrix construction
        with different type of matrices in SAGE.
        
    """
        
    def __init__(self, M = []):
        """
        Initiates a block diagonal matrix based on input
        
        Argument:
            'M':
                A python list of square matrices as blocks.
                Matrices can be in different formats:
                    Sage matrix;
                    Numpy matrix;
                    Double indexed python list;
                    CvxOpt matrix.
        """
        
        if type(M) is not list:
            ErrorString = "The input should be a list of matrices or numbers."
            print ErrorString
            return
        self.Blocks = []
        self.Type = []
        for A in M:
            self.add_block(A)
        
    def add_block(self, M):
        """
        Appends a square block 'M' at the end of current block matrix.
        
        Argument:
            'M':
                A square matrix of the following types is acceptable:
                    Sage matrix;
                    Numpy matrix;
                    Double indexed python list;
                    CvxOpt matrix.
        """
        
        if is_constant(M):
            self.Blocks.append([M])
            self.Type.append(1)
        else:
            B = matrix_converter(M, 'numpy')
            if B.shape[0] != B.shape[1]:
                ErrorString = "All blocks should be square."
                print ErrorString
                return
            self.Blocks.append(B)
            self.Type.append(B.shape[0])
        self.TotalSize = sum(self.Type)
        self.NumBlocks = len(self.Type)
    
    def Matrix2Block(self, M, typ):
        """
        Split a matrix into blocks of square matrices.
        
        Arguments:
            'M':
                A block diagonal matrix of the following type:
                    Sage matrix;
                    Numpy matrix;
                    Double indexed python list;
                    CvxOpt matrix.
            'typ':
                A list on positive integers which determines the size of each blocks.
        """
        
        self.Type = typ
        self.TotalSize = sum(typ)
        self.NumBlocks = len(typ)
        begin = 0
        for end in typ:
            self.Blocks.append(numpy.matrix(M[begin:begin+end, begin:begin+end]))
            begin += end
    
    def Block2Matrix(self):
        """
        Constructs a single block diagonal matrix in 'numpy' format from class data.
          
        """
        
        tmp = numpy.matrix(block_diag(*self.Blocks))
        return tmp
    
########################## Semidefinite Program class ##########################

class SemidefiniteProgram:
    """
    The 'SemidefiniteProgram' class provides a user friendly interface to solve
    primal or dual semidefinite programming (sdp) problems within Sage.
    
    A typical semidefinite program in primal form is the following:  
            
            minimize 
                a^t*y
            subject to
                y_1 A_{11} + ... + y_n A_{1n} - C_1 = Z_1
                    .    .
                    .    .
                    .    .
                y_1 A_{m1} + ... + y_n A_{mn} - C_m = Z_m
            where
                Z = [Z_1 0 ... 0
                     0 Z_2 0...0
                     ...   ...
                     ...   ...
                     0 ... 0 Z_m]
                                 is psd.

    A typical semidefinite program in dual form is the following:
                
            maximize 
                tr(C*X)
            subject to
                tr(A_1 * X) = a_1
                    .    .
                    .    .
                    .    .
                tr(A_m * X) = a_m
            where
                    X is  psd.
            and each A_i is a symmetric block diagonal matrix 

    Currently supports 4 solvers subject to availability:
    'CvxOpt':
        The default solver which comes with Sage.
        see http://cvxopt.org/
    'csdp':
        If CSDP is installed for Sage.
        see https://projects.coin-or.org/Csdp/
    'SeDuMi':
        Available through MATLAB, if MATLAB is installed and it 
        is added to system PATH.
        see http://sedumi.ie.lehigh.edu/
    'SDPNAL'
        Available through MATLAB, if MATLAB is installed and it 
        is added to system PATH.
        see http://www.math.nus.edu.sg/~mattohkc/SDPNAL.html
    """
    
    def __init__(self, primal = True, solver = 'cvxopt'):
        """
        Constructor for 'SemidefiniteProgram' class.
        
        Arguments:
            'primal':
                A boolean value, True or False to indicate type of the sdp. 
            
            'solver':
                Currently supports 4 solvers subject to availability:
                'CvxOpt':
                    The default solver which comes with Sage.
                    see http://cvxopt.org/
                'csdp':
                    If CSDP is installed for Sage.
                    see https://projects.coin-or.org/Csdp/
                'SeDuMi':
                    Available through MATLAB, if MATLAB is installed and it 
                    is added to system PATH.
                    see http://sedumi.ie.lehigh.edu/
                'SDPNAL'
                    Available through MATLAB, if MATLAB is installed and it 
                    is added to system PATH.
                    see http://www.math.nus.edu.sg/~mattohkc/SDPNAL.html
        """
        
        self.Primal = primal
        self.Solver = solver.lower()
        self.Interface = sdp(self.Solver)
        self.ConstraintsSize = []
        self.NumConstraints = 0
        self.HasBlockMat = False
        self.BlockMatType = []
        self.a = []
        self.C = BlockMat()
        self.A = []
        
    def set_objective(self, M):
        """
        Sets the objective of the 'SemidefiniteProgram'.
        
        Argument:
            'M':
                Primal:
                    A vector of real numbers ('a').
                Dual:
                    A sparse square matrix ('C') in one of the following types:
                        BlockMat;
                        Sage matrix;
                        Numpy matrix;
                        Double indexed python list;
                        CvxOpt matrix.
        """
        
        ## Primal case:
        if self.Primal:
            converted_M = matrix_converter(M, 'numpy')
            self.ObjectiveSize = converted_M.shape
            if (self.ObjectiveSize[0] != 1):
                ErrorString = "A semidefinite program in primal form requires a vector for objective function:\n"
                ErrorString += "Primal objective function format: a^T * X"
                raise ValueError(ErrorString)
            self.NumConstraints = self.ObjectiveSize[1]
            self.a = converted_M[0].tolist()[0]
        ## Dual case:
        else:
            if isinstance(M, BlockMat):
                self.C = M
                self.HasBlockMat = True
                self.BlockMatType = M.Type
            else:
                converted_M = matrix_converter(M, 'numpy')
                self.ObjectiveSize = converted_M.shape
                if (self.ObjectiveSize[0] != self.ObjectiveSize[1]):
                    ErrorString = "A semidefinite program in dual format requires a square matrix to form the objective function:\n"
                    ErrorString += "Dual objective funation format: tr(C * X)"
                    raise ValueError(ErrorString)
                if self.HasBlockMat:
                    C.Matrix2Block(converted_M, self.BlockMatType)
                else:
                    self.C = [converted_M]
        ##############################################################
        
    def add_constraint(self, M1, M2):
        """
        Adds a constraint to the 'SemidefiniteProgram'.
        
        Arguments:
            Primal:
                'M1':
                    A python list of equi-size square matrices ('A_{i1},..., A_{in}') 
                    in either of the following formats:
                        Sage matrix;
                        Numpy matrix;
                        Double indexed python list;
                        CvxOpt matrix.
                'M2':
                    A square matrix ('C_i') in a format similar to 'M1'.
                    
            Dual:
                'M1':
                    A square or diagonal block matrix ('A_i') in one of the following types:
                        BlockMat;
                        Sage matrix;
                        Numpy matrix;
                        Double indexed python list;
                        CvxOpt matrix.
                'M2':
                    A constant real number ('a_i').
                    
        """
        
        ## SDP in primal form
        if self.Primal:
            ## Check the input
            ### First argument
            if not isinstance(M1, list):
                ErrorString = "A constraint in primal semidefinite program is of the form y1.A1 + ... + yn.An >= C\n"
                ErrorString += "The first argument should be a list of square matrices."
                raise ValueError(ErrorString)
            ### Second argument
            if is_constant(M2):
                ErrorString = "A constraint in primal semidefinite program is of the form y1.A1 + ... + yn.An >= C\n"
                ErrorString += "The second argument should be a square matrix."
                raise ValueError(ErrorString)
            N = len(M1)
            if self.A == []:
                for _ in range(N):
                    self.A.append(BlockMat())
            shp = None
            for i in range(N):
                M = M1[i]
                converted_M = matrix_converter(M, 'numpy')
                if shp == None:
                    shp = converted_M.shape
                    if shp[0] != shp[1]:
                        ErrorString = "Non-square matrix is given."
                        raise ValueError(ErrorString)
                else:
                    if shp != converted_M.shape:
                        ErrorString = "All matrices in the first argument should of the same size."
                        raise ValueError(ErrorString)
                self.A[i].add_block(converted_M)
            
            converted_M2 = matrix_converter(M2, 'numpy')
            if shp != converted_M2.shape:
                ErrorString = "The size of second argument is different from the first."
                raise ValueError(ErrorString)
                
            self.ConstraintsSize.append(converted_M.shape[0])
            self.C.add_block(converted_M2)
        
        ## SDP in dual form
        else:
            if not is_constant(M2):
                ErrorString = "A constraint in dual semidefinite program is of the form tr(A*X) = a\n"
                ErrorString += "The second argument should be a real number."
                raise ValueError(ErrorString)
            if isinstance(M1, BlockMat):
                if not self.HasBlockMat:
                    self.HasBlockMat = True
                    self.BlockMatType = M1.Type
                    for i in range(len(self.A)):
                        temp_block_mat = BlockMat()
                        temp_block_mat.Matrix2Block(self.A[i][0], self.BlockMatType)
                        self.A[i] = temp_block_mat
                    if not isinstance(self.C, BlockMat):
                        temp_block_mat = BlockMat()
                        temp_block_mat.Matrix2Block(self.C[0], self.BlockMatType)
                        self.C = temp_block_mat
                self.A.append(M1)
            else:
                converted_M1 = matrix_converter(M1, 'numpy')
                if converted_M1.shape[0] != converted_M1.shape[1]:
                    ErrorString = "A constraint in dual semidefinite program is of the form tr(A*X) = a\n"
                    ErrorString += "'A' should be either a square matrix or a block matrix."
                    raise ValueError(ErrorString)
                if self.HasBlockMat:
                    temp_block_mat = BlockMat()
                    temp_block_mat.Matrix2Block(converted_M1, self.BlockMatType)
                    self.A.append(temp_block_mat)
                else:
                    self.A.append([converted_M1])
                    self.ConstraintsSize.append(converted_M1.shape[0])
            self.a.append(M2)
        return
    
    def solve(self):
        """
        Solves the 'SemidefiniteProgram'.
        
        Argument:
            None.
        
        Output:
            Returns a dictionary 'SemidefiniteProgram.Result' with following keys:
                
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
        
        A_Blocks = []
        
        ## SDP in primal form:
        if self.Primal:
            for Obj in self.A:
                A_Blocks.append(Obj.Blocks)
            self.Interface.solve(self.C.Blocks, self.a, A_Blocks)
        
        ## SDP in dual form:
        else:
            if self.HasBlockMat:
                for Obj in self.A:
                    A_Blocks.append(Obj.Blocks)
                self.Interface.solve(self.C.Blocks, self.a, A_Blocks)
            else:
                self.Interface.solve(self.C, self.a, self.A)
        
        self.Result = self.Interface.Info
        return self.Result