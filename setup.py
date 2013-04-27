from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os

sourcefiles = ['pycsdp.pyx']
ext_modules = [Extension("pycsdp", 
				version='1.0.0',
				author='Mehdi Ghasemi',
				author_email='mghasemi@ntu.edu.sg',
				sourcefiles,
				library_dirs = ["./lib"],
				libraries = ["sdp", "lapack", "blas", "gfortran", "m"])]

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
    )
