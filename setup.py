from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import os

SAGE_LIB = os.environ['SAGE_LOCAL']+'/lib'
CSDP_INCLUDES = os.environ['SAGE_ROOT'] + '/spkg/build/csdp-6.2/src/include/csdp'

sourcefiles = ['pycsdp.pyx']
setup(
    name = 'SDP',
    version = '1.0.0',
    author = 'Mehdi Ghasemi',
    author_email = 'mehdi.ghasemi@gmail.com',
    packages = ['SDP'],
    url = 'https://github.com/mghasemi/SDP.git',
    license = 'GNU GNU Public License (GPL)',
    description = 'A generic python wraper for SDP solvers.',
    long_description = open('README.txt').read(),
)

ext_modules = [Extension("pycsdp", sourcefiles,
				library_dirs = [SAGE_LIB],
				define_macros = [('NOSHORTS',None)],
				include_dirs = [CSDP_INCLUDES],
				libraries = ["sdp", "lapack", "blas", "gfortran", "m"])]

setup(
	name = 'pycsdp',
	version = '1.0.0',
	author = 'Mehdi Ghasemi',
	author_email = 'mehdi.ghasemi@gmail.com',
	url = 'https://github.com/mghasemi/pycsdp.git',
	license = 'GNU GNU Public License (GPL)',
	description = 'Cython CSDP wrapper for Python and Sage.',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
    )
