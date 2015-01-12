from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy
import os

SAGE_LIB = os.environ['SAGE_LOCAL']+'/lib'
CSDP_INCLUDES = os.environ['SAGE_LOCAL'] + '/include/csdp'

sourcefiles = ['cpycsdp.pyx']
setup(
    name = 'SDP',
    version = '1.0.0',
    author = 'Mehdi Ghasemi',
    author_email = 'mehdi.ghasemi@gmail.com',
    packages = ['SDP'],
    url = 'https://github.com/mghasemi/SDP.git',
    license = 'GNU GNU Public License (GPL)',
    description = 'A generic python wraper for SDP solvers.',
    long_description = open('README.md').read(),
)

ext_modules = [Extension("cpycsdp", sourcefiles,
				library_dirs = [SAGE_LIB],
				define_macros = [('NOSHORTS',None)],
				include_dirs = [CSDP_INCLUDES, numpy.get_include()],
				libraries = ["sdp", "lapack", "blas", "gfortran", "m"])]

# For Sage 6.5+ use libraries = ["sdp", "atlas", "gfortran", "m"]) instead.

setup(
	name = 'cpycsdp',
	version = '1.0.0',
	author = ['Dmitrii V Pasechnik', 'Mehdi Ghasemi'],
	author_email = ['dima@ntu.edu.sg', 'mehdi.ghasemi@gmail.com'],
	url = 'https://github.com/mghasemi/pycsdp.git',
	license = 'GNU GNU Public License (GPL)',
	description = 'Cython CSDP wrapper for Python and Sage.',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
    )
