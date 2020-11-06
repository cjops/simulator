from distutils.core import setup, Extension
from numpy.distutils.misc_util import get_numpy_include_dirs

module1 = Extension('simulator',
                    include_dirs = get_numpy_include_dirs(),
                    sources = ['py.cpp', 'simulator.cpp'],
                    extra_compile_args=['-std=c++11'])

setup (name = 'simulator',
       version = '1.0',
       description = 'This is a simulation package',
       ext_modules = [module1])
