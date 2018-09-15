from distutils.core import setup, Extension

module1 = Extension('simulator',
                    include_dirs = ['C:\\Users\\jops\\AppData\\Local\\Programs\\Python\\Python36\\Lib\\site-packages\\numpy\\core\\include'],
                    sources = ['main.cpp', 'Simulator.cpp'])

setup (name = 'Simulator',
       version = '1.0',
       description = 'This is a simulation package',
       ext_modules = [module1])
