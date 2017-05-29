from distutils.core import setup, Extension

module1 = Extension('ccmap',
                    libraries = ['m'],
                    include_dirs = ['./include'],
                    sources = ['ccmapmodule.c', './src/mesh.c', './src/python_logging.c'])

setup (name = 'ccmapModule',
       version = '1.0',
       description = 'This is the C implementation of the mesh based contact map',
       ext_modules = [module1])