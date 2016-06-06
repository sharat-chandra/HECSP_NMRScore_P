
from distutils.core import setup
import os
import sys

# Packages in MSMT toolkit
packages = ['lib', 'coodtrans', 'moltools', 'ringpath']

## Modules
#modules = ['pymsmtexp']

# Scripts
scripts = ['calHcsp/hecsp.py']

## See if our Python version will support OpenMM. Of the ParmEd-supported
## Pythons, only 2.4 and 2.5 do not work with OpenMM
#major, minor = sys.version_info[:2]

if __name__ == '__main__':

    try:
        from distutils.command.build_py import build_py_2to3 as build_py
        from distutils.command.build_scripts import build_scripts_2to3 as build_scripts
        PY3 = True
    except ImportError:
        from distutils.command.build_py import build_py
        from distutils.command.build_scripts import build_scripts
        PY3 = False

    setup(name='HECSP',
          version='1.0', # For AmberTools 15
          description='Proton Chemical shift Perturbation calculation',
          author='Zhuoqin Yu',
          author_email='zhuoqiny -at- gmail.com',
          license='GPL v2 or later',
          packages=packages,
#          py_modules=modules,
          cmdclass={'build_py':build_py, 'build_scripts':build_scripts},
          scripts=scripts)