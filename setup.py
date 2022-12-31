###############################################################################
#
# SWIT v1.1: Seismic Waveform Inversion Toolbox
#
#   A Python package for seismic waveform inversion
#   By Haipeng Li at USTC & Stanford
#   Email: haipengl@mail.ustc.edu.cn, haipeng@stanford.edu 
#
#   Setup script for installing SWIT
#
###############################################################################

import os
import sys
import subprocess

from setuptools import setup
from setuptools.command.install import install

def compile_solver():
    ''' Compile the solver implemented in Fortran 90.
    '''

 
    # check if gfortran and make is installed
    try:
        print('\nInstallation: checking gfortran, make, mpif90, mpirun are installed...\n')
        subprocess.check_call(['which', 'gfortran'])
        subprocess.check_call(['which', 'make'])
    except:
        print('\n On Ubuntu, you can install "gfortran" and "make" by running: \n')
        print('     $ sudo apt-get install build-essential')
        print('     $ sudo apt install gfortran')
        print('\n On MacOS, you can install "gfortran" and "make" by running: \n')
        print('     $ brew install gfortran')
        print('     $ brew install make')
        print(' For more information about brew, please visit https://brew.sh/ \n')
        print('\n')
        sys.exit(1)
    

    # check if mpi is installed
    try:
        subprocess.check_call(['which', 'mpif90'])
        subprocess.check_call(['which', 'mpirun'])
        print('\nInstallation: All required compilers are installed, ready to compile the solver...')
    except:
        print('\nNo mpif90 or mpirun is found. Please install openmpi first. On Ubuntu or Mac, you can install openmpi by running: \n')
        print('     $ wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.1.tar.gz')
        print('     $ tar -xvf openmpi-4.1.1.tar.gz')
        print('     $ cd openmpi-4.1.1')
        print('     $ ./configure --prefix=/usr/local/openmpi CC=gcc FC=gfortran')
        print('     $ make')
        print('     $ sudo make install')
        print('\nIf no sudo permission, set the path to the local directory, e.g., "/home/username/software/openmpi"')
        print('Then, run the configure command with new local directory, make, and make install without sudo.\n')
        
        print('Next, add the following lines to your ~/.bashrc file:\n')
        print('     $ echo export PATH=/usr/local/openmpi/bin:$PATH')
        print('     $ echo export LD_LIBRARY_PATH=/usr/local/openmpi/lib:$LD_LIBRARY_PATH')

        print('\nFinally, source the ~/.bashrc file and you are good to go.\n')
        print('     $ source ~/.bashrc')
        print('\n')
        sys.exit(1)


    # compile the solver
    solver_dir = os.path.join(os.path.dirname(__file__), 'solver')
    try:
        # compile the solver
        print('Installation: compiling the solver...')
        subprocess.check_call(['make', 'clean'], cwd=solver_dir)
        subprocess.check_call(['make'], cwd=solver_dir)
        print('Installation: the solver is compiled successfully.')

        # add the solver to the path
        cmd = 'echo export PATH={}:$PATH >> ~/.bashrc'.format(os.path.join(os.path.dirname(__file__), 'bin'))
        os.system(cmd)
        os.system('. ~/.bashrc')
        print('\nInstallation: the solver is added to the path.\n')
    except:
        print('\nInstallation: failed to compile the solver. Please check the error message above.\n')
        sys.exit(1)


def read(fname):
    ''' Read the README file.
    '''
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


class CustomInstall(install):
    ''' Custom install class to compile the solver.
    '''

    def run(self):
        compile_solver()
        super().run()


setup(
    name='SWIT',
    version='1.1.0',
    description='Seismic Waveform Inversion Toolbox',
    author='Haipeng Li',
    author_email='haipeng@stanford.edu',
    url='https://github.com/Haipeng-ustc/SWIT-1.0',
    license='GNU',
    install_requires=[
        'numpy',
        'scipy',
        'obspy',
        'matplotlib',
        'multiprocess',
        # 'pandas==0.23.3',
        # 'numpy>=1.14.5'
    ],

    long_description=read('README.md'),
    packages=['toolbox-dev', ],

    # compile the solver
    cmdclass={'install': CustomInstall})
