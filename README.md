# SWIT-1.0

Seismic Waveform Inversion Toolbox 1.0

### Install SWIT 

#### Step 1: Install gcc and gfortran

```bash
#Install gcc and gfortran
sudo apt-get install build-essential
sudo apt install gfortran
```

#### Step 2 : Install OpenMPI

```bash
#Download the latest OpenMPI package, or go to  http://www.open-mpi.org/software/ompi to download the desired version
wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.1.tar.gz 
tar xvfz openmpi-4.1.1.tar.gz
cd openmpi-4.1.1

#Configure the installation files and install OpenMPI (this would takes a while)
./configure --prefix=/usr/local/openmpi CC=gcc FC=gfortran
make
sudo make install

#Add env path 
export PATH=/usr/local/openmpi/bin:$PATH
source ~/.bashrc

#Check OpenMPI is successfully installed
mpirun --version
```

#### Step 3 : Install Anaconda Environment  

```bash
#Anaconda is recommended. For installing Anaconda, refers to https://docs.anaconda.com/anaconda/install/linux/
#1. download package from: https://www.anaconda.com/products/individual/download-success
#2. bash ~/your_Anaconda_package

#Once the Anaconda is installed, create the conda environment for SWIT
conda create --name SWIT python=3
conda activate SWIT
#Install dependencies
pip install numpy obspy scipy matplotlib
pip install multiprocess PySimpleGUI psutil Pillow
```

#### Step 4 : Install SWIT  

```bash
#Download from: https://github.com/Haipeng-ustc/SWIT-1.0
git clone https://github.com/Haipeng-ustc/SWIT-1.0

#Complie the fd2dmpi forward solver (Fortran version)
#Edit the Makefile.config file, make sure FCC is correct in line 18
cd ~/SWIT-1.0/fd2dmpi/
rm *.mod
make clean   
make

#Add fd2dmpi to the env path
export PATH=~/SWIT-1.0/bin:$PATH
source ~/.bashrc

#Install SWIT package
cd ~/SWIT-1.0/toolbox/
python setup.py install
```

### SWIT Examples 

#### Case 1: Overthrust Model Inversion   (GUI)

```bash
~~~
```

#### Case 2: Marmousi Model Inversion   (Python script)

```bash
~~~
```

