# Seismic Waveform Imaging Toolbox (SWIT-1.2)

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
#Extract the file
tar xvfz openmpi-4.0.2.tar.gz
cd openmpi-openmpi-4.0.2
#Configure the installation files 
./configure --prefix=/usr/local/openmpi CC=gcc FC=gfortran
#Make and install OpenMPI, this would takes a while
make
sudo make install
#Add env paths 
export PATH=/usr/local/openmpi/bin:$PATH
source ~/.bashrc
#Check that OpenMPI is successfully installed
mpirun --version
```

#### Step 3 : Install Python Environment  

```bash
#Anaconda is recommended. For installing Anaconda, refers to https://docs.anaconda.com/anaconda/install/linux/
#Create the conda environment for SWIT
conda create --name swit python=3
conda activate swit
```

#### Step 4 : Install SWIT  

```bash
# complie the finite-difference forward solver (Fortran version)
cd ~/swit/fd2dmpi/
# edit the Makefile.config file, make sure FCC is correct in line 18
make clean; make
export PATH=~/swit/bin:$PATH

# Install SWIT package
cd ~/swit/toolbox/
python setup.py install

```
