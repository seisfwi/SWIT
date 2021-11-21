# **S**eismic **W**aveform **I**nversion **T**oolbox   (SWIT-1.0)

By Haipeng Li @ USTC

Contact: haipengl@mail.ustc.edu.cn

### First look at SWIT

<img src="./doc/SWIT-GUI.png" style="zoom:25%;" />

### Contents of SWIT



<img src="./doc/SWIT-Contents.png" style="zoom:25%;" />

### Workflow of SWIT 

<img src="./doc/SWIT-Workflow.png" style="zoom:25%;" />

### Waveform Selection in FWI

<img src="./doc/SWIT-Waveform-selection.png" style="zoom:25%;" />



## SWIT Installation 

#### Step 1: Install  gfortran

```bash
sudo apt-get install build-essential
sudo apt install gfortran
```

#### Step 2 : Install OpenMPI

```bash
# Download the latest OpenMPI package, 
# or go to: http://www.open-mpi.org/software/ompi to download the desired version
wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.1.tar.gz 
tar xvfz openmpi-4.1.1.tar.gz
cd openmpi-4.1.1

# Configure & install OpenMPI (this would take quite a while)
./configure --prefix=/usr/local/openmpi CC=gcc FC=gfortran
make; sudo make install

# Add env path at your ~/.bashrc
vim ~/.bashrc
export PATH=/usr/local/openmpi/bin:$PATH
source ~/.bashrc

# Check
which mpirun
```

#### Step 3 : Install Python Dependencies

```bash
# Anaconda is recommended. For installing Anaconda: https://docs.anaconda.com/anaconda/install/linux/
# create a new envirment for SWIT
# conda create --name SWIT python=3.7.5
# conda activate SWIT

# Install dependencies using USTC mirrors (whether use Anaconda or not)
pip install numpy obspy scipy matplotlib multiprocess PySimpleGUI psutil Pillow -i https://pypi.mirrors.ustc.edu.cn/simple/
```

#### Step 4 : Install SWIT  

```bash
# Complie the fd2dmpi forward solver with the default fortran compiler and mpif90.
# Edit ~/SWIT-1.0/fd2dmpi/Makefile.config (line 18) to change.
cd /your/own/path/to/SWIT-1.0/fd2dmpi/
rm *.mod; make clean; make

# Add fd2dmpi and Python toolbox to the env path at your ~/.bashrc 
vim ~/.bashrc 
export PATH=/your/own/path/to/SWIT-1.0/bin:$PATH
export PYTHONPATH=/your/own/path/to/SWIT-1.0/toolbox
source ~/.bashrc
```

#### Step 5: Run SWIT  

```bash
# Option 1. Run SWIT via GUI
cd /your/own/path/to/SWIT-1.0/toolbox/
python runswit_Linux.py    
# python runswit_MacOS.py 

# Option 2. Run SWIT via the Python script
cd /your/own/path/to/SWIT-1.0/example/some_case/
./run_workflow     

# You need to modify all the paths in the Python script before running
```

## FWI examples (keep updating)

| No.  | Acquisition |   Model    |         Misfit         |       Features        | Optimization |     Size      |
| :--: | :---------: | :--------: | :--------------------: | :-------------------: | :----------: | :-----------: |
|  1   |    Land     |  Marmousi  |        Waveform        |           -           |     NLCG     | 481x121, 25 m |
|  2   |    Land     | Overthrust |        Waveform        |           -           |     NLCG     | 401x101, 25 m |
|  3   |   Marine    |  Marmousi  |        Waveform        |           -           |     NLCG     | 481x141, 25 m |
|  4   |   Marine    | Overthrust |        Waveform        |           -           |     NLCG     | 401x121, 25 m |
|  5   |    Land     |  Marmousi  | Traveltime &  Waveform |   1D initial model    |     NLCG     | 401x121, 25 m |
|  6   |    Land     | Overthrust |        Waveform        | Multi-scale Inversion |     NLCG     | 401x101, 25 m |
|  7   |    Some     |  scripts   |       to convert       |       your data       |     for      |   SWIT-1.0    |

## Citations :   

```
If you find SWIT is useful, please cite the following reference:

1. Li, H., Li, J., Liu, B., & Huang, X. (2021). Application of full-waveform tomography on deep seismic profiling data set for tectonic fault characterization. In First International Meeting for Applied Geoscience & Energy Expanded Abstracts (pp. 657–661). https://doi.org/10.1190/segam2021-3583190.1

2. Schuster, G. T. (2017). Seismic inversion. Society of Exploration Geophysicists. https://library.seg.org/doi/book/10.1190/1.9781560803423
```
## Few more words:
1. Simplicity is the greatest virtue ever.

2. The seismic WIT always lies within.
