###############################################################################
#
# SWIT v1.1: Seismic Waveform Inversion Toolbox
#
#   A Python package for seismic waveform inversion
#   Developed by Haipeng Li at USTC, updated on 2022-12-21 at Stanford
#   haipengl@mail.ustc.edu.cn, haipeng@stanford.edu
#
#   Data preprocessing module
#
###############################################################################


# Use one thread in calling scipy to do the filtering
os.environ["OMP_NUM_THREADS"] = "1" # export OMP_NUM_THREADS=1
