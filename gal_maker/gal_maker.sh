#!/bin/bash -l
#SBATCH -q debug
#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -t 00:25:00
#SBATCH -J mjwilson-gmaker
#SBATCH -o mjwilson-gmaker-stdout 
#SBATCH -e mjwilson-gmaker-err


export  BEAST=/global/homes/m/mjwilson/desi/BEAST
source $BEAST/env.sh

##  INCLUDED IN DESI ENVIRONMENT. 
##  export HDF5_USE_FILE_LOCKING=FALSE

##  DO NOT USE.
##  http://www.nersc.gov/users/data-analytics/data-management/i-o-libraries/hdf5-2/h5py/

##  srun -n 1 python py/gal_maker.py
srun -n 1 python py/lbg_maker.py
