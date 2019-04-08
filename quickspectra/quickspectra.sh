#!/bin/bash -l                                                                                                                                             
#SBATCH -q debug                                                                                                                                          
#SBATCH -N 1                                                                                                                                              
#SBATCH -C haswell                                                                                                                                        
#SBATCH -t 00:30:00                                                                                                                                       


export  BEAST=/global/homes/m/mjwilson/desi/BEAST
source $BEAST/env.sh

##  INCLUDED IN DESI ENVIRONMENT.                                                                                                                         
##  export HDF5_USE_FILE_LOCKING=FALSE

##  DO NOT USE.                                                                                                                                           
##  http://www.nersc.gov/users/data-analytics/data-management/i-o-libraries/hdf5-2/h5py/                                                                  
##                                                                                                                                                       
##  module load h5py-parallel                                                                                                                            
##  module add hdf5-parallel/1.10.1                                                                                                                        

##  SPECTRA=spec-elg-o2flux-8e-17
##  SPECTRA=spec-qso-z2.4-rmag23.00

##  SURVEY='desi'
##  TARGET='lbg'
##  PROGRAM='DARK'  ## 'DARK' / 'BRIGHT'
##  REPEAT=2

echo 'SOLVING FOR '$EXPOSURE' '$SURVEY' '$TARGET' '$PROGRAM' '$REPEAT

##  --moonfrac 0.9 --moonalt 70 --moonsep 20
##  srun -n 1
quickspectra -i $INPUT -o $OUTPUT --repeat $REPEAT --exptime $EXPOSURE --seeing 0.8 --specsimconfig $SURVEY --program $PROGRAM --fiberloss table

