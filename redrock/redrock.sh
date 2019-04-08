#!/bin/bash -l
#SBATCH -q regular
#SBATCH -N 9
#SBATCH -t 00:25:00
#SBATCH -C haswell

##  Note:  SBATCH -N should match the number of exposures.  
##
##  If this script goes to completeness on a worker, and you try to rerun,                                                                                  
##  will fail do to empty NERSC_HOST.  So set up to be later unset.
export  NERSC_HOST=cori
export  BEAST=/global/homes/m/mjwilson/desi/BEAST
source $BEAST/env.sh

##  To equip custom install.                                                                                                                                
source /project/projectdirs/desi/software/desi_environment.sh  ##  18.7                                                                                     

##  $DESIMODULES/install_jupyter_kernel.sh 18.7                                                                                                              
##  JUPYTER @ NERSC                                                                                                                                        
##  https://jupyter-dev.nersc.gov/                                                                                                                           

module unload desimodel
module unload desisim
module unload desispec
module unload desitarget
module unload desiutil
module unload redrock
module unload redrock-achetypes
module unload redrock-templates
module unload specsim
module unload specter
module unload tutorials

export  ROOT=/global/homes/m/mjwilson/desi/BEAST/
export  BEAST=/global/homes/m/mjwilson/desi/BEAST
export  SIMSPEC=$CSCRATCH/desi/simspec/

export  DESIHUB=/global/homes/m/mjwilson/desi/BEAST/desihub/
export  DESIMODEL=$DESIHUB/desimodel/

export  DESI_ROOT=/project/projectdirs/desi/

export  DESI_BASIS_TEMPLATES=$DESIHUB/basis_templates/v3.0/

export  PYTHONPATH=/global/homes/m/mjwilson/desi/BEAST/gal_maker/py/:$PYTHONPATH
export  PYTHONPATH=$BEAST/gal_maker/py/:$PYTHONPATH

export  PYTHONPATH=$DESIHUB/desimodel/py:$PYTHONPATH
export  PATH=$DESIHUB/desimodel/bin:$PATH

export  PYTHONPATH=$DESIHUB/desispec/py:$PYTHONPATH
export  PATH=$DESIHUB/desispec/bin:$PATH

export  PYTHONPATH=$DESIHUB/desiutil/py:$PYTHONPATH
export  PATH=$DESIHUB/desiutil/bin:$PATH

export  PYTHONPATH=$DESIHUB/redrock/py:$PYTHONPATH
export  PATH=$DESIHUB/redrock/bin:$PATH

export  PYTHONPATH=$DESIHUB/specter/py:$PYTHONPATH
export  PATH=$DESIHUB/specter/bin:$PATH

export  PYTHONPATH=$DESIHUB/desisim/py:$PYTHONPATH
export  PATH=$DESIHUB/desisim/bin:$PATH

export  PYTHONPATH=$DESIHUB/desitarget/py:$PYTHONPATH
export  PATH=$DESIHUB/desitarget/bin:$PATH

export  PYTHONPATH=$DESIHUB/specsim:$DESIHUB/specsim/py:$PYTHONPATH
export  PATH=$DESIHUB/specsim/bin:$PATH

export  PYTHONPATH=/global/homes/m/mjwilson/desi/BEAST/desihub/empca/:$PYTHONPATH

##
export  OMP_NUM_THREADS=1

##  INCLUDED IN DESI ENVIRONMENT. 
##  export  HDF5_USE_FILE_LOCKING=FALSE

##
echo   $NERSC_HOST
unset   NERSC_HOST
echo   $NERSC_HOST

##  export SURVEY='desi'
##  export TARGET='lbg'

for COUNT in 1 2 3 4 5 6 7 8 9
do
    export EXPOSURE="$((1500 * $COUNT))"

    echo 'SOLVING FOR '$EXPOSURE' '$SURVEY' '$TARGET  

    export ODIR=/global/cscratch1/sd/mjwilson/desi/runs/$RUN/
    
    if [ $SURVEY == 'pfs' ]; then
	export INPUT=$ODIR/$TARGET-$SURVEY-spectra-nearir-exp-$EXPOSURE.fits
	export RRH5=$ODIR/$TARGET-$SURVEY-spectra-nearir-exp-$EXPOSURE-rr.h5
	export ZBEST=$ODIR/$TARGET-$SURVEY-spectra-nearir-exp-$EXPOSURE-zbest.fits
	export TEMPLATEDIR=$BEAST/desihub/redrock/py/redrock/templates-nearir/
    
    else
	export INPUT=$ODIR/$TARGET-$SURVEY-spectra-exp-$EXPOSURE.fits
	export RRH5=$ODIR/$TARGET-$SURVEY-spectra-exp-$EXPOSURE-rr.h5
	export ZBEST=$ODIR/$TARGET-$SURVEY-spectra-exp-$EXPOSURE-zbest.fits
	export TEMPLATEDIR=$BEAST/desihub/redrock/py/redrock/templates/
    fi
    if [ $SBATCH == 'False' ]; then
      ##  --dependency=afterok:11254323  
      rrdesi $INPUT --output $RRH5 --zbest $ZBEST --mp 4 --templates $TEMPLATEDIR
    else
      ##  Memory overload on edison with n = 24. 
      srun -N 1 -n 24 rrdesi_mpi $INPUT --output $RRH5 --zbest $ZBEST --templates $TEMPLATEDIR & 
    fi
done
wait
