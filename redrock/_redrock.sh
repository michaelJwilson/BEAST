## If this script goes to completeness on a worker, and you try to rerun, 
## will fail due to empty NERSC_HOST.  So set up to be later unset. 
export  NERSC_HOST=cori

export  BEAST=/global/homes/m/mjwilson/desi/BEAST
source $BEAST/env.sh

##  To equip custom install.                                                                                                                                  
source /project/projectdirs/desi/software/desi_environment.sh ## 18.7                                                                                       

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

##  mpi4py
##  /global/common/software/desi/cori/desiconda/20180709-1.2.6-spec/aux/lib/python3.6/site-packages/mpi4py/MPI.cpython-36m-x86_64-linux-gnu.so

##  DO NOT USE.
##  http://www.nersc.gov/users/data-analytics/data-management/i-o-libraries/hdf5-2/h5py/
##  
##  module load h5py-parallel
##  module add hdf5-parallel/1.10.1

##  export SPECTRA=spec-qso-z2.4-rmag23.00

export SURVEY='desi'
export TARGET='elg'

for COUNT in {2..2..1};
do
    export EXPOSURE="$((1000 * $COUNT))"

    echo 'SOLVING FOR '$EXPOSURE' '$SURVEY' '$TARGET  

    ##  export INPUT=$BEAST/quickspectra/dat/$SPECTRA-exp-$EXPOSURE.fits
    export INPUT=/global/cscratch1/sd/mjwilson/desi/simspec/$TARGET-$SURVEY-spectra-exp-$EXPOSURE.fits

    ##  export RRH5=$BEAST/redrock/dat/rrh5/$SPECTRA-exp-$EXPOSURE-rr.h5
    ##  export ZBEST=$BEAST/redrock/dat/zbest/$SPECTRA-exp-$EXPOSURE-zbest.h5

    export RRH5=/global/cscratch1/sd/mjwilson/desi/simspec/$TARGET-$SURVEY-spectra-exp-$EXPOSURE-rr.h5
    export ZBEST=/global/cscratch1/sd/mjwilson/desi/simspec/$TARGET-$SURVEY-spectra-exp-$EXPOSURE-zbest.fits

    ##  --dependency=afterok:11254323
    ##  rrdesi $INPUT --output $RRH5  --zbest $ZBEST --mp 1

    ##  Memory overload on edison with n = 24. 
    srun -N 2 -n 48 rrdesi_mpi $INPUT --output $RRH5 --zbest $ZBEST
done

echo 'DONE.'
