##  To equip custom install.
source /project/projectdirs/desi/software/desi_environment.sh 18.7
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

## export  PYTHONPATH=$DESIHUB/specsim:$DESIHUB/specsim/py:$PYTHONPATH
## export  PATH=$DESIHUB/specsim/bin:$PATH

export  PYTHONPATH=$BEAST/sandbox/specsim:$BEAST/sandbox/specsim/py:$PYTHONPATH                                                                          
export  PATH=$BEAST/sandbox/specsim/bin:$PATH

export  AUX=/global/homes/m/mjwilson/desi/BEAST/sandbox/specsim/aux

export  PYTHONPATH=/global/homes/m/mjwilson/desi/BEAST/desihub/empca/:$PYTHONPATH

echo 'NEW ENVIRONMENT SET.'
