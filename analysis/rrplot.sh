##  https://github.com/desihub/tutorials/blob/master/redrock/RedrockPlotSpec.md
ROOT=$CSCRATCH/desi/simspec/

## Input files
## export rrfile=$ROOT'/precious/lbg-beast-spectra-exp-8000-rr.h5'
## export specfile=$ROOT'/lbg-beast-spectra-exp-8000.fits'

export rrfile='/global/cscratch1/sd/mjwilson/desi/simspec/lbg-desi-spectra-exp-7000-rr.h5'
export specfile='/global/cscratch1/sd/mjwilson/desi/simspec/lbg-desi-spectra-exp-7000.fits'

##  Now actually run it.
rrplot --specfile $specfile --rrfile $rrfile
