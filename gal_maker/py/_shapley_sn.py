import os
import numpy        as     np
import pylab        as     pl
import matplotlib   as     mpl

from   lbg_maker    import lbg_maker
from   prep_filters import prep_filters
from   madau        import lephare_madau
from   extinction   import calzetti00
from   extinction   import apply          as ext_apply
from   scipy.signal import medfilt


root              =  os.environ['BEAST']

redshifts         =   1.5 + np.linspace(0., 4.0, 4)
magnitudes        =  22.5 * np.ones_like(redshifts)

flux, wave, meta  =  lbg_maker(ngal=5, restframe=False, printit=False, test=False, seed=314, redshifts=redshifts, magnitudes=magnitudes, hiwave=1.e4)

filters           =  prep_filters(['LSST', 'VIDEO'], normed=True)

##  Subsample to 1A. 
wave              =  wave[::5]
flux              =  flux[:, ::5]

smooth_flux       =  medfilt(flux[1,:], 11)

pl.plot(wave, flux[1,:] - smooth_flux)
pl.xlim(3550., 1.e4)
pl.show()
