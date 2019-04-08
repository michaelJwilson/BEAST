import os 
import sys 
import numpy             as      np
import pylab             as      pl

import astropy.io.fits   as      fits
import astropy.units     as      u
import matplotlib.pyplot as      plt

from   scipy.signal      import  medfilt
from   lya               import  Gauss
from   scipy.interpolate import  interp1d
from   utils             import  latexify


depths = [25.2, 24.0]     ## At z=2.5
depths = [24.8, 23.6]     ## At z=5.0   

zs     = np.arange(2.5, 7.5, 2.5)
wave   = 1216. * (1. + zs) * u.AA                     ##  Angstrom.

for mag, EW in zip(depths, [52.63, 11.0, -1.1, -14.92]):
  flux  = 10. ** (-(mag + 48.60) / 2.5)               ##  ergs/s/cm2/Hz.
  flux  = flux * u.erg / u.s / u.cm / u.cm  / u.Hz

  flux  = flux.to(u.erg / u.s / u.cm / u.cm / u.AA, equivalencies=u.spectral_density(wave)) 
  flux *= 1.e17                                       ##  10^-17 erg/s/cm2/Angstrom.

  print(EW, flux * EW)
    
print('\n\nDone.\n\n')
