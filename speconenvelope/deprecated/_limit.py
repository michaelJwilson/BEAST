import os, sys 
import numpy             as      np
import pylab             as      pl

import astropy.io.fits   as      fits
import astropy.units     as      u

from   scipy.signal      import  medfilt
from   lya               import  Gauss
from   scipy.interpolate import  interp1d


if len(sys.argv) == 1:
  EW        =   52.63                  ##  [52.63, 11.00, -1.1, -14.92]
  bexposure =    1500                  ##  base exposure time. 
  nexp      =       2                  ##  Scaling of base exposure time. 
  nsig      =     0.3                  ##  Number of sigma over continuum for successful redshift. 
  survey    =   'psf'
  run       =  'four'
  lw        =   1000

else:
  EW        =  np.float(sys.argv[1])   ##  [52.63, 11.00, -1.1, -14.92]                                                                                     
  bexposure =  np.int(sys.argv[2])     ##  base exposure time.                                                                                      
  nexp      =  np.int(sys.argv[3])     ##  Scaling of base exposure time.                                                                                    
  nsig      =  np.int(sys.argv[4])     ##  Number of sigma over continuum for successful redshift.                                                           
  survey    =  sys.argv[5]
  run       =  sys.argv[6]
  lw        =  np.int(sys.argv[7])

##  print(EW, bexposure, nexp, nsig, survey, run)
print('Solving for: EW = %f, bexposure = %d, nexp = %d, nsig = %d, survey=%s, run=%s, lw=%d' % (EW, bexposure, nexp, nsig, survey, run, lw))

mags        =  np.arange(23., 26., 0.5)

## 
'''
for mag in mags:
  ##  mag  = -2.5 np.log10(flux) - 48.60
  flux     = 10. ** (-(mag + 48.60) / 2.5)               ##  ergs/s/cm2/Hz.
  flux     = flux * u.erg / u.s / u.cm / u.cm  / u.Hz

  wave     = np.arange(3.e3, 1.e4, 1.) * u.AA            ##  Angstrom.

  flux     = flux 
  flux     = flux.to(u.erg / u.s / u.cm / u.cm / u.AA, equivalencies=u.spectral_density(wave)) 
  flux    *= 1.e17                                       ##  10^-17 erg/s/cm2/Angstrom.

  pl.plot(wave, flux, label=mag)
'''
##
base      = True  ##  Utilise baseline 1500s exposure and scale with nexp.
fs        = []

exposure  = nexp * bexposure

if survey == 'pfs':
    nearir  = '-nearir'

else:
    nearir  = ''

## 
if base:
  fname   = os.environ['CSCRATCH'] + '/desi/runs/%s/lbg-%s-spectra%s-exp-%d.fits' % (run, survey, nearir,     1500)
  dat     = fits.open(fname)

else:
  fname   = os.environ['CSCRATCH'] + '/desi/runs/%s/lbg-%s-spectra%s-exp-%d.fits' % (run, survey, nearir, exposure)
  dat     = fits.open(fname)

##  
for ARM in ['B', 'R', 'Z']:
  wave   = dat['%s_WAVELENGTH' % ARM].data
  flux   = dat['%s_FLUX'       % ARM].data               ##  10^-17 erg/s/cm2/Angstrom.'

  dwave  = wave[1] - wave[0] 

  ##  assert  np.allclose(dwave, np.diff(wave)[1:])

  ivar   = dat['%s_IVAR' % ARM].data
  sig    = 1. / np.sqrt(ivar)

  zs     = wave / 1215.24 - 1.0
  Gs     = []

  ##  Normalised Gaussian max. 
  for z in zs:
    lambda0  = 1215.24 * (1. + z)
    sigma    = lambda0 * (1. * lw / 2.9979e5)            ##  E.g.  600 km/s line profile.

    if survey == 'pfs':
      S2N2     = 1.5 ** 2.

    else:
      S2N2     = 2.0 ** 2.

    ##  Empirical for -1.1A suggests nsig=1 sufficient for redshift.
    proxy    = 1.4 * S2N2 + np.sum((np.abs(EW) * Gauss(wave, lambda0, sigma)) ** 2.) * dwave
    nsig     = np.sqrt(S2N2 / proxy)

    Gs.append(nsig)

  Gs      = np.array(Gs)
  back    = sig[0,:] * 1.e-17 * u.erg / u.s / u.cm / u.cm / u.AA
  
  if base:
    back *= np.sqrt(1. / nexp) 

  Flim    = Gs * back

  Flim    = Flim.to(u.erg / u.s / u.cm / u.cm/ u.Hz, equivalencies=u.spectral_density(wave * u.AA))
  Flim    = -2.5 * np.log10(Flim.value) - 48.60
  
  ##  pl.plot(zs, medfilt(Flim, 13))
    
  interim = interp1d(zs, medfilt(Flim, 13), kind='linear', copy=True, bounds_error=False, fill_value=0.0, assume_sorted=False)
  fs.append(interim)

wave   = np.arange(3.e3, 1.e4, 1.0) 
zs     = wave / 1215.24 - 1.0

result = lambda x: np.array([f(x) for f in fs]).max() 
##  pl.plot(zs, np.array([result(z) for z in zs]), 'k-')

##  pl.xlabel('$z$')
##  pl.ylabel(r'$10^{-17}$ erg/$s$/cm$^2$/$\AA$')

##  pl.xlim(1.5,  5.5)
##  pl.ylim(22., 26.5)

##  pl.legend(ncol=3)
##  pl.title(survey.upper() + ' in %d s' % exposure)

##  pl.show()
##  pl.savefig('sn_%d.pdf' % exposure)

##  Save. 
np.savetxt('dat/limit_%s_%d_%s_%d.dat' % (survey, exposure, EW, lw), np.c_[zs, np.array([result(z) for z in zs])])

print('\n\nDone.\n\n')
