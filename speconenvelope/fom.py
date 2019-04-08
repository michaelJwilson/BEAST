import pandas as pd
import numpy  as np

from   astropy.cosmology import FlatLambdaCDM


##
h               = 0.673
H0              = 100. * h

cosmo           = FlatLambdaCDM(H0=H0, Om0=0.3)

area            = 14.e3

##  fNL.
##  log10nbar   = -4

##  RSD.
log10nbar       = np.log10(3.e-4)

##                                                                                                                                                         
print('\n\nWelcome to figure-of-merit (n = %le) .\n\n' % log10nbar)

zmin            =   2.0
zmax            =   5.0

volume          = (cosmo.comoving_volume(zmax).value - cosmo.comoving_volume(zmin).value) 
volume         *=  h **3.

fsky            = area / 41252.96
ngal            = 10. ** log10nbar * fsky * volume

dat             = pd.read_csv('surveys.txt', sep='\s+', comment='#', names=['Survey', 'Diameter', 'R', 'Multiplex', 'FOV'])
dat.set_index(dat['Survey'], drop=True, append=False, inplace=True, verify_integrity=True)

dat['FOV']     /= 60. ##  arcmin to deg. 
dat['FOV']      = np.pi * (dat['FOV'] / 2.) ** 2. 

dat['nFov']       = np.ceil(area / dat['FOV'])
dat['nMultiplex'] = np.ceil(ngal / dat['Multiplex'])

dat['npoint']     = np.ceil(np.maximum(ngal / dat['Multiplex'], area / dat['FOV']))
dat['exposure']   = 3000. * np.ones_like(dat['npoint']) / 60.                                         ##  Minutes.

for survey in dat['Survey']:
  dat.at[survey, 'exposure'] *= (dat.at['BEAST', 'Diameter']**2.   / dat.at[survey, 'Diameter']**2.)
  dat.at[survey, 'exposure'] *= (dat.at['BEAST', 'R'] / dat.at[survey, 'R'])

dat['FOM%.1le' % log10nbar]  =  1. / (dat['exposure'] * dat['npoint'] / 60. / 24. / 30. / 12. / 10.)  ##  Exposure to decades. 

print(dat)

print('\n\nDone.\n\n')
