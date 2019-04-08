import astropy.io.fits as fits
import numpy as np
import os
import pylab as pl
import pandas as pd


cols     = ['ID', 'RA', 'DEC', 'TYPE', 'ZSPEC', 'LYAFLUX', 'FLUXERR', 'EW', 'EWPERR', 'EWNERR', 'FWHM', 'EFWHM', 'FLAG', 'SKEW'] ## 'SKEWERR'            
Mallery  = pd.read_csv(os.environ['BEAST']+'/COSMOS/DEIMOS/Mallery12/Table2.dat', engine='python', error_bad_lines=False, delim_whitespace=True, names=cols)


## cosmos = fits.open('asu.fits')
## ras    = cosmos[1].data['RAJ2000']
## decs   = cosmos[1].data['DEJ2000']

cosmos    = np.genfromtxt('cosmos_phot_20060103.tbl', dtype=str) ##  max_rows=10  
ras       = cosmos[:,2].astype(np.float)
decs      = cosmos[:,3].astype(np.float)

arcsec    = 1. / 60. / 60.

precision = [5, 2, 1]
matches   = []

dtheta    = 6.5 / 60. / 60.
count     = 0
degen     = 0

for i, x in enumerate(Mallery['RA']):
  indx    = np.where(((x - dtheta) <= ras) & (ras <= (x + dtheta)))[0] 

  if len(indx) == 0:
    print('Warning:  no match for %.3lf' % x)
  
  else:
    iindx      = np.where(((Mallery['DEC'][i] - dtheta) <= decs[indx]) & (decs[indx] <= (Mallery['DEC'][i] + dtheta)))[0]
    indx       = np.array(np.array(indx)[iindx])

    if len(indx) > 0:
      count += 1

    if len(indx) > 0:
      degen += len(indx)

    out      = '\n(%.6lf, %.6lf),\n' % (x, Mallery['DEC'][i])
    out     += '  '.join('(%.6lf, %.6lf), ' % (ras[ind], decs[ind]) for ind in indx)

    print(out)

print(len(Mallery['RA']), count, degen)
