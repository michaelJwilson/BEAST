import os
import numpy as np
import pylab as pl

root  = os.environ['BEAST']

iacs  = np.loadtxt(root + '/filters/hst/acs_f775w.pb')
zacs  = np.loadtxt(root + '/filters/hst/acs_f850lp.pb')

ulsst = np.loadtxt(root + '/filters/lsst/total_i.dat')
glsst = np.loadtxt(root + '/filters/lsst/total_i.dat')
rlsst = np.loadtxt(root + '/filters/lsst/total_i.dat')
ilsst = np.loadtxt(root + '/filters/lsst/total_i.dat')
zlsst = np.loadtxt(root + '/filters/lsst/total_z.dat')

USTD  = np.loadtxt(root + '/filters/steidel/U.pb')
GSTD  = np.loadtxt(root + '/filters/steidel/G.pb')
RSTD  = np.loadtxt(root + '/filters/steidel/R.pb')

filters    = [iacs, zacs, ulsst, glsst, rlsst, ilsst, zlsst, USTD, GSTD, RSTD]
wavescales = [1., 1., 10., 10., 10., 10., 10., 1., 1., 1.],
labels     = [r'$i$-ACS', r'$z$-ACS', r'$u$-LSST', r'$g$-LSST', r'$r$-LSST', r'$i$-LSST', r'$z$-LSST', r'$U$-Steidel', r'$G$-Steidel', r'$R$-Steidel'] 

for x, wavescale, label in zip(filters, wavescales, labels):
  pl.plot(x[:,0], x[:,1], label=label)

pl.legend(ncol=3)

pl.xlabel(r'$\lambda$')
pl.ylabel(r'$T$')

## pl.xlim(6.e3, 11.5e3)

pl.savefig('plots/stark.pdf')
