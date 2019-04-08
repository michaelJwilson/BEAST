import numpy as np
import pylab as pl


dat       = np.loadtxt('../dat/du_ew.dat')
ngal      = dat[:,1].sum()

dat[:,1] *= 100. / ngal
dat[:,1]  = np.cumsum(dat[:,1])

print(dat)
