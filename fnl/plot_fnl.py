import numpy as np
import pylab as pl
import matplotlib.pyplot as plt

from   scipy.interpolate import interp1d


cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
sigs  = []

for kk, (target, nbar, label) in enumerate(zip(['bmbx', 'u', 'g', 'r'], [9.8e-4, 1.2e-4, 1.0e-4, 0.4e-4], ['2', '3', '4', '5'])): 
  dat    = np.loadtxt('dat/fnl_%sdrop.dat' % target)
  pl.loglog(dat[:,0], dat[:,1], label=r'$z \simeq $' + label, c=cycle[kk])
  pl.axvline(nbar, ymin=0., ymax=1., c=cycle[kk], linestyle='--')

  fnlsig = interp1d(dat[:,0], dat[:,1], kind='linear', copy=True, bounds_error=True, assume_sorted=False)
  fnlsig = fnlsig(nbar)
  
  sigs.append(fnlsig)

sigs  = np.array(sigs)
sigs  = 1. / sigs ** 2.
sigs  = sigs.sum()

final = 1. / np.sqrt(sigs)

pl.axhline(y=final, xmin=0., xmax=1., label=r'BEAST: $f_{NL}=%.2lf$' % final, c='k')

pl.xlim(5.e-3, 1.e-5)
pl.ylim(0.5,      10.)

pl.ylabel(r'$\sigma_{f_{NL}}$')
pl.xlabel(r'$\bar n \ [(h^{-1} \rm{Mpc})^3]$')

pl.legend(loc=3, ncol=3)

pl.savefig('beast_fnl.pdf')
