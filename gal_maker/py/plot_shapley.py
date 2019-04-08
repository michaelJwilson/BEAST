import os
import numpy              as     np
import pylab              as     pl
import matplotlib         as     mpl
import matplotlib.pyplot  as     plt

from   lbg_maker          import lbg_maker
from   prep_filters       import prep_filters
from   madau              import lephare_madau
from   extinction         import calzetti00
from   extinction         import apply          as ext_apply
from   utils              import latexify 
 

latexify(columns=1, equal=True, fontsize=10, ratio=None, ggplot=True, usetex=True)

root              =  os.environ['BEAST']

redshifts         =   1.5 + np.linspace(0., 4.0, 4)
magnitudes        =  22.5 * np.ones_like(redshifts)

flux, wave, meta  =  lbg_maker(ngal=None, restframe=False, printit=False, test=True, seed=314, redshifts=None, magnitudes=None, hiwave=2.e4)

filters           =  prep_filters(['LSST', 'VIDEO'], normed=True)

for i, band in enumerate(list(filters.keys())):
    pl.fill_between(filters[band]['ls'] / 1.e4, 2. * filters[band]['Ts'], alpha=0.4, label='')  ##  label=filters[band]['ppkey']

for i, x in enumerate(flux):
    rwave = wave / (1. + meta['REDSHIFT'][i])
    x     = ext_apply(calzetti00(rwave, a_v=0.2, r_v=4.05, unit='aa'), x)

    ##  Madau extinction. 
    x    *= lephare_madau(rwave, meta['REDSHIFT'][i]) 

    pl.semilogy(wave / 1.e4, x, label=r'$%+ .1lf\AA, z=%.1lf$' % (meta['Lya-EW'][i], meta['REDSHIFT'][i]))

pl.xlabel(r'$\lambda \ [\mu m]$')
pl.ylabel(r'$F_\lambda \ [10^{-17} \ \rm{ergs} / s / \rm{cm}^2 / \AA]$')

pl.xlim(0.3, 1.9)
pl.ylim(0.1,   30.1)

pl.legend(ncol=1, loc=1, frameon=False)
plt.tight_layout()

##  pl.show()
pl.savefig(root + '/gal_maker/plots/shapley.pdf')
