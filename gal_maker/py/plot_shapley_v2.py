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
 

latexify(columns=2, equal=False, fontsize=10, ratio=None, ggplot=True, usetex=True)

root              =  os.environ['BEAST']

redshifts         =   1.5 + np.linspace(0., 4.0, 4)
magnitudes        =  22.5 * np.ones_like(redshifts)

flux, wave, meta  =  lbg_maker(ngal=None, restframe=False, printit=False, test=True, seed=314, redshifts=None, magnitudes=None, hiwave=2.e4)

filters           =  prep_filters(['LSST', 'EUCLID-NIR'], normed=True)

rcol = 'coral'
bcol = 'lightblue'

for band, norm, color in zip(['H', 'J', 'Y', 'u', 'g', 'r', 'i', 'z', 'y'], 3.e2 * np.array([1.2, 1.1, 1.0] + [1.0] * 6), [rcol, rcol, rcol, bcol, bcol, bcol, bcol, bcol, bcol]):
    pl.fill_between(np.log10(filters[band]['ls'] / 1.e4), norm * filters[band]['Ts'], alpha=0.6, label='', color=color)  ##  label=filters[band]['ppkey']

for i, x in enumerate(flux):
    rwave = wave / (1. + meta['REDSHIFT'][i])
    x     = ext_apply(calzetti00(rwave, a_v=0.2, r_v=4.05, unit='aa'), x)

    ##  Madau extinction. 
    x    *= lephare_madau(rwave, meta['REDSHIFT'][i]) 

    pl.semilogy(np.log10(wave / 1.e4), wave * x, label=r'$%+ .1lf\AA, z=%.1lf$' % (meta['Lya-EW'][i], meta['REDSHIFT'][i]))

pl.xlabel(r'$\log_{\rm{10}}|\lambda / \mu m|$')
pl.ylabel(r'$\lambda \cdot F_\lambda \ [10^{-17} \ \rm{ergs} / s / \rm{cm}^2 / \AA]$')

pl.xlim(-0.5, 0.32)
pl.ylim(1.e2, 1.e5)

pl.legend(ncol=2, loc=1, frameon=False)
plt.tight_layout()

##  pl.show()
pl.savefig(root + '/gal_maker/plots/shapley_v2.pdf')
