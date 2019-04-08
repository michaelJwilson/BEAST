import  os 
import  sys 
import  numpy                   as      np
import  pylab                   as      pl
 
import  astropy.io.fits         as      fits
import  astropy.units           as      u
import  matplotlib.pyplot       as      plt

from    scipy.signal            import  medfilt
from    lya                     import  Gauss
from    scipy.interpolate       import  interp1d
from    utils                   import  latexify
from    scipy.signal            import  medfilt
from    lbg_maker               import  lbg_maker 
from    prep_filters            import  prep_filters 
from    desispec.interpolation  import  resample_flux


latexify(columns=2, equal=True, fontsize=10, ratio=None, ggplot=True, usetex=False)

root              =  os.environ['BEAST']

ngal              =  50

redshifts         =   1.5 + np.linspace(0.5, 4.0, 4)
magnitudes        =  24.3 * np.ones_like(redshifts)

##  Flux [1e-17 ergs/s/cm2/AA].
flux, wave, meta  =  lbg_maker(ngal=ngal, restframe=False, printit=False, test=False, seed=314, redshifts=redshifts, magnitudes=magnitudes, hiwave=2.e4)

##  
print(meta)

##  Subsample to 1A intervals in wavelength.                                                                                                                                         
swave             =  wave[::5]
flux              =  flux[:,::5]

##  Choose the nth SED and calculate residual from 100A smoothed version. 
nth               =  20
sflux             =  flux[nth,:] - medfilt(flux[nth,:], 101)

pl.clf()

fig               =  pl.gcf()

ax1               =  fig.add_subplot(211)
ax2               =  fig.add_subplot(212)

## 
axarr             = [ax1, ax2]

##                                                                                                                                                      
fnames            = [os.environ['CSCRATCH'] + '/desi/runs/three/lbg-beast-spectra-exp-3000.fits',\
                     os.environ['CSCRATCH'] + '/desi/runs/four/lbg-pfs-spectra-nearir-exp-3000.fits']

colors            =  plt.rcParams['axes.prop_cycle'].by_key()['color']

for exp, ax, fname in zip(['BEAST', 'PFS'], axarr, fnames):
  ## 
  dat      = fits.open(fname)
  tSN2     = 0.0

  for i, ARM in enumerate(['B', 'R', 'Z']):
    wave   = dat['%s_WAVELENGTH' % ARM].data
    flux   = dat['%s_FLUX'       % ARM].data               ##  10^-17 erg/s/cm2/Angstrom.'
      
    ivar   = dat['%s_IVAR' % ARM].data
    sig    = 1. / np.sqrt(ivar)
      
    back   = sig[0,:] * u.erg / u.s / u.cm / u.cm / u.AA
   
    SN2    = resample_flux(wave, swave, sflux) ** 2. / sig[0,:] ** 2.
    tSN2  += np.sum(SN2)

    ax.plot(wave / 1.e3, SN2, c=colors[i], label = str(meta['Lya-EW'][nth]) + ' ' + str(meta['REDSHIFT'][nth]))

  print(exp, np.sqrt(tSN2))

  ##  ax.set_ylim(0., 0.6)
  ##  ax.set_ylabel(r'$10^{-17}$ erg/$s$/cm$^2$/$\AA$')

  ax.set_axis_on()

  ax.spines['bottom'].set_color('black')
  ax.spines['top'].set_color('black')
  ax.spines['left'].set_color('black')
  ax.spines['right'].set_color('black')

  ax.set_xlim(3.55, 13.) 

pl.xlabel('Wavelength $[\mu \ m]$')
pl.legend()
plt.tight_layout()
pl.show()
##  pl.savefig('plots/back.pdf')

print('\n\nDone.\n\n')
