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
from    scipy.ndimage.filters   import  gaussian_filter1d


latexify(columns=1, equal=True, fontsize=8, ratio=None, ggplot=True, usetex=True)

root              =  os.environ['BEAST']

THRESHOLD         =   5.00
ngal              =   1000

redshifts         =   1.5 + np.linspace(0.5, 4.0, 500)
magnitudes        =  24.3 * np.ones_like(redshifts)

##  Flux [1e-17 ergs/s/cm2/AA].
flux, wave, meta  =  lbg_maker(ngal=ngal, restframe=False, printit=False, test=False, seed=314, redshifts=redshifts, magnitudes=magnitudes, hiwave=2.e4)

##  
print(meta)

##  Subsample to 1A intervals in wavelength.                                                                                                                                         
swave             =    wave[::5]
flux              =  flux[:,::5]

##                                                                                                                                                                                                         
fnames            = [os.environ['CSCRATCH'] + '/desi/runs/three/lbg-beast-spectra-exp-3000.fits',\
                     os.environ['CSCRATCH'] + '/desi/runs/four/lbg-pfs-spectra-nearir-exp-3000.fits']

colors            =  plt.rcParams['axes.prop_cycle'].by_key()['color']

for exp, fname in zip(['BEAST', 'PFS'], fnames):
  pl.clf()
  
  ##  Choose the nth SED and calculate residual from 100A smoothed version.                                                                                                                               
  ax          = pl.gca()
 
  dat         = fits.open(fname)

  result      = {'Q0': [], 'Q2': [], 'Q3': [], 'Q4': []}

  for nth in np.arange(ngal):
    ##  sflux = flux[nth,:] - medfilt(flux[nth,:], 101)
    gflux     = gaussian_filter1d(flux[nth,:], 100, output=None, mode='reflect') 

    sflux     = flux[nth,:] - gflux

    ## 
    tSN2      = 0.0
    
    for i, ARM in enumerate(['B', 'R', 'Z']):
      _wave   = dat['%s_WAVELENGTH' % ARM].data      
      _ivar   = dat['%s_IVAR' % ARM].data
      _sig    = 1. / np.sqrt(_ivar)
      
      _back   = _sig[0,:] * u.erg / u.s / u.cm / u.cm / u.AA
   
      SN2     = resample_flux(_wave, swave, sflux) ** 2. / _sig[0,:] ** 2.
      tSN2   += np.sum(SN2)

    print(meta['SUBTYPE'][nth], meta['REDSHIFT'][nth], meta['AB mag'][nth], tSN2)
    
    T   =  3000. * (THRESHOLD / np.sqrt(tSN2)) ** 2.
    T  /=  1500.  ## [1500s]

    result[meta['SUBTYPE'][nth]].append([meta['REDSHIFT'][nth], meta['AB mag'][nth], T])

  for i, key in enumerate(result):
    result[key] = np.array(result[key])
    ind         = np.argsort(result[key][:,0])  ##  sort by redshift. 
      
    pl.semilogy(result[key][ind,0], result[key][ind,2], label=key, c=colors[i])

    np.savetxt('Tz/Tz_%s_%s_%d.txt' % (exp, key, THRESHOLD), np.c_[result[key][ind,0], result[key][ind,2]])
  
  pl.xlabel(r'$z$')
  pl.ylabel(r'$T \ / \ 10^{-0.8(24.3 - m)} \ / \ (\theta/\theta_f)^2 \ / \ 1500s$')

  ##  pl.scatter(result[:,0], result[:,1], c=result[:,2], s=8)

  pl.ylim(5.e-2, 15.)
  ##  pl.colorbar()
  ##  plt.clim(1.0, 12.0)
  
  ax.set_axis_on()

  ax.spines['bottom'].set_color('black')
  ax.spines['top'].set_color('black')
  ax.spines['left'].set_color('black')
  ax.spines['right'].set_color('black')

  plt.tight_layout()
  pl.legend(ncol=2, frameon=False, loc=1)
  ##  pl.show()
  pl.savefig('plots/Tz_%s.pdf' % exp)

print('\n\nDone.\n\n')
