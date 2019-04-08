import  os 
import  sys 
import  numpy              as      np
import  pylab              as      pl

import  astropy.io.fits    as      fits
import  astropy.units      as      u
import  matplotlib.pyplot  as      plt

from    scipy.signal       import  medfilt
from    lya                import  Gauss
from    scipy.interpolate  import  interp1d
from    utils              import  latexify


latexify(columns=2, equal=False, fontsize=10, ggplot=True, usetex=True)

## 
mags      = np.arange(23., 26., 0.5)

pl.clf()
fig       = pl.gcf()

ax1       = fig.add_subplot(211)
ax2       = fig.add_subplot(212)

## 
axarr     = [ax1, ax2]

##                                                                                                                                                      
fnames    = [os.environ['CSCRATCH'] + '/desi/runs/three/lbg-beast-spectra-exp-3000.fits',\
             os.environ['CSCRATCH'] + '/desi/runs/four/lbg-pfs-spectra-nearir-exp-3000.fits']

colors    =  plt.rcParams['axes.prop_cycle'].by_key()['color']

for ax, fname in zip(axarr, fnames):
  for mag in mags:
    flux     = 10. ** (-(mag + 48.60) / 2.5)               ##  ergs/s/cm2/Hz.
    flux     = flux * u.erg / u.s / u.cm / u.cm  / u.Hz

    wave     = np.arange(3.e3, 1.3e4, 1.) * u.AA           ##  Angstrom.

    print(wave)
    
    flux     = flux 
    flux     = flux.to(u.erg / u.s / u.cm / u.cm / u.AA, equivalencies=u.spectral_density(wave)) 
    flux    *= 1.e17                                       ##  10^-17 erg/s/cm2/Angstrom.
    
    ax.plot(wave / 1.e4, flux, label=mag, c='k', lw=1., alpha=0.5)
    
  ## 
  dat      = fits.open(fname)

  for i, ARM in enumerate(['B', 'R', 'Z']):
    wave   = dat['%s_WAVELENGTH' % ARM].data
    flux   = dat['%s_FLUX'       % ARM].data               ##  10^-17 erg/s/cm2/Angstrom.'
      
    ivar   = dat['%s_IVAR' % ARM].data
    sig    = 1. / np.sqrt(ivar)
      
    back   = sig[0,:] * u.erg / u.s / u.cm / u.cm / u.AA
   
    ax.plot(wave / 1.e4, back, c=colors[i], lw=1.)

  ax.set_ylim(0., 0.6)
  ax.set_ylabel(r'$B \ [10^{-17}$ erg/$s$/cm$^2$/$\AA$]')

  ax.set_axis_on()

  ax.spines['bottom'].set_color('black')
  ax.spines['top'].set_color('black')
  ax.spines['left'].set_color('black')
  ax.spines['right'].set_color('black')

  ax.set_xlim(0.355, 1.3) 

ax1.tick_params(
  axis='x',          # changes apply to the x-axis
  which='both',      # both major and minor ticks are affected
  bottom=False,      # ticks along the bottom edge are off
  top=False,         # ticks along the top edge are off
  labelbottom=False) # labels along the bottom edge are off

pl.xlabel('Wavelength $[\mu \ m]$')
plt.tight_layout()
##  pl.show()
pl.savefig('plots/back.pdf')

print('\n\nDone.\n\n')
