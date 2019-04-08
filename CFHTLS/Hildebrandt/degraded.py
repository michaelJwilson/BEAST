import  os
import  numpy               as      np
import  pylab               as      pl 
import  astropy.io.fits     as      fits
import  matplotlib.pyplot   as      plt

from    magABsource         import  magAB
from    HilderandtSN        import  merr, ferr
from    luptitudes          import  luptitude
from    utils               import  latexify
from    depths              import  get_depths
from    colourcut_dNdz      import  colourcut
from    app_mags            import  get_colors


'''
Hildebrandt Selection criteria:

--  Colour / colour.
    E.g.  'u_g_ISO_cor2' and 'g_r_ISO_cor2'

--  dat[1].data['MASK_u']      == 0
--  dat[1].data['MASK_g']      == 0
--  dat[1].data['MASK_r']      == 0
--  dat[1].data['MASK_i']      == 0
--  dat[1].data['MASK_z']      == 0
--  dat[1].data['MASK_STARS']  == 0 
--  dat[1].data['masksa']      == 0
--  dat[1].data['maskgscw2']   == 0
--  dat[1].data['CLASS_STAR'] < 0.9
'''

latexify(columns=1, equal=True, fontsize=12, ggplot=True, usetex=True)

##  Selection type. 
ttype   =        'u'  ##  ['u', 'g']
FULL    =      True   ##
DODEPTH =     'Full'  ##  ['Full', 'Degraded']

np.random.seed(seed=314)

##  Load Hendrik's catalogues. 
root    = os.environ['CSCRATCH']
dat     = fits.open(root + '/Hildebrandt/CFHTLS/D1.fits')

##  Good imaging cut. 
good    = (dat[1].data['MASK_u'] == 0) & (dat[1].data['MASK_g'] == 0) & (dat[1].data['MASK_r'] == 0) &\
          (dat[1].data['MASK_i'] == 0) & (dat[1].data['MASK_z'] == 0) & (dat[1].data['MASK_STARS'] == 0) &\
          (dat[1].data['masksa'] == 0) & (dat[1].data['maskgscw2'] == 0) & (dat[1].data['CLASS_STAR'] < 0.9)

##  dat[1].header

umag    = dat[1].data['MAG_ISO_u']
gmag    = dat[1].data['MAG_ISO_g']
rmag    = dat[1].data['MAG_ISO_r']
imag    = dat[1].data['MAG_ISO_i']
zmag    = dat[1].data['MAG_ISO_z']

umg     = dat[1].data['u_g_ISO_cor2'] 
gmr     = dat[1].data['g_r_ISO_cor2']
rmi     = dat[1].data['r_i_ISO_cor2']
imz     = dat[1].data['i_z_ISO_cor2']

##  Redshifts. 
HypZCWW = dat[1].data['Z_PHOT']          ##  HyperZ run with the CWW.                                                                                                                                                                                                                                           
HypZBC3 = dat[1].data['Z_PHOT_BC']       ##  HyperZ run with the BC03.                                                                                                                                                                                                                                        
BPZ     = dat[1].data['Z_B_V3']          ##  BpZ.   

##  Degraded depths. 
degraded_depths  =  get_depths()

ngal             =  len(umag)

if not FULL:
   nrows            =  500
   magdict          =  dict(zip(['u', 'g', 'r', 'i', 'z'], [umag[good][:nrows], gmag[good][:nrows], rmag[good][:nrows], imag[good][:nrows], zmag[good][:nrows]]))
   photozs          =  HypZBC3[good][:nrows]

else:
   magdict = dict(zip(['u', 'g', 'r', 'i', 'z'], [umag[good], gmag[good], rmag[good], imag[good], zmag[good]]))
   photozs = HypZBC3[good]

for band in magdict:
    Flux          =  10. ** (-(magdict[band] + 48.60) / 2.5)
    SigF          =  ferr(magdict[band], degraded_depths[band], estar=0.2, alphab=-0.25, alphaf=0.22, lim_snr=None)
    
    mSigF         =  np.ma.masked_invalid(SigF, copy=True)     
    Noise         =  np.random.normal(loc=np.zeros_like(mSigF), scale=mSigF)

    ##  See also eqn. (10.2) of Chromey, Introduction to Observational Astronomy.                                                                                                                                                                                                                               
    ##  dmag      =  -2.5 * np.log10(Flux + Noise) - 48.60
    lup           =  luptitude(Flux + Noise, SigF)
        
    ##  print('%s \t %.6lf \t %.6lf' % (band, magdict[band], lup))

    if DODEPTH == 'Degraded':
      magdict[band] = lup

## 
detected  = 0
udetected = 0

## 
if   ttype == 'u':
    bcol = 'u-g'
    rcol = 'g-r'

    bmin = -0.5
    bmax =  2.5

    rmin = -0.3
    rmax =  1.2

    infrac  = 1.00
    outfrac = 0.01
    
elif ttype == 'g':
    bcol = 'g-r'
    rcol = 'r-i'
    
    bmin = -0.5
    bmax =  2.5

    rmin = -0.3
    rmax =  1.2
    
    infrac  = 1.00
    outfrac = 0.01

else:
    raise  UserWarning()

##  dN/dz binning.                 
dz        = 0.25                                                                                                                                                                                                                             
zbins     = np.arange(0.0, 6.0 + dz, dz)
midz      = zbins[:-1] + dz/2.

drop_zs   = []

for i, row in enumerate(magdict['u']):
   mags   = dict(zip(magdict.keys(), [magdict[x][i] for x in magdict.keys()]))
   colors = get_colors(mags, get_colors=['g-r', 'r-i', 'u-g', 'g-r', 'u-z', 'g-i'], fname = None)

   if ttype     == 'u':
       is_detected  = (mags['r'] < degraded_depths['r'])

   elif ttype == 'g':
       is_detected  = (mags['i'] < degraded_depths['i'])

   else:
       is_detected  = False

   ##  if DODEPTH == 'Full':
   ##    is_detected  =  True

   if is_detected:
       detected  += 1
       is_lbg     = colourcut(mags, dropband=ttype, good=True, fourthlimit=False, BzK_type='all')

       draw       = np.random.uniform(0.0, 1.0, 1)

       if draw <= 0.01:
         print('%d/%d\t%+.4lf \t %+.4lf \t %+.3lf \t %+.3lf \t %s \t %+.3lf \t %+.3lf' % (i, ngal, mags['u'], mags['g'], mags['r'],\
                                                                                          photozs[i], str(is_lbg), degraded_depths['r'],\
                                                                                          degraded_depths['i']))

       if is_lbg:
          drop_zs.append(photozs[i])

          ##  cax = plt.scatter(colors[rcol], colors[bcol], c=photozs[i], s=10, vmin=0.0, vmax=5.0, rasterized=True, alpha=0.8)

       else:
          draw       = np.random.uniform(0.0, 1.0, 1)
 
          if draw <= outfrac:
             pass
             ##  plt.scatter(colors[rcol], colors[bcol], c=photozs[i], marker='x', alpha=0.2, s=10, vmin=0.0, vmax=5.0, rasterized=True)
  
   else:
       udetected += 1

print('Detection rates: %d \t %d' % (detected, udetected))

##  Histogram of color-cut redshifts.                                                                                                                                                                                                                         
(dNdz, bins) = np.histogram(drop_zs, bins = zbins)
dNdz         = dNdz.astype(np.float)

output       = np.c_[midz, dNdz]
np.savetxt('dNdz/%s_%sdrops_dz_%.2lf.txt' % (DODEPTH, ttype, dz), output, fmt='%.6le')

exit(1)

ax  = pl.gca()
fig = pl.gcf()

ax.spines['bottom'].set_color('black')
ax.spines['top'].set_color('black')
ax.spines['left'].set_color('black')
ax.spines['right'].set_color('black')

pl.xlabel(r'$%s$' % rcol)
pl.ylabel(r'$%s$' % bcol)

stretch = 0.0

pl.xlim(rmin - stretch, rmax + stretch)
pl.ylim(bmin - stretch, bmax + stretch)

##  [left, bottom, width, height].                                                                                                                                                                                                                                                                                
cbaxes = fig.add_axes()   ## [0.8, 0.1, 0.03, 0.8]                                                                                                                                                                                                                                                                
cb     = plt.colorbar(ax=ax, cax=cbaxes, label=r'redshift')

cb.set_alpha(1.0)
cb.draw_all()

if DODEPTH == 'Full':
    cb.set_alpha(0.0)
    cb.set_label('')
    cb.set_ticklabels([])
    cb.set_ticks([])
    cb.draw_all()
    
plt.tight_layout()

pl.show()

print('\n\nDone.\n\n')
