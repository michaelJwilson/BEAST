import  os
import  pylab              as      pl
import  numpy              as      np
import  desisim
import  desisim.templates
import  desispec.io
import  desisim.io
import  desisim.templates

from    astropy.io         import  fits
from    astropy.table      import  Table, Column
from    matplotlib         import  colors         as mcolors
from    matplotlib.lines   import  Line2D
from    utils              import  latexify


def type2color(types):
    unique  = list(set(types))

    nlen    = len(unique)
    ucolors = list(mcolors.BASE_COLORS.keys())[:len(unique)]

    result  = types.copy()

    for i, x in enumerate(result):        
        result[i] = ucolors[unique.index(result[i])]

    return  np.array(result)

def get_table(target, repeat=1, path=None):
  if path is None:
      path =  os.environ['BEAST'] + '/gal_maker/dat/Tables/galmaker-%s-meta.txt' % target

  t        =  Table.read(path, format='ascii')
  zs       =  np.tile(t['REDSHIFT'], 2)
  types    =  np.tile(np.array(t['OBJTYPE']), 2) 
  mags     =  np.tile(t['r'], 2)
  umags    =  np.unique(mags)[:-1]
  ngal     =  repeat * len(mags)

  markers  =  Line2D.filled_markers[1:]

  return  t, zs, types, mags, umags, ngal, markers

if __name__ == '__main__':
  latexify(fig_width=15., fig_height=8., columns=2, equal=False, fontsize=10, ratio=None, ggplot=True, usetex=False)

  ##  y = x.                                                                                                                                                            
  pl.plot(np.arange(0.0, 10., 0.01), np.arange(0.0, 10., 0.01), 'k-', alpha=0.3)

  repeat     =       2

  for exposure in [1e3, 4e3, 7e3]: 
    for target in ['elg', 'qso', 'lbg']:
      pl.clf()

      print('%s \t %d' % (target, exposure))
     
      for survey, color in zip(['desi', 'beast'], ['r', 'b']):
        ##  Reset zs, mags after binning-subtraction.
        t, zs, types, mags, umags, ngal, markers = get_table(target=target)

        maglims    = np.arange(np.floor(mags[np.isfinite(mags)].min()), np.ceil(mags[np.isfinite(mags)].max()), 0.5)[::-1]
  
        zbest      = Table.read('/global/cscratch1/sd/mjwilson/desi/simspec/safe/%s-%s-spectra-exp-%d-zbest.h5' % (target, survey, exposure), 'ZBEST')
        rr_zs      = np.array(zbest['Z'])
 
        zmin, zmax = np.floor(zs.min()), np.ceil(zs.max())

        for i, (maglim, marker) in enumerate(zip(maglims, markers)):
          pl.scatter(zs[mags >= maglim], rr_zs[mags >= maglim], c=color, label=survey.upper() + ' ' + str(maglim), s=20., marker=marker)

          zs       = zs[mags < maglim]
          rr_zs    = rr_zs[mags < maglim]
          mags     = mags[mags < maglim]
          markers  = markers[1:]
  
      pl.xlabel(r'$z_{\rm{True}}$')
      pl.ylabel(r'$z_{\rm{redrock}}$')

      pl.xlim(zmin, zmax)
      pl.ylim(zmin, zmax)

      pl.legend(ncol=2)
      pl.title('TARGETS:  %s for EXPOSURE = %.2lf min' % (target.upper(), exposure / 60.))

      pl.savefig(os.environ['BEAST'] + '/analysis/plots/zz-%s-exposure-%d.pdf' % (target, exposure), bbox_inches='tight')
  
  print('\n\nDone.\n\n')
