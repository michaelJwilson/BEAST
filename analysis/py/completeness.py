import  os
import  pylab              as      pl
import  numpy              as      np
import  desisim
import  desisim.templates
import  desispec.io
import  desisim.io
import  desisim.templates
import  matplotlib.pyplot  as      plt

from    astropy.io         import  fits
from    astropy.table      import  Table, Column
from    matplotlib         import  colors         as mcolors
from    matplotlib.lines   import  Line2D
from    utils              import  latexify
from    zz                 import  get_table


latexify(columns=2, equal=False, fontsize=10, ratio=None, ggplot=True, usetex=False)

target   =     'lbg'
survey   =   'beast'
exposure =     4.e3
repeat   =        2

gal_table, zs, types, mags, umags, ngal, markers = get_table(target, repeat=repeat,\
                                                             path=os.environ['BEAST'] + '/gal_maker/dat/Tables/galmaker-%s-meta.txt' % target)

##  Find the min. / max. non-zero. mag.
finite   = np.isfinite(mags)
print(mags[finite].min(), mags[finite].max())

zbest    = Table.read('/global/cscratch1/sd/mjwilson/desi/simspec/%s-%s-spectra-exp-%d-zbest.fits' % (target, survey, exposure), 'ZBEST')
rrzs     = np.array(zbest['Z'])

good_z                      = np.ones_like(rrzs)
good_z[zbest['ZWARN'] != 0] = 0.0             

EWs               = np.unique(gal_table['Lya-EW'])
maglims           = np.arange(20., 26.1, 0.2)

for EW in EWs:
    cut           = np.where(gal_table['Lya-EW'] == EW)
    EWCUT         = good_z[cut]
    
    completeness  = []
    
    for maglim in maglims:
        magcut    = EWCUT[(maglim <= mags[cut]) & (mags[cut] < (maglim + 0.2))]
        completeness.append(100.00 * np.sum(magcut) / len(magcut))

    completeness  = np.array(completeness)

    pl.plot(maglims, completeness, label=str(EW))

pl.xlabel(r'$r$')
pl.ylabel('Completeness [%]')
pl.title('%s for %s with %ds exposure' % (target.upper(), survey.upper(), exposure))
pl.legend()
plt.savefig('../plots/zgood_%s_%s_exp_%d.pdf' % (target, survey, exposure), bbox_inches='tight')
