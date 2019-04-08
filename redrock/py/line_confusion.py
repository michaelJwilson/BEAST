import  os
import  pickle
import  pylab                    as      pl
import  numpy                    as      np
import  matplotlib.pyplot        as      plt
import  matplotlib               as      mpl
import  matplotlib.colors        as      colors

from    astropy.table            import  Table, vstack
from    utils                    import  latexify
from    redrock.results          import  read_zscan
from    mpl_toolkits.axes_grid1  import  make_axes_locatable
from    collections              import  OrderedDict


run        =   'three'          ##  four for pfs, three otherwise. 
target     =     'lbg'
survey     =   'beast'
base_exp   =     1500.          ##  [Seconds].                                                                                                               
repeat     =        1
min_dChi2  =        9

if survey == 'pfs':
  nearir   = '-nearir' 

else:
  nearir   = ''

root       = os.environ['BEAST']

##
latexify(columns=1, equal=True, ggplot=True, usetex=True, fontsize=10)

pl.clf()

##  Truth values.                                                                                                                                
if target == 'lbg':
  names    = ['OBJTYPE', 'SUBTYPE', 'REDSHIFT', 'Lya-EW', '"AB mag"', 'dbands']

else:
  names    = ['OBJTYPE', 'SUBTYPE', 'REDSHIFT', 'AGE', 'METALLICITY', 'IMF', 'Tau', 'CalzEBV', 'dband', '"AB mag"']

##  
Truth      =  os.environ['CSCRATCH'] + '/desi/runs/%s/galmaker-%s-meta%s.txt' % (run, target, nearir)
Truth      =  Table.read(Truth, format='ascii', names=names)
Truth      =  vstack([Truth] * repeat)

if target == 'BC03':
  ##  Null Ly-EW col.                                                                                                                                    
  Truth['Lya-EW'] = np.ones_like(Truth['REDSHIFT'].quantity) * 1.e99

TRUEZS     =  np.array(Truth['REDSHIFT'].quantity).astype(np.float)
TRUEEWS    =  np.array(Truth['Lya-EW'].quantity).astype(str)
TRUEMAGS   =  np.array(Truth['"AB mag"'].quantity).astype(np.float)

EWs        =  sorted(list(set([str(x) for x in Truth['Lya-EW']])))
lws        =  {'-14.92': 2400, '-1.1': 600, '11.0': 600, '52.63': 600}

lw         =  0.4

zs         =  np.arange(1.30, 1.99, 0.01)  ##  No Ly-a cut.  

##  Emission lines
dat        =  np.loadtxt('../dat/sdss7_emlines.dat', usecols=(0))

##  Shapley absorption.
shapa      =  np.loadtxt('../dat/shapley_abslines.dat', usecols=(0))

##  Shapley emission.                                                                                                                                                            
shape      =  np.loadtxt('../dat/shapley_emlines.dat', usecols=(0))

alphas     =  np.arange(0.1, 1.1, 0.2)[::-1]

for ii, true_line in enumerate(np.concatenate([shapa[:2], shape[:2]])):
  for VAC in dat:
    AIR  = VAC / (1.0 + 2.735182e-4 + 131.4182 / VAC**2. + 2.76249e8 / VAC ** 4.)

    zp   = true_line * (1. + zs) / AIR - 1.
    pl.plot(zs, zp, 'k-', lw=lw, alpha=alphas[ii])

'''
##  Absorption lines. 
dat        = np.loadtxt('../dat/sdss7_abslines.dat', usecols=(0))

for row in dat:
  zp   = true_line * (1. + zs) / row - 1.
  pl.plot(zs, zp, 'c--', lw=lw)
'''
'''
## 
dat        = np.loadtxt('../dat/sdss7_skylines.dat', usecols=(0))

for row in dat:
  zp   = true_line * (1. + zs) / row - 1.
  pl.plot(zs, zp, 'm--', lw=lw)
'''

##  Define true line for confusion.                                                                                                                                                
true_line  =  dat[0]
zs         =  np.arange(1.98, 5.5, 0.01)  ##  Ly-a cut.

for VAC in dat:
  AIR  = VAC / (1.0 + 2.735182e-4 + 131.4182 / VAC**2. + 2.76249e8 / VAC ** 4.)
  
  zp   = true_line * (1. + zs) / AIR - 1.
  pl.plot(zs, zp, 'k-', lw=lw, alpha=alphas[ii])

##  H alpha.                                                                                                                                              
VAC  = 6564.61
AIR  = VAC / (1.0 + 2.735182e-4 + 131.4182 / VAC**2. + 2.76249e8 / VAC ** 4.)

zp   = true_line * (1. + zs) / AIR - 1.
pl.plot(zs, zp, 'r', lw=2.*lw, label=r'H$\alpha$')
    
##  OIII                 
VAC  = 5008.240
AIR  = VAC / (1.0 + 2.735182e-4 + 131.4182 / VAC**2. + 2.76249e8 / VAC ** 4.) 
                                                                                                                          
zp   = true_line * (1. + zs) / AIR - 1.
pl.plot(zs, zp, 'b', lw=2.*lw, label=r'OIII')

##  Hd                          
VAC  = 4102.89
AIR  = VAC / (1.0 + 2.735182e-4 + 131.4182 / VAC**2. + 2.76249e8 / VAC ** 4.)                                                                             
zp   = true_line * (1. + zs) / AIR - 1.

##  pl.plot(zs, zp, 'c', lw=2.*lw, label=r'H$\delta$')

##  OIIa                         
VAC  = 3727.092
AIR  = VAC / (1.0 + 2.735182e-4 + 131.4182 / VAC**2. + 2.76249e8 / VAC ** 4.)
                                                                                                                         
pl.plot(zs, true_line * (1. + zs) / AIR - 1., 'g', lw=2.*lw, label='OII')

##  OIIb
VAC  = 3729.875
AIR  = VAC / (1.0 + 2.735182e-4 + 131.4182 / VAC**2. + 2.76249e8 / VAC ** 4.)

pl.plot(zs, true_line * (1. + zs) / AIR - 1., 'g', lw=2.*lw)

##  MgII
VAC  = 2799.117
AIR  = VAC / (1.0 + 2.735182e-4 + 131.4182 / VAC**2. + 2.76249e8 / VAC ** 4.)

pl.plot(zs, true_line * (1. + zs) / AIR - 1., 'orange', lw=2.*lw, label='MgII')

print(EWs)

## 
for marker, EW in zip(['^',  '*', 'p', 'D'], EWs):
    print('\n\nSolving for %s' % EW)

    ##  Pickle it.                                                                                                                                      
    es     =  pickle.load(open(root + '/redrock/pickle/%s/%s/exptimes_%d.pkl' % (survey, target, min_dChi2), 'rb'))
    es     =  np.array(es).astype(np.float)
    es    /=  1500.                           ##  Seconds to exposures.  

    zbests =  pickle.load(open(root + '/redrock/pickle/%s/%s/zbests_%d.pkl' % (survey, target, min_dChi2), 'rb'))

    ##  Find the redshift confusion.                                                                                                                       
    cut    =  np.where((TRUEEWS == EW) & (es < 2.e4))

    args   =  np.argsort(TRUEZS[cut])

    xx     =  zbests[cut][args]
    yy     =  TRUEZS[cut][args]

    pl.scatter(yy + np.random.uniform(-5.e-2, 5.e-2), xx, s=16, c='k', alpha=0.6, marker=marker, label=EW + r'$\AA$')

##  pl.axvline(x=1.98, ymin=0., ymax=1., c='k')
##  pl.legend(ncol=4, loc=1)

ax = pl.gca()

ax.spines['bottom'].set_color('black')
ax.spines['top'].set_color('black')
ax.spines['left'].set_color('black')
ax.spines['right'].set_color('black')

pl.xlim(1.96, 5.0)
pl.ylim(0.0,  2.0)

pl.xlabel(r'$z$')
pl.ylabel(r'$\hat z$')

plt.tight_layout()      

##  pl.show()
pl.savefig('../plots/lines.pdf') 
