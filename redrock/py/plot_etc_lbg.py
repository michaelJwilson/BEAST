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
from    scipy.signal             import  medfilt
from    scipy.interpolate        import  interp1d


run        =     'four'          ##  four for pfs, three otherwise. 
target     =      'lbg'
survey     =      'pfs'
base_exp   =      1500.          ##  [Seconds].                                                                                                               
repeat     =         1
min_dChi2  =        25

if survey == 'pfs':
  nearir   = '-nearir' 

else:
  nearir   = ''

root       = os.environ['BEAST']

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

EWs        =  list(set([str(x) for x in Truth['Lya-EW']]))

for EW in EWs:
    print('\n\nSolving for %s' % EW)

    ##  Pickle it.                                                                                                                                      
    es     =  pickle.load(open(root + '/redrock/pickle/%s/%s/exptimes_%d.pkl' % (survey, target, min_dChi2), 'rb'))
    es     =  np.array(es).astype(np.float)
    es    /=  1500.                           ##  Seconds to exposures.  

    zbests =  pickle.load(open(root + '/redrock/pickle/%s/%s/zbests_%d.pkl' % (survey, target, min_dChi2), 'rb'))

    ##  And plot.                                                                                                                                      
    pl.clf()

    latexify(fig_width=None, fig_height=None, columns=1, equal=True, ggplot=True, fontsize=10)
    
    ax = pl.gca()
    ax.set_axis_on()

    ax.spines['bottom'].set_color('black')
    ax.spines['top'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.spines['right'].set_color('black')

    ## 
    fig      = pl.gcf()

    cmap     = plt.cm.tab20c
    cmap     = colors.ListedColormap([cmap(i) for i in np.array([1, 4, 5, 8, 9, 12, 13, 17, 18])[::1]])  ##  [2, 16, 13]                              

    ##  Set bad to be white.                                                                                                                              
    cmap.set_bad(color = 'white', alpha = 1.)
    
    ##  
    QS         = {'-14.92': 'Q0', '-1.1': 'Q2', '11.0': 'Q3', '52.63': 'Q4'}  
    CONFIDENCE = 5

    fname      = root + '/speconenvelope/Tz/Tz_%s_%s_%d.txt' % (survey.upper(), QS[EW], CONFIDENCE)
    _zs, _Ts   = np.loadtxt(fname, unpack=True)
    _Ts        = interp1d(_zs, _Ts, kind='linear', copy=True, bounds_error=True, assume_sorted=False)

    zs         = np.arange(2.1, 5.0, 0.01)
    Ts         = _Ts(zs) * (6. / CONFIDENCE)**2.
 
    for ii in np.arange(1, 3, 1):
      ms       = 24.3 + np.log10(ii / Ts) / 0.8

      pl.plot(zs, ms, alpha=0.45, lw=0.5, c=cmap(ii-1))

    ## 
    bounds   = np.arange(1, 8, 1)
    norm     = colors.BoundaryNorm(bounds, cmap.N, clip=False)

    ##  print(type(EW), np.array(Truth['Lya-EW']).astype(str), np.array(Truth['Lya-EW'])[np.array(Truth['Lya-EW']).astype(str) == EW])

    cut      = np.where((TRUEEWS == EW) & (es > 2.e4))
    ##  pl.scatter(TRUEZS[cut], TRUEMAGS[cut], c='k', marker='x', alpha=0.2)

    print('Number with no redshifts in maximum exposure: %d' % len(es[cut]))

    cut      = np.where((TRUEEWS == EW) & (es <= 6))
    ##  plt.pcolormesh(xv, yv, zv, cmap=cmap, norm=norm) 
    pl.scatter(TRUEZS[cut], TRUEMAGS[cut], c=es[cut].astype(np.float), marker='^', cmap=cmap, norm=norm, s=15)

    ##  Good z dot. 
    pl.plot(zbests[cut], TRUEMAGS[cut], 'ko', markersize=1)

    if target != 'BC03':
      ##  4.15
      plt.text(1.65, 26.7, r'$%s\AA$' % EW, bbox=dict(facecolor='none', edgecolor='none'),\
               horizontalalignment='left', verticalalignment='bottom', fontsize=8)

    pl.xlim(0.98 * np.unique(TRUEZS)[0], 1.02 * np.unique(TRUEZS)[-1])

    ##  [0.98 * np.unique(TRUEMAGS[np.isfinite(TRUEMAGS)])[0] + 2.0, 1.02 * np.unique(TRUEMAGS[np.isfinite(TRUEMAGS)])[-1] -1.0]
    pl.ylim(23.0, 27.0)

    pl.xlabel(r'$z$')
    pl.ylabel(r'$m_{UV}$')
    
    ##  pl.title(survey)

    ax      = pl.gca()
    divider = make_axes_locatable(ax)
    cax     = divider.append_axes("right", size="5%", pad=0.05)

    cb      = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, spacing='proportional', ticks=bounds+0.5,\
                                        boundaries=bounds, format='%.2lf')

    cb.set_ticklabels(bounds, update_ticks=True)
    cb.set_label(r'$\rm{Minimum \ exposure \ time} \ [1500 s]$', rotation=270, labelpad=15)
    
    plt.tight_layout()

    ofname  = os.environ['BEAST'] + '/redrock/plots/%s/%s/expgrid_%s_%s_%d.pdf' % (survey, target, survey, EW.replace('.', 'p'), min_dChi2)
    
    ##  pl.show()
    pl.savefig(ofname, bbox_inches='tight')

print('\n\nDone.\n\n')
