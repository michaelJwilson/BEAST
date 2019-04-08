import  os
import  itertools
import  numpy                        as      np
import  _pickle                      as      cPickle
import  astropy.units                as      u
import  pylab                        as      pl

from    utils                        import  latexify
from    prep_filters                 import  prep_filters
from    BC03_maker                   import  BC03_maker
from    app_mags                     import  get_appmags
from    sklearn                      import  linear_model
from    sklearn.metrics              import  mean_squared_error, r2_score
from    isubaru                      import  bands                          as  isubaru
from    matplotlib                   import  pyplot                         as   plt
from    matplotlib.patches           import  Ellipse

from    astroML.utils                import  pickle_results
from    astroML.density_estimation   import  XDGMM
from    astroML.plotting.tools       import  draw_ellipse


np.random.seed(314)

##  latexify(fig_width=None, fig_height=None, columns=2, equal=False)                                                                                                                                                                                                                                                          
@np.vectorize
def merr(m, mstar, estar=0.2, alphab=-0.25, alphaf=0.22):
  ##  Magnitude error,                                                                                                                                                                                                                                                                                                         
  if m < mstar:
    return  estar * 10. ** (0.4 * (alphab + 1.) * (m - mstar))

  else:
    return  estar * np.exp(10. ** (alphaf * (m - mstar))) / 2.72

def ferr(m, mstar, estar=0.2, alphab=-0.25, alphaf=0.22, lim_snr=None):
  flux = 10. ** -((m + 48.60) / 2.5)  ##  Assumes AB mag.  [erg/s/cm2/Hz].                                                                                                                                                                                                                                                     

  sigm = merr(m, mstar, estar=estar, alphab=alphab, alphaf=alphaf)

  ##  E.g. eqn. (6) of https://arxiv.org/pdf/1509.00870.pdf                                                                                                                                                                                                                                                                    
  sigf = flux * sigm * np.log(10.) / 2.5

  if lim_snr is not None:
    sigf = np.sqrt(sigf**2. + (flux / lim_snr) ** 2.)

  return  sigf


if __name__ == '__main__':
  print('\n\nWelcome to (BC03) LBG maker.')

  ngal              =    1000  ##  Over written if rest-frame.                                                                                                                                                                                                                                                               
  test              =   False
  target_type       =   'BC03'
  save              =   False
  restframe         =   False

  redshifts         =   3.5 + np.linspace(0.,   1.0,    20)
  magnitudes        =  20.0 +   np.arange(0.,   7.0,   0.1)

  flux, wave, meta  =  BC03_maker(ngal=ngal, restframe=restframe, printit=False, test=test, redshifts=redshifts, magnitudes=magnitudes,\
                                  alliseds=True, calzetti=True, madau=True)
  
  uwave             =  wave * u.AA
  vs                =  uwave.to(u.Hz, equivalencies=u.spectral())

  print
  print(wave.shape)
  print(flux.shape)
  print(meta)

  ##                                                                                                                                                                                                                                                                                                                       
  filters     =  prep_filters(['LSST'])

  wave        =  wave * u.AA
  vs          =  wave.to(u.Hz, equivalencies=u.spectral())

  Z           =  []

  gunn        =  ['u', 'g', 'r', 'i', 'z', 'y']
  classic     =  ['U', 'G', 'R', 'B', 'V', 'I']
  
  for i, x in enumerate(flux):
    interim   =  x * u.erg / u.s / u.cm / u.cm / u.AA
    Fv        =  interim.to(u.erg / u.s / u.cm / u.cm / u.Hz, equivalencies = u.spectral_density(wave))
    Fv       *=  1.e-17
    
    mags      =  get_appmags(vs.value[::-1], Fv.value[::-1], filters, printit = False)

    if meta['AGE'][i] < 0.3:
      Z.append([mags[i] for i in ['g', 'r', 'i', 'z']])  ##  gunn

  Z           =  np.array(Z)
  Zerr        =  []

  for row in Z:
    row_errs  = []
    
    for band in row:
      _merr   = merr(25.0, 26.0, estar=0.2, alphab=-0.25, alphaf=0.22)
      row_errs.append(_merr)

    Zerr.append(np.diag(np.array([x ** 2. for x in row_errs])))

  Zerr = np.array(Zerr)

  print(Z)
  print(Zerr)

  ##  Compute and save the results.
  ##  @pickle_results("XD_toy.pkl")
  def compute_XD_results(n_components=10, max_iter=500):
    clf = XDGMM(n_components, max_iter=max_iter, tol=1e-03, verbose=False, random_state=None)
    clf.fit(Z, Zerr)

    return  clf

  clf    = compute_XD_results(25, 100)
  sample = clf.sample(ngal)

  ##  Plot the results
  fig    = plt.figure()

  fig.subplots_adjust(left=0.1, right=0.95,
                      bottom=0.1, top=0.95,
                      wspace=0.02, hspace=0.02)

  ax1 = fig.add_subplot(221)
  ax1.scatter(Z[:,0], Z[:,1], s=4, lw=0, c='k')

  ax2 = fig.add_subplot(222)
  ##  ax2.scatter(x, y, s=4, lw=0, c='k')
  
  ax3 = fig.add_subplot(223)
  ax3.scatter(sample[:,0], sample[:,1], s=4, lw=0, c='k')

  ax4 = fig.add_subplot(224)
  ##  ax4.scatter(sample[:,0], sample[:,1], s=4, lw=0, c='k', alpha=0.1)

  print('\n\nFitted coefficients:\n\n')

  for i in range(clf.n_components):
    print(clf.mu[i], clf.V[i])

    draw_ellipse(clf.mu[i], clf.V[i],  ax=ax4, scales=[2],
                 ec='k', fc='gray', alpha=0.2)

  titles = ['True Distribution', 'Noisy Distribution',
            'Extreme Deconvolution\n  resampling',
            'Extreme Deconvolution\n  cluster locations']

  ax     = [ax1, ax2, ax3, ax4]
  
  for i in range(4):
    ax[i].xaxis.set_major_locator(plt.MultipleLocator(4))
    ax[i].yaxis.set_major_locator(plt.MultipleLocator(5))

    ax[i].text(0.05, 0.95, titles[i], ha='left', va='top', transform=ax[i].transAxes)

    if i in (0, 1):
        ax[i].xaxis.set_major_formatter(plt.NullFormatter())

    else:
        ax[i].set_xlabel('x')

    if i in (1, 3):
        ax[i].yaxis.set_major_formatter(plt.NullFormatter())

    else:
        ax[i].set_ylabel('y')
  
  plt.show()
