import  os
import  numpy           as      np
import  astropy.units   as      u
import  pylab           as      pl
import  itertools

from    utils           import  latexify
from    prep_filters    import  prep_filters
from    BC03_maker      import  BC03_maker
from    app_mags        import  get_appmags
from    sklearn         import  linear_model
from    sklearn.metrics import  mean_squared_error, r2_score
from    isubaru         import  bands                          as  isubaru


##  latexify(fig_width=None, fig_height=None, columns=2, equal=False)                                                                                     

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
filters     =  prep_filters(['LSST', 'STEIDEL', 'SUBARU', 'JKC', 'HUBBLE', 'ISUBARU'])

wave        =  wave * u.AA
vs          =  wave.to(u.Hz, equivalencies=u.spectral())

Z           = []

gunn        =  ['u', 'g', 'r', 'i', 'z', 'y']
classic     =  ['U', 'G', 'R', 'B', 'V', 'I']

for i, x in enumerate(flux):    
  interim   =  x * u.erg / u.s / u.cm / u.cm / u.AA
  Fv        =  interim.to(u.erg / u.s / u.cm / u.cm / u.Hz, equivalencies = u.spectral_density(wave))
  Fv       *=  1.e-17

  mags      =  get_appmags(vs.value[::-1], Fv.value[::-1], filters, printit = False)

  if meta['AGE'][i] < 0.3:
    ##  Z.append([mags[i] for i in gunn + isubaru])
    Z.append([mags[i] for i in gunn + classic])

##
Z           = np.array(Z)

print(Z)
print('\n\n')

xs          = np.arange(20., 30., 0.01)
pl.plot(xs, xs, 'k-', lw=0.1)

for kk, band in enumerate(['g']):
  '''
  Y         = Z[:, kk].reshape(-1, 1)    ##  e.g. U-band.
  
  ## -- Classic -- ##  
  ##  Fig. 1 of https://www.aanda.org/articles/aa/pdf/2006/46/aa6082-06.pdf
  if band in  ['U', 'B']:
    inbands = ['u', 'g']
    X       = Z[:, 6:8]                  ##  LSST u and g.

  if band in ['V']:
    inbands = ['g', 'r']
    X       = Z[:, 7:9]                  ##  LSST g and r.

  if band in ['R']:
    inbands = ['r', 'i']
    X       = Z[:, 8:10]                 ##  LSST r and i. 

  if band in ['I']:
    inbands = ['i', 'z']
    X       = Z[:, 9:11]                 ##  LSST i and z. 
  '''
  
  if band in ['g']:
    ##  g-proxy based on Subaru intermediate bands. 
    '''
    ##  No 484 or 527 in COSMOS.
    inbands = ['I427', 'I464', 'I505', 'I574']

    Y       = Z[:, 1].reshape(-1, 1)  
    X       = Z[:, [6, 7, 9, 11]]        ##  I:  427, 464, 484, 505, 527, 574.
    '''

    ##  g-proxy based on (B-V).
    inbands = ['B', 'V']

    Y       = Z[:, 1].reshape(-1, 1)                                                                                                                         
    X       = Z[:, [9, 10]]              ##  Change Z.append above.  

  '''
  ##  Test set. 
  X         = np.c_[np.arange(10), np.arange(10) ** 2., np.arange(10) ** 3., np.random.uniform(-0.5, 0.5, 10)]
  Y         = X[:,1] + 2. * X[:,0] + 3. * X[:,2] + np.random.uniform(-0.1, 0.1, 10)
  '''

  ##  Linear regression fit.
  reg       = linear_model.LinearRegression(fit_intercept=True).fit(X, Y)

  print('\n\n%s-band Score: %.3lf' % (band, reg.score(X, Y)))
  print('Regression:  ' + '  '.join('%.6lf * %s' % (x, inbands[ii]) for ii, x in enumerate(reg.coef_[0]))\
                        + ' + %.6lf (%.6lf)' % (reg.intercept_, np.std(Y - reg.predict(X))))

  ##  print('Regression intercept:  %.6le' % reg.intercept_)
  ##  print('RMS residual: %.6le' % np.std(Y - reg.predict(X)))

  prediction  = reg.predict(X) 

  ##  for i, row in enumerate(X):
  ##    print(Y[i], prediction[i])

  pl.plot(Y, prediction, 'o', label=band, markersize=2)

pl.xlim(22.,  27.)
pl.ylim(22.,  27.)

##  pl.xlabel(r'Classic')
##  pl.ylabel(r'|Classic - LSST|')

pl.legend(loc=3, ncol=3)

pl.show(block=True)

print('\n\nDone.\n\n')
