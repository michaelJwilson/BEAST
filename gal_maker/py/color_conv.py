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


##  latexify(fig_width=None, fig_height=None, columns=2, equal=False)                                                                                     

print('\n\nWelcome to (BC03) LBG maker.')

ngal              =     400  ##  Over written if rest-frame.                                                                                           
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
filters     =  prep_filters(['LSST', 'STEIDEL', 'SUBARU', 'JKC', 'HUBBLE'])
allmags     =  ['U', 'G', 'R', 'B', 'V', 'I', 'u', 'g', 'r', 'i', 'z', 'y']

wave        =  wave * u.AA
vs          =  wave.to(u.Hz, equivalencies=u.spectral())

Z           = []

for i, x in enumerate(flux):    
  interim   =  x * u.erg / u.s / u.cm / u.cm / u.AA
  Fv        =  interim.to(u.erg / u.s / u.cm / u.cm / u.Hz, equivalencies = u.spectral_density(wave))
  Fv       *=  1.e-17

  mags      =  get_appmags(vs.value[::-1], Fv.value[::-1], filters, printit = False)

  if meta['AGE'][i] < 0.3:
    ##  Will cull templates to less than ngal.
    Z.append([mags[i] for i in allmags])

##
Z           = np.array(Z)

print(Z.shape)
print('\n\n')

## 
class_cols  =  ['R-I', 'G-R', 'R-z', 'U-G', ]

for kk, col in enumerate(class_cols):
  bands     = []

  print('Solving for: %s \t %s' % (col[0], col[2]))

  for band in [col[0], col[2]]:
    if band in ['U', 'B']:
      bands.append('u')
      bands.append('g')
    
    if band in ['V']:
      bands.append('g')
      bands.append('r')
      
    if band in ['R']:
      bands.append('r')
      bands.append('i')

    if band in ['I']:
      bands.append('i')
      bands.append('z')

    if band in ['z']:
      bands.append('z')   

  ##  Only need unique bands.
  bands    = list(set(bands))

  ##  And unique colors. 
  new_cols = list(set(itertools.combinations(bands, 2)))

  ##
  mmap     = dict(zip(allmags, np.arange(len(allmags))))

  print(new_cols)

  ##  Data vector to fit:  classic colour. 
  ##  Fig. 1 of https://www.aanda.org/articles/aa/pdf/2006/46/aa6082-06.pdf 
  Y    =  Z[:, mmap[col[0]]]  - Z[:, mmap[col[2]]]
  X    = [Z[:, mmap[ xx[0]]]  - Z[:, mmap[ xx[1]]] for xx in new_cols]

  X    = np.array(X).T

  print(X.shape)
  print(Z.shape)
  print(Y.shape)

  ##  Linear regression fit.
  reg = linear_model.LinearRegression(fit_intercept=True).fit(X, Y)

  print('\n\n%s-band Score: %.3lf'                             % (band, reg.score(X, Y)))
  print('Regression coefficients:  ' + '  '.join('%.6lf %s-%s' % (x, new_cols[ii][0], new_cols[ii][1]) for ii, x in enumerate(reg.coef_)))
  print('Regression intercept:  %.6le'                         % reg.intercept_)
  print('RMS residual: %.6le'                                  % np.std(Y - reg.predict(X)))

  prediction  = reg.predict(X) 

  print('\n\n')

  for i, row in enumerate(X):
    print(Y[i], prediction[i])

  print('\n\n')
    
  pl.plot(Y,          Y, 'k-', label=band, markersize=2, alpha=0.4)
  pl.plot(Y, prediction,  'o', label=band, markersize=2)
  
  ##  pl.xlim(22.,  27.)
  ##  pl.ylim(1e-3, 1e0)
  
  pl.xlabel(col)
  pl.ylabel(r'$\hat{%s}$' % col)

  pl.legend(loc=3, ncol=3)
  
  pl.show(block=True)

print('\n\nDone.\n\n')
