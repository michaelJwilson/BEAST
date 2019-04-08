import  os
import  numpy                        as      np
import  astropy.units                as      u
import  pylab                        as      pl
import  itertools

from    utils                        import  latexify
from    prep_filters                 import  prep_filters
from    BC03_maker                   import  BC03_maker
from    app_mags                     import  get_appmags
from    sklearn                      import  linear_model
from    sklearn.metrics              import  mean_squared_error, r2_score
from    isubaru                      import  bands                          as  isubaru
from    sompy.sompy                  import  SOMFactory
from    sompy.visualization.mapview  import  View2D
from    sompy.visualization.bmuhits  import  BmuHitsView


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
    ##  Z.append([mags[i] for i in gunn + isubaru])
    Z.append([mags[i] for i in gunn])

##
Z           = np.array(Z)

print(Z)
print('\n\n')

sm    = SOMFactory().build(Z, normalization = 'var', initialization='random', component_names=gunn)
sm.train(n_job=1, verbose=False, train_rough_len=2, train_finetune_len=5)

topographic_error  = sm.calculate_topographic_error()
quantization_error = np.mean(sm._bmu[1])

print ("Topographic error = %s; Quantization error = %s" % (topographic_error, quantization_error))

vhts    = BmuHitsView(10, 10, 'Hits Map', text_size=7)
vhts.show(sm, anotate=True, onlyzeros=False, labelsize=12, cmap='Greys', logaritmic=False)

pl.show()

view2D  = View2D(10, 10, 'rand data', text_size=10)
view2D.show(sm, col_sz=4, which_dim='all', denormalize=True)

pl.show()
