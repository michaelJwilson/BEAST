import   numpy                  as       np
import   pylab                  as       pl

from     astropy.cosmology      import   FlatLambdaCDM


##  Default.
h      =   0.7
Om     =   0.3
Tcmb   =   2.725

'''
##  Last column of https://arxiv.org/pdf/1807.06209.pdf
h      =   0.6766
Om     =   0.3111

##  https://en.wikipedia.org/wiki/Cosmic_microwave_background
Tcmb   =   2.7255
'''
'''
##  BX 'drops'.
ttype  = 'BX'
ngal   =  2000   ## [per sq. deg.]
 
loz    =  2.00
hiz    =  2.50
'''
'''
##  u-drops
ttype  = 'u'
ngal   =  500.
loz    =  2.5
hiz    =  3.5
'''
'''
##  g-drops.
ttype  =  'g'
ngal   =   300.  ## [per sq. deg.]
loz    =   3.5
hiz    =  4.25
'''

##  r-drops
ttype  =   'r'
ngal   =   100.  ## [per sq. deg.]                                                                                                             
loz    =   4.5                                                                                                                                  
hiz    =   5.2    


print('\n\nSolving for %s:' % ttype)

cosmo  =  FlatLambdaCDM(H0=100. * h, Om0=Om, Tcmb0=Tcmb)

hichi  =  h * cosmo.comoving_distance(hiz)   ## [Mpc / h]
lochi  =  h * cosmo.comoving_distance(loz)   ## [Mpc / h]

fsky   =  14000. / 41252.961

##  number of galaxies per sq. deg. to number in footprint. 
ngal  *=  14000.                             

vol    =  4. * np.pi * (hichi ** 3. - lochi ** 3.) / 3. 
vol   *=  fsky

print('\n\nvolume \t\t ngal [M] \t ngal per cubic Mpc/h \t ngal per sq deg. to saturate fNL (10^-4).')
print('%.4le \t %d \t\t %.4le \t\t %.4lf' % (vol.value, ngal / 1.e6, ngal / vol.value, 1.e-4 * vol.value / 14000.))

print('\n\nDone.\n\n')
