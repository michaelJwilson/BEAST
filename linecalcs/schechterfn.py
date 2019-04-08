import  pylab  as  pl
import  numpy  as  np

from    scipy.special      import  gamma, gammaincc
from    astropy.cosmology  import  FlatLambdaCDM


cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)

'''                                                                                                                                                                                                                                
See https://pythonhosted.org/Astropysics/coremods/models.html                       
                                                                                                                                                                                                                                  
https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gammaincc.html#scipy.special.gammaincc
https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.special.gamma.html
'''

def SchechterLfn(L, phi_star, L_star, alpha):                                                                                                                                                                             
  ##  Schechter function in luminosity.                                                                                                                                                                                                    

  LL = L / L_star

  return phi_star * (LL ** alpha) * np.exp(-LL) / L_star               ## Phi(L) wrt dL.

def Sobral_Lya(z, Lmin=None, printit=False, h=0.7):
    ##  Table 6 of https://arxiv.org/pdf/1712.04451.pdf
    ##  redshift, alpha, log10(L*) [erg/s] log10(Phi*) [Mpc ** -3]
    fits       =  np.array([[2.2, -2.0, 42.82, -3.59], [2.5, -1.72, 42.71, -3.10], [3.1, -1.63, 42.77, -3.06], [3.9, -2.26, 42.93, -3.66], [4.7, -2.35, 43.28, -4.25], [5.4, -1.98, 43.28, -3.83]])
    ind        =  np.where(np.abs(fits[:,0] - z) == np.abs(fits[:,0] - z).min())

    alpha      =         fits[ind,1]
    L_star     =  10. ** fits[ind,2]
    phi_star   =  10. ** fits[ind,3] / h ** 3.

    if Lmin == None:
      Lmin     =  L_star

    else:
      Lmin     =  10. ** Lmin

    zmin       =  np.log10(Lmin / L_star)

    zs         =  np.arange(zmin, zmin + 4.e0, 0.1)
    xs         =  zs * np.log(10)
    ys         =  10. ** zs

    dz         =  zs[1] - zs[0]

    nbar       =  phi_star * np.exp(-ys) * 10. **((1. + alpha) * zs)
    nbar       =  np.sum(nbar) * dz * np.log(10)

    meanlum    =  phi_star * np.exp(-ys) * 10. **((1. + alpha) * zs) * ys * L_star
    meanlum    =  np.sum(meanlum) * dz * np.log(10)
    meanlum   /=  nbar
    meanlum    =  np.log10(meanlum)
    
    if printit:
      print('%.3lf  %.3lf  %.3lf  %.3le  %.3le  %.3le  %.3le' % (z, fits[ind,0], alpha, L_star, phi_star, nbar, meanlum))

    return  [nbar, meanlum]

def dVols(zs, cosmo, h=0.7):
  Vs    = cosmo.comoving_volume(zs).value                                           ##  Get volume to each redshift slice.                                                                                                                                                                                              
  dVs   = Vs - np.roll(Vs, 1)
  dVs  *= h **3.                                                                    ## [h^-1 Mpc]^3                                                                                                                                                                                                                   

  zs    = zs[:-1]
  dVs   = dVs[1:]

  return zs, dVs
   
def Mhalo(Lya, z):
  ## eqn. (12) of https://arxiv.org/pdf/1811.00556.pdf
  ## implicit dependence on redshift via L*

  fits       =  np.array([[2.2, -2.0, 42.82, -3.59], [2.5, -1.72, 42.71, -3.10], [3.1, -1.63, 42.77, -3.06], [3.9, -2.26, 42.93, -3.66], [4.7, -2.35, 43.28, -4.25], [5.4, -1.98, 43.28, -3.83]])
  ind        =  np.where(np.abs(fits[:,0] - z) == np.abs(fits[:,0] - z).min())

  L_star     =  10. ** fits[ind,2]
  result     =  10. ** 11.91 * ((10. ** Lya) / L_star) ** 1.44 

  return  np.log10(result)

def tinker_bias(nu, Delta=200):
    ## Tinker bias fn. (https://arxiv.org/pdf/1001.3162.pdf)                                                                                                                                                                                                                                                            
    y = np.log10(Delta)

    A = 1.0 + 0.24 * y * np.exp(- (4/y)**4.)
    a = 0.44 * y - 0.88
    B = 0.183
    b = 1.5
    C = 0.019 + 0.107 * y + 0.19 * np.exp(-(4./y)**4.)
    c = 2.4

    return  1. - A * nu ** a /(nu ** a + 1.686 ** a) + B * nu ** b + C * nu **c


if __name__ == '__main__':
    print('\n\nWelcome to Schecter.\n\n')

    zs       = np.arange(2., 6.,  0.1)
    dVs      = dVols(zs, cosmo, h=0.7)

    results  = [] 

    ## nus          =  1.686 / sigmas                                                                                                                                                                                                                                                                               
    ## tinker_bs    =  tinker_bias(nus, Delta=200)  

    for z in zs:
      nbar, meanlum = Sobral_Lya(z, 42.5, printit=True)
      mhalo         = Mhalo(meanlum, z)

      results.append([np.log10(nbar), meanlum, mhalo])

    results = np.array(results)
      
    pl.plot(zs, results[:,0], label='density [$(h^{-1} \rm{Mpc})^3$]')
    ## pl.plot(zs, results[:,1], label='<L> [ergs/s]')
    ## pl.plot(zs, results[:,2], label=r'Mhalo [M$_\odot / h$]') 

    pl.xlabel(r'$z$')
    ## pl.ylabel(r'$\bar n(z, L > L_*)$')
    pl.yscale('linear')
    pl.legend()
    pl.show()

    print('\n\nDone.\n\n')
