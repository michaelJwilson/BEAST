import   numpy                  as       np
import   pylab                  as       pl
from     astropy.cosmology      import   FlatLambdaCDM


h     = 0.7
z     = 3.8

##  Simple app. mag. selection of ly-a emitters from g-drops. 
cosmo = FlatLambdaCDM(H0=100. * h, Om0=0.3, Tcmb0=2.725)

##  Ono eqn. (13) https://arxiv.org/pdf/1704.06004.pdf
def MUV(m, z):
  DL    = cosmo.luminosity_distance(z).value  ## Mpc                                                                                                                                       
  DL   *= 1.e6                                ##  pc
  DL   *= h                                   ## [h^-1 pc]      
  
  return  m + 2.5 * np.log10(1. + z) - 5. * np.log10(DL / 10.)


if __name__ == '__main__':
  print('\n\nWelcome.\n\n')

  for m in np.arange(23., 26.5, 0.25):
    print('%.4lf \t %.4lf' % (m, MUV(m,z)))

  print('\n\nDone.\n\n')
    
