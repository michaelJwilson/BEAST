import numpy as np
import pylab as pl


def Gauss(x, mu, sigma):
    return  np.exp(-(x-mu)**2. / 2. / sigma / sigma) / np.sqrt(2. * np.pi) / sigma 


if __name__ == '__main__':
    lambdas  = np.arange(3.e3, 1.e4, 1.)

    zs       = np.arange(1.5, 5., 1.0)

    LINEFLUX = 1.0 ##  1.55E-17 

    for z in zs:
      lambda0  = 1215.24 * (1. + z)
      sigma    = lambda0 * (600. / 2.9979e5)

      pl.plot(lambdas, LINEFLUX * Gauss(lambdas, lambda0, sigma))

      print(Gauss(lambdas, lambda0, sigma).max())

    pl.ylabel('erg/s/cm$^2$')
    pl.show()

    print('\n\nDone.\n\n')
