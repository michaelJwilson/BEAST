import numpy             as np
import pylab             as pl
import matplotlib.pyplot as plt

from   scipy.optimize    import curve_fit

def func(x, a, b, c, d):
  return  a + b * np.exp((x - c)**2. / d) 

dat        = np.loadtxt('result_Q0.txt')


popt, pcov = curve_fit(func, dat[:,0], dat[:,1])

diff       = dat[:,1] / func(dat[:,0], *popt)

##  plt.plot(dat[:,0], dat[:,1], 'r-')
##  plt.plot(dat[:,0], func(dat[:,0], *popt), 'c-')

##  pl.show()

FT         = np.fft.rfft(diff[:-1])

MOD2       = FT * FT
BIG        = MOD2.argsort()[-15:][::-1]
SMALL      = np.arange(len(MOD2))

SMALL      = [x for x in SMALL if x not in BIG]

np.put(FT, SMALL, 0.0)

iFT        = np.fft.irfft(FT)

plt.plot(dat[:,0],  diff, 'g-')
plt.plot(dat[:-1,0], iFT, 'c-')

pl.show()
