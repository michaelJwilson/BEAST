import os
import numpy             as np
import pylab             as pl 
import astropy.io.fits   as fits
import matplotlib.pyplot as plt


'''
Hildebrandt Selection criteria:

--  Colour / colour.
    E.g.  'u_g_ISO_cor2' and 'g_r_ISO_cor2'

--  dat[1].data['MASK_u']      == 0
--  dat[1].data['MASK_g']      == 0
--  dat[1].data['MASK_r']      == 0
--  dat[1].data['MASK_i']      == 0
--  dat[1].data['MASK_z']      == 0
--  dat[1].data['MASK_STARS']  == 0 
--  dat[1].data['masksa']      == 0
--  dat[1].data['maskgscw2']   == 0
--  dat[1].data['CLASS_STAR'] < 0.9
'''

def isCooke_udrop(u, g, r, i, z):
  ##  Page 4 of https://arxiv.org/pdf/1305.0562.pdf.                                                                                                    
  ##  Eqns. (1) - (6).

  for cut in [((u-g) > 0.7), ((u-g) > 1.2 * (g-r) + 0.9), ((-1.  < (g-r)), ((g-r) < 1.0)), ((u-g) > (g-i) + 0.7),\
              ((-1.  < (g-i)), ((g-i) < 1.3)), ((r-i  < 0.4))]: 

    if not cut:
      return  False
    
    else:
      ##  Passed this cut, try the next. 
      pass

  ##  All cuts passed, return True. 
  return  True


##  Load Hendrik's catalogues. 
root   = os.environ['CSCRATCH']
dat    = fits.open(root + '/Hildebrandt/CFHTLS/D1.fits')

##  Good imaging cut. 
good   = (dat[1].data['MASK_u'] == 0) & (dat[1].data['MASK_g'] == 0) & (dat[1].data['MASK_r'] == 0) &\
         (dat[1].data['MASK_i'] == 0) & (dat[1].data['MASK_z'] == 0) & (dat[1].data['MASK_STARS'] == 0) &\
         (dat[1].data['masksa'] == 0) & (dat[1].data['maskgscw2'] == 0) & (dat[1].data['CLASS_STAR'] < 0.9)

##  dat[1].header

umag   = dat[1].data['MAG_ISOCOR_u']
gmag   = dat[1].data['MAG_ISOCOR_g']
rmag   = dat[1].data['MAG_ISOCOR_r']
imag   = dat[1].data['MAG_ISOCOR_i']
zmag   = dat[1].data['MAG_ISOCOR_z']

umg    = dat[1].data['u_g_ISO_cor2'] 
gmr    = dat[1].data['g_r_ISO_cor2']
rmi    = dat[1].data['r_i_ISO_cor2']
imz    = dat[1].data['i_z_ISO_cor2']

'''
plt.scatter(umg[good], umag[good] - gmag[good], c='k', marker='o', rasterized=True)
pl.xlim(-5., 5.)
pl.ylim(-5., 5.)
pl.show()
'''

##  Defines Henrik's dropout samples. 
udrop  = (1.5 < umg) & (-1.0 < gmr) & (gmr < 1.2) & (1.5 * gmr < umg -0.75)
ucat   = good & udrop

gdrop  = (1.0 < gmr) & (-1.0 < rmi) & (rmi < 1.0) & (1.5 * rmi < gmr -0.80)
gcat   = good & gdrop

rdrop  = (1.2 < rmi) & (-1.0 < imz) & (imz < 0.7) & (1.5 * imz < rmi -1.00)
rcat   = good & rdrop

##  Cooke u-dropout sample.
cudrop = good & [isCooke_udrop(umag[i], gmag[i], rmag[i], imag[i], zmag[i]) for i, u in enumerate(umag)]

##  i-lim
ilimd  = (imag <= 26.4)

print(len(umag), len(umag[ucat]), len(umag[cudrop]), len(umag[ucat & ilimd] == True), len(umag[cudrop & ilimd] == True))

##  u-drops.
##  plt.scatter(gmr[ucat],   umg[ucat],   c='k', marker='x', rasterized=True) 
plt.scatter(gmr[cudrop], umg[cudrop], c='k', marker='o', rasterized=True)

pl.xlabel(r'$(g-r)$')
pl.ylabel(r'$(u-g)$')

pl.xlim(-0.5, 4.0)
pl.ylim(-0.5, 4.0)

pl.savefig('udrop_comp.pdf')

pl.clf()

exit(0)

'''
##  g-drops
pl.plot(rmi[gcat], gmr[gcat], 'kx')

pl.xlabel(r'$(r-i)$')
pl.ylabel(r'$(g-r)$')

pl.xlim(-0.5, 4.0)
pl.ylim(-0.5, 4.0)

pl.show()


##  r-drops                                                                                                                                                                                                                                                                            
pl.plot(imz[rcat], rmi[rcat], 'kx')

pl.xlabel(r'$(i-z)$')
pl.ylabel(r'$(r-i)$')

pl.xlim(-0.5, 4.0)
pl.ylim(-0.5, 4.0)

pl.show()
'''
'''
pl.clf()

##  Photometric redshifts.                                                                                                                              
HypZCWW = dat[1].data['Z_PHOT']          ##  HyperZ run with the CWW.
HypZBC3 = dat[1].data['Z_PHOT_BC']       ##  HyperZ run with the BC03.  
BPZ     = dat[1].data['Z_B_V3']          ##  BpZ.

for cat, band in zip([ucat, gcat, rcat], ['u', 'g', 'r']): 
  pl.clf()
  pl.hist(dat[1].data['Z_PHOT'][cat],    bins=25, label='%s-drop HPZ-CWW' % band, alpha=0.4)
  pl.hist(dat[1].data['Z_PHOT_BC'][cat], bins=25, alpha=0.4, label='HPZ-BC03')
  pl.hist(dat[1].data['Z_B_V3'][cat],    bins=25, alpha=0.4, label='BPZ')

  pl.xlabel(r'$z$')
  pl.ylabel(r'$dN/dz$')

  pl.legend()
  pl.show()
'''

sige   = 0.25
siga   = 0.23

lsstu  = 26.1 

ilims  = np.arange(22., 28., 0.25)

whatu  = udrop
##  whatu  = cudrop 

result = []

for ilim in ilims:
  nwhatu     =  len(np.where((whatu == True) & (imag <= ilim))[0])

  ##  Cooke++ Ly-alpha cut.                                                                               
  aLBG       =  whatu & ((gmag-imag) >= (0.38 * imag - 8.9 + 2.0 * siga))
  eLBG       =  whatu & ((gmag-imag) <= (0.38 * imag - 8.9 - 2.0 * sige))
  
  limd_aLBG  =  aLBG[imag <= ilim]
  limd_eLBG  =  eLBG[imag <= ilim]

  nlimd_aLBG =  len(np.where(limd_aLBG == True)[0])
  nlimd_eLBG =  len(np.where(limd_eLBG == True)[0])

  aeff       =  100. * nlimd_aLBG / nwhatu
  eeff       =  100. * nlimd_eLBG / nwhatu

  result.append([aeff, eeff])

  print('%.2lf \t %.1lf \t %.1lf' % (ilim, aeff, eeff))

## 
result = np.array(result)

pl.plot(ilims, result[:,0], label='ALBG frac.')
pl.plot(ilims, result[:,1], label='ELBG frac.')
pl.xlabel(r'$i_{\rm{lim}}$')
pl.ylabel(r'Fraction of $u$-dropouts($i$ < $i_{\rm{lim}}$) $\in$ [ALBG, ELBG]')
pl.legend()

##  pl.show()
##  pl.savefig('Cooke.pdf')

print('\n\nDone.\n\n')
