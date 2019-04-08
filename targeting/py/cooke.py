import  astropy.io.fits    as  fits
import  numpy              as  np
import  pylab              as  pl
import  matplotlib.pyplot  as  plt


def isCooke_udrop(u, g, r, i, z):
  ##  Page 4 of https://arxiv.org/pdf/1305.0562.pdf.                                                                                                        
  ##  Eqns. (1) - (6).                                                                                                                                      
  for cut in [(u-g) > 0.7,\
              (u-g) > 1.2 * (g-r) + 0.9,\
               -1.  < (g-r),\
              (g-r) < 1.0,\
              (u-g) > (g-i) + 0.7,\
               -1.  < (g-i),\
              (g-i) < 1.3,\
              (r-i) < 0.4]:

    if not cut:
      return  False

  ##  All cuts passed, return True.                                                                                                                         
  return  True

'''
--  Sextractor binary flags --
     1 = object has neighbors or bad pixels that bias the photometry
     2 = object originally blended with another one
     4 = at least one pixel of the object is saturated
     8 = object truncated (too close to an image boundary)
    16 = aperture data are incomplete or corrupted
    32 = isophotal data are incomplete or corrupted

--  TERAPIX binary flags --:
     1 = star
     2 = saturated star in one of the filters
     4 = masked region in z-band
     8 = masked region in y-band
    16 = masked region in i-band
    32 = masked region in r-band
    64 = masked region in g-band
   128 = masked region in u*-band
'''

if __name__ == '__main__':
  print('\n\nWelcome to a Cooke calculator.\n\n')

  data             = fits.open('/global/cscratch1/sd/mjwilson/CFHTLS/asu.fit')  ##  data[1].header

  u, g, r, i, z, y = data[1].data['umag'], data[1].data['gmag'], data[1].data['rmag'],\
                     data[1].data['imag'], data[1].data['zmag'], data[1].data['ymag'] 

  ##  i-band limited catalogue. 
  imax             = 24.5
  nsig             =  1.5

  ##  Total deep field area [deg].                                                                                                                       
  Area             =  4.0

  ##  Position 
  ra               = data[1].data['RAJ2000']
  dec              = data[1].data['DEJ2000']

  ##  Star-galaxy cut (SEXTRACTOR stellarity).
  sg               = data[1].data['gcl']

  ## 
  sfl              = data[1].data['sfl']

  ## 
  Tfl              = data[1].data['Tfl']

  ##  Maximum extinction of 0.033
  EBV              = data[1].data['E_B-V_']

  ##  Quality cut flag.
  good             = (sfl == 0.0) & (Tfl == 0.0)

  ##  Stellarity
  stellarity       = sg < 0.9

  ##  Mag. lim. cut.                                                                                                                                        
  ilim             = good & (i <= imax)

  is_udrop         = good & [isCooke_udrop(u[k], g[k], r[k], i[k], z[k]) for k, _ in enumerate(u)]
  
  ##  Cooke++ Ly-alpha EW cut; Table 2 of https://arxiv.org/pdf/1305.0562.pdf; Data rows for i < 25.5 (first panel). 
  sige             = 0.21
  siga             = 0.25

  aLBG             = good & ((g-i) >= (0.38 * i - 8.9 + nsig * sige))
  eLBG             = good & ((g-i) <= (0.38 * i - 8.9 - nsig * siga))

  ##  --  LSST Y-10 cut -- 
  ## 
  ##  cuts.append(u < 26.1)
  ##  cuts.append(g < 26.8)
  ##  cuts.append(r < 27.2)
  ##  cuts.append(i < 27.0)
  ##  cuts.append(z < 25.7)

  N                = len(data[1].data['umag'])                                                                                      
  Ni               = len(data[1].data['umag'][ilim])
  Nu               = len(data[1].data['umag'][is_udrop])
  Niu              = len(data[1].data['umag'][is_udrop & ilim])
  Na               = len(data[1].data['umag'][is_udrop & aLBG])
  Ne               = len(data[1].data['umag'][is_udrop & eLBG])
  Nia              = len(data[1].data['umag'][is_udrop & ilim & aLBG])
  Nie              = len(data[1].data['umag'][is_udrop & ilim & eLBG])

  print('\n\nNumber of galaxies in CFHTLS Deep: %d'           % N)
  print('Number of iAB < %.2lf galaxies in CFHTLS Deep: %d'   % (imax, Ni))
  print('Number of u-dropouts in CFHTLS Deep: %d'             % Nu)
  print('Number of iAB < %.2lf u-dropouts in CFHTLS Deep: %d' % (imax, Niu))                     
  print('Number of aLBG in CFHTLS Deep: %d'                   % Na)
  print('Number of eLBG in CFHTLS Deep: %d'                   % Ne)
  print('Number of iAB < %.2lf aLBGs in CFHTLS Deep: %d'      % (imax, Nia))
  print('Number of iAB < %.2lf eLBGs in CFHTLS Deep: %d'      % (imax, Nie))

  ##  Scatter..
  pl.scatter(i[is_udrop], g[is_udrop] - i[is_udrop], rasterized=True, c='k', marker='^', s=4, alpha=0.3)

  ##  Colour cut. 
  ii = np.arange(22., 29., 0.1)

  pl.plot(ii, 0.38 * ii - 8.9)

  pl.plot(ii, 0.38 * ii - 8.9 + nsig * sige, label='%.1lf aLBG per sq. deg.' % (Nia / Area))
  pl.plot(ii, 0.38 * ii - 8.9 - nsig * siga, label='%.1lf eLBG per sq. deg.' % (Nie / Area))

  pl.axvline(x=imax, ymin=0., ymax=1., c='k')

  pl.xlim( 23.,  27.)
  pl.ylim(-0.8, 1.5)

  pl.xlabel(r'$i$')
  pl.ylabel(r'$(g-i)$')

  pl.legend(loc=4)

  pl.show()

print('\n\nDone.\n\n')
