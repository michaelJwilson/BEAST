#!/usr/bin/env python
# By Simone Ferraro (Princeton) - simo.ferraro@gmail.com

# NOTE FOR CAMB: the linear matter power spectrum and the zeta power are related by 
# Delta^2matter(kvec) = littleh**4*kvec**4*Tcamb(kvec)**2*Delta2zeta(kvec)

# NOTE ON UNITS: k is in h/Mpc, masses are in M_sun/h

import sys
from   numpy import *
import numpy as np
import scipy
from   scipy import interpolate
from   scipy import integrate
import pylab as pl


# *********** SURVEY PARAMS **********

##  BM/BX
drop_type   = 'bmbx'
bg_L        = 1.5                # Lagrangian bias
V           = 2.9 * 1.e10        # Survey Volume in (Mpc/h)^3
zm          = 2.2                # Median redshift
deltaP0oP0  = 0.2                # Fractional uncertainty on constant shotnoise component to Power Spectrum.
kmax        = 0.1
'''
##  u-dropouts
drop_type   = 'u'
bg_L        = 3.0                # Lagrangian bias                                                                                                       
V           = 5.7 * 1.e10        # Survey Volume in (Mpc/h)^3                                                                                              
zm          = 3.0                # Median redshift                                                                                                         
deltaP0oP0  = 0.2                # Fractional uncertainty on constant shotnoise component to Power Spectrum.                                               
kmax        = 0.1
'''
'''
##  g-dropouts
drop_type   = 'g'
bg_L        = 2.50               # Lagrangian bias                                                                                                          
V           = 3.95 * 1.e10       # Survey Volume in (Mpc/h)^3                                                                                               
zm          = 3.75               # Median redshift                                                                                                          
deltaP0oP0  = 0.2                # Fractional uncertainty on constant shotnoise component to Power Spectrum.                                                
kmax        = 0.1
'''
'''
##  r-dropouts
drop_type   = 'r'
bg_L        = 4.5                # Lagrangian bias                                                                                                           
V           = 3.3 * 1.e10        # Survey Volume in (Mpc/h)^3                                                                                               
zm          = 4.75               # Median redshift                                                                                                          
deltaP0oP0  = 0.2                # Fractional uncertainty on constant shotnoise component to Power Spectrum.                                                
kmax        = 0.1
'''

# *********** Change below at your own risk. ********** 
bg_E        = 1. + bg_L               # Eulerian bias
am          = 1./(1 + zm)
kmin        = 2. * np.pi / V**(1./3)

# *********** COSMO PARAMS ***********
(kminc, kmaxc) = (0.0001, 4.)

Aphi    = 1.56e-8   # Check!
littleh = 0.67
k0      = 0.05 / littleh
ns      = 0.9619
zeq     = 3403
Omegam  = 0.3183
Omegar  = Omegam / (1. + zeq)
cc      = 2.99792458e5   # Speed of light in Km / s 
Omegal  = 1. - Omegam - Omegar
rhobar  = 2.78e11 * Omegam # in M_sun/h * h^3/Mpc^3
aeq     = 1. / (1. + zeq)
deltac  = 1.42

# ************************************

n       = np.logspace
a       = np.loadtxt('PLANCK_matterpower0.dat')
splinea = interpolate.splrep(np.log(a[:,0]), np.log(a[:,1]), s=0)

b       = np.loadtxt('PLANCK_transfer0.dat')
splineb = interpolate.splrep(np.log(b[:,0]), np.log(b[:,1]), s=0)

def growth(z): # Normalized such that D(a)=a in matter domination
    a = 1. / (1. + z)
    (amin, amax) = (0., a)
    return 2.5 * Omegam * E(z) * integrate.quad(lambda aact: pow(aact * E(1. / aact - 1), -3.), amin, amax, epsabs = 0, epsrel = 1.0e-4)[0]

def Pm(k, z):
    return np.exp(interpolate.splev(np.log(k), splinea, der=0)) * (growth(z) / growth(0.))**2

def chi(z): # Comoving distance to redshift z (chi = 0 at z = 0). In Mpc/h
    (zmin, zmax) = (0., z)
    return (cc / 100.) * integrate.quad(lambda zz: 1. / E(zz) , zmin, zmax, epsabs=0, epsrel=1.0e-5)[0]

def Tcamb(k):
    return np.exp(interpolate.splev(np.log(k), splineb, der=0))
    
def T(k):
    return Tcamb(k) / Tcamb(k=0.0001)

def E(z):   # H(z) / H0. Normalized such that E(z=0) = 1
	return np.sqrt(Omegam * (1. + z)**3 + Omegar * (1. + z)**4 + Omegal)

def alpha(k, z):
	return 5./3. * littleh**2 * k**2 * Tcamb(k) * growth(z)/growth(0.)


# *************** NEW WAY *****************

keq = aeq * (100. / cc) * E(zeq)    # In h / Mpc

def R0(z):
    return np.sqrt(2./3. * growth(z) / Omegam) * (cc / 100.)

def n0(z):
    return 1. / (Aphi * R0(z)**3)

def n1(z):
    return n0(z) / (keq * R0(z))
    
@np.vectorize
def k1(z, nbar):
    return  keq * (bg_E**2 * nbar) / n1(z)

@np.vectorize
def gamma(z, nbar):
    return (bg_E**2 * nbar) / n1(z)

def fa(z):
    return [4. * deltac * bg_E * bg_L * Aphi * R0(z)**2, 2. * bg_E * Aphi * R0(z)**4, 1.]  

def ga(k):
    return [1./k * T(k) * (k/k0)**(ns-1.), k * T(k)**2 * (k/k0)**(ns-1.), 1]

FprimeMat = np.matrix(np.zeros([3, 3]))
priorMat  = np.matrix(np.zeros([3, 3]))
FMat      = np.matrix(np.zeros([3, 3]))


'''
# VOLUME DEPENDENCE 
Vvec = np.logspace(8, 11, 15)
kminvec = 2. * np.pi / Vvec**(1./3)

FprimeInvfnl = np.zeros(len(kminvec))

# FIT fnl as a function of k1
K = 0
for kmin in kminvec:
    for I in range(3):
        for J in range(3):
            FprimeMat[J, I] = 1./(2. * np.pi**2) * integrate.quad(lambda kk: kk**2 * ga(kk)[I] * ga(kk)[J] / (kk * T(kk)**2 * (kk/k0)**(ns-1.) + keq**2/k1(zm, nbar))**2, kmin, kmax, epsabs = 0, epsrel = 1.0e-5)[0]
    priorMat[2,2] = 1. / (deltaP0oP0 * 1./nbar)**2
    #FprimeMat = FprimeMat + priorMat
    FprimeInv = FprimeMat.I
    print FprimeInv
    FprimeInvfnl[K] = FprimeInv[0,0]  # Marginalization here
    K += 1

# Fit relation between V and F'
splineV = interpolate.splrep(np.log(Vvec), np.log(FprimeInvfnl), s=0)

V0 = 5.0e9
def logFprimeoflogV(logV):
    return interpolate.splev(logV, splineV, der=0)

logh = np.log(1.0e-1)
epsilon = (logFprimeoflogV(np.log(V0) + logh) - logFprimeoflogV(np.log(V0) - logh)) / (2. * logh)
run = (logFprimeoflogV(np.log(V0) + logh) + logFprimeoflogV(np.log(V0) - logh) - 2. * logFprimeoflogV(np.log(V0))) / logh**2

print epsilon, run

pl.loglog(Vvec / 1.0e9, FprimeInvfnl)
pl.loglog(Vvec / 1.0e9, np.exp(logFprimeoflogV(np.log(Vvec))))
pl.loglog(Vvec / 1.0e9, np.exp(logFprimeoflogV(np.log(V0))) * (Vvec / V0)**(epsilon + 0.5 * run * np.log(Vvec / V0)))
pl.xlabel(r'$V$   [Gpc/$h$]$^3$')
pl.ylabel(r'Fprime')
pl.show()
'''

V0 = 5.0e9

# FIXING THE n DEPENDENCE 

epsSV = -0.330044769793
runSV = 0.0102514428279
epsP  = -0.204984698966
runP  = 0.0731393778569

nbarfitvec = np.logspace(-9, -1, 20)
gammavec   = gamma(zm, nbarfitvec)

FprimeInvfnl = np.zeros(len(nbarfitvec))
Sigmafnl     = np.zeros(len(nbarfitvec))

K = 0

for nbar in nbarfitvec:
    for I in range(3):
        for J in range(3):
            FprimeMat[J, I] = 1./(2. * np.pi**2) * integrate.quad(lambda kk: kk**2 * ga(kk)[I] * ga(kk)[J] / (kk * T(kk)**2 * (kk/k0)**(ns-1.) + keq**2/k1(zm, nbar))**2, kmin, kmax, epsabs = 0, epsrel = 1.0e-4)[0]
            FMat[J, I] = 0.5 * V * fa(zm)[I] * fa(zm)[J] / (bg_E**4 * Aphi**2 * R0(zm)**8) * FprimeMat[J, I]
    priorMat[2,2] = 1. / (deltaP0oP0 * 1./nbar)**2
    FprimeInv = FprimeMat.I
    print FprimeInv
    FMat = FMat + priorMat 
    Sigmafnl[K] = np.sqrt(FMat.I[0,0])
    FprimeInvfnl[K] = FprimeInv[0,0]
    K += 1

# np.savetxt('fnl_marg.txt', (nbarfitvec, Sigmafnl))

# Excellent fit to F' in both cases!!
'''
AP  = 6.70597135e-01
ASV = 6.35745798e-02
V0  = 5.0e9

FitPvec  = AP * gammavec**(-2) * (V / V0)**(epsP + 0.5 * runP * np.log(V / V0))
FitSVvec = ASV * np.ones(len(gammavec))  * (V / V0)**(epsSV + 0.5 * runSV * np.log(V / V0))


pl.loglog(gammavec, FprimeInvfnl)
pl.loglog(gammavec, FitPvec)
pl.loglog(gammavec, FitSVvec)
pl.show()
'''


# Wrapping up:  Sigma(fnl)

def unorm_fit_sigma_SV(z, V, gamma):
    return growth(z) * bg_E / (bg_E - 1)  # FOR V = V0 !!!!!!!!!
     
A_fit_SV = Sigmafnl[-1] / unorm_fit_sigma_SV(zm, V0, gammavec[-1]) 
print A_fit_SV


def unorm_fit_sigma_P(z, V, gamma):
    return bg_E / (bg_E - 1) * growth(z) / gamma# FOR V = V0 !!!!!!!!!
     
A_fit_P = Sigmafnl[0] / unorm_fit_sigma_P(zm, V0, gammavec[0]) 
print A_fit_P
'''
'''
A_fit_SV = 15.9095384565
A_fit_P = 54.2249401469

fid_SV = -2./3
fid_P = -1./2

epsprimeSV = 0.5*(epsSV - 1.) - fid_SV
epsprimeP = 0.5*(epsP - 1.) - fid_P
runprimeSV = runSV /2
runprimeP = runP /2

@np.vectorize
def fit_sigma_SV(z, V, gamma):
    return A_fit_SV * growth(z) * bg_E / (bg_E - 1) * (V / V0)**(fid_SV + epsprimeSV + 0.5 * runprimeSV * np.log(V/V0))

def fit_sigma_P(z, V, gamma):
    return A_fit_P * growth(z) * bg_E / (bg_E - 1) * (V / V0)**(fid_P + epsprimeP + 0.5 * runprimeP * np.log(V/V0)) / gamma

##  print epsprimeSV, epsprimeP, runprimeSV, runprimeP

_nbars = gammavec * n1(zm) / bg_E**2.

pl.loglog(_nbars, Sigmafnl)
##  pl.loglog(gammavec, fit_sigma_SV(zm, V, gamma(zm, nbarfitvec)) + fit_sigma_P(zm, V, gamma(zm, nbarfitvec)))
##  pl.loglog(gammavec, fit_sigma_P(zm, V, gamma(zm, nbarfitvec)))

pl.xlim(1.2e-2, 5.e-5)
pl.ylim(0.1, 10.)

pl.ylabel(r'$\sigma_{f_{NL}}$')
pl.xlabel(r'$\bar n \ [(h^{-1} \rm{Mpc})^3]$')

##  pl.savefig('plots/fnl.pdf')

##  pl.semilogx(gammavec, (fit_sigma_SV(zm, V, gamma(zm, nbarfitvec))**0.9 + fit_sigma_P(zm, V, gamma(zm, nbarfitvec))**0.9)**(1./0.9)/Sigmafnl -1.)
pl.show()

np.savetxt('dat/fnl_%sdrop.dat' % drop_type, np.c_[_nbars, Sigmafnl], fmt='%.6le')
