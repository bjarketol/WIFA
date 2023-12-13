#import fonctions_universelles as fct
import os
from os import sep
import numpy as np
import scipy.integrate as integrate
from cmath import *

#Universal functions
def phim_hogstrom_s(z, dlmo):
    zeta = z * dlmo
    if zeta<0.5:
        return 1+4.8*zeta
    elif zeta<10:
        return 7.9 - 4.25/zeta + (1/zeta)**2
    else:
        return 0.7485*zeta
    
def phih_hogstrom_s(z, dlmo):
    zeta = z * dlmo
    return 0.95+7.8*zeta

#Nieuwstadt variable profiles
def psim_vplus_nieuwstadt(fonction, z, zi, dlmo, z0, alpha):
    f = lambda x: fonction(x, dlmo)*pow(1-x/zi, alpha/2-1)*np.sin(np.sqrt(alpha*(alpha+2))/2*np.log(1-x/zi))/x
    return integrate.quad(f, z0, z)[0]

def psim_uplus_nieuwstadt(fonction, z, zi, dlmo, z0, alpha):
    f = lambda x: fonction(x, dlmo)*pow(1-x/zi, alpha/2-1)*np.cos(np.sqrt(alpha*(alpha+2))/2*np.log(1-x/zi))/x

    return integrate.quad(f, z0, z)[0]

def psih_nieuwstadt(fonction, z, zi, dlmo, z0, alpha):
    f = lambda x: fonction(x, dlmo)*pow(1 - x/zi, alpha-2)/x
    return integrate.quad(f, z0, z)[0]

def tke_nieuwstadt(fonctionh, fonctionm, z, zi, dlmo, alpha):
    d = 9.7
    c = 6
    c2 = 0.1
    c3 = -0.8
    ctheta = 1.4
    w0 = 1/3 - 2/(3*c) + 2/c + 4*c3/(3*c)
    w1 = 2*c2/(3*c) - 2/c - 4*c3/(3*c)
    zeta = dlmo*z
    rich = fonctionh(z, dlmo)/pow(fonctionm(z, dlmo), 2)*zeta
    a1 = 0.5 + 1.5*pow(rich, 2) - pow(rich, 3)
    tke = 0.5*np.sqrt(np.abs(d*(fonctionm(z, dlmo) - zeta)/(fonctionh(z, dlmo)*(w0+w1*fonctionm(z, dlmo)/(fonctionm(z, dlmo)-zeta)-(1-a1)*zeta/(ctheta*(fonctionm(z, dlmo)-zeta))))))*pow(1-z/zi, alpha/2+1)
    tke = np.nan_to_num(tke)
    return tke

def epsilon_nieuwstadt(fonction, z, zi, dlmo, alpha):
    zeta = dlmo*z
    epsilon = (fonction(z, dlmo)-zeta)*pow(1-z/zi, alpha)/z
    return epsilon

#Nieuwstadt profiles
def Nieuwstadt_profile(altitude=np.linspace(0.0001,500,1000),Lmo=0.0, z0=0.0001, T0=293.15, dh=100.,
                       dtheta=5.0, gamma=0.004, G=10.0, lat=55., hub_height=150.0, hub_dir=270.0):

    alphai = 1 #alpha parameter in Nieuwstadt formulation

    kappa = 0.41
    dlmo = 1/Lmo
    zi = 1000
    omegat = 7.292115e-5
    f = 2*omegat*np.sin(2*np.pi*lat/360)
    eps = 10
    count = 0

    uadim = np.empty(len(altitude))
    vadim = np.empty(len(altitude))
    theta_adim = np.zeros(len(altitude))
    tke_adim = np.zeros(len(altitude))
    epsilon_adim = np.zeros(len(altitude))
    while eps > 1e-4 and count < 10:
        for i, z in enumerate(altitude):
            if z < zi:
                uadim[i] = psim_uplus_nieuwstadt(phim_hogstrom_s, z, zi, dlmo, z0, alphai)/kappa
                vadim[i] = psim_vplus_nieuwstadt(phim_hogstrom_s, z, zi, dlmo, z0, alphai)/kappa
                theta_adim[i] = psih_nieuwstadt(phih_hogstrom_s, z, zi, dlmo, z0, alphai)/kappa
                tke_adim[i] = tke_nieuwstadt(phih_hogstrom_s, phim_hogstrom_s, z, zi, dlmo, alphai)
                epsilon_adim[i] = epsilon_nieuwstadt(phim_hogstrom_s, z, zi, dlmo, alphai)/kappa
            else:
                uadim[i] = uadim[i-1]
                vadim[i] = vadim[i-1]
                tke_adim[i] = 0
                epsilon_adim[i] = 0

        idx_zalpha = np.argmin(np.abs(altitude - hub_height))
        rotation_angle = np.radians(270.0-hub_dir) - np.arctan2(vadim[idx_zalpha], uadim[idx_zalpha])
        rotation_matrix = np.array([[np.cos(rotation_angle), -np.sin(rotation_angle)],
                            [np.sin(rotation_angle), np.cos(rotation_angle)]])
        uadim_rot, vadim_rot = np.dot(rotation_matrix, np.array([uadim, vadim]))

        ustar = G/uadim_rot[-1]
        thetastar = ustar**2/(kappa*9.81/T0*Lmo)
        uprof = ustar*uadim_rot
        vprof = ustar*vadim_rot
        thetaprof = thetastar*theta_adim + T0
        tke = pow(ustar, 2)*tke_adim
        epsilon = pow(ustar, 3)*epsilon_adim
        zi_temp = zi
        zi = 0.4*np.sqrt(ustar*Lmo/f)
        eps = np.abs((zi-zi_temp)/zi)
        #print("eps = "+str(eps))
        count += 1

    iz = np.max(np.where(altitude <= zi))
    # Capping inversion
    ksi = 1.5
    c = 1.0/(2*ksi)
    l = zi + dh/2
    eta = (altitude-l)/(c*dh)
    f = (np.tanh(eta)+1.0)/2.0
    with np.errstate(over='ignore'):
        g = (np.log(2*np.cosh(eta))+eta)/2.0
    for j in np.where(np.isinf(g)):
        g[j] = (np.abs(eta[j])+eta[j])/2.0
    b = gamma*dh/3
    a = dtheta + b
    th = (thetaprof[iz-1] + a*f + b*g)
    thetaprof[iz-1:] = th[iz-1:]

    return thetaprof, uprof, vprof,  tke, epsilon
