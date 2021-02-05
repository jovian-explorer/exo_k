# -*- coding: utf-8 -*-
"""
Created on Jun 20 2019

@author: jeremy leconte

Library of useful functions for radiative transfer calculations
"""
import numpy as np
from scipy.linalg import solve_banded
import numba

def solve_2stream_nu(NG, source_nu, dtau_nu, omega0_nu=0., g_asym_nu=0.,
                mu0=0.5, flux_top_dw_nu=0., alb_surf=0., verbose=False):
    """Deals with the spectral axis
    """
    NLEV=source_nu.shape[0]
    NW = dtau_nu.shape[1]
    if NG is None:
        flux_up=np.zeros((NLEV, NW))
        flux_dw=np.zeros((NLEV, NW))
        flux_net=np.zeros((NLEV, NW))
        for iW in range(NW):
            flux_up[:,iW], flux_dw[:,iW], flux_net[:,iW] = \
                solve_2stream(source_nu[:,iW], dtau_nu[:,iW],
                    omega0_nu[:,iW], g_asym_nu[:,iW],
                    mu0=mu0, flux_top_dw=flux_top_dw_nu,
                    alb_surf=alb_surf, verbose=verbose)
    else:
        flux_up=np.zeros((NLEV, NW, NG))
        flux_dw=np.zeros((NLEV, NW, NG))
        flux_net=np.zeros((NLEV, NW, NG))
        #source_nu2=np.zeros((NLEV, NW, NG), dtype=float)
        #print(flux_up.shape, flux_dw.shape, source_nu.shape, dtau_nu.shape,
        #  omega0_nu.shape, g_asym_nu.shape)
        for iW in range(NW):
            for iG in range(NG):
                #print(iG,NG)
                #source_nu2[:,iW,iG]=source_nu[:,iW]
                flux_up[:,iW, iG], flux_dw[:,iW, iG], flux_net[:,iW, iG] = \
                    solve_2stream(source_nu[:,iW], dtau_nu[:,iW,iG], 
                        omega0_nu[:,iW,iG], g_asym_nu[:,iW,iG],
                        mu0=mu0, flux_top_dw=flux_top_dw_nu,
                        alb_surf=alb_surf, verbose=verbose)
    return flux_up, flux_dw, flux_up-flux_dw

#@numba.jit(nopython=True,fastmath=True)
def solve_2stream(source, dtau, omega0=0., g_asym=0.,
                mu0=0.5, flux_top_dw=0., alb_surf=0., verbose=False):
    """After Toon et al. (1989)
    
    emis_surf=1.-alb_surf
    
    As we only consider hemispheric mean or quadrature, mu1==m0

    Parameters
    ----------
        source: array
            pi*B (Planck function) at each of the Nlay+1 level interfaces of the model
            (top to bottom). Thanks to the pi factor it can be compared to the incoming
            flux.
        dtau: array
            optical depth of each of the Nlay model layers.
        omega0: float or array
            single scattering albedo of each layer
        g_asym: float or array
            asymmetry factor (g in Toon et al.)
        mu0: float
            1./2. yields hemisperic mean approx
            1./sqrt(3) yields quadrature
        flux_top_dw: float
            Incoming difuse flux at the upper boundary
        alb_surf: float
            Surface albedo. Emissivity is assumed to be 1.-alb_surf

    """
    Nlay=dtau.size
    gam_1, gam_2=_gammas_toon(omega0, g_asym, mu0=mu0)
    res=matrix_toon(Nlay, source, dtau, gam_1, gam_2, mu0,
                flux_top_dw, alb_surf, verbose=verbose)
    flux_net = res[0]
    J4pimu = res[1]
    return 0.5*(J4pimu+flux_net), 0.5*(J4pimu-flux_net), flux_net

@numba.jit(nopython=True,fastmath=True)
def _gammas_toon(omega0, g_asym, mu0=0.5):
    """mu0=0.5 yields hemispheric mean
    
    Parameters
    ----------
        omega0: float or array
            single scattering albedo of each layer
        g_asym: float or array
            asymmetry factor (g in Toon et al.)
        mu0: float
            1./2. yields hemisperic mean approx
            1./sqrt(3) yields quadrature

    Returns
    -------
        gamma1, gamma2: floats or arrays
            gammas defined in Table 1 of Toon et al. 
    """
    return (2.-omega0*(1.+g_asym))/(2.*mu0), omega0*(1.-g_asym)/(2.*mu0)

#@numba.jit(nopython=True,fastmath=True)
def matrix_toon(Nlay, source, dtau, gam_1, gam_2, mu1, flux_top_dw, alb_surf, verbose=False):
    """
    Returns
    -------
        flux_net: array
            Net flux at the bottom of the Nlay layers.
    """
    e_1, e_2, e_3, e_4 = e_i_toon(dtau, gam_1, gam_2)
    c_up_top, c_dw_top, c_up_bot, c_dw_bot = c_planck(source, dtau, gam_1, gam_2, mu1)
    A=np.empty((2*Nlay))
    B=np.empty((2*Nlay))
    D=np.empty((2*Nlay))
    E=np.empty((2*Nlay))
    # upper boundary
    # A[0]=0. #no need because of the way diagonals are treated
    B[0]=e_1[0]
    D[0]=0.
    D[1]=-e_2[0]
    E[0]=flux_top_dw-c_dw_top[0]
    #even 
    A[:-2:2]=e_1[:-1]*e_2[1:]-e_3[:-1]*e_4[1:]
    B[1:-1:2]=e_2[:-1]*e_2[1:]-e_4[:-1]*e_4[1:]
    D[2:-1:2]=e_1[1:]*e_4[1:]-e_2[1:]*e_3[1:]
    E[1:-1:2]=(c_up_top[1:]-c_up_bot[:-1])*e_2[1:]-(c_dw_top[1:]-c_dw_bot[:-1])*e_4[1:]
    # middle sign above different in my calculations and toon (+ in Toon)
    #odd
    A[1:-2:2]=e_2[:-1]*e_3[:-1]-e_4[:-1]*e_1[:-1]
    B[2:-1:2]=e_1[:-1]*e_1[1:]-e_3[:-1]*e_3[1:]
    D[3::2]=e_3[:-1]*e_4[1:]-e_1[:-1]*e_2[1:]
    E[2:-1:2]=(c_up_top[1:]-c_up_bot[:-1])*e_3[:-1]-(c_dw_top[1:]-c_dw_bot[:-1])*e_1[:-1]
    #surface
    A[-2]=e_1[-1]-alb_surf*e_3[-1]
    B[-1]=e_2[-1]-alb_surf*e_4[-1]
    #D[-1]=0.
    E[-1]=(1.-alb_surf)*source[-1]-c_up_bot[-1]+alb_surf*c_dw_bot[-1]
    #return mat,E
    Y=solve_banded((1,1),[D,B,A],E)
    flux_net = np.empty((Nlay+1))
    J4pimu   = np.empty((Nlay+1))
    flux_net[1:] = Y[::2]*(e_1-e_3)+Y[1::2]*(e_2-e_4)+c_up_bot-c_dw_bot
    J4pimu[1:]   = Y[::2]*(e_1+e_3)+Y[1::2]*(e_2+e_4)+c_up_bot+c_dw_bot
    flux_net[0]  = Y[0]*e_3[0]-Y[1]*e_4[0]+c_up_top[0]
    J4pimu[0]    = flux_net[0] + flux_top_dw
    flux_net[0] -= flux_top_dw
    if not verbose:
        return flux_net, J4pimu
    else:
        return flux_net, J4pimu, e_1, e_2, e_3, e_4, c_up_top, \
            c_dw_top, c_up_bot, c_dw_bot, A, B, D, E, Y


@numba.jit(nopython=True,fastmath=True)
def lambda_toon(gam_1, gam_2):
    """lambda from eq 21 of Toon et al.
    """
    return np.sqrt(gam_1*gam_1-gam_2*gam_2)

@numba.jit(nopython=True,fastmath=True)
def lambda_GAMMA(gam_1, gam_2):
    """lambda and GAMMA from eq 21 and 22 of Toon et al.
    For GAMMA, the positive alterative is used (center in eq 22)
    """
    lamb=lambda_toon(gam_1, gam_2)
    GAMMA=gam_2/(lamb+gam_1)
    return lamb,GAMMA

@numba.jit(nopython=True,fastmath=True)
def c_planck(source, dtau, gam_1, gam_2, mu1):
    """c_up/dw is for c+/- whithout direct beam scattering. 
    _top is for tau equal 0 (top of the layer)
    _bot is for tau=dtau (bottom of the layer)
    """
    cst=2*mu1
    #print(gam_1+gam_2)
    inv_dtaugam=1./(dtau*(gam_1+gam_2))
    #print(1./(dtau*(gam_1+gam_2)), 1./dtau*(gam_1+gam_2))
    c_up_top=cst*( source[:-1]*(1.-inv_dtaugam)+source[1:]*    inv_dtaugam )
    c_dw_top=cst*( source[:-1]*(1.+inv_dtaugam)-source[1:]*    inv_dtaugam )
    c_up_bot=cst*(-source[:-1]*    inv_dtaugam +source[1:]*(1.+inv_dtaugam))
    c_dw_bot=cst*( source[:-1]*    inv_dtaugam +source[1:]*(1.-inv_dtaugam))
    #print(c_up_top, c_dw_top, c_up_bot, c_dw_bot)
    return c_up_top, c_dw_top, c_up_bot, c_dw_bot

@numba.jit(nopython=True,fastmath=True)
def e_i_toon(dtau, gam_1, gam_2):
    """e_i factors defined 
    """
    lamb,GAMMA=lambda_GAMMA(gam_1, gam_2)
    expdtau=np.exp(-lamb*dtau)
    #print('lamb,GAMMA,expdtau',lamb,GAMMA,expdtau)
    e_1=1.+GAMMA*expdtau
    e_2=1.-GAMMA*expdtau
    e_3=GAMMA+expdtau
    e_4=GAMMA-expdtau
    return e_1, e_2, e_3, e_4



