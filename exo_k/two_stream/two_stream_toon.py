# -*- coding: utf-8 -*-
"""
Created in Jan 2021

@author: jeremy leconte
"""
import numpy as np
from scipy.linalg import solve_banded
import numba

@numba.jit(nopython=True, fastmath=True, cache=True)
def solve_2stream_nu_xsec(source_nu, dtau_nu, omega0_nu, g_asym_nu,
                flux_top_dw_nu, mu0=0.5, alb_surf=0., flux_at_level=False):
    """Deals with the spectral axis
    """
    NLEV, NW = source_nu.shape
    flux_up=np.empty((NLEV, NW))
    flux_dw=np.empty((NLEV, NW))
    flux_net=np.empty((NLEV, NW))
    for iW in range(NW):
        flux_up[:,iW], flux_dw[:,iW], flux_net[:,iW] = \
            solve_2stream(source_nu[:,iW], dtau_nu[:,iW],
                omega0_nu[:,iW], g_asym_nu[:,iW],
                mu0=mu0, flux_top_dw=flux_top_dw_nu[iW],
                alb_surf=alb_surf, flux_at_level=flux_at_level)
    return flux_up, flux_dw, flux_net

@numba.jit(nopython=True, fastmath=True, cache=True)
def solve_2stream_nu_corrk(source_nu, dtau_nu, omega0_nu, g_asym_nu,
                flux_top_dw_nu, mu0=0.5, alb_surf=0., flux_at_level=False):
    """Deals with the spectral axis
    """
    NLEV, NW = source_nu.shape
    NG = dtau_nu.shape[-1]
    flux_up=np.empty((NLEV, NW, NG))
    flux_dw=np.empty((NLEV, NW, NG))
    flux_net=np.empty((NLEV, NW, NG))
    for iW in range(NW):
        for iG in range(NG):
            flux_up[:,iW, iG], flux_dw[:,iW, iG], flux_net[:,iW, iG] = \
                solve_2stream(source_nu[:,iW], dtau_nu[:,iW,iG], 
                    omega0_nu[:,iW,iG], g_asym_nu[:,iW,iG],
                    mu0=mu0, flux_top_dw=flux_top_dw_nu[iW],
                    alb_surf=alb_surf, flux_at_level=flux_at_level)
    return flux_up, flux_dw, flux_net

@numba.jit(nopython=True, fastmath=True, cache=True)
def solve_2stream(source, dtau, omega0, g_asym,
                mu0=0.5, flux_top_dw=0., alb_surf=0., flux_at_level=False):
    """After Toon et al. (JGR, 1989). Equation number refer to this article.
    
    emis_surf=1.-alb_surf
    
    As we only consider hemispheric mean or quadrature, mu1==mu0

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
        flux_at_level: bool, optional
            if flux_at_level is True, fluxes are calculated at the level surfaces.
            If False, fluxes are computed at the middle of the layers.
            The top of atmosphere flux is always computed at the top of the uppermost layer
            (1st level).


    """
    Nlay=dtau.size
    gam_1, gam_2=_gammas_toon(omega0, g_asym, mu0=mu0)
    flux_net, J4pimu=matrix_toon_tridiag(Nlay, source, dtau, gam_1, gam_2, mu0,
                flux_top_dw, alb_surf, flux_at_level)
    return 0.5*(J4pimu+flux_net), 0.5*(J4pimu-flux_net), flux_net

@numba.jit(nopython=True, fastmath=True, cache=True)
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


@numba.jit(nopython=True, fastmath=True, cache=True)
def matrix_toon_tridiag(Nlay, source, dtau, gam_1, gam_2, mu1,
        flux_top_dw, alb_surf, flux_at_level):
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
    A[0]=0.
    D[0]=-e_2[0]
    E[0]=flux_top_dw-c_dw_top[0]
    B[0]=e_1[0]
    #even in Toon numerotation. (remember python arrays start at 0.)
    A[1:-1:2]=e_1[:-1]*e_2[1:]-e_3[:-1]*e_4[1:]
    B[1:-1:2]=e_2[:-1]*e_2[1:]-e_4[:-1]*e_4[1:]
    D[1:-1:2]=e_1[1:]*e_4[1:]-e_2[1:]*e_3[1:]
    E[1:-1:2]=(c_up_top[1:]-c_up_bot[:-1])*e_2[1:]-(c_dw_top[1:]-c_dw_bot[:-1])*e_4[1:]
    # middle sign above different in my calculations and toon (+ in Toon)
    #odd
    A[2:-1:2]=e_2[:-1]*e_3[:-1]-e_4[:-1]*e_1[:-1]
    B[2:-1:2]=e_1[:-1]*e_1[1:]-e_3[:-1]*e_3[1:]
    D[2:-1:2]=e_3[:-1]*e_4[1:]-e_1[:-1]*e_2[1:]
    E[2:-1:2]=(c_up_top[1:]-c_up_bot[:-1])*e_3[:-1]-(c_dw_top[1:]-c_dw_bot[:-1])*e_1[:-1]
    #surface
    A[-1]=e_1[-1]-alb_surf*e_3[-1]
    B[-1]=e_2[-1]-alb_surf*e_4[-1]
    D[-1]=0.
    E[-1]=(1.-alb_surf)*source[-1]-c_up_bot[-1]+alb_surf*c_dw_bot[-1]
    #return mat,E
    Y=DTRIDGL(2*Nlay,A,B,D,E)
    flux_net = np.empty((Nlay+1))
    J4pimu   = np.empty((Nlay+1))
    if flux_at_level:
        c_up_mid, c_dw_mid = c_planck_mid(source, dtau, gam_1, gam_2, mu1)
        #print('mid:',np.amax(np.abs((c_up_mid-c_up_top)/c_up_top)))
        #print('bot:',np.amax(np.abs((c_up_bot-c_up_top)/c_up_top)))
        mid_factor1, mid_factor2 = mid_factor_toon(dtau, gam_1, gam_2)
        flux_net[1:] = Y[1::2]*mid_factor2+c_up_mid-c_dw_mid
        J4pimu[1:]   = Y[::2]*mid_factor1+c_up_mid+c_dw_mid
        #flux_net[-1] = Y[-2]*(e_1[-1]-e_3[-1])+Y[-1]*(e_2[-1]-e_4[-1])+c_up_bot[-1]-c_dw_bot[-1]
        #J4pimu[-1]   = Y[-2]*(e_1[-1]+e_3[-1])+Y[-1]*(e_2[-1]+e_4[-1])+c_up_bot[-1]+c_dw_bot[-1]
    else:
        flux_net[1:] = Y[::2]*(e_1-e_3)+Y[1::2]*(e_2-e_4)+c_up_bot-c_dw_bot
        J4pimu[1:]   = Y[::2]*(e_1+e_3)+Y[1::2]*(e_2+e_4)+c_up_bot+c_dw_bot
    flux_net[0]  = Y[0]*e_3[0]-Y[1]*e_4[0]+c_up_top[0]
    J4pimu[0]    = flux_net[0] + flux_top_dw
    flux_net[0] -= flux_top_dw
    return flux_net, J4pimu

@numba.jit(nopython=True, fastmath=True, cache=True)
def lambda_toon(gam_1, gam_2):
    """lambda from eq 21 of Toon et al.
    """
    return np.sqrt(gam_1*gam_1-gam_2*gam_2)

@numba.jit(nopython=True, fastmath=True, cache=True)
def lambda_GAMMA(gam_1, gam_2):
    """lambda and GAMMA from eq 21 and 22 of Toon et al.
    For GAMMA, the positive alterative is used (center in eq 22)
    """
    lamb=lambda_toon(gam_1, gam_2)
    GAMMA=gam_2/(lamb+gam_1)
    #GAMMA=(gam_1-lamb)/(gam_2)
    return lamb,GAMMA

@numba.jit(nopython=True, fastmath=True, cache=True)
def c_planck(source, dtau, gam_1, gam_2, mu1):
    """c_up/dw is for c+/- whithout direct beam scattering. 
    _top is for tau equal 0 (top of the layer)
    _bot is for tau=dtau (bottom of the layer)
    removed a pi factor because source is pi*B
    """
    #cst=2*mu1
    #cst=1. # this factor seems to avoid emissivity greater than one. 
           # but what is the analytical basis for this ???
    #print(gam_1+gam_2)
    inv_dtaugam=1./(dtau*(gam_1+gam_2))
    #print(1./(dtau*(gam_1+gam_2)), 1./dtau*(gam_1+gam_2))
    c_up_top=( source[:-1]*(1.-inv_dtaugam)+source[1:]*    inv_dtaugam )#*cst
    c_dw_top=( source[:-1]*(1.+inv_dtaugam)-source[1:]*    inv_dtaugam )#*cst
    c_up_bot=(-source[:-1]*    inv_dtaugam +source[1:]*(1.+inv_dtaugam))#*cst
    c_dw_bot=( source[:-1]*    inv_dtaugam +source[1:]*(1.-inv_dtaugam))#*cst
    #print(c_up_top, c_dw_top, c_up_bot, c_dw_bot)
    return c_up_top, c_dw_top, c_up_bot, c_dw_bot

@numba.jit(nopython=True, fastmath=True, cache=True)
def c_planck_mid(source, dtau, gam_1, gam_2, mu1):
    """c_up/dw is for c+/- whithout direct beam scattering.
    These are computed at the middle of the layer, i.e. at tau=dtau/2.
    removed a pi factor because source is pi*B
    """
    cst=1. # this factor seems to avoid emissivity greater than one. 
           # but what is the analytical basis for this ???
    inv_dtaugam=1./(dtau*(gam_1+gam_2))
    c_up_mid=cst*( source[:-1]*(.5-inv_dtaugam)+source[1:]*(0.5+inv_dtaugam))
    c_dw_mid=cst*( source[:-1]*(.5+inv_dtaugam)+source[1:]*(0.5-inv_dtaugam))
    return c_up_mid, c_dw_mid

@numba.jit(nopython=True, fastmath=True, cache=True)
def e_i_toon(dtau, gam_1, gam_2):
    """e_i factors defined in eq 44.
    """
    lamb,GAMMA=lambda_GAMMA(gam_1, gam_2)
    expdtau=np.exp(-lamb*dtau)
    #print('lamb,GAMMA,expdtau',lamb,GAMMA,expdtau)
    e_1=1.+GAMMA*expdtau
    e_2=1.-GAMMA*expdtau
    e_3=GAMMA+expdtau
    e_4=GAMMA-expdtau
    return e_1, e_2, e_3, e_4

@numba.jit(nopython=True, fastmath=True, cache=True)
def mid_factor_toon(dtau, gam_1, gam_2):
    """Factors to recover the flux at mid layer. 
    """
    lamb,GAMMA=lambda_GAMMA(gam_1, gam_2)
    expdtaumid=np.exp(-lamb*dtau*0.5)
    mid_fac1=2.*(1.+GAMMA)*expdtaumid
    mid_fac2=2.*(1.-GAMMA)*expdtaumid
    return mid_fac1, mid_fac2

@numba.jit(nopython=True, fastmath=True, cache=True)
def DTRIDGL(L,AF,BF,CF,DF):
    """
          DIMENSION AF(L),BF(L),CF(L),DF(L),XK(L)
          DIMENSION AS(2*L),DS(2*L)

    !*    THIS SUBROUTINE SOLVES A SYSTEM OF TRIDIAGIONAL MATRIX
    !*    EQUATIONS. THE FORM OF THE EQUATIONS ARE:
    !*    A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = D(I)
    !======================================================================!
    """
    AS=np.empty_like(AF)
    DS=np.empty_like(AF)
    XK=np.empty_like(AF)
    AS[-1] = AF[-1]/BF[-1]
    DS[-1] = DF[-1]/BF[-1]

    for I in range(1,L):
        X         = 1./(BF[L+1-I-2] - CF[L+1-I-2]*AS[L+2-I-2])
        AS[L+1-I-2] = AF[L+1-I-2]*X
        DS[L+1-I-2] = (DF[L+1-I-2]-CF[L+1-I-2]*DS[L+2-I-2])*X
 
    XK[0]=DS[0]
    for I in range(1,L):
        XKB   = XK[I-1]
        XK[I] = DS[I]-AS[I]*XKB
    return XK
    
#@numba.jit(nopython=True,fastmath=True, cache=True)
def matrix_toon(Nlay, source, dtau, gam_1, gam_2, mu1, flux_top_dw, alb_surf):
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
    # A[0]=0. #no need because of the way diagonals are treated (see solve_banded)
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
    return flux_net, J4pimu
