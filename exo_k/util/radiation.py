# -*- coding: utf-8 -*-
"""
Created on Jun 20 2019

@author: jeremy leconte

Library of useful functions for radiative transfer calculations
"""
import numpy as np
import numba
import scipy.integrate as integrate
from exo_k.util.cst import PI,PLANCK,C_LUM,KBOLTZ,SIG_SB

PLANCK_CST1MKS=2.*PLANCK*C_LUM**2
PLANCK_CST1=1.e2*PLANCK_CST1MKS
PLANCK_CST2=PLANCK*C_LUM/(KBOLTZ)
PLANCK_CST1_lamb=1.e-6*2.*PLANCK*C_LUM**2




@numba.jit(nopython=True,fastmath=True)
#@numba.vectorize([float64(float64,float64)],fastmath=True)
def Bnu(nu, T):
    """Computes the Planck law in wavenumber domain.

    Parameters
    ----------
        nu: float
            The wavenumber in cm^-1. 
        T: float
            Temperature (K).

    Returns
    -------
        float
            Bnu in W/m^2/str/cm^-1 (not SI units).

    sigma*T^4=Int(Bnu dnu) where nu is expressed in cm^-1
    """
    sigma=nu*1.e2
    return PLANCK_CST1*sigma**3/(np.exp(PLANCK_CST2*sigma/T)-1.)


def Bnu_integral(nu_edges, T):
    """Computes the integral of the Planck function in wavenumber bins.
    
    Parameters
    ----------
        nu_edges: array
            Edges of the wavenumber bins in cm^-1.
        T: float
            Temperature (K).

    Returns
    -------
        array
            Bnu_integral in W/m^2/str.

    sigma*T^4/PI=sum of Bnu_integral
    """
    nu_edges=np.array(nu_edges)
    res=np.empty(nu_edges.size-1)
    for ii in range(nu_edges.size-1):
        res[ii]=integrate.quad(lambda nu: Bnu(nu,T),nu_edges[ii],nu_edges[ii+1])[0]
    return res

@numba.jit(nopython=True, fastmath=True)
def Bnu_integral_num(nu_edges, T, n=30):
    """Computes the integral of the Planck function in wavenumber bins.
    
    Parameters
    ----------
        nu_edges: array
            Edges of the wavenumber bins in cm^-1.
        T: float
            Temperature (K).
        n: int
            number of terms to take into account in the black body calculation.
    
    Returns
    -------
        Array of shape (T_array.size,nu_edges.size-1)
            Bnu_integral_num: Integral of the source function at 
            temperatures T_array inside the bins in W/m^2/str
    
    sigma*T^4=sum of Bnu_integral_num
    """
    kp=PLANCK_CST2/T
    res=np.zeros(nu_edges.size)
    Nn=nu_edges.size
    for ii in range(Nn):
        kpnu=kp*nu_edges[ii]*1.e2
        for iteration in range(1,n+1):
            kpnun=kpnu*iteration
            #res[ii]+=np.exp(-kpnun)*(6+6*kpnun+3*kpnun**2+kpnun**3)/iteration**4
            res[ii]+=np.exp(-kpnun)*(6.+6.*kpnun+3.*kpnun**2+kpnun**3)/iteration**4
    for ii in range(Nn-1):
        res[ii]-=res[ii+1]

    return res[:-1]*PLANCK_CST1MKS/kp**4

@numba.jit(nopython=True, fastmath=True)
def Bnu_integral_array(nu_edges, T_array, Nw, Nt, n=30):
    """Computes the integral of the Planck function in wavenumber bins.
    
    Parameters
    ----------
        nu_edges: array
            Edges of the wavenumber bins in cm^-1.
        T: array
            Array of temperature (K).
        Nw: int
            Number of spectral bins.
        Nt: int
            Number of temperatures.
        n: int, optional
            number of terms to take into account in the black body calculation.
    
    Returns
    -------
        Array of shape (T_array.size,nu_edges.size-1)
            Bnu_integral_num: Integral of the source function at
            temperatures T_array inside the bins in W/m^2/str
    
    sigma*T^4=sum of Bnu_integral_num
    """
    res=np.empty((Nt,Nw))
    for iT in range(Nt):
        kp=PLANCK_CST2/T_array[iT]
        mult=PLANCK_CST1MKS/kp**4
        kpnu=kp*nu_edges[0]*1.e2
        tmp=0.
        for iteration in range(1,n+1):
            kpnun=kpnu*iteration
#            tmp+=np.exp(-kpnun)*(6+6*kpnun+3*kpnun**2+kpnun**3)/iteration**4
            tmp+=np.exp(-kpnun)*(6.+6.*kpnun+3.*kpnun**2+kpnun**3)/iteration**4
        for ii in range(1,Nw+1):
            tmp2=0.
            kpnu=kp*nu_edges[ii]*1.e2
            for iteration in range(1,n+1):
                kpnun=kpnu*iteration
#                tmp2+=np.exp(-kpnun)*(6+6*kpnun+3*kpnun**2+kpnun**3)/iteration**4
                tmp2+=np.exp(-kpnun)*(6.+6.*kpnun+3.*kpnun**2+kpnun**3)/iteration**4
            res[iT,ii-1]=(tmp-tmp2)*mult
            tmp=tmp2
    return res

def Bnu_integral_old(nu_edges, T):
    """Deprecated.
    
    Computes the integral of the Planck function in wavenumber bins.
    nu_edges is an array of the edges of the wavenumber bins in cm^-1.
    
    sigma*T^4=sum of Bnu_integral

    Bnu_integral is in W/m^2/str
    """
    a1=1.e8*2.*PLANCK*C_LUM**2
    a2=1.e2*PLANCK*C_LUM/(KBOLTZ*T)
    res=np.empty(nu_edges.size-1)
    for ii in range(nu_edges.size-1):
        res[ii]=integrate.quad(lambda nu: nu**3/(np.exp(nu*a2)-1.),nu_edges[ii],nu_edges[ii+1])[0]
    return res*a1


@numba.jit(nopython=True)
def Bmicron(lamb, T):
    """Computes the Planck law in wavelength domain.

    Parameters
    ----------
        lamb: float
            The wavelength in microns. 
        T: float
            Temperature (K).
 
    Returns:
        float
            Bmicron in W/m^2/str/micron (not SI units).

    sigma*T^4=Int(Bmicron d lamb) where lamb is expressed in microns
    """
    lambda_si=lamb*1.e-6
    return PLANCK_CST1_lamb/(lambda_si**5   \
        *(np.exp(PLANCK_CST2/(lambda_si*T))-1.))


@numba.jit(nopython=True)
def BBnu_Stellar_Spectrum(nu, T, Flux):
    """Computes the outgoing spectral flux of a black body at temperature T so that its 
    bolometric flux is equal to the input 'Flux'

    Parameters
    ----------
        nu: float
            The wavenumber in cm^-1. 
        T: float
            Temperature (K).
        Flux: float
            bolometric flux for rescaling.
    
    Returns
    -------
        float
            BBnu_Stellar_Spectrum is in W/m^2/cm^-1 (not SI units).
    
    sigma*T^4=Int(Bnu dnu) where nu is expressed in cm^-1
    """
    Scaling_Factor=PI/(SIG_SB*T**4)*Flux
    return Scaling_Factor*Bnu(nu,T)
    
@numba.jit(nopython=True)
def BBmicron_Stellar_Spectrum(lamb, T, Flux):
    """Computes the outgoing spectral flux of a black body at temperature T so that its 
    bolometric flux is equal to the input 'Flux'

    Parameters
    ----------
        lamb: float
            The wavelength in microns. 
        T: float
            Temperature (K).
        Flux: float
            bolometric flux for rescaling.
    
    Returns
    -------
        float
            BBmicron_Stellar_Spectrum is in W/m^2/micron (not SI units).
            
    sigma*T^4=Int(Bmicron d lamb) where lamb is expressed in microns.
    
    """
    Scaling_Factor=PI/(SIG_SB*T**4)*Flux
    return Scaling_Factor*Bmicron(lamb,T)

def wavelength_grid_R(lambda_min, lambda_max, R):
    """Creates a wavelength grid starting from lambda_min and ending at lambda_max
    with a resolution of R (roughly).

    Parameters
    ----------
        lambda_min/lambda_max: int or float
            Min/max wavelength (in micron).
        R: int or float
            Resolution.

    Returns
    -------
        array
            Grid of wavelength.
    """
    return np.append(np.exp(np.arange(np.log(lambda_min),np.log(lambda_max),1./R)),lambda_max)

def wavenumber_grid_R(nu_min, nu_max, R):
    """Creates a wavenumber grid starting from nu_min and ending at nu_max
    with a resolution of R (roughly).

    Parameters
    ----------
        nu_min/nu_max: int or float
            Min/max wavenumber (in cm^-1).
        R: int or float
            Resolution.

    Returns
    -------
        array
            Grid of wavenumber.
    """
    return np.append(np.exp(np.arange(np.log(nu_min),np.log(nu_max),1./R)),nu_max)

@numba.njit
def rad_prop_corrk(dcol_density, opacity_prof, mu0):
    """Computes the optical depth of each of the layers
    for the opacity given for every wavelength and g point.

    Parameters
    ----------
        dcol_density: array
            Column density per layer (molecules/m^2/layer).
        opacity_prof: array
            effective cross section per molecule (m^2/molecule).
        mu0: cosine of the equivalent zenith angle for the rays.
            e.g. Cos(60deg)=0.5 is chosen (See, for example, chap 4 of Pierrehumbert 2010)

    Returns
    -------
        tau: Array
            cumulative optical depth from the top (with zeros at the top of the first layer)
        dtau: Array
            optical depth of each layer for each wavenumber.
   """
    Nlev,Nw,Ng=opacity_prof.shape
    OvMu=1./mu0
    dtau=np.empty((Nlev,Nw,Ng))
    tau=np.zeros((Nlev+1,Nw,Ng))
    for lev in range(Nlev):
        for iW in range(Nw):
            for ig in range(Ng):
                dtau_tmp=dcol_density[lev]*opacity_prof[lev,iW,ig]*OvMu
                dtau[lev,iW,ig]=dtau_tmp
                tau[lev+1,iW,ig]=tau[lev,iW,ig]+dtau_tmp
    return tau, dtau

@numba.njit
def rad_prop_xsec(dcol_density, opacity_prof, mu0):
    """Computes the optical depth of each of the layers
    for the opacity given for every wavelength.

    Parameters
    ----------
        dcol_density: array
            Column density per layer (molecules/m^2/layer).
        opacity_prof: array
            effective cross section per molecule (m^2/molecule).
        mu0: cosine of the equivalent zenith angle for the rays.
            e.g. Cos(60deg)=0.5 is chosen (See, for example, chap 4 of Pierrehumbert 2010)

    Returns
    -------
        tau: Array
            cumulative optical depth from the top (with zeros at the top of the first layer)
        dtau: Array
            optical depth of each layer for each wavenumber.
    """
    Nlev,Nw=opacity_prof.shape
    OvMu=1./mu0
    dtau=np.empty((Nlev,Nw))
    tau=np.zeros((Nlev+1,Nw))
    for lev in range(Nlev):
        for iW in range(Nw):
            dtau_tmp=dcol_density[lev]*opacity_prof[lev,iW]*OvMu
            dtau[lev,iW]=dtau_tmp
            tau[lev+1,iW]=tau[lev,iW]+dtau_tmp
    return tau,dtau

@numba.njit(nogil=True,fastmath=True)
def path_integral_corrk(Nlay, Nw, Ng, tangent_path, density_prof, opacity_prof, weights):
    """Computes the transmittance for each layer of the atmosphere for k-coefficients.

    Parameters
    ----------
        Nlay: int
            Number of layers.
        Nw: int
            Number of wavenumber bins.
        Ng: int
            Number of gauss points.
        tangent_path: List of arrays
            Triangular array of tangent path for each layer.
        density_prof: array
            Array with the number density of molecules in the atmosphere.
        opacity_prof: array
            Effective cross section of the atmsophere for each (layer,wavenumber,g-point).
        weights: array
            Weights for the quadrature in g-space.

    Returns
    -------
        transmittance: array
            Transmittance for each layer and wavenumber bin (Exp(-tau_sigma)).
    """
    exp_min_tau=np.zeros((Nlay,Nw,Ng))
    transmittance=np.zeros((Nlay,Nw))
    for ilay in range(Nlay):
        for jlay in range(ilay+1):
            dm=tangent_path[ilay][jlay]*density_prof[jlay]
            for iW in range(Nw):
                for ig in range(Ng):
                    exp_min_tau[ilay,iW,ig]+=dm*opacity_prof[jlay,iW,ig]
    exp_min_tau=np.exp(-exp_min_tau)
    for ilay in range(Nlay):
        for iW in range(Nw):
            trans=0.
            for ig in range(Ng):
                trans+=exp_min_tau[ilay,iW,ig]*weights[ig]
            transmittance[ilay,iW]=trans
    return transmittance


@numba.njit(nogil=True,fastmath=True)
def path_integral_xsec(Nlay, Nw, tangent_path, density_prof, opacity_prof):
    """Computes the transmittance for each layer of the atmosphere from cross sections.

    Parameters
    ----------
        Nlay: int
            Number of layers.
        Nw: int
            Number of wavenumber bins.
        tangent_path: List of arrays
            Triangular array of tangent path for each layer.
        density_prof: array
            Array with the number density of molecules in the atmosphere.
        opacity_prof: array
            Effective cross section of the atmsophere for each (layer,wavenumber).

    Returns
    -------
        transmittance: array
            Transmittance for each layer and wavenumber bin (Exp(-tau_sigma)).
    """
    transmittance=np.zeros((Nlay,Nw))
    for ilay in range(Nlay):
        for jlay in range(ilay+1):
            dm=tangent_path[ilay][jlay]*density_prof[jlay]
            for iW in range(Nw):
                transmittance[ilay,iW]+=dm*opacity_prof[jlay,iW]
    transmittance=np.exp(-transmittance)
    return transmittance

