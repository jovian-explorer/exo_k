
# -*- coding: utf-8 -*-
"""
@author: jeremy leconte
"""
import numpy as np
import numba
import exo_k.util.cst as cst


@numba.jit(nopython=True, fastmath=True, cache=True)
def dry_convective_adjustment(timestep, Nlay, t_lay, exner, dmass,
    tracer_array, Mgas, verbose = False):
    r"""Computes the heating rates needed to adjust unstable regions 
    of a given atmosphere to a convectively neutral T profile on
    a given timestep.
    
    Parameters
    ----------
        timestep: float
            Duration of the adjustment.
            If given in seconds, H=dT/timestep is in K/s.
            In the current model, timestep is in second over cp.
            This ensures that H=dT/timestep is in W/kg.
        Nlay: int
            Number of atmospheric layers
        t_lay: array
            Temperatures of the atmospheric layers (K)
        exner: array
            Exner function computed at the layer centers ((p/psurf)**rcp)

            .. math::
              \Pi=(p / p_{s})^{R/c_p}

        dmass: array
            Mass of gas in each layer (kg/m^2)

    Returns
    -------
        array
            Heating rate in each atmospheric layer (W/kg).  
    """
    theta_lev=t_lay/exner
    new_theta_lev = np.copy(theta_lev)
    new_Mgas = np.copy(Mgas)
    new_tracers = tracer_array.copy()
    exner_dmass=dmass*exner
    H_conv=np.zeros(Nlay)
    n_iter=0
    if verbose: print('enter convection')
#    if verbose: print(new_theta_lev, Mgas, theta_ov_mu)
    while True:
        theta_ov_mu = new_theta_lev/new_Mgas
        #conv=np.nonzero(new_theta_lev[:-1]-new_theta_lev[1:]<epsilon)[0]
        conv=np.nonzero(theta_ov_mu[:-1]<theta_ov_mu[1:])[0]
        # find convective layers
        if verbose: print(conv)
        #if verbose:
        #    print('start at, end at:',conv[0]-4,conv[-1]+3)
        #    print(theta_ov_mu[conv[0]-4:conv[-1]+4])
        N_conv=conv.size
        if N_conv==0: # no more convective regions, normal exit
            if verbose: print(conv)
            return H_conv, new_tracers
        i_conv=0
        i_top=conv[i_conv] #upper unstable layer
        while i_conv<N_conv-1: #search from the top of the 1st unstable layer for its bottom
            if conv[i_conv+1]==conv[i_conv]+1:
                i_conv+=1
                continue
            else:
                break
        i_bot=conv[i_conv]+1
        if verbose: print('it,ib,:',i_top,i_bot,i_conv)
        mass_conv=0.
        intexner=0.
        theta_mean=0.
        Mmean=0.
        for ii in range(i_top,i_bot+1): # compute new mean potential temperature
            intexner += exner_dmass[ii]
            mass_conv += dmass[ii]
            theta_mean += exner_dmass[ii] * (new_theta_lev[ii] - theta_mean) / intexner
            Mmean += dmass[ii] * (new_Mgas[ii] - Mmean) / mass_conv
        #if verbose: print('theta_mean,Mmean,:',theta_mean,Mmean,theta_mean/Mmean)
        i_top_last, ibot_last = -1, -1
        while (i_top != i_top_last) or (i_bot != ibot_last):
            i_top_last, ibot_last = i_top, i_bot
            while i_top>0:#look for newly unstable layers above
                if theta_ov_mu[i_top-1]<theta_mean/Mmean:
                    i_top -= 1
                    intexner += exner_dmass[i_top]
                    mass_conv += dmass[i_top]
                    theta_mean += exner_dmass[i_top] * (new_theta_lev[i_top] - theta_mean) / intexner
                    Mmean += dmass[i_top] * (new_Mgas[i_top] - Mmean) / mass_conv
                else:
                    break
            while i_bot<Nlay-1: #look for newly unstable layers below
                if theta_ov_mu[i_bot+1]>theta_mean/Mmean:
                    i_bot+=1
                    intexner+=exner_dmass[i_bot]
                    mass_conv+=dmass[i_bot]
                    theta_mean+=exner_dmass[i_bot] * (new_theta_lev[i_bot] - theta_mean) / intexner
                    Mmean += dmass[i_bot] * (new_Mgas[i_bot] - Mmean) / mass_conv
                else:
                    break
        #if verbose: print('it,ib, mconv1,2, M, th:',
        #    i_top,i_bot, mass_conv, np.sum(dmass[i_top:i_bot+1]), Mmean, theta_mean, theta_mean/Mmean)
        # compute heating and adjust before looking for a new potential unstable layer
        H_conv[i_top:i_bot+1] += (theta_mean-new_theta_lev[i_top:i_bot+1]) \
            *exner[i_top:i_bot+1]/timestep
        new_theta_lev[i_top:i_bot+1] = theta_mean
        new_Mgas[i_top:i_bot+1] = Mmean
        # mix tracers
        for q in new_tracers:
            q[i_top:i_bot+1]=np.sum(q[i_top:i_bot+1]*dmass[i_top:i_bot+1])/mass_conv
        n_iter+=1
        if n_iter>Nlay+1:
            if verbose : print('oops, went crazy in convadj')
            break
    return H_conv, new_tracers # we exit through here only when we exceed the max number of iteration

@numba.jit(nopython=True, fastmath=True, cache=True)
def turbulent_diffusion(timestep, Nlay, p_lay, p_lev, dmass,
        t_lay_ov_mu, g, Kzz, tracer_array, verbose = False):
    r"""Solves turbulent diffusion equation:

    .. math::
      \rho frac{\partial q}{\partial t} = \frac{\partial F_{diff}}{\partial z}
    
    with a diffusive flux given by 

    .. math::
      F_{diff} = - \rho K_{zz} \frac{\partial q}{\partial z}

    The equation is solved with an implicit scheme assuming that
    there is no flux at the top and bottom boundaries
    (evaporation must be treated separately for now).

    Parameters
    ----------
        timestep: float
            Time step in seconds.
        Nlay: int
            Number of atmospheric layers
        t_lay_ov_mu: array
            Temperatures of the atmospheric layers divided by the molar_mass in kg/mol
        p_lay: array
            Pressure at the layer centers (Pa)
        p_lev: array
            Presure at the Nlay+1 level boundaries (Pa)
        dmass: array
            Mass of gas in each layer (kg/m^2)
        g: float
            Gravity (m/s^2)
        Kzz: float
            Eddy mixing coefficient (m^2/s)
        tracer_array: array (Ntrac, Nlay)
            Array containing the mass mixing ratio of all tracers at each layer
            before the mixing

    Returns
    -------
        new_tracers: array (Ntrac, Nlay)
            Array containing the mass mixing ratio of all tracers at each layer
            after the mixing
    """
    mid_density = p_lev[1:-1]*2./(cst.RGP*(t_lay_ov_mu[1:]+t_lay_ov_mu[:-1]))
    mid_factor = - g * g * timestep * mid_density**2 / np.diff(p_lay) * 0.5*(Kzz[1:]+Kzz[:-1])
    if verbose:
        print(mid_factor)
        print(dmass)
    A = np.zeros(Nlay)
    B = np.copy(dmass)
    C = np.zeros(Nlay)
    A[1:] = mid_factor
    C[:-1] = mid_factor
    B += - C - A
    D = np.zeros(Nlay)
    new_tracers = tracer_array.copy()
    for i_q, q in enumerate(new_tracers):
        D = dmass * q
        new_q = DTRIDGL(Nlay,A,B,C,D)
        new_tracers[i_q] = new_q
        #mix_rate[name] = (new_q-q)/timestep
    return new_tracers

@numba.jit(nopython=True, fastmath=True, cache=True)
def molecular_diffusion(timestep, Nlay, p_lay, p_lev, dmass,
        tlay, mu, g, Dmol, verbose = False):
    r"""Solves turbulent diffusion equation:

    .. math::
      \rho frac{\partialT}{\partial t} = \frac{\partial F_{diff}}{\partial z}
    
    with a diffusive flux given by 

    .. math::
      F_{diff} = - \rho D_{mol} \frac{\partial T}{\partial z}

    The equation is solved with an implicit scheme assuming that
    there is no flux at the top and bottom boundaries
    (evaporation must be treated separately for now).

    Parameters
    ----------
        timestep: float
            Time step in seconds.
        Nlay: int
            Number of atmospheric layers
        t_lay_ov_mu: array
            Temperatures of the atmospheric layers divided by the molar_mass in kg/mol
        p_lay: array
            Pressure at the layer centers (Pa)
        p_lev: array
            Presure at the Nlay+1 level boundaries (Pa)
        dmass: array
            Mass of gas in each layer (kg/m^2)
        g: float
            Gravity (m/s^2)
        Dmol: float
            molecular diffusion coefficient (m^2/s)

    Returns
    -------
        new_tlay: array (Nlay)
            Array containing the temperature at each layer
            after the mixing
    """
    mid_density = p_lev[1:-1]*2.*(mu[1:]+mu[:-1])/(cst.RGP*(tlay[1:]+tlay[:-1]))
    mid_factor = - g * g * timestep * mid_density**2 / np.diff(p_lay) * Dmol
    if verbose:
        print(mid_factor)
        print(dmass)
    A = np.zeros(Nlay)
    B = np.copy(dmass)
    C = np.zeros(Nlay)
    A[1:] = mid_factor
    C[:-1] = mid_factor
    B += - C - A
    D = dmass * tlay
    new_tlay = DTRIDGL(Nlay,A,B,C,D)
    H_diff = (new_tlay-tlay)/timestep
    return H_diff

@numba.jit(nopython=True, fastmath=True, cache=True)
def DTRIDGL(L,AF,BF,CF,DF):
    """
    !  GCM2.0  Feb 2003

    !     DOUBLE PRECISION VERSION OF TRIDGL

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
    