
# -*- coding: utf-8 -*-
"""
@author: jeremy leconte
"""
import numpy as np
import numba
import exo_k.util.cst as cst
from .condensation import moist_adiabat, compute_condensation_parameters, Tsat_P

@numba.jit(nopython=True, fastmath=True, cache=True)
def dry_convective_adjustment_numba(timestep, Nlay, t_lay, exner, dmass,
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
def convective_acceleration_numba(timestep, Nlay, H_rad, rad_layers, tau_rad, tau_rads, dmass,
    convective_acceleration_mode = 0, verbose = False):
    r"""Computes a heating rate for the whole convective region to accelerate convergence
    
    Parameters
    ----------
        timestep: float
            Duration of the adjustment.
            If given in seconds, H=dT/timestep is in K/s.
            In the current model, timestep is in second over cp.
            This ensures that H=dT/timestep is in W/kg.
        Nlay: int
            Number of atmospheric layers
        H_rad: array
            Radiative heating rate
        rad_layers: array of bool
            Elements in the array are true if layer is purely radiative
        tau_rad: array
            Baseline radiative timescale for the atmosphere. (e.g. the min of tau_rads)
            Should use the same units as timestep.
        tau_rads: array
            Radiative timescale for each layer. Should use the same units as timestep.
        dmass: array
            Mass of gas in each layer (kg/m^2)
        convective_acceleration_mode: int
            convective_acceleration_mode=0 is the default mode
            =1 emulates an old (and bad) behavior.

    Returns
    -------
        array
            Heating rate in each atmospheric layer (W/kg).  
    """
    H_acc = np.zeros(Nlay)
    n_iter=0
    # find convective layers
    conv=np.nonzero(np.invert(rad_layers))[0]
    if verbose: print('in conv acc:',conv)
    N_conv=conv.size
    i_conv=-1
    while True:
        #if verbose:
        #    print('start at, end at:',conv[0]-4,conv[-1]+3)
        #    print(theta_ov_mu[conv[0]-4:conv[-1]+4])
        i_conv += 1
        if i_conv == N_conv: # no more convective regions, normal exit
            if verbose: print('normal exit')
            return H_acc
        i_top=conv[i_conv] #upper unstable layer
        while i_conv<N_conv-1: #search from the top of the 1st unstable layer for its bottom
            if conv[i_conv+1]==conv[i_conv]+1:
                i_conv+=1
                continue
            else:
                break
        i_bot=conv[i_conv]
        if verbose: print('it,ib,:',i_top,i_bot,i_conv)
        mass_conv=0.
        H_rad_mean=0.
        for ii in range(i_top,i_bot+1): # computeaverage radiative heating
            mass_conv += dmass[ii]
            H_rad_mean += dmass[ii] * (H_rad[ii] - H_rad_mean) / mass_conv
        tau = np.amin(tau_rads[i_top:i_bot+1])
        # compute heating and adjust before looking for a new potential unstable layer
        if convective_acceleration_mode == 0:
            H_acc[i_top:i_bot+1] += H_rad_mean * tau / tau_rad
        else:
            H_acc[i_top:i_bot+1] += H_rad_mean * tau / timestep
        n_iter+=1
        if n_iter>Nlay+1:
            if verbose : print('oops, went crazy in convective_acceleration')
            break
    return H_acc # we exit through here only when we exceed the max number of iteration

@numba.jit(nopython=True, fastmath=True, cache=True)
def turbulent_diffusion_numba(timestep, Nlay, p_lay, p_lev, dmass,
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
def compute_condensation_numba(timestep, Nlay, tlay, play, cp, Mgas, qarray,
        idx_vap, idx_cond, thermo_parameters, latent_heating = True,
        condensation_timestep_reducer = None,
        verbose = False):
    r"""Computes the heating rates needed to bring sursaturated 
    regions back to saturation

    Parameters
    ----------
        timestep
        Nlay: float
            Number of layers
        tlay: array
            Layer temperatures
        cp: float
            Specific heat capacity
        Mgas: array
            Layer molar mass
        qarray: array
            Array of tracer specific concentrations
        idx_vap, idx_cond: int
            Indices for condensing species in tracer array
        thermo_parameters: array
            Array of thermodynamical paramaters for condensing species.
            See `:class:`~exo_k.atm_evolution.condensation.Condensation_Thermodynamical_Parameters`
        latent_heating: bool
            Whether latent heating should be taken into account

    Returns
    -------
        H_cond: array
            Heating rate in W/kg
    """
    itermax = 1000
    if condensation_timestep_reducer is None:
        alpha = 1.
    else:
        alpha = condensation_timestep_reducer
    H_cond = np.zeros(Nlay)
    if latent_heating:
        for ii in range(Nlay):
            i_iter = 0
            T_tmp = tlay[ii]
            qvap = qarray[idx_vap,ii]
            qcond = qarray[idx_cond,ii]
            if qvap == 0.: continue
            while i_iter < itermax:
                Lvap, qsat, dqsat_dt = compute_condensation_parameters(T_tmp, play[ii], Mgas[ii], 
                        thermo_parameters[1], thermo_parameters[2], thermo_parameters[3], thermo_parameters[4],
                        thermo_parameters[5], thermo_parameters[6], thermo_parameters[7], thermo_parameters[8],
                        thermo_parameters[9])
                dqvap = alpha * (qsat - qvap) / (1. + Lvap * dqsat_dt / cp)
                if dqvap > qcond:
                    dqvap = qcond
                qcond -= dqvap
                qvap += dqvap
                T_tmp -= Lvap * dqvap / cp
                if np.abs(dqvap/(qvap*alpha)) <= 1.e-9: break
                i_iter += 1
            H_cond[ii] = (T_tmp - tlay[ii]) / timestep
            qarray[idx_vap,ii] = qvap
            qarray[idx_cond,ii] = qcond
            if verbose: print('in cond, i, iter:', ii, i_iter+1, dqvap/qvap)
            if verbose: print('in cond, RH, T:', qvap/qsat, T_tmp, tlay[ii], dqsat_dt)

    else:
        # here we treat whole arrays at once.
        Lvap, qsat, dqsat_dt = compute_condensation_parameters(tlay, play, Mgas, 
                thermo_parameters[1], thermo_parameters[2], thermo_parameters[3], thermo_parameters[4],
                thermo_parameters[5], thermo_parameters[6], thermo_parameters[7], thermo_parameters[8],
                thermo_parameters[9])
        dqvap= np.where(qsat <= qvap, qsat-qvap, 0.)
        qarray[idx_vap] += dqvap
        qarray[idx_cond] -= dqvap
        if verbose: print('in cond, RH, T:', qarray[idx_vap]/qsat, tlay, dqsat_dt)
    return H_cond

@numba.jit(nopython=True, fastmath=True, cache=True)
def moist_convective_adjustment_numba(timestep, Nlay, tlay, play, dmass, cp, Mgas, q_array,
        i_vap, i_cond, thermo_parameters, 
        moist_inhibition = True, verbose = False):
    r"""Computes the heating rates needed to adjust unstable regions 
    of a given atmosphere to a moist adiabat.

    Parameters
    ----------
        timestep
        Nlay: float
            Number of layers
        tlay: array
            Layer temperatures
        play:array
            Pressure at layer centers
        dmass: array
            mass of layers in kg/m^2
        cp: float
            specific heat capacity at constant pressure
        q_array: array
            mass mixing ratio of tracers
        i_vap: int
            index of vapor tracer in qarray
        i_cond: int
            index of condensate tracer in qarray
        qsat: array
            Saturation mmr for each layer
        dqsat_dt: array
            d qsat / dT in each layer
        Lvap: array
            Latent heat of vaporization (can have different values
            in each layer if Lvap=f(T))
        dlnt_dlnp_moist: array
            threshold thermal gradient (d ln T / d ln P) for a moist atmosphere
            computed at the layer centers.
        q_crit: array
            Critical mass mixing ratio for the inhibition of moist convection
            (Eq. 17 of Leconte et al. 2017)

    Returns
    -------
        H_madj: array
            Heating rate in each atmospheric layer (W/kg). 
        new_q: array
            tracer mmr array after adjustment.
        new_t: array
            Temperature of layers after adjustment. 
    """
    if verbose: print('enter moist adj')
    H_madj=np.zeros(Nlay)
    new_q = q_array.copy()
    new_t = tlay.copy()
    dp = np.diff(play)

    dlnt_dlnp_moist, Lvap, psat, qsat, dqsat_dt, q_crit = \
        moist_adiabat(new_t, play, cp, Mgas, thermo_parameters[0],
            thermo_parameters[1], thermo_parameters[2], thermo_parameters[3], thermo_parameters[4],
            thermo_parameters[5], thermo_parameters[6], thermo_parameters[7], thermo_parameters[8],
            thermo_parameters[9])

    nabla_interlayer = tlay * dlnt_dlnp_moist /play
    nabla_interlayer = 0.5*(nabla_interlayer[:-1]+nabla_interlayer[1:])
    dTmoist_array = nabla_interlayer * dp
    dT_inter_lay = np.diff(tlay)
    qvap = new_q[i_vap]
    mvap_sursat_array = (qvap-qsat) * dmass
    if moist_inhibition:
        q_crit_criterion = qvap/q_crit < 1. # convection possible if True. qcri can be negative.
    else:
        q_crit_criterion = qvap<2. #should always be true
    #print('dT:', dT_inter_lay)
    #print('dTmoist:', dTmoist_array)
    #dT_unstab = np.nonzero(dT_inter_lay>dTmoist_array)[0]
    #saturated = np.nonzero(mvap_sursat_array>0.)[0]
    conv = np.nonzero((dT_inter_lay>dTmoist_array)*(mvap_sursat_array[:-1]>0.) \
            *q_crit_criterion[:-1])[0]# find convective layers
    if verbose: 
        print(conv)
        print(np.nonzero(dT_inter_lay>dTmoist_array)[0])
        print(np.nonzero(mvap_sursat_array[:-1]>0.)[0])
        print(np.nonzero(q_crit_criterion)[0])
    N_conv=conv.size
    if N_conv==0: # no more convective regions, can exit
        return H_madj, new_q, new_t
    i_top=conv[0] #upper unstable layer
    T_top = new_t[i_top]
    mvap_sursat = mvap_sursat_array[i_top]
    dqsdm = dqsat_dt[i_top]*dmass[i_top]
    int_dqsdm = dqsdm
    C = cp*dmass[i_top] + Lvap[i_top]*dqsdm
    B = C * new_t[i_top] + Lvap[i_top] * mvap_sursat_array[i_top]
    dT_moist = 0.
    int_m_cond = mvap_sursat_array[i_top] + dqsdm*(new_t[i_top] - dT_moist)
    i_bot=i_top+1
    while i_bot<Nlay: #search for the bottom of the 1st unstable layer from its top
        tmp_sursat = mvap_sursat + mvap_sursat_array[i_bot]
        tmp_dT_moist = dT_moist + dTmoist_array[i_bot-1]
        dqsdm = dqsat_dt[i_bot] * dmass[i_bot]
        tmp_int_dqsdm = int_dqsdm + dqsdm
        tmp_int_m_cond = int_m_cond + mvap_sursat_array[i_bot] + dqsdm * (new_t[i_bot] - tmp_dT_moist)
        tmp = cp *dmass[i_bot] + Lvap[i_bot]* dqsdm
        tmp_C = C + tmp
        tmp_B = B + tmp * (new_t[i_bot]-tmp_dT_moist) + Lvap[i_bot] * mvap_sursat_array[i_bot]
        tmp_new_Ttop = tmp_B / tmp_C
        tmp_m_cond = tmp_int_m_cond - tmp_int_dqsdm * tmp_new_Ttop
        if tmp_sursat>0. and tmp_dT_moist<new_t[i_bot]-T_top and q_crit_criterion[i_bot] and tmp_m_cond>0.:
            dT_moist = tmp_dT_moist
            mvap_sursat = tmp_sursat
            int_dqsdm = tmp_int_dqsdm
            int_m_cond = tmp_int_m_cond
            C = tmp_C
            B = tmp_B
            m_cond = tmp_m_cond
            i_bot += 1
            continue
        else:
            i_bot -= 1
            break
    if verbose: print('it,ib=', i_top, i_bot)
    if i_top == i_bot: # need at least 2 layers to convect, so exit
        return H_madj, new_q, new_t
    new_Ttop = B / C
    if verbose: print(new_Ttop, m_cond, dT_moist)
    dT = new_Ttop - new_t[i_top]
    qvap[i_top] = qsat[i_top] + dqsat_dt[i_top] * dT
    #if verbose: print('top i, dT, qv, qs, dqs, qf', i_top, dT, q_array[i_vap, i_top], qsat[i_top], dqsat_dt[i_top], qvap[i_top])
    new_t[i_top] = new_Ttop
    H_madj[i_top] = dT / timestep
    for ii in range(i_top+1, i_bot+1):
        dT = new_t[ii-1] + dTmoist_array[ii-1] - new_t[ii]
        #print(ii, new_t[ii-1], dTmoist_array[ii-1],  new_t[ii], new_t[ii-1] + dTmoist_array[ii-1] - new_t[ii])
        qvap[ii] = qsat[ii] + dqsat_dt[ii] * dT
        new_t[ii] += dT
        # compute heating and adjust before looking for a new potential unstable layer
        H_madj[ii] = dT / timestep
        #if verbose: print('i, dT, qv, qs, dqs, qf', ii, dT, q_array[i_vap, ii], qsat[ii], dqsat_dt[ii], qvap[ii])
    # put ice
    m_cond_2 = np.sum((q_array[i_vap, i_top:i_bot+1]-qvap[i_top:i_bot+1])*dmass[i_top:i_bot+1])
    dTmoist_array[i_top-1]=0.
    if m_cond<0.:
        print('Negative condensates in moist adj, i:', i_top, i_bot)
        print(q_array[i_vap, i_top:i_bot+1], qvap[i_top:i_bot+1], q_array[i_vap, i_top:i_bot+1]-qvap[i_top:i_bot+1])
    m_conv = np.sum(dmass[i_top:i_bot+1])
    new_q[i_cond, i_top:i_bot+1] += m_cond / m_conv
    if verbose: 
        print('m_cond, m_conv, m_cond2', m_cond, m_conv, m_cond_2)
    return H_madj, new_q, new_t

@numba.jit(nopython=True, fastmath=True, cache=True)
def compute_rainout_numba(timestep, Nlay, tlay, play, dmass, cp, Mgas, qarray,
        idx_vap, idx_cond, thermo_parameters, evap_coeff, qvap_deep,
        verbose = False):
    r"""Computes the heating rates needed to adjust unstable regions 
    of a given atmosphere to a moist adiabat.

    Parameters
    ----------
        timestep
        Nlay: float
            Number of layers
        tlay: array
            Layer temperatures
    """
    H_rain=np.zeros(Nlay)
    Lvap, qsat, dqsat_dt = compute_condensation_parameters(tlay, play, Mgas, 
            thermo_parameters[1], thermo_parameters[2], thermo_parameters[3], thermo_parameters[4],
            thermo_parameters[5], thermo_parameters[6], thermo_parameters[7], thermo_parameters[8],
            thermo_parameters[9])

    Tsat_p = Tsat_P(play, thermo_parameters[5], thermo_parameters[8], thermo_parameters[9])

    if verbose: print('in rainout, RH, T:', qarray[idx_vap]/qsat, tlay)
    mass_cond = 0.
    for i_lay in range(Nlay):
        #  if evap_coeff =1, rain vaporisation in an undersaturated layer can fill the layer up to the (estimated) saturation
        #  if 0 < evap_coeff < 1, rain vaporisation in one layer is limited to a fraction of the amount that would saturate the layer
        #  This allows not to exceed saturation, to spread rain vaporization in more and denser layers
        qvap = qarray[idx_vap,i_lay]
        if (tlay[i_lay] >= Tsat_p[i_lay]): # above boiling temperature, try to evaporate everything
            dqvap = (qsat[i_lay] - qvap) / (1. + Lvap[i_lay]*dqsat_dt[i_lay]/cp)
        else:
            dqvap = evap_coeff * (qsat[i_lay] - qvap) / (1. + Lvap[i_lay]*dqsat_dt[i_lay]/cp)
        dqvap = np.core.umath.minimum(dqvap, 1.- qvap)
        mass_cond += qarray[idx_cond,i_lay]*dmass[i_lay]
        qarray[idx_cond,i_lay] = 0.
        # dqvap is the change in vapor content to reach saturation and accounting for the temperature change.
        # dqvap < 0 implies condensation, meaning that there is a remaining excess of vapor after the previous 
        # condensation step. In such case we apply this new change in vapor content and temperature and
        # increase the amount of falling condensed species (mass_cond).
        # dqvap > 0 implies evaporation. Here there are two possibilities:
        #   - the amount of condensates is lower than dqvap. All condensates are vaporised in the layer
        #    and the tempertaure change is -Ldqv/cp where dqv is the actual change in vapor.
        #   - the amount of condensates is larger than dqvap. We then apply dqvap and the corresponding 
        #    change in temperature and transfer the remaining condensate in the falling rain reservoir.
        if dqvap < 0: # more condensation in the layer
            qarray[idx_vap][i_lay] += dqvap
            H_rain[i_lay] = -Lvap[i_lay]*dqvap/(cp*timestep)
            mass_cond -= dqvap*dmass[i_lay]               
            if verbose: print('dqvap < 0:',i_lay, dqvap, qvap, mass_cond, mass_cond/dmass[i_lay])
        else: # evaporation of rain
            mass_dvap = dqvap*dmass[i_lay]     
            if mass_dvap > mass_cond: # evaporate everything
                if verbose: print('mass_dvap > mass_cond:', i_lay, dqvap, qvap, mass_cond, mass_cond/dmass[i_lay], mass_dvap, dmass[i_lay],(mass_dvap > mass_cond), (tlay[i_lay] >= Tsat_p[i_lay]))
                qarray[idx_vap,i_lay] += mass_cond/dmass[i_lay]
                H_rain[i_lay] = - Lvap[i_lay] * mass_cond / (dmass[i_lay]*cp*timestep)
                mass_cond = 0.
            else:
                if verbose: print('mass_dvap < mass_cond:', i_lay, dqvap, qvap, mass_cond, mass_cond/dmass[i_lay])
                qarray[idx_vap,i_lay] += dqvap
                H_rain[i_lay] = -Lvap[i_lay]*dqvap/(cp*timestep)
                mass_cond -= dqvap*dmass[i_lay]
    if mass_cond != 0.:
        qarray[idx_cond,-1]+= mass_cond/dmass[-1]
        # Issue : qarray[idx_cond,-1] sometimes can become large (when starting with a hot atmosphere dominated with H2O that cools to saturation)
    if qvap_deep>=0.:
        qarray[idx_vap,-1] = qvap_deep
    return H_rain



@numba.jit(nopython=True, fastmath=True, cache=True)
def molecular_diffusion_numba(timestep, Nlay, p_lay, p_lev, dmass,
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
    