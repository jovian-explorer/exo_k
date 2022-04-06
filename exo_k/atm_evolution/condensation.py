# -*- coding: utf-8 -*-
"""
@author: jeremy leconte
"""
import numpy as np
import numba
import exo_k.util.cst as cst


class Condensing_species(object):

    def __init__(self, Latent_heat_vaporization = 0., cp_vap = 0., Mvap = 0.,
            T_ref = 0., Psat_ref = 0., cp_cond = None):
        """

        Parameters
        ----------
            Latent_heat_vaporization: float
                specific Latent heat vaporization (J/kg)
            cp_vap: float
                Specific heat capacity of the vapor (J/kg/K)
            Mvap: float
                Molar mass of vapor (kg/mol)
            T_ref: float
                Reference temperature
            Psat_ref: float
                Saturation vapor pressure at the reference temperature (Pa)
            cp_cond: float (optional)
                Specific heat capacity of the condensate (J/kg/K).
                Assumed equal to cp_vap if not provided.

        """
        self.Latent_heat_vaporization = Latent_heat_vaporization
        self.cp_vap = cp_vap
        self.T_ref = T_ref
        self.Psat_ref = Psat_ref
        self.Mvap = Mvap
        self.Rvap =  cst.RGP / self.Mvap
        if cp_cond is None:
            self.cp_cond = cp_vap
        else:
            self.cp_cond = cp_cond
        self.delta_cp = self.cp_vap - self.cp_cond 
        self.delta_cp_R = self.delta_cp / self.Rvap
        self.LovR = self.Latent_heat_vaporization / self.Rvap
        self.c1 = self.LovR / self.T_ref - self.delta_cp_R
        self.c2 = self.delta_cp_R * self.T_ref - self.LovR 
 
    def Lvap(self, T):
        """Latent heat at temperature T

        Parameters
        ----------
            T: array
                Temperature in layers (K)
        """
        return self.Latent_heat_vaporization + self.delta_cp * (T-self.T_ref)

    def Psat(self, T):
        """Saturation vapor pressure for the condensing species

        Parameters
        ----------
            T: array
                Temperature in layers (K)
        """
        #return self.Psat_ref*np.exp( - self.LovR * (1./T - 1/self.T_ref))
        return self.Psat_ref*np.exp( self.c1 + self.c2/T + self.delta_cp_R*np.log(T/self.T_ref))

    def Tsat(self, P):
        """Boiling temperature for the condensing species

        NOT EXACT IF DELTA CP != 0.

        Parameters
        ----------
            P: array
                Pressure in layers (Pa)
        """
        return self.c2/(np.log(P/self.Psat_ref)-self.c1)


    def qsat(self, psat, p, epsilon):
        """Saturation vapor mass mixing ratio for the condensing species

        Parameters
        ----------
            psat: array
                saturation vapor pressure
            p: array
                pressure at the layer center
            epsilon: float or array
                Ratio of the molar mass of the vapor over the background molar mass
        """
        psat_tmp = np.core.umath.minimum(psat,p)
        fac =  p + (epsilon -1.) * psat_tmp
        qsat = epsilon * psat_tmp / fac
        return qsat

    def dPsat_dT(self, T):
        """Saturation vapor pressure derivative for the condensing species 

        Parameters
        ----------
            T: array
                Temperature in layers (K)
        """
        RT2 = self.Rvap * T * T
        Lvap = self.Lvap(T)
        psat = self.Psat(T)
        return psat * Lvap / RT2, psat, Lvap

    def dlnPsat_dlnT(self, T):
        """Saturation vapor pressure for the condensing species and its derivative 

        Parameters
        ----------
            T: array
                Temperature in layers (K)
        """
        RT = self.Rvap * T
        Lvap = self.Lvap(T)
        return Lvap / RT, Lvap

    def moist_adiabat(self, T, P, cp, Mgas):
        """Computes the threshold thermal gradient (d ln T / d ln P) for a moist atmosphere.

        Parameters
        ----------
            T: array
                Temperature in layers (K)
            P: array
                pressure at the layer center
            cp: float
                specific heat capacity at constant pressure of the background gas
            Mgas: float or array
                Molar mass of the background atmosphere
        
        Returns
        -------
            array
                Moist adiabat lapse rate
            Lvap: array
                Latent heat at temperature T
            psat: array
                Saturation vapor pressure for the condensing species (Pa)
            qsat: array
                Saturation vapor mass mixing ratio for the condensing species
            dqsat_dt: array
                Derivative of qsat with respect to temperature at fixed pressure
            q_crit: array
                Critical mass mixing ratio for the inhibition of moist convection
                (Eq. 17 of Leconte et al. 2017)

        """
        epsilon = self.Mvap / Mgas
        psat = self.Psat(T)
        qsat = self.qsat(psat, P, epsilon)
        dlnPsat_dlnT, Lvap = self.dlnPsat_dlnT(T)
        dqsat_dt = qsat**2 * P * dlnPsat_dlnT / (epsilon*psat*T)
        qa=1.-qsat
        qLvt = qsat * Lvap / T 
        fac = qsat * self.cp_vap + qa * cp + qLvt * dlnPsat_dlnT
        q_crit = epsilon / ((epsilon - 1.) * dlnPsat_dlnT)
        return ( (1.-qsat) * cst.RGP / Mgas + qLvt ) / fac, Lvap, psat, qsat, dqsat_dt, q_crit # is missing p/p_a terms
    
    def compute_condensation_parameters(self, T, P, Mgas):
        """Computes necessary quantities to compute
        large scale condensation.

        Parameters
        ----------
            T: array
                Temperature in layers (K)
            P: array
                pressure at the layer center
            Mgas: float or array
                Molar mass of the background atmosphere
        """
        epsilon = self.Mvap / Mgas
        psat = self.Psat(T)
        qsat = self.qsat(psat, P, epsilon)
        dlnPsat_dlnT, Lvap = self.dlnPsat_dlnT(T)
        dqsat_dt = qsat**2 * P * dlnPsat_dlnT / (epsilon*psat*T)
        return Lvap, qsat, dqsat_dt



    def __repr__(self):
        """Method to output parameters
        """
        output="""
        Lvap         : {lvap}
        cp (vap)     : {cp}
        T0        (K): {t0}
        Psat(T0) (Pa): {p0}""".format(lvap=self.Latent_heat_vaporization,
            cp=self.cp_vap, t0=self.T_ref, p0=self.Psat_ref)
        return output

class Condensation_Thermodynamical_Parameters(object):

    def __init__(self, Latent_heat_vaporization = 0., cp_vap = 0., Mvap = 0.,
            T_ref = 0., Psat_ref = 0., cp_cond = None):
        """

        Parameters
        ----------
            Latent_heat_vaporization: float
                specific Latent heat vaporization (J/kg)
            cp_vap: float
                Specific heat capacity of the vapor (J/kg/K)
            Mvap: float
                Molar mass of vapor (kg/mol)
            T_ref: float
                Reference temperature
            Psat_ref: float
                Saturation vapor pressure at the reference temperature (Pa)
            cp_cond: float (optional)
                Specific heat capacity of the condensate (J/kg/K).
                Assumed equal to cp_vap if not provided.

        """
        Rvap =  cst.RGP / Mvap
        if cp_cond is None:
            cp_cond = cp_vap
        else:
            cp_cond = cp_cond
        delta_cp = cp_vap - cp_cond 
        delta_cp_R = delta_cp / Rvap
        LovR = Latent_heat_vaporization / Rvap
        c1 = LovR / T_ref - delta_cp_R
        c2 = delta_cp_R * T_ref - LovR 
        self.th_params = np.array([cp_vap, Mvap, Rvap,
                Latent_heat_vaporization, T_ref, Psat_ref,
                delta_cp, delta_cp_R, c1, c2])

@numba.jit(nopython=True, fastmath=True, cache=True)
def Lvap_T(T, Latent_heat_vaporization, T_ref, delta_cp):
    """Latent heat at temperature T

    Parameters
    ----------
        T: array
            Temperature in layers (K)
    """
    return Latent_heat_vaporization + delta_cp * (T-T_ref)

@numba.jit(nopython=True, fastmath=True, cache=True)
def Psat_T(T, T_ref, Psat_ref, c1, c2, delta_cp_R):
    """Saturation vapor pressure for the condensing species

    Parameters
    ----------
        T: array
            Temperature in layers (K)
    """
    return Psat_ref*np.exp( c1 + c2/T + delta_cp_R*np.log(T/T_ref))

@numba.jit(nopython=True, fastmath=True, cache=True)
def Tsat_P(P, Psat_ref, c1, c2):
    """Boiling temperature for the condensing species

    NOT EXACT IF DELTA CP != 0.
    
    Parameters
    ----------
        P: array
            Pressure in layers (Pa)
    """
    return c2/(np.log(P/Psat_ref)-c1)

@numba.jit(nopython=True, fastmath=True, cache=True)
def Qsat(psat, p, epsilon):
    """Saturation vapor mass mixing ratio for the condensing species

    Parameters
    ----------
        psat: array
            saturation vapor pressure
        p: array
            pressure at the layer center
        epsilon: float or array
            Ratio of the molar mass of the vapor over the background molar mass
    """
    psat_tmp = np.core.umath.minimum(psat,p)
    fac =  p + (epsilon -1.) * psat_tmp
    qsat = epsilon * psat_tmp / fac
    return qsat

@numba.jit(nopython=True, fastmath=True, cache=True)
def dPsat_dT(T, Latent_heat_vaporization, T_ref, Psat_ref, Rvap, delta_cp, delta_cp_R, c1, c2):
    """Saturation vapor pressure derivative for the condensing species 

    Parameters
    ----------
        T: array
            Temperature in layers (K)
    """
    RT2 = Rvap * T * T
    Lvap = Lvap_T(T, Latent_heat_vaporization, T_ref, delta_cp)
    psat = Psat_T(T, T_ref, Psat_ref, c1, c2, delta_cp_R)
    return psat * Lvap / RT2, psat, Lvap

@numba.jit(nopython=True, fastmath=True, cache=True)
def dlnPsat_dlnT(T, Latent_heat_vaporization, T_ref, delta_cp, Rvap):
    """Saturation vapor pressure for the condensing species and its derivative

    Parameters
    ----------
        T: array
            Temperature in layers (K)
    """
    RT = Rvap * T
    Lvap = Lvap_T(T, Latent_heat_vaporization, T_ref, delta_cp)
    return Lvap / RT, Lvap

@numba.jit(nopython=True, fastmath=True, cache=True)
def moist_adiabat(T, P, cp, Mgas, cp_vap, Mvap, Rvap, Latent_heat_vaporization, T_ref, Psat_ref, delta_cp, delta_cp_R, c1, c2):
    """Computes the threshold thermal gradient (d ln T / d ln P) for a moist atmosphere.

    Parameters
    ----------
        T: array
            Temperature in layers (K)
        P: array
            pressure at the layer center
        cp: float
            specific heat capacity at constant pressure of the background gas
        Mgas: float or array
            Molar mass of the background atmosphere
    
    Returns
    -------
        array
            Moist adiabat lapse rate
        Lvap: array
            Latent heat at temperature T
        psat: array
            Saturation vapor pressure for the condensing species (Pa)
        qsat: array
            Saturation vapor mass mixing ratio for the condensing species
        dqsat_dt: array
            Derivative of qsat with respect to temperature at fixed pressure
        q_crit: array
            Critical mass mixing ratio for the inhibition of moist convection
            (Eq. 17 of Leconte et al. 2017)
    """
    epsilon = Mvap / Mgas
    psat = Psat_T(T, T_ref, Psat_ref, c1, c2, delta_cp_R)
    qsat = Qsat(psat, P, epsilon)
    dlpsat_dlT, Lvap = dlnPsat_dlnT(T, Latent_heat_vaporization, T_ref, delta_cp, Rvap)
    dqsat_dt = qsat**2 * P * dlpsat_dlT / (epsilon*psat*T)
    qa=1.-qsat
    qLvt = qsat * Lvap / T 
    fac = qsat * cp_vap + qa * cp + qLvt * dlpsat_dlT
    q_crit = epsilon / ((epsilon - 1.) * dlpsat_dlT)
    return ( (1.-qsat) * cst.RGP / Mgas + qLvt ) / fac, Lvap, psat, qsat, dqsat_dt, q_crit # is missing p/p_a terms

@numba.jit(nopython=True, fastmath=True, cache=True)
def compute_condensation_parameters(T, P, Mgas, Mvap, Rvap, Latent_heat_vaporization, T_ref, Psat_ref, delta_cp, delta_cp_R, c1, c2):
    """Computes necessary quantities to compute
    large scale condensation.

    Parameters
    ----------
        T: array
            Temperature in layers (K)
        P: array
            pressure at the layer center
        Mgas: float or array
            Molar mass of the background atmosphere
    """
    epsilon = Mvap / Mgas
    psat = Psat_T(T, T_ref, Psat_ref, c1, c2, delta_cp_R)
    qsat = Qsat(psat, P, epsilon)
    dlpsat_dlT, Lvap = dlnPsat_dlnT(T, Latent_heat_vaporization, T_ref, delta_cp, Rvap)
    dqsat_dt = qsat**2 * P * dlpsat_dlT / (epsilon*psat*T)
    return Lvap, qsat, dqsat_dt