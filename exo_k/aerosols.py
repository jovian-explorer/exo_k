# -*- coding: utf-8 -*-
"""
@author: jeremy leconte
"""
import numpy as np
from .util.spectral_object import Spectral_object

class Aerosols(Spectral_object):
    """Dict-like class to handle gas composition (with background gas) and molar mass.

    If `logp_array`, `t_array`, and radiative databases are provided, :any:`cross_section`
    can be used to compute the opacity of the gas
    """

    def __init__(self, aer_reffs_densities={}, a_database=None):
        """__init__ Instantiates
        an aerosols object.
        """
        self.set_a_database(a_database=a_database)
        self.set_aer_reffs_densities(aer_reffs_densities=aer_reffs_densities)
        self._wn_range=None

    def set_aer_reffs_densities(self, aer_reffs_densities={}):
        self.aer_reffs_densities = aer_reffs_densities

    def set_a_database(self, a_database=None):
        """Change the radiative database attached to the current instance of aerosols

        Parameters
        ----------
            a_database: :class:`~exo_k.kdatabase.Kdatabase` object
                New Adatabase to use.
        """
        self.adatabase=a_database
        if self.adatabase is not None and self.adatabase.r_eff_unit != 'm' :
            print("""
            You're being Bad!!! You are trying *NOT* to use MKS units
            for the effective radii in the aerosol database!!!""")
            raise RuntimeError("Bad units in the Adatabase used with aerosols.")

    def _compute_spectral_range(self, wn_range=None, wl_range=None):
        """Converts an unordered spectral range in either wavenumber or wavelength
        in an ordered wavenumber range.

        Parameters
        ----------
            wn_range: list or array of size 2
                Minimum and maximum wavenumber (in cm^-1).
            wl_range: list or array of size 2
                Minimum and maximum wavelength (in micron)
        """
        if wl_range is not None:
            if wn_range is not None:
                print('Cannot specify both wl and wn range!')
                raise RuntimeError()
            else:
                _wn_range=np.sort(10000./np.array(wl_range))
        else:
            if wn_range is not None:
                _wn_range=np.sort(np.array(wn_range))
            else:
                _wn_range=self._wn_range
        return _wn_range

    def optical_properties(self, aer_reffs_densities=None, wl_range=None, wn_range=None,
            log_interp=True, compute_all_opt_prop=True, **kwargs):
        """bla
        """
        if self.adatabase is None: raise RuntimeError("""
            a_database not provided. 
            Use the a_database keyword during initialization or use the set_a_database method.""")
        if self.adatabase.wns is None: raise RuntimeError("""
            All tables in the Adatabase should have the same wavenumber grid to proceed.
            You should probably use sample().""")
        if aer_reffs_densities is not None:
            self.set_aer_reffs_densities(aer_reffs_densities=aer_reffs_densities)

        local_wn_range=self._compute_spectral_range(wl_range=wl_range, wn_range=wn_range)

        [k, k_scat, g] = self.adatabase.optical_properties(self.aer_reffs_densities,
            wngrid_limit=local_wn_range, log_interp=log_interp,
            compute_all_opt_prop=compute_all_opt_prop)
        return [k, k_scat, g]
    
    def absorption_coefficient(self, aer_reffs_densities=None, wl_range=None, wn_range=None,
            log_interp=True, **kwargs):
        """bla
        """
        if self.adatabase is None: raise RuntimeError("""
            a_database not provided. 
            Use the a_database keyword during initialization or use the set_a_database method.""")
        if self.adatabase.wns is None: raise RuntimeError("""
            All tables in the Adatabase should have the same wavenumber grid to proceed.
            You should probably use sample().""")
        if aer_reffs_densities is not None:
            self.set_aer_reffs_densities(aer_reffs_densities=aer_reffs_densities)

        local_wn_range=self._compute_spectral_range(wl_range=wl_range, wn_range=wn_range)

        return self.adatabase.absorption_coefficient(self.aer_reffs_densities,
            wngrid_limit=local_wn_range, log_interp=log_interp)

