# -*- coding: utf-8 -*-
"""
@author: jeremy leconte
Class for Rayleigh opacties. Only for H2 and He at the moment. 
"""
import numpy as np
from .util.cst import PI, KBOLTZ
from .util.singleton import Singleton

class Rayleigh(Singleton):
    """Class to access Rayleigh opacities
    """

    def init(self, *args, **kwds):
        """Initializes various parameters for Rayleigh computations
        """
        N_std=1.e5/(KBOLTZ*273.15)
        self._mult_factor=24.*PI**3/(N_std)**2
        self._mult_factor=self._mult_factor*100.**4
        # last 100.**4 is because we are using wavenumbers in cm^-1
        # instead of wavelength in m (see eq 12 of Caldas 2019)
        self._ABC={'He':np.array([0.2283])}

    def sigma(self,wns,abundances_dict):
        """Computes the Rayleigh cross section for the gas.
        Parameters:
            wns: array
                array of wavenumbers
            abundances_dict: dict
                Keys are molecule names. Values are the volume mixing ratios
        Output:
            res: array
                Rayleigh cross section for the whole gas in m^2/molecule
        """
        res=np.zeros(wns.size)
        wave = 1e8/wns
        for mol,x in abundances_dict.items():
            if mol is 'H2':
                res+=x*((8.14E-13)*(wave**(-4.))* \
                    (1+(1.572E6)*(wave**(-2.))+(1.981E12)*(wave**(-4.))))*1E-4
                #res+=x*(8.14e-33+1.28e-42*wns**2+1.61e-51*wns**4)*wns**4/1.e16
            if mol is 'He':
                res+=x*((5.484E-14)*(wave**(-4.))*(1+(2.44E5)*(wave**(-2.))))*1E-4
            else:
                pass

        return res
