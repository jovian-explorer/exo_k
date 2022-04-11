# -*- coding: utf-8 -*-
"""
@author: jeremy leconte

Class for Rayleigh opacties. 
"""
import numpy as np
from .util.cst import PI, KBOLTZ
from .util.singleton import Singleton

class Rayleigh(Singleton):
    """Class to compute Rayleigh opacities
    """

    def init(self, *args, **kwds):
        """Initializes various parameters for Rayleigh computations
        """
        N_std=1.e5/(KBOLTZ*273.15)
        self._mult_factor=24.*PI**3/(N_std)**2
        self._mult_factor=self._mult_factor*100.**4
        # last 100.**4 is because we are using wavenumbers in cm^-1
        # instead of wavelength in m (see eq 12 of Caldas 2019)

    def sigma(self, wns, vmr):
        """Computes the Rayleigh cross section for the gas.
        This one is faster than sigma_array, but can be used only
        when vmr values are constants.

        Parameters
        ----------
            wns: array
                array of wavenumbers

            vmr: dict of arrays
                Keys are molecule names. Values are the volume mixing ratios.
                For speedup, only the first value will be used because we assume
                that the vmr arrays are constant

        Returns
        -------
            array of shape (wns.size)
                Rayleigh cross section for the whole gas in m^2/molecule
        """
        res=np.zeros(wns.size)
        wn2 = wns*wns
        wn4 = wn2*wn2
        for mol, x in vmr.items():
            to_add, tmp = self.sigma_mol(mol, wn2, wn4)
            if to_add: res+=x[0]*tmp

        return res

    def sigma_array(self, wns, vmr):
        """Computes the Rayleigh cross section for the gas.

        Parameters
        ----------
            wns: array
                array of wavenumbers

            vmr: dict of arrays
                Keys are molecule names. Values are arrays the volume mixing ratios

        Returns
        -------
            array of shape (vmr.values.size, wns.size)
                Rayleigh cross section for the whole gas in m^2/molecule
        """
        first_mol=True
        wn2 = wns*wns
        wn4 = wn2*wn2
        for mol, x in vmr.items():
            x_array=np.array(x)
            if first_mol:
                res=np.zeros((x_array.size,wns.size))
                first_mol=False
            to_add, tmp = self.sigma_mol(mol, wn2, wn4)
            if to_add: res+=x_array[:,None]*tmp

        return res

    def sigma_mol(self, mol, wn2, wn4):
        """Intermediary function to compute rayleigh for each molecule.

        Parameters
        ----------
            mol: str
                Molecule name.
            wn2, wn4: arrays
                Array of the wavenumber (in cm^-1) to the 2nd and 4th power.
                (To avoid recomputing it each time).
        
        Returns
        -------
            to_add: bool
                Says whether the molecule has been found
                and the contribution needs to be added.
            tmp: array of size self.wns or None
                The cross section for the molecule as a function of wns.
                None if the molecule has not
                been found.

        """
        to_add=True
        tmp=None
        if mol == 'H2':
            #tmp=((8.14E-13)*(wave**(-4.))* \
            #    (1+(1.572E6)*(wave**(-2.))+(1.981E12)*(wave**(-4.))))*1E-4
            tmp=(8.14e-49+1.28e-58*wn2+1.61e-67*wn4)*wn4
        elif mol == 'He':
            #tmp=((5.484E-14)*(wave**(-4.))*(1+(2.44E5)*(wave**(-2.))))*1E-4
            tmp=(5.484e-50+1.338e-60*wn2)*wn4
        elif mol=='N2':
            tmp=self._mult_factor * wn4 * (1.034 + 3.17e-12*wn2) * \
                (6.4982e-5 + 3.0743305e6/(1.44e10-wn2))**2 * 4./9.
                # 4./9. is approximately ((n+1)/(n**2+2))**2
                # from eq 12 of caldas et al. 2019
    #    elif mol=='O2':
    #        tmp=mult_factor * wn4 * (1.096 + 1.385e-11*wn2 + 1.448e-20*wn4) * \
    #           (2.1351e-4 + 0.218567e-6/(0.409e10-wn2))**2 * 4./9.
        elif mol=='O2':
            tmp=self._mult_factor * wn4 * (1.096 + 1.385e-11*wn2 + 1.448e-20*wn4) * \
                (2.05648e-4 + 2.480899e5/(0.409e10-wn2))**2 * 4./9.
        elif mol=='CO2':
            # Error corrected "by hand" as there was a 1e6 departure between formula and measured value.
            # The error seems to be in the ref used by Caldas et al.
            # Also Caldas et al. mistakenly used the (n**2 -1) formula as the formula for
            # ((n**2 -1)/(n**2 +2))**2 hence a 4./9. error
            tmp=self._mult_factor * wn4 * (1.1364+2.53e-11*wn2)* 1.e-6*\
                (1.1427e6 * (5.799e3/(1.661750e10-wn2) + 1.2005e2/(7.960886e9-wn2) \
                + 5.3334/(5.630626e9-wn2)  + 4.3244/(4.601954e9-wn2)  + 1.218145e-5/(5.847382e6-wn2)))**2 * 4./9.   
        elif mol=='H2O':
            # this formula is only valid for wavelengths > 0.23 microns
            # Caldas+2019 includes another formula for shorter wavelengths
            tmp=self._mult_factor * wn4 * (4.92303e6/(2.380185e10-wn2) + 1.42723e5/(5.73262e9-wn2))**2 * 4./9.
        elif mol=='CH4':
            tmp=self._mult_factor * wn4 * (4.6662e-4+4.02e-14*wn2)**2 *  4./9.
        elif mol=='CO':
            # The present formula does not reproduce measured values by a factor of ~2
            # The reason is not yet understood.
            # All formula (except for H2 and He) come from Sneep and Ubachs JQSRT 2005.
            tmp=self._mult_factor * wn4 * (2.2851e-4 + 0.456e4/(5.101816e9 - wn2))**2 *  4./9.
        else:
            to_add=False
        return to_add, tmp
