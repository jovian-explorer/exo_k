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
        This one is faster than sigma_array, but can be used only when vmr values are floats.

        Parameters
        ----------
            wns: array
                array of wavenumbers

            vmr: dict of floats
                Keys are molecule names. Values are the volume mixing ratios

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
            if to_add: res+=x*tmp

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
                (6.4982e-5 + 3.0743305e-10/(1.44e10-wn2))**2 * 4./9.
                # 4./9. is approximately ((n+1)/(n**2+2))**2
                # from eq 12 of caldas et al. 2019
        elif mol=='O2':
            tmp=self._mult_factor * wn4 * (1.096 + 1.385e-11*wn2 + 1.448e-20*wn4) * \
                (2.1351e-4 + 0.218567e-10/(0.409e10-wn2))**2 * 4./9.
        elif mol=='CO2':
            tmp=self._mult_factor * wn4 * (1.1364 + 2.53e-11*wn2) * \
                (1.1427e-2 * (5.799e3/(1.6617e10-wn2) + 1.2005e2/(7.960e9-wn2) \
                    + 5.3334/(5.630e9-wn2)))**2 * 4./9.
        else:
            to_add=False
        return to_add, tmp
