# -*- coding: utf-8 -*-
"""
@author: jeremy leconte
"""
from exo_k.util.singleton import Singleton
from .molar_mass import Molar_mass
from .cst import RGP

class Molar_heat_capacity(Singleton):
    """A class to compute the molar heat capacity at constant pressure
    (in J/mol/K) of regular molecules.
    This class can also store the molar mass of custom gases with arbitrary names
    (for example: My_gas, earth_background).
    """
    
    def init(self, *args, **kwds):
        """Initializes empty dictionary of custom heat capacities masses.
        """
        self._custom_molar_cp={}

    def add_species(self, species_dict):
        """Add one or several species to the database. 

        Parameters
        ----------
            species_dict: dict
                Keys are gases names (they do not have to be real molecules).
                Values are molar cp.
        """
        self._custom_molar_cp.update(species_dict)

    def __repr__(self):
        """Print the currently known species in the database. 
        """
        return self._custom_molar_cp.__repr__()

    def fetch(self, molecule_name):
        """yields the molar cp of a molecule
        
        Parameters:
            molecule_name: str
                Name of the molecule.

        Returns:
            float:
                Molar cp in J/mol/K
        """
        if molecule_name in self._custom_molar_cp.keys():
            return self._custom_molar_cp[molecule_name]
        if molecule_name in _specific_cp_at_300k.keys():
            Mmol = Molar_mass().fetch(molecule_name)
            molar_cp = _specific_cp_at_300k[molecule_name] * Mmol
            self._custom_molar_cp[molecule_name] = molar_cp
            return molar_cp
        else: # if we do not know, assume monoatomic gas
            return 1.5*RGP

# standard data at 300K
_specific_cp_at_300k = { \
'H2': 14310.,
'CO2': 744.,
'N2': 1040.,
'He': 5193.,
'H2O': 1864.,
'SO2': 640.,
'H2S': 1003.,
'CH4': 2226.,
'NH3': 2175,
'CO': 1040.,
'HCN': 1328.,
'O2': 918.,
'C2H6': 1763.,
'C2H2': 1575.,
}