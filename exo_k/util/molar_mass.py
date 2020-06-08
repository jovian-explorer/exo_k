# -*- coding: utf-8 -*-
"""
@author: jeremy leconte
"""
from re import findall
from exo_k.util.singleton import Singleton

class Molar_mass(Singleton):
    """A class to compute the molar mass (in kg/mol) of regular molecules
    with a name written in a regular way (e.g. CO2, H2O, etc.).
    This class can also store the molar mass of custom gases with arbitrary names
    (for example: My_gas, earth_background).
    """
    
    def init(self, *args, **kwds):
        """Initializes empty dictionary of custom molecular masses.
        """
        self._custom_mol_mass={}

    def add_species(self, species_dict):
        """Add one or several species weights to the database. 

        Parameters
        ----------
            species_dict: dict
                Keys are gases names (they do not have to be real molecules).
                Values are molecular weight in kg/mol.
        """
        self._custom_mol_mass.update(species_dict)

    def __repr__(self):
        """Print the currently known species in the database. 
        """
        return self._custom_mol_mass.__repr__()

    def fetch(self, molecule_name):
        """Computes the molar mass of a molecule in kg/mol
        
        Parameters:
            molecule_name: str
                Name of the molecule.

        Returns:
            float:
                Molar mass in kg/mol.
        """
        if molecule_name in self._custom_mol_mass.keys():
            return self._custom_mol_mass[molecule_name]
        #s = re.findall('([A-Z][a-z]?)([0-9]*)', molecule_name)
        s = findall('([A-Za-z][a-z_]*)([0-9]*)', molecule_name)
        molecule_mass = 0
        for element, count in s:
            count = int(count or '1')
            try:
                molecule_mass += _atomic_mass_in_amu[element] * count
            except KeyError:
                print("""The following molecule name was not recognized: {mol}
                    If this is a custom gas
                    add its molecular weight to the library by running:
                    Molar_mass().add_species({{\'{mol}\':weight in kg/mol}})
                    """.format(mol=molecule_name) )
                raise RuntimeError()
        return molecule_mass*1.e-3 # conversion to kg/mol

_atomic_mass_in_amu = {
  "H":	1.00794,
  "He":	4.002602,
  "Li":	6.941,
  "Be":	9.012182,
  "B":	10.811,
  "C":	12.011,
  "N":	14.00674,
  "O":	15.9994,
  "F":	18.9984032,
  "Ne":	20.1797,
  "Na":	22.989768,
  "Mg":	24.3050,
  "Al":	26.981539,
  "Si":	28.0855,
  "P":	30.973762,
  "S":	32.066,
  "Cl":	35.4527,
  "Ar":	39.948,
  "K":	39.0983,
  "Ca":	40.078,
  "Sc":	44.955910,
  "Ti":	47.88,
  "V":	50.9415,
  "Cr":	51.9961,
  "Mn":	54.93805,
  "Fe":	55.847,
  "Co":	58.93320,
  "Ni":	58.6934,
  "Cu":	63.546,
  "Zn":	65.39,
  "Ga":	69.723,
  "Ge":	72.61,
  "As":	74.92159,
  "Se":	78.96,
  "Br":	79.904,
  "Kr":	83.80,
  "Rb":	85.4678,
  "Sr":	87.62,
  "Y":	88.90585,
  "Zr":	91.224,
  "Nb":	92.90638,
  "Mo":	95.94,
  "Tc":	98,
  "Ru":	101.07,
  "Rh":	102.90550,
  "Pd":	106.42,
  "Ag":	107.8682,
  "Cd":	112.411,
  "In":	114.82,
  "Sn":	118.710,
  "Sb":	121.757,
  "Te":	127.60,
  "I":	126.90447,
  "Xe":	131.29,
  "Cs":	132.90543,
  "Ba":	137.327,
  "La":	138.9055,
  "Ce":	140.115,
  "Pr":	140.90765,
  "Nd":	144.24,
  "Pm":	145,
  "Sm":	150.36,
  "Eu":	151.965,
  "Gd":	157.25,
  "Tb":	158.92534,
  "Dy":	162.50,
  "Ho":	164.93032,
  "Er":	167.26,
  "Tm":	168.93421,
  "Yb":	173.04,
  "Lu":	174.967,
  "Hf":	178.49,
  "Ta":	180.9479,
  "W":	183.85,
  "Re":	186.207,
  "Os":	190.2,
  "Ir":	192.22,
  "Pt":	195.08,
  "Au":	196.96654,
  "Hg":	200.59,
  "Tl":	204.3833,
  "Pb":	207.2,
  "Bi":	208.98037,
  "Po":	209,
  "At":	210,
  "Rn":	222,
  "Fr":	223,
  "Ra":	226.0254,
  "Ac":	227,
  "Th":	232.0381,
  "Pa":	213.0359,
  "U":	238.0289,
  "Np":	237.0482,
  "Pu":	244,
  "Am":	243,
  "Cm":	247,
  "Bk":	247,
  "Cf":	251,
  "Es":	252,
  "Fm":	257,
  "Md":	258,
  "No":	259,
  "Lr":	260,
  "Rf":	261,
  "Db":	262,
  "Sg":	263,
  "Bh":	262,
  "Hs":	265,
  "Mt":	266,
}

