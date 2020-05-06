# -*- coding: utf-8 -*-
"""
@author: jeremy leconte

Module containing a class to handle a database of CIA tables to compute opacities with it.
"""
import numpy as np

from .cia_table import Cia_table

class CIAdatabase(object):
    """Class to group CIA tables and combine them in radiative transfer
    """

    def __init__(self, filenames, *str_filters, remove_zeros=True, **kwargs):
        """Initializes cia tables and supporting data from a list of filenames.

        Parameters
        ----------
            filenames: list
                List of names (not full path) of the input cia files. 
                The files must be in the global search path.
                
        A local search path can be specified with 'search_path='

        See the options of Cia_table.__init__()
        """
        self.cia_tables={}
        self.wns=None
        self.Nw=None
        self.abs_coeff_unit=None
        for filename in filenames:
            tmp_cia_table=Cia_table(*([filename]+list(str_filters)),
                remove_zeros=remove_zeros, **kwargs)
            self.add_cia_tables(tmp_cia_table)
    
    def add_cia_tables(self, *cia_tables):
        """Adds news cia tables to a CIA database.

        Parameters
        ----------
            cia_tables: :class:`CIA_table`
                As many cia tables as you want.
        """
        for cia_table in cia_tables:
            if self.abs_coeff_unit is None:
                self.abs_coeff_unit=cia_table.abs_coeff_unit
            elif self.abs_coeff_unit!=cia_table.abs_coeff_unit:
                print('CIAdatabase units for cia coeff: {k}'.format(k=self.abs_coeff_unit))
                print('{mol1}-{mol2} units for cia coeff: {k}'.format(\
                    mol1=cia_table.mol1, mol2=cia_table.mol2, k=cia_table.abs_coeff_unit))
                raise RuntimeError("""You naughty:
                all CIAdatabase in a database must have the same coeff units""")

            if cia_table.mol1 in self.cia_tables:
                if cia_table.mol2 in self.cia_tables[cia_table.mol1]:
                    continue
                else:
                    self.cia_tables[cia_table.mol1][cia_table.mol2]=cia_table
            else:
                self.cia_tables[cia_table.mol1]={cia_table.mol2:cia_table}

    def __getitem__(self, molecule):
        """Overrides getitem so that CIAdatabase['mol'] directly accesses 
        the database for that molecule.
        """
        if molecule not in self.cia_tables:
            raise KeyError('The requested molecule is not available.')
        return self.cia_tables[molecule]

    def sample(self, wngrid, remove_zeros=False):
        """Samples all the cia_table in the database on the same wavenumber grid
        to be able to use them in radiative transfer modules.

        Parameters
        ----------
            wngrid : array
                new wavenumber grid (cm-1)
        """
        for mol1 in self.cia_tables.values():
            for cia_table in mol1.values():
                cia_table.sample(wngrid,remove_zeros)
        self.wns=np.array(wngrid)
        self.Nw=self.wns.size

    def convert_to_mks(self):
        """Converts units of all Cia_tables to MKS.
        """
        first=True
        for mol1 in self.cia_tables.values():
            for cia_table in mol1.values():
                cia_table.convert_to_mks()
                if first:
                    self.abs_coeff_unit=cia_table.abs_coeff_unit
        

    def cia_cross_section(self, logP_array, T_array, gas_comp, wngrid_limit=None):
        """Computes the absorption coefficient in m^-1 for the whole mix specified 
        (assumes data in MKS).

        Parameters
        ----------
            logP_array: array

            T_array   : array
                log10 Pressure and temperature profiles

            gas_comp  : :class:`chemistry.Gas_mix` object
                behaves like a dict with mol names as keys and vmr as values

            wngrid_limit: array, optional
                Smaller and bigger wavenumbers inside which to perform the calculation

        Returns:
            array:
                The cia effective cross section coefficient profile for the whole gas (in m^2).
        """
        logP_array=np.array(logP_array)
        T_array=np.array(T_array)
        first=True
        if self.wns is None: raise RuntimeError("""cia tables must be sampled
            on the same grid before calling cia_cross_section.
            Please use the sample() method.""")
        for mol,x1 in gas_comp.items():
            if mol in self.cia_tables.keys():
                for mol2 in self.cia_tables[mol].keys():
                    if mol2 in gas_comp.keys():
                        x2=gas_comp[mol2]
                        if first:
                            res=self.cia_tables[mol][mol2].effective_cross_section( \
                                logP_array,T_array,x1,x2,wngrid_limit=wngrid_limit)
                            first=False
                        else:
                            res+=self.cia_tables[mol][mol2].effective_cross_section( \
                                logP_array,T_array,x1,x2,wngrid_limit=wngrid_limit)
        return res

        

