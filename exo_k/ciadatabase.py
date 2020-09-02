# -*- coding: utf-8 -*-
"""
@author: jeremy leconte

Module containing a class to handle a database of CIA tables to compute opacities with it.
"""
import numpy as np

from .cia_table import Cia_table
from .settings import NoFileFoundError

class CIAdatabase(object):
    """Class to group :class:`~exo_k.cia_table.Cia_table` objects
    and combine them in radiative transfer
    """

    def __init__(self, *str_filters, filenames=None, molecule_pairs=None,
        molecules=None, remove_zeros=True, **kwargs):
        """Initializes cia tables and supporting data.
        The files to load can be specified either using (in order of precedence)
        a list of filenames, molecule pairs, or molecules. See below for details.

        Any number of filters (strings) can be provided to refine the search to
        files whose names contain the filters. 

        Parameters
        ----------
            filenames: list, optional
                List of names (not full path) of the input cia files. 
                The files must be in the global search path.
            molecule_pairs: list of size 2 lists, optional
                List of the molecule pairs we want to consider, 
                specified as an array with two strings (like ['H2','H2'] or ['N2','H2O']).
                The order of the molecules in the pair is irrelevant. 
            molecules: list of strings, optional
                A list of all the molecules we want to consider. `Exo_k` will look
                for all the possible pairs for which a cia file is found. 

        The default path searched is set with :func:`exo_k.settings.Settings.set_cia_search_path`
        or :func:`~exo_k.settings.Settings.add_cia_search_path`.
        A local search path can be specified with the `search_path` keyword.

        See the options of :class:`~exo_k.cia_table.Cia_table` __init__ method.
        """
        self.cia_tables={}
        self.wns=None
        self.Nw=None
        self.abs_coeff_unit=None
        if filenames is not None:
            for filename in filenames:
                tmp_cia_table=Cia_table(*([filename]+list(str_filters)),
                    remove_zeros=remove_zeros, **kwargs)
                self.add_cia_tables(tmp_cia_table)
        elif molecule_pairs is not None:
            for mol_pair in molecule_pairs:
                tmp_cia_table=Cia_table(*str_filters, molecule_pair=mol_pair,
                    remove_zeros=remove_zeros, **kwargs)
                self.add_cia_tables(tmp_cia_table)
        elif molecules is not None:
            for mol1 in molecules:
                for mol2 in molecules:
                    try:
                        tmp_cia_table=Cia_table(*str_filters, molecule_pair=[mol1,mol2],
                            remove_zeros=remove_zeros, **kwargs)
                        self.add_cia_tables(tmp_cia_table)
                    except NoFileFoundError:
                        pass


    
    def add_cia_tables(self, *cia_tables):
        """Adds new :class:`~exo_k.cia_table.Cia_table` objects to a CIA database.

        Parameters
        ----------
            cia_tables: :class:`~exo_k.cia_table.Cia_table`
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
                all Cia_tables in a database must have the same coeff units""")

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

    def copy(self):
        """Creates a new instance of :class:`CIAdatabase` object and (deep) copies data into it
        """
        res=CIAdatabase(None)
        for dic_cia_tabs in self.cia_tables.values():
            for cia_tab in dic_cia_tabs.values():
                res.add_cia_tables(cia_tab.copy())
        return res

    def __repr__(self):
        """Method to output
        """
        output='The available molecule pairs are: \n'
        for mol1, dico in self.cia_tables.items():
            for mol2, cia_tab in dico.items():
                output+=mol1+'-'+mol2+'->'+cia_tab.filename+'\n'
        if self.wns is not None:
            output+='All tables share a common spectral grid\n'
        else:
            output+='All tables do NOT have common spectral grid\n'
            output+='You will need to run sample before using the database\n'

        return output

    def sample(self, wngrid, remove_zeros=False, use_grid_filter=False):
        """Samples all the cia_table in the database on the same wavenumber grid
        to be able to use them in radiative transfer modules.

        Parameters
        ----------
            wngrid : array
                new wavenumber grid (cm-1)

        .. seealso::
            See :func:`exo_k.cia_table.Cia_table.sample for further details on options.

        """
        for mol1 in self.cia_tables.values():
            for cia_table in mol1.values():
                cia_table.sample(wngrid, remove_zeros=remove_zeros,
                    use_grid_filter=use_grid_filter)
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

            T_array: array
                log10 Pressure and temperature profiles

            gas_comp: :class:`~exo_k.gas_mix.Gas_mix` object
                behaves like a dict with mol names as keys and vmr as values.

            wngrid_limit: array, optional
                Smaller and bigger wavenumbers inside which to perform the calculation.

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
                                logP_array, T_array, x1, x2, wngrid_limit=wngrid_limit)
                            first=False
                        else:
                            res+=self.cia_tables[mol][mol2].effective_cross_section( \
                                logP_array, T_array, x1, x2, wngrid_limit=wngrid_limit)
        if first: # means that no molecule was in the database, we need to initialize res
            res=np.zeros((T_array.size,self.wns.size))
        return res

        

