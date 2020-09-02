# -*- coding: utf-8 -*-
"""
@author: jeremy leconte

Module containing a class to handle a database of CIA tables to compute opacities with it.
"""
import numpy as np

from .atable import Atable
from .util.spectral_object import Spectral_object

class Adatabase(Spectral_object):
    """Class to group :class:`~exo_k.atable.Atable` objects
    and combine them in radiative transfer
    """

    def __init__(self, *str_filters, filenames=None, remove_zeros=False, **kwargs):
        """Initializes aerosol tables and supporting data from a list of filenames.

        Parameters
        ----------
            filenames: list
                List of names (not full path) of the input cia files. 
                The files must be in the global search path.
                
        A local search path can be specified with 'search_path='

        See the options of :class:`~exo_k.atable.Atable` __init__ method.
        """
        self.atables={}
        self.r_eff_unit=None
        self.wns=None
        self.Nw=None
        if filenames is not None:
            for filename in filenames:
                tmp_atable=Atable(*([filename]+list(str_filters)),
                    remove_zeros=remove_zeros, **kwargs)
                self.add_atables(tmp_atable)
    
    def add_atables(self, *atables):
        """Adds new :class:`~exo_k.atable.Atable` objects to a Adatabase.

        Parameters
        ----------
            atables: :class:`~exo_k.atable.Atable`
                As many Atables as you want.
        """
        for atable in atables:
            if self.r_eff_unit is None:
                self.r_eff_unit=atable.r_eff_unit
            elif self.r_eff_unit!=atable.r_eff_unit:
                print('Adatabase units for r_eff: {k}'.format(k=self.r_eff_unit))
                print('{aer} units for cia coeff: {k}'.format(\
                    aer=atable.aerosol_name, k=atable.r_eff_unit))
                raise RuntimeError("""You naughty:
                all Atables in a database must have the same r_eff units""")
            # should ultimately check for wn_unit as well !

            self.atables[atable.aerosol_name]=atable

    @property
    def names(self):
        """Gives the names of the aerosols in the database
        """
        return list(self.atables.keys())

    def __getitem__(self, name):
        """Overrides getitem so that Adatabase['name'] directly accesses 
        the database for that aerosol name.
        """
        if name not in self.atables:
            raise KeyError('The requested aerosol is not available.')
        return self.atables[name]

    def copy(self):
        """Creates a new instance of :class:`CIAdatabase` object and (deep) copies data into it
        """
        res=Adatabase(None)
        for atab in self.atables.values():
            res.add_atables(atab.copy())
        return res

    def __repr__(self):
        """Method to output
        """
        output='The available aerosols are: \n'
        for aer,atab in self.atables.items():
            output+=aer+'->'+atab.filename+'\n'
        if self.wns is not None:
            output+='All tables share a common spectral grid\n'
        else:
            output+='All tables do NOT have common spectral grid\n'
            output+='You will need to run sample before using the database\n'

        return output

    def sample(self, wngrid, remove_zeros=True, use_grid_filter=True,
            sample_all_vars=False):
        """Samples all the Atables in the database on the same wavenumber grid
        to be able to use them in radiative transfer modules.

        .. important::
            For this method, the default options are remove_zeros=True and use_grid_filter=True.
            This ensures that you do not count several contributinns twice when
            you use two different aerosol tables to represent
            the same aerosol type over two different
            spectral window (typically VI and IR)

        Parameters
        ----------
            wngrid : array
                new wavenumber grid (cm-1)

        .. seealso::
            See :func:`exo_k.atable.Atable.sample for further details on options.
        """
        for atable in self.atables.values():
            atable.sample(wngrid, remove_zeros=remove_zeros, use_grid_filter=use_grid_filter,
                sample_all_vars=sample_all_vars)
        self.wns=np.array(wngrid)
        self.Nw=self.wns.size

    def convert_to_mks(self):
        """Converts units of all Cia_tables to MKS.
        """
        first=True
        for atable in self.atables.values():
            atable.convert_to_mks()
            if first:
                self.r_eff_unit=atable.r_eff_unit
        

    def absorption_coefficient(self, aer_reffs_densities,
            wngrid_limit=None, log_interp=None):
        """Computes the absorption coefficient in m^-1 for the whole mix specified 
        (assumes data in MKS).

        Parameters
        ----------
            aer_reffs_densities: dict
                A dictionary with aerosol names as keys and lists containing 2
                floats (or arrays) as values. The values are the particle effective radii
                and number densities.

            wngrid_limit: array, optional
                Smaller and bigger wavenumbers inside which to perform the calculation.

        Returns:
            array:
                The effective cross section coefficient profile for the aerosols (in m^2).
        """
        first=True
        if self.wns is None: raise RuntimeError("""Atables must be sampled
            on the same grid before calling cross_section.
            Please use the sample() method.""")
        for aer, values in aer_reffs_densities.items():
            if aer in self.atables.keys():
                if first:
                    res=self.atables[aer].absorption_coefficient(values[0], values[1],
                        wngrid_limit=wngrid_limit, log_interp=log_interp)
                    first=False
                else:
                    res+=self.atables[aer].absorption_coefficient(values[0], values[1],
                        wngrid_limit=wngrid_limit, log_interp=log_interp)
        if first: # means that no molecule was in the database, we need to initialize res
            res=np.zeros((self.wns.size))
        return res

        

