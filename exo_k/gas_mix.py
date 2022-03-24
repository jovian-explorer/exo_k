# -*- coding: utf-8 -*-
"""
@author: jeremy leconte
"""
import numpy as np
from .util.molar_mass import Molar_mass
from .util.molar_heat_capacity import Molar_heat_capacity
from .util.interp import RandOverlap_2_kdata_prof, rm_molec
from .rayleigh import Rayleigh
from .util.spectral_object import Spectral_object
from .util.singleton import Singleton

class Gas_mix(Spectral_object):
    """Dict-like class to handle gas composition (with background gas) and molar mass.

    If `logp_array`, `t_array`, and radiative databases are provided, :any:`cross_section`
    can be used to compute the opacity of the gas
    """

    def __init__(self, composition=None, logp_array=None, t_array=None,
        k_database=None, cia_database=None):
        """__init_ Instantiates
        a Gas_mix object and computes the vmr of the 'background' gas.
        """
        if composition is None: composition=dict()
        self.Narray=None
        self._wn_range=None
        self.iw_min=None
        self.iw_max=None
        self.wns=None
        self.wnedges=None
        self.Nw=None
        self.need_to_recompute_vmr_array = True

        self.set_composition(composition)
        self.set_logPT(logp_array=logp_array, t_array=t_array)
        self.set_k_database(k_database)
        self.set_cia_database(cia_database)


    def set_composition(self, composition):
        """Reset composition and computes the vmr of the 'background' gas.

        Parameters
        ----------
            composition: dict
                Keys are molecule names. Values are vmr or arrays of vmr.
            bg_gas: str
                Name of the background molecule. 
                If None, it is inferred from the molecule for which vmr='background'.
        """
        self.composition = dict()
        self.bg_gas = None
        for mol,vmr in composition.items():
            if isinstance(vmr, list):
                self.composition[mol]=np.array(vmr)
            elif isinstance(vmr, str):
                self.bg_gas = mol
            else:
                self.composition[mol] = vmr
        if (self.bg_gas is not None):
            self.composition['inactive_gas'] = 0.
        else:
            self.composition['inactive_gas'] = 'background'
            self.bg_gas = 'inactive_gas'
        self.get_background_vmr()
        self.need_to_recompute_vmr_array = True

    def set_logPT(self, logp_array=None, t_array=None):
        """Sets the pressure (in Pa) and temperature fields.
        """
        if logp_array is not None:
            if hasattr(logp_array, "__len__"):
                self.logp_array=np.array(logp_array, dtype=float)
            else:
                self.logp_array=np.array([logp_array], dtype=float)

            if hasattr(t_array, "__len__"):
                self.t_array=np.array(t_array)
            else:
                self.t_array=np.array([t_array])

            if self.logp_array.size != self.t_array.size:
                raise TypeError(\
                    'logp_array and t_array should be 1d lists or arrays of the same size')
            self.Narray=self.logp_array.size

    def get_background_vmr(self):
        """Computes the volume mixing ratio of the background gas in a mix
        Uses the fact that self.composition is a dictionary
        with the name of the molecules as keys and the vmr as values.
        vmr are either float or arrays.
        
        At this stage, the background gas should be identified by `self.bg_gas`,
        and its Vol. Mix. Ratio will be updated in the composition dict. 
        """
        other_vmr=0.
        for mol,vmr in self.composition.items():
            if mol==self.bg_gas:
                continue
            try:
                other_vmr+=vmr
            except ValueError:
                raise TypeError('Incompatible shapes in Gas_mix arrays.')
        #if np.amax(other_vmr)>1.:
        #    print("""Careful: the sum of the vmr of your gas components is > 1.
        #    If there is a background gas, its vmr will become negative.
        #    I hope you know what you are doing.""")
        self.composition[self.bg_gas]=1.-other_vmr

    def normalize(self):
        """Renormalizes the vmr of all the gases so that
        the total be equal to 1.

        This works only if the vmr of 'inactive_gas' is not 0.
        """
        vmr_inac = self.composition['inactive_gas']
        for mol,vmr in self.composition.items():
            if mol!='inactive_gas':
                self.composition[mol]=vmr/(1.-vmr_inac)
            else:
                self.composition[mol]=0.

    def molar_mass(self):
        """Computes and returns the molar mass of a mix of gases

        Returns
        -------
            float or array:
                Molar mass of the active gases in kg/mol
        """
        mol_mass_active_gases=0.
        vmr_active_gases=0.

        for mol,vmr in self.composition.items():
            if mol!='inactive_gas':
                Mmol=Molar_mass().fetch(mol)
                mol_mass_active_gases+=vmr*Mmol
                vmr_active_gases+=vmr
        mol_mass=mol_mass_active_gases/vmr_active_gases     
        return mol_mass

    def cp(self):
        """Computes and returns the specific heat capacity (cp)
        of a mix of gases

        Returns
        -------
            float or array:
                Specific cp of the gas mix
        """
        mol_cp_active_gases=0.
        vmr_active_gases=0.

        for mol,vmr in self.composition.items():
            if mol != 'inactive_gas':
                cp_mol = Molar_heat_capacity().fetch(mol)
                mol_cp_active_gases += vmr * cp_mol
                vmr_active_gases += vmr
        mol_mass = self.molar_mass()
        return mol_cp_active_gases/vmr_active_gases/mol_mass    

    def get_vmr_array(self, sh=None):
        """Returns a dictionary with an array of vol. mix. ratios for each species. 

        Parameters
        ----------
            sh: set or list
                shape of the array wanted if all the vmr are floats.
                If some are already arrays, check whether the shape is the correct one. 

        Returns
        -------
            vmr_array: dict
                A dictionary with the an array of vmr per species.
            cst_array: boolean
                Is True if all the values in the arrays are constant.
        """
        vmr_array=dict()
        cst_array=True
        for mol,vmr in self.composition.items():
            if isinstance(vmr,(float,int)):
                vmr_array[mol]=np.ones(sh)*vmr
            else:
                cst_array=False
                vmr_array[mol]=np.array(vmr) # this np.array could probably go because vmr should be an array at this stage. 
                if not np.array_equal(vmr_array[mol].shape, sh):
                    print('molecule:',mol)
                    print('requested shape:',sh,', molecule shape:',vmr_array[mol].shape)
                    raise RuntimeError('Wrong shape in get_vmr_array')
        return vmr_array, cst_array

    def get_q_array(self, sh=None):
        """Returns a dictionary with an array of specific concentrations for each species. 

        Parameters
        ----------
            sh: set or list
                shape of the array wanted if all the vmr are floats.
                If some are already arrays, check whether the shape is the correct one. 

        Returns
        -------
            q_array: dict
                A dictionary with the an array of specific concentration per species.
            cst_array: boolean
                Is True if all the values in the arrays are constant.
        """
        q_array=dict()
        cst_array=True
        Mg = self.molar_mass()
        for mol,vmr in self.composition.items():
            if mol!='inactive_gas':
                mmol=Molar_mass().fetch(mol)
            else:
                mmol=0.
            if isinstance(vmr,(float,int)):
                q_array[mol]=np.ones(sh)*vmr*mmol/Mg
            else:
                cst_array=False
                q_array[mol]=np.array(vmr)*mmol/Mg # this np.array could probably go because vmr should be an array at this stage. 
                if not np.array_equal(q_array[mol].shape, sh):
                    print('molecule:',mol)
                    print('requested shape:',sh,', molecule shape:',q_array[mol].shape)
                    raise RuntimeError('Wrong shape in get_vmr_array')
        return q_array, cst_array


    def get_vmr_array_basic_molecules(self, sh=None):
        """Returns a dictionary with an array of vol. mix. ratios for each and decomposes
        gas mixtures into basic molecules. 

        Parameters
        ----------

        Returns
        -------
            vmr_array_basic_mol: dictionary of arrays
                A dictionary with the an array of vmr per species.
                Each species is considered a basic molecule.
            vmr_array: dictionary of arrays
                A dictionary with the an array of vmr per species.
                Species in this dictionary can be a gas mixture. 
            cst_array: boolean
                Is True if all the values in the arrays are constant.
        """
        vmr_array_basic_mol=dict()
        cst_array=True
        vmr_array, tmp_cst_array = self.get_vmr_array(sh=sh)
        if not tmp_cst_array: cst_array=False
        for mol,vmr in vmr_array.items():
            if mol in Known_composite_species().keys():
                tmp, _, tmp_cst_array= \
                    Known_composite_species()[mol].get_vmr_array_basic_molecules(sh=sh)
                if not tmp_cst_array: cst_array=False
                for mol2,vmr2 in tmp.items():
                    if mol2 in vmr_array_basic_mol.keys():
                        vmr_array_basic_mol[mol2]+=vmr*vmr2
                    else:
                        vmr_array_basic_mol[mol2]=vmr*vmr2
            else:
                if mol in vmr_array_basic_mol.keys():
                    vmr_array_basic_mol[mol]+=vmr
                else:
                    vmr_array_basic_mol[mol]=vmr
        return vmr_array_basic_mol, vmr_array, cst_array


    def set_k_database(self, k_database=None):
        """Change the radiative database attached to the current instance of Gas_mix

        Parameters
        ----------
            k_database: :class:`~exo_k.kdatabase.Kdatabase` object
                New Kdatabase to use.
        """
        self.kdatabase=k_database
        if self.kdatabase is None:
            self.Ng=None
        else:
            self.Ng=self.kdatabase.Ng
            if self.kdatabase.p_unit != 'Pa' or rm_molec(self.kdatabase.kdata_unit) != 'm^2':
                print("""You're being Bad!!! You are trying *NOT* to use MKS units!!!
                You can convert to mks using convert_to_mks on your Kdatabase.
                More generally, you can specify exo_k.Settings().set_mks(True) 
                to set MKS system as default for all newly loaded data,
                but beware that this global setting is overridden by local options
                specified during the loading.
                You will have to reload all your data though.
                (A good thing it does not take so long). """)
                raise RuntimeError("Bad units in the Kdatabase used with Gas_mix.")
            if (not self.kdatabase.consolidated_p_unit) \
                or (not self.kdatabase.consolidated_kdata_unit):
                raise RuntimeError( \
                    """All tables in the database should have the same units to proceed.
                    You should probably use convert_to_mks().""")

    def set_cia_database(self, cia_database=None):
        """Changes the CIA database attached to the current instance of Gas_mix

        Parameters
        ----------
            cia_database: :class:`~exo_k.ciadatabase.CIAdatabase` object
                New CIAdatabase to use.
        """
        self.cia_database=cia_database
        if self.cia_database is not None:
            if self.cia_database.abs_coeff_unit != 'm^5':
                print("""You're being Bad!!! You are trying *NOT* to use MKS units!!!
                You can convert to mks using convert_to_mks on your CIAdatabase.
                More generally, you can specify exo_k.Settings().set_mks(True) 
                to set MKS system as default for all newly loaded data,
                but beware that this global setting is overridden by local options
                specified during the loading.
                You will have to reload all your data though.
                (A good thing it does not take so long). """)

                raise RuntimeError("Bad units in the CIAdatabase used with Gas_mix.")

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

    def set_spectral_range(self, wn_range=None, wl_range=None):
        """Sets the default spectral range in which computations will be done by specifying
        either the wavenumber or the wavelength range.

        Parameters
        ----------
            wn_range: list or array of size 2
                Minimum and maximum wavenumber (in cm^-1).
            wl_range: list or array of size 2
                Minimum and maximum wavelength (in micron)
        """
        self._wn_range=self._compute_spectral_range(wn_range=wn_range, wl_range=wl_range)

    def _compute_wn_range_indices(self, wn_range=None):
        """Compute the starting and ending indices to be used for current wn_range
        """
        if wn_range is None:
            local_wn_range=self._wn_range
        else:
            local_wn_range=wn_range
        if local_wn_range is None:
            self.iw_min=0
            self.iw_max=self.kdatabase.Nw
        else:
            self.iw_min, self.iw_max = np.where((self.kdatabase.wnedges > local_wn_range[0]) \
                & (self.kdatabase.wnedges <= local_wn_range[1]))[0][[0,-1]]
            # to be consistent with interpolate_kdata

    def cross_section(self, composition=None, logp_array=None, t_array=None,
            wl_range=None, wn_range=None, rayleigh=True,
            write=0, random_overlap=False, logp_interp=True, use_basic_molecules=False, **kwargs):
        """Computes the cross section (m^2/total number of molecule) for the mix
        at each of the logPT points as a function of wavenumber (and possibly g point).

        Parameters
        ----------
            wl_range: array or list of two values, optional
                Wavelength range to cover.
            wn_range: array or list of two values, optional
                Wavenumber range to cover.
            rayleigh: boolean, optional
                Whether to compute rayleigh scattering.
            random_overlap: boolean, optional
                Whether Ktable opacities are added linearly (False),
                or using random overlap method (True).

        Returns
        -------
            kdata_array: array
                Cross section array of shape (layer number, Nw (, Ng if corrk)).

        After every computation, the following variables are updated to account for any possible
        change in spectral range:

          * self.Nw, Number of wavenumber bins
          * self.wns, Wavenumber array
          * self.wnedges, Wavenumber of the edges of the bins
        """
        if self.kdatabase is None: raise RuntimeError("""k_database not provided. 
        Use the k_database keyword during initialization or use the set_k_database method.""")
        if not self.kdatabase.consolidated_wn_grid: raise RuntimeError("""
            All tables in the database should have the same wavenumber grid to proceed.
            You should probably use bin_down().""")
        if self.cia_database is not None and \
            (not np.array_equal(self.cia_database.wns,self.kdatabase.wns)):
            raise RuntimeError("""CIAdatabase not sampled on the right wavenumber grid.
            You should probably run something like CIAdatabase.sample(Kdatabase.wns).""")
        
        self.set_logPT(logp_array=logp_array, t_array=t_array)
        if self.Narray is None:
            raise RuntimeError('You must prescribe logP (in Pa) and T arrays first.')

        if composition is not None:
            self.set_composition(composition)
        if self.need_to_recompute_vmr_array:
            self.vmr_arr_base_mol, self.vmr_arr, self.cst_array = \
                self.get_vmr_array_basic_molecules((self.Narray,))
            self.need_to_recompute_vmr_array = False
            if write>6 :
                print('recomputed vmr array:', self.vmr_arr_base_mol, self.vmr_arr) 
        if use_basic_molecules:
            molecs=self.vmr_arr_base_mol.keys()
        else:
            molecs=self.vmr_arr.keys()
        if write>6 :
            print('molecs:', molecs) 
        mol_to_be_done=list(set(molecs).intersection(self.kdatabase.molecules))
        if not mol_to_be_done:
            if 'total' in self.kdatabase.molecules:
                #special case where you have one ktable that describes the whole gas.
                mol_to_be_done=['total']
                use_basic_molecules=False
                self.vmr_arr['total']=np.ones((self.Narray))
            else:
                raise RuntimeError("""The k_database you provided 
                    should contain at least one molecule in your atm,
                    or a ktable with the 'total' key.""")
        if all(elem in self.kdatabase.molecules for elem in molecs):
            if write>3 : print("""I have all the molecules present in the atmosphere
              in ktables provided:""")
        else:
            if write>3 : print("""Some missing molecules in my database,
             I ll compute opacites with the available ones:""")
        if write>3 : print(mol_to_be_done)

        # does what needs to be done to reduce the spectral range
        local_wn_range=self._compute_spectral_range(wl_range=wl_range, wn_range=wn_range)
        self._compute_wn_range_indices(wn_range=local_wn_range)
        self.wnedges=np.copy(self.kdatabase.wnedges[self.iw_min:self.iw_max+1])
        self.dwnedges=np.diff(self.wnedges)
        self.wns=np.copy(self.kdatabase.wns[self.iw_min:self.iw_max]) # 2021: are these copy needed?
        self.Nw=self.wns.size

        first_mol=True
        if use_basic_molecules:
            vmr_arr_to_use = self.vmr_arr_base_mol
        else:
            vmr_arr_to_use = self.vmr_arr
        for mol in mol_to_be_done:
            tmp_kdata=self.kdatabase[mol].interpolate_kdata(logp_array=self.logp_array,
                t_array=self.t_array, x_array=vmr_arr_to_use[mol], wngrid_limit=local_wn_range,
                logp_interp=logp_interp,
                **kwargs)
            if first_mol:
                kdata_array=tmp_kdata
                first_mol=False
            else:
                if random_overlap and (self.kdatabase.Ng is not None):
                    kdata_array=RandOverlap_2_kdata_prof(self.Narray,
                        self.Nw, self.kdatabase.Ng, 
                        kdata_array,tmp_kdata, self.kdatabase.weights,
                        self.kdatabase.ggrid)
                else:
                    kdata_array+=tmp_kdata

        if rayleigh or (self.cia_database is not None):
            # continua always use base molecules
            cont_sig=np.zeros((self.Narray,self.Nw))
            if rayleigh:
                if self.cst_array:
                    cont_sig+=Rayleigh().sigma(self.wns, self.vmr_arr_base_mol)
                else:
                    cont_sig+=Rayleigh().sigma_array(self.wns, self.vmr_arr_base_mol)
                self.kdata_scat=np.copy(cont_sig)
            if self.cia_database is not None:
                cont_sig+=self.cia_database.cia_cross_section(self.logp_array,
                    self.t_array, self.vmr_arr_base_mol, wngrid_limit=local_wn_range, Nw=self.Nw)
            if self.kdatabase.Ng is None:
                kdata_array+=cont_sig
            else:
                kdata_array+=cont_sig[:,:,None]

        return kdata_array
    
    def __getitem__(self, molecule):
        """Overrides getitem
        """
        return self.composition[molecule]

    def __setitem__(self, molecule, vmr):
        """Overrides setitem and recomputes the background gas mixing ratio if needed
        """
        if isinstance(vmr, str):
            self.bg_gas=molecule
            self.composition['inactive_gas']=0.
        else:
            if isinstance(vmr, list):
                self.composition[molecule]=np.array(vmr)
            else:
                self.composition[molecule]=vmr
            if molecule==self.bg_gas:
                self.bg_gas='inactive_gas'
        self.get_background_vmr()
    
    def items(self):
        """Emulates dict.items() method
        """
        return self.composition.items()

    def values(self):
        """Emulates dict.values() method
        """
        return self.composition.values()

    def keys(self):
        """Emulates dict.keys() method
        """
        return self.composition.keys()

    def __repr__(self):
        """Method to output
        """
        output='Volume mixing ratios in Gas_mix: \n'
        for mol, vmr in self.composition.items():
            output+=mol+'->'+str(vmr)+'\n'
        return output

    def copy(self):
        """Deep copy of the dict and arrays. 
        The databases are not deep copied. 
        """
        res=Gas_mix(self.composition, bg_gas=self.bg_gas,
            k_database=self.kdatabase, cia_database=self.cia_database)
        res.logp_array=np.copy(self.logp_array)
        res.t_array=np.copy(self.t_array)
        res.Narray=self.Narray
        res._wn_range=np.copy(self._wn_range)
        res.iw_min=self.iw_min
        res.iw_max=self.iw_max

    def mix_with(self, other_gas, vmr_other_gas):
        """Mix with other Gas_mix.
        """
        raise NotImplementedError()

    def __mul__(self, vmr):    
        """Defines multiplication
        """
        raise NotImplementedError()
        #composition=dict()
        #for mol,mol_vmr in self.composition.items():
        #    composition[mol]=vmr*mol_vmr
        #res=Gas_mix(composition, bg_gas=self.bg_gas)
        #return res

    __rmul__ = __mul__


class Known_composite_species(Singleton):
    """A class to store composite gas mixes that can be used in other gas mixes
    This class can also store the molar mass of custom gases with arbitrary names
    (for example: My_gas, earth_background).
    """
    
    def init(self, *args, **kwds):
        """Initializes empty dictionary of custom molecular masses.
        """
        self._known_composite_species={}

    def add_species(self, species_dict):
        """Add one or several composite species to the database. 

        Parameters
        ----------
            species_dict: dict
                Keys are gases names (they do not have to be real molecules).
                Values are a gas_mix object. 
        """
        self._known_composite_species.update(species_dict)
        for species, gas_mix in species_dict.items():
            Molar_mass().add_species({species: gas_mix.molar_mass()})

    def items(self):
        """Emulates dict.items() method
        """
        return self._known_composite_species.items()

    def values(self):
        """Emulates dict.values() method
        """
        return self._known_composite_species.values()

    def keys(self):
        """Emulates dict.keys() method
        """
        return self._known_composite_species.keys()

    def __getitem__(self, molecule):
        """Overrides getitem
        """
        return self._known_composite_species[molecule]

    def __setitem__(self, molecule, mix):
        """Overrides setitem
        """
        self._known_composite_species[molecule]=mix
