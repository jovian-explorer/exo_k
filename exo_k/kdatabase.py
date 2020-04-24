# -*- coding: utf-8 -*-
"""
Created on May 2 2019

@author: jeremy leconte
Doc can be created with "pydoc -w corrk_lib"
"""
import numpy as np
from .ktable import Ktable
from .xtable import Xtable
from .chemistry import gas_mix
from .settings import Settings

class Kdatabase(object):
    """This object contains mainly a dictionary of individual Ktable objects for each molecule. 
    In addition, the informations about the P,T,Wn,g grids
    are reloaded as atributes of the Kdatabase object.
    """

    def __init__(self, molecules, *str_filters, search_path=None,
        remove_zeros=True, **kwargs):
        """Initializes k coeff tables and supporting data from a list of molecules
        Parameters:
            molecules: list or dict
                If a list of molecules is provided,
                    the file starting with the molecule name and containing all the str_filters
                    are searched in the Settings()._search_path
                If a dictionary is provided, the keys are the molecules to load,
                    and the values are the path to the corresponding file.
                    If a None value is given, the str_filters will be used as above.
        Options: 
            See the options of Ktable.__init__()
            search_path: str
                If search_path is provided, it locally overrides the global _search_path settings
                and only files in search_path are returned.            
        """
        self.ktables={}
        self._settings=Settings()
        self.consolidated_wn_grid=True
        self.consolidated_PT_grid=True

        delim=self._settings._delimiter

        self.molecules=None
        if molecules is None:
            return
        if isinstance(molecules,list):
            for mol in molecules:
                try:
                    tmp_ktable=Ktable(*([mol+delim]+list(str_filters)), mol=mol,
                        remove_zeros=remove_zeros, search_path=search_path, **kwargs)
                except:
                    tmp_ktable=Xtable(*([mol+delim]+list(str_filters)), mol=mol,
                        remove_zeros=remove_zeros, search_path=search_path, **kwargs)
                self.add_ktables(tmp_ktable)
        else:
            for mol,filename in molecules.items():
                try:
                    # below, we still provide  *([mol+delim]+list(str_filters)) 
                    # as an input in case filename is None
                    tmp_ktable=Ktable(*([mol+delim]+list(str_filters)), filename=filename, mol=mol,
                        remove_zeros=remove_zeros, search_path=search_path, **kwargs)
                except:
                    tmp_ktable=Xtable(*([mol+delim]+list(str_filters)), filename=filename, mol=mol,
                        remove_zeros=remove_zeros, search_path=search_path, **kwargs)
                self.add_ktables(tmp_ktable)
        

    def add_ktables(self, *ktables):
        """Adds as many Ktables to the database as you want.
        """
        for tmp_ktable in ktables:
            if self.molecules is None:
                self.ktables[tmp_ktable.mol]=tmp_ktable
                self.pgrid=tmp_ktable.pgrid
                self.logpgrid=tmp_ktable.logpgrid
                self.tgrid=tmp_ktable.tgrid
                self.wns=tmp_ktable.wns
                self.wnedges=tmp_ktable.wnedges
                self.Np=tmp_ktable.Np
                self.Nt=tmp_ktable.Nt
                self.Nw=tmp_ktable.Nw
                self.Ng=tmp_ktable.Ng
                if tmp_ktable.Ng is not None:
                    self.Ng=tmp_ktable.Ng
                    self.weights=tmp_ktable.weights
                    self.ggrid=tmp_ktable.ggrid
                    self.gedges=tmp_ktable.gedges
            else:
                if (self.Ng is None)^(tmp_ktable.Ng is None): raise RuntimeError( \
                    'All elements in a database must have the same type (Ktable or Xtable).')
                if (self.Ng is not None) and (not np.array_equal(tmp_ktable.ggrid,self.ggrid)):
                    raise RuntimeError('All Ktables in a database must have the same g grid.')
                self.ktables[tmp_ktable.mol]=tmp_ktable
                if not np.array_equal(tmp_ktable.wns,self.wns):
                    self.consolidated_wn_grid=False
                    print("""Carefull, not all tables have the same wevelength grid.
                        You'll probably need to run bin_down()""")
                    self.wns    = None
                    self.wnedges= None
                    self.Nw     = None
                if not (np.array_equal(tmp_ktable.pgrid,self.pgrid) \
                    and np.array_equal(tmp_ktable.tgrid,self.tgrid)) :
                    self.consolidated_PT_grid=False
                    self.pgrid   = None
                    self.logpgrid= None
                    self.tgrid   = None
                    self.Np      = None
                    self.Nt      = None
                    print("""Carefull, not all tables have the same PT grid.
                        You'll probably need to run remap_logPT()""")
            self.molecules=list(self.ktables.keys())

    def __getitem__(self,molecule):
        """Overrides getitem so as to access directly a Ktable with Kdatabase['mol']
        """
        if molecule not in self.ktables.keys():
            raise KeyError('The requested molecule is not available.')
        return self.ktables[molecule]

    @property
    def shape(self):
        """Returns the shape of self.kdata
        """
        return np.array([self.Np,self.Nt,self.Nw,self.Ng])

    @property
    def wls(self):
        """Returns the wavelength array for the bin centers
        """
        return 10000./self.wns

    @property
    def wledges(self):
        """Returns the wavelength array for the bin edges
        """
        return 10000./self.wnedges

    def remap_logPT(self,logp_array=None,t_array=None):
        """Applies the bin_down method to all the tables in the database. This can be used 
        to put all the tables onthe same PT grid.
        See data_table.remap_logPT() for details.
        """
        for mol in self.molecules:
            self.ktables[mol].remap_logPT(logp_array=logp_array,t_array=t_array)
        self.logpgrid=np.array(logp_array)
        self.pgrid   =10**self.logpgrid
        self.tgrid   =np.array(t_array)
        self.Np      =self.logpgrid.size
        self.Nt      =self.tgrid.size
        self.consolidated_PT_grid=True

    def bin_down(self,wnedges=None,**kwargs):
        """Applies the bin_down method to all the tables in the database. This can be used 
        to put all the tables onthe same wavenumber grid.
        See Ktable.bin_down() or Xtable.bin_down() for details.
        """
        first=True
        for mol in self.molecules:
            self.ktables[mol].bin_down(wnedges=wnedges,**kwargs)
            if first:
                self.wns=self.ktables[mol].wns
                self.wnedges=self.ktables[mol].wnedges
                self.Nw=self.ktables[mol].Nw
                self.ggrid=self.ktables[mol].ggrid
                self.weights=self.ktables[mol].weights
                self.Ng=self.ktables[mol].Ng
                self.consolidated_wn_grid=True

    def sample(self,wngrid,**kwargs):
        """Applies the bin_down method to all the tables in the database. This can be used 
        to put all the tables onthe same wavenumber grid.
        See Ktable.bin_down() or Xtable.bin_down() for details.
        """
        first=True
        if self.Ng is not None: raise RuntimeError('sample is only available for Xtable objects.')
        for mol in self.molecules:
            self.ktables[mol].sample(wngrid,**kwargs)
            if first:
                self.wns=self.ktables[mol].wns
                self.wnedges=self.ktables[mol].wnedges
                self.Nw=self.ktables[mol].Nw
                self.consolidated_wn_grid=True


    def create_mix_ktable(self, composition, inactive_species=[]):
        """creates the kdata table for a mix of molecules
        Parameters:
            composition: dict
                Keys are the molecule names (they must match the names in the database).
                Values are either numbers or arrays of volume mixing ratios
                with shape (pgrid.size,tgrid.size).
                This composition will instantiate a gas_mix object.
                In particular, if a value is 'background', this gas will
                be used to fill up to sum(vmr)=1 (See chemistry.gas_type for details).
            
                For each (P,T) point, the sum of all the mixing ratios
                should be lower or equal to 1.
                If it is lower, it is assumed that the rest of the gas is transparent.

            inactive_species: list, optional
                List the gases that are in composition but for which we do not want the 
                opacity to be accounted for. 

        Returns:
            res: Ktable object
                A new ktable for the mix
        """
        if self.Ng is None: raise RuntimeError("""
           Create_mix_ktable cannot work with a database of cross sections.
           Please load correlated-k data.""")
        if not self.consolidated_wn_grid: raise RuntimeError("""
           All tables in the database should have the same wavenumber grid to proceed.
           You should probably use bin_down().""")
        if not self.consolidated_PT_grid: raise RuntimeError("""
           All tables in the database should have the same PT grid to proceed.
           You should probably use remap_logPT().""")
        mol_to_be_done=set(composition.keys())
        mol_to_be_done=mol_to_be_done-set(inactive_species)
        if not mol_to_be_done:
            print("""You are creating a mix without any active gas:
                This will be awfully transparent""")
            res=self[self.molecules[0]].copy(cp_kdata=False)
            res.kdata=np.zeros(res.shape)
            return res
        gas_mixture=gas_mix(composition)
        if all(elem in self.molecules for elem in mol_to_be_done):
            print('I have all the requested molecules in my database')
            print(mol_to_be_done)
        else:
            print('Do not have all the molecules in my database')
            mol_to_be_done=mol_to_be_done.intersection(set(self.molecules))
            print('Molecules to be treated: ',mol_to_be_done)
        first_mol=True
        for mol in mol_to_be_done:
            if first_mol:
                res=self.ktables[mol].copy(cp_kdata=True)
                try:
                    res.kdata=res.VolMixRatioNormalize(gas_mixture[mol])
                except TypeError:
                    print('gave bad mixing ratio format to VolMixRatioNormalize')
                    raise TypeError('bad mixing ratio type')            
                if len(mol_to_be_done)==1:
                    print('only 1 molecule:',mol)
                    return res
                first_mol=False
            else:
              print('treating molecule ',mol)
              res.kdata=res.RandOverlap(self.ktables[mol],None,gas_mixture[mol])
              # no need to re normalize with respect to 
              # the abundances of the molecules already done.
        return res  

    def create_mix_ktable5d(self, bg_comp={}, vgas_comp={}, x_array=None,
            bg_inac_species=[], vgas_inac_species=[], **kwargs):
        """Creates a Ktable5d for a mix of molecules with a variable gas.
        In essence, the regular create_mix_ktable is called to create
        two mixes:
            - the background mix specified by bg_comp={}, bg_inac_species=[]
            - the variable gas specified by vgas_comp={}, vgas_inac_species=[]
        See create_mix_ktable for details.

        These two gases are then mixed together for an array of vmr x_array where
        var_gas has a vmr of x_array and the background gas has a vmr of 1.-x_array

        Returns:
            res: Ktable5D object
                A new ktable for the mix
        """
        if x_array is None:
            raise RuntimeError('x_array is None: pas bien!!!')
        background_mix=self.create_mix_ktable(bg_comp,inactive_species=bg_inac_species)
        var_gas_mix=self.create_mix_ktable(vgas_comp,inactive_species=vgas_inac_species)
        ktab5d=var_gas_mix.copy(ktab5d=True)
        ktab5d.xgrid=np.array(x_array)
        ktab5d.Nx=ktab5d.xgrid.size
        print(ktab5d.shape)
        new_kdata=np.zeros(ktab5d.shape)
        for iX, vmr in enumerate(ktab5d.xgrid):
            new_kdata[:,:,iX,:,:]=var_gas_mix.RandOverlap(background_mix,vmr,1.-vmr, **kwargs)
        ktab5d.set_kdata(new_kdata)
        return ktab5d
