# -*- coding: utf-8 -*-
"""
@author: jeremy leconte
"""
import numpy as np
from scipy.interpolate import RectBivariateSpline
from .util.molar_mass import Molar_mass


class EquChemTable(object):
    """Class to load and interpolate chemistry data.
    """
    
    def __init__(self, filename=None, remove_zeros=False):
        """Initializes chemical composition tables from various init files

        Parameters
        ----------
            filename : str
                Name of the input file
        """
        self.filename=filename
        self.tgrid=None
        self.pgrid=None
        self.logpgrid=None
        self.p_unit=None
        self.Nt=None
        self.Np=None
        self.molecules=None
        self.vol_mix_ratio={}
        self.logx_interp={}
        if self.filename.lower().endswith('.in'):
            self.read_composition_in(filename=self.filename, remove_zeros=remove_zeros)
        elif self.filename.lower().endswith('.dat'):
            self.read_composition_dat(filename=self.filename, remove_zeros=remove_zeros)
        self.setup_interpolation()

    def read_composition_in(self, filename=None, remove_zeros=False, skiprows=7):
        """Initializes chemical composition tables from composition.in files

        Parameters
        ----------
            filename : str
                Name of the input file
            skiprows : int, optional
                Number of lines to skip in the file
        """
        data = np.loadtxt(filename,skiprows=skiprows,unpack=True)
        self.Nt=np.where(data[2,1:]==data[2,0])[0][0]+1
        self.tgrid=data[2,0:self.Nt]
        #print('T grid:',self.tgrid)
        self.Np=data[0].size//self.Nt
        self.pgrid=data[1,::self.Nt]/0.986923267*1.e5
        self.p_unit='Pa'
        #print('P grid:',self.pgrid)
        #print('Np, Nt:',self.Np, self.Nt)
        self.logpgrid=np.log10(self.pgrid)
        with open(filename, 'r') as file:
            last_line=file.readline()
            while True:
                new_line=file.readline()
                if new_line.split()[0]=='z':
                    break
                last_line=new_line
        molecules=last_line.split()
        #print(molecules)
        for ii, mol in enumerate(molecules):
            self.vol_mix_ratio[mol]=data[ii+3].reshape((self.Np,self.Nt))
        if remove_zeros:
            self.remove_zeros()

    def read_composition_dat(self, filename=None, remove_zeros=False, skiprows=1):
        """Initializes chemical composition tables from composition.in files

        Parameters
        ----------
            filename : str
                Name of the input file
            skiprows : int, optional
                Number of lines to skip in the file
        """
        with open(filename, 'r') as file:
            molecules=file.readline().split()[2:]
        molecules[0]=molecules[0].replace('[mol/mol]','')
        data = np.loadtxt(filename,skiprows=skiprows,unpack=True)
        self.pgrid=np.sort(np.array(list(set(data[0]))))
        self.logpgrid=np.log10(self.pgrid)
        self.tgrid=np.sort(np.array(list(set(data[1]))))
        self.Np=self.pgrid.size
        self.Nt=self.tgrid.size
        for ii, mol in enumerate(molecules):
            self.vol_mix_ratio[mol]=data[ii+2].reshape((self.Np,self.Nt))[:,-1::-1]
        if remove_zeros:
            self.remove_zeros()

    def remove_zeros(self,deltalog_min_value=30.):
        """Finds zeros in the chem data and set them to (10.^-deltalog_min_value)
        times the minimum positive value
        in the table. This is to be able to work in logspace. 
        """
        for mol,x_ratio in self.vol_mix_ratio.items():
            mask = np.zeros((self.Np,self.Nt),dtype=bool)
            mask[np.nonzero(x_ratio)] = True
            minvalue=np.amin(x_ratio[mask])
            self.vol_mix_ratio[mol][~mask]=minvalue/(10.**deltalog_min_value)

    def setup_interpolation(self):
        """Creates interpolating functions to be called later on.
        """
        for mol,x_ratio in self.vol_mix_ratio.items():
            self.logx_interp[mol]=RectBivariateSpline(self.logpgrid,self.tgrid,np.log(x_ratio))

    def __getitem__(self, molecule):
        """Overrides getitem so that EquChemTable['molecule'] directly accesses 
        the database for that molecule.
        """
        if molecule not in self.vol_mix_ratio:
            raise KeyError('The requested molecule is not available.')
        return self.vol_mix_ratio[molecule]

    def vmr(self, logP, T, mol):
        """Interpolates a single molecule on a logP-T value
        """
        return np.exp(self.logx_interp[mol](logP,T,grid=False)[0])

    def interpolate_vmr(self, logp_array=None,t_array=None, mols=None,grid=False):
        """Interpolates all molecules in mols on a logP-T grid
        """
        res={}
        for mol in mols:
            res[mol]=np.exp(self.logx_interp[mol](logp_array, t_array, grid=grid))
        return res

#    def vmr_logPT(self, logPgrid, tgrid, mols, is_sorted=False):
#        """Interpolates all molecules in mols on a logP-T grid
#        """
#        res={}
#        #lP,T=np.meshgrid(logPgrid,tgrid)
#        for mol in mols:
#            #res[mol]=10**self.logx_interp[mol](logPgrid.reshape((self.Np,1)),tgrid)
#            res[mol]=np.exp(self.logx_interp[mol](logPgrid, tgrid, \
#                assume_sorted=is_sorted).transpose())
#        return res

class gas_mix(object):
    """Dict-like class to handle gas composition (with background gas_mix) and molar mass.
    """

    def __init__(self, composition, bg_gas=None):
        """Instantiates a gas_mix object and computes the vmr of the 'background' gas.
        """
        self.composition=composition.copy()
        self.bg_gas=bg_gas
        if (True in [isinstance(val,str) for val in self.composition.values()]) \
                or (self.bg_gas is not None):
            self.composition['inactive_gas']=0.
        else:
            self.composition['inactive_gas']='background'
            self.bg_gas='inactive_gas'
        self.get_background_vmr()

    def get_background_vmr(self):
        """Computes the volume mixing ratio of the background gas in a mix
        Uses the fact that self.composition is a dictionary
        with the name of the molecules as keys and the vmr as values.
        vmr are either float or arrays.
        The background gas must have a vmr='background'
        
        Returns
        -------
            float or array:
                Vol. Mix. Ratio of the background gas
        """
        other_vmr=0.
        if self.bg_gas is None:
            for mol,vmr in self.composition.items():
                if isinstance(vmr,str):
                    self.bg_gas=mol
                    continue
                other_vmr+=vmr
        else:
            for mol,vmr in self.composition.items():
                if mol==self.bg_gas:
                    continue
                other_vmr+=vmr
        if np.amax(other_vmr)>1.:
            print("""Carefull: the sum of the vmr of your gas components is > 1.
            If there is a background gas, its vmr will become negative.
            I hope you know what you are doing.""")
        self.composition[self.bg_gas]=1.-other_vmr

    def molar_mass(self):
        """Computes the molar mass of a mix of gases

        Parameters
        ----------
            composition: dict
                Dictionary with the name of the molecules as keys
                and the vmr as values. vmr are either float or arrays.
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

    def get_vmr_array(self, sh=None):
        """Returns a dictionary with an array of vol. mix. ratios for each species. 

        Parameters
        ----------
            sh: set or list
                shape of the array wanted if all the vmr are floats.
                If some are already arrays, check whether the shape is the correct one. 
        """
        res=dict()
        cst_array=True
        for mol,vmr in self.composition.items():
            if isinstance(vmr,(float,int)):
                res[mol]=np.ones(sh)*vmr
            else:
                cst_array=False
                res[mol]=np.array(vmr)
                if not np.array_equal(res[mol].shape, sh):
                    raise RuntimeError('Wrong shape in get_gas_array')
        return res, cst_array
    
    def mix_with(self, other_gas, vmr_other_gas):
        """Mix with other gas_mix.
        """

    def __mul__(self, vmr):    
        """Defines multiplication
        """
        composition=dict()
        for mol,mol_vmr in self.composition.items():
            composition[mol]=vmr*mol_vmr
        res=gas_mix(composition, bg_gas=self.bg_gas)
        return res

    __rmul__ = __mul__

    def __getitem__(self, molecule):
        """Overrides getitem
        """
        return self.composition[molecule]

    def __setitem__(self, molecule, vmr):
        """Overrides setitem and recomputes the background gas mixing ratio if needed
        """
        if isinstance(vmr,str):
            self.bg_gas=molecule
            self.composition['inactive_gas']=0.
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

