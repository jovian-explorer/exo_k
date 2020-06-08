# -*- coding: utf-8 -*-
"""
@author: jeremy leconte
"""
import numpy as np
from scipy.interpolate import RectBivariateSpline

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

