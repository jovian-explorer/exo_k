# -*- coding: utf-8 -*-
"""
@author: jeremy leconte
"""
import numpy as np
import h5py

class Kspectrum(object):
    """A class defining a Kspectrum object.
    """

    def __init__(self, filename, **kwargs):
        
        self.pressure=None
        self.p_unit=None
        self.temperature=None
        self.kdata=None
        self.kdata_unit=None
        self.wns=None
        self.is_xsec=None

        if filename.lower().endswith(('.hdf5', '.h5')):
            self.read_hdf5(filename, **kwargs)
        else:
            self.read_ascii(filename, **kwargs)


    def read_ascii(self, filename, skiprows=0):
        """Read native kspectrum format

        Parameters
        ----------
            filename: str
                Initial kspectrum filename.
            skiprows: int, optional
                Number of header lines to skip. For the latest Kspectrum format,
                the header is skipped automatically. 
        """
        with open(filename, 'r') as file:
            tmp = file.readline().split()
            if tmp[0]=='Pressure': #new kspectrum format
                self.is_xsec=True
                self.pressure = float(file.readline().split()[0])
                self.pressure = self.pressure/9.8692e-6 # conversion from atm to Pa
                tmp = file.readline().split()
                self.temperature = float(file.readline().split()[0])
                tmp = file.readline().split()
                nb_mol = int(file.readline().split()[0])
                skiprows=skiprows+9+5*nb_mol
            else:
                self.is_xsec=False
        raw=np.genfromtxt(filename,skip_header=skiprows,usecols=(0,1),names=('wns','sigma')) 
        self.kdata=raw['sigma']
        self.wns=raw['wns']
        if self.is_xsec:
            self.p_unit='Pa'
            self.kdata_unit='m^2/molec'
            self.kdata=self.kdata*1.e-4 # conversion from cm^2 in file to m^2
        else:
            self.kdata_unit='m^-1'

    def write_hdf5(self, filename):
        """Writes kspectrum file to hdf5
        """
        f = h5py.File(filename, 'w')
        f.attrs["p"] = self.pressure
        if self.is_xsec:
            f.attrs["data_type"] = 'xsec'
        else:
            f.attrs["data_type"] = 'k'
        f.attrs["p_unit"] = 'Pa'
        f.attrs["t"] = self.temperature
        f.create_dataset("wns", data=self.wns,compression="gzip")
        f.create_dataset("kdata", data=self.kdata,compression="gzip")
        f["kdata"].attrs["units"]=self.kdata_unit
        f.close()    

    def read_hdf5(self, filename):
        """Reads kspectrum file from hdf5
        """
        f = h5py.File(filename, 'r')
        data_type=f.attrs['data_type']
        self.is_xsec=(data_type=='xsec')
        self.wns=f['wns'][...]
        self.kdata=f['kdata'][...]
        if self.is_xsec:
            self.pressure=f.attrs['p']
            self.p_unit=f.attrs['p_unit']
            self.temperature=f.attrs['t']
            self.kdata_unit=f['kdata'].attrs['units']

        f.close()


