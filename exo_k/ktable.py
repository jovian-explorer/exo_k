# -*- coding: utf-8 -*-
"""
@author: jeremy leconte
Doc can be created with "pydoc -w corrk_lib"
"""
import os
from math import log10
import time
import pickle
import h5py
import numpy as np
import astropy.units as u
from .data_table import Data_table
from .ktable5d import Ktable5d
from .util.interp import rm_molec, interp_ind_weights, rebin_ind_weights, rebin, is_sorted, \
        gauss_legendre, spectrum_to_kdist, kdata_conv_loop, bin_down_corrk_numba
from .hires_spectrum import Hires_spectrum
from .util.filenames import create_fname_grid_Kspectrum_LMDZ, select_kwargs, read_nemesis_binary
from .util.cst import KBOLTZ, nemesis_hitran_id_numbers


class Ktable(Data_table):
    """A class that handles 4D tables of k-coefficients.
    Based on the Data_table class that handles basic operations common to Xtable.
    """
        
    def __init__(self, *filename_filters, filename=None, xtable=None, path=None,
        p_unit='unspecified', file_p_unit='unspecified',
        kdata_unit='unspecified', file_kdata_unit='unspecified',
        remove_zeros=False, search_path=None, mol=None,
        **kwargs):
        """Initializes k coeff table and supporting data from various sources
        (see below by order of precedence)

        Parameters
        ----------
            filename: str, optional
                Relative or absolute name of the file to be loaded. 
            filename_filters: sequence of string
                As many strings as necessary to uniquely define
                a file in the global search path defined in
                :class:`~exo_k.settings.Settings`.
                This path will be searched for a file
                with all the filename_filters in the name.
                The filename_filters can contain '*'.
            xtable: Xtable object
                If no filename nor filename_filters are provided, this xtable object will be used to
                create a ktable. In this case, wavenumber bins must be given
                with the wnedges keyword.
            path: str
                If none of the above is specifed, path can point to
                a directory with a LMDZ type k coeff table.
                In this case, see read_LMDZ for the keywords to specify.

        If none of the parameters above is specified,
        just creates an empty object to be filled later.

        Other Parameters
        ----------------
            p_unit: str, optional
                String identifying the pressure units to convert to (e.g. 'bar', 'Pa', 'mbar', 
                or any pressure unit recognized by the astropy.units library).
                If ='unspecified', no conversion is done.
            file_p_unit: str, optional
                String to specify the current pressure unit if it is unspecified or if 
                you have reasons to believe it is wrong (e.g. you just read a file where
                you know that the pressure grid and the pressure unit do not correspond)
            kdata_unit: str, optional
                String to identify the unit to convert to.
                Accepts 'cm^2', 'm^2'
                or any surface unit recognized by the astropy.units library.
                If ='unspecified', no conversion is done.
                In general, kdata should be kept in 'per number' or 'per volume'
                units (as opposed to 'per mass' units) as composition will
                always be assumed to be a number or volume mixing ratio.
                Opacities per unit mass are not supported yet.
                Note that you do not need to specify the '/molec' or '/molecule' in the unit.
            file_kdata_unit: str, optional
                String to specify the current kdata unit if it is unspecified or if 
                you have reasons to believe it is wrong (e.g. you just read a file where
                you know that the kdata grid and the kdata unit do not correspond)
            remove_zeros: boolean, optional
                If True, the zeros in the kdata table are replaced by
                a value 10 orders of magnitude smaller than the smallest positive value
            mol: str, optional
                The name of the gas or molecule described by the :class:`Ktable`
            search_path: str, optional
                If search_path is provided,
                it locally overrides the global _search_path
                in :class:`~exo_k.settings.Settings`
                and only files in search_path are returned.
        """
        super().__init__()

        if filename is not None:
            self.filename=filename
        elif filename_filters or (mol is not None and (xtable is None and path is None)): 
        # a none empty sequence returns a True in a conditional statement
            self.filename=self._settings.list_files(*filename_filters, molecule = mol,
                only_one=True, search_path=search_path)[0]

        if self.filename is not None:
            if self.filename.lower().endswith('pickle'):
                self.read_pickle(filename=self.filename, mol=mol)
            elif self.filename.lower().endswith(('.hdf5', '.h5')):
                self.read_hdf5(filename=self.filename, mol=mol)
            elif self.filename.lower().endswith('.kta'):
                self.read_nemesis(filename=self.filename, mol=mol)
            else:
                raise NotImplementedError("""Requested format not recognized.
            Should end with .pickle, .hdf5, .h5, or .kta""")
        elif xtable is not None:
            self.xtable_to_ktable(xtable=xtable, **kwargs)
        elif path is not None:
            self.read_LMDZ(path=path, mol=mol, **kwargs)
        else:                  #if there is no input file, just create an empty object 
            self.wnedges=None
            self.weights=None
            self.ggrid=None
            self.gedges=None

        super().finalize_init(p_unit=p_unit, file_p_unit=file_p_unit,
            kdata_unit=kdata_unit, file_kdata_unit=file_kdata_unit,
            remove_zeros=remove_zeros)


    @property
    def shape(self):
        """Returns the shape of self.kdata
        """
        return np.array([self.Np,self.Nt,self.Nw,self.Ng])

    def read_pickle(self, filename=None, mol=None):
        """Initializes k coeff table and supporting data from an Exomol pickle file

        Parameters
        ----------
            filename: str
                Name of the input pickle file
            mol: str, optional
                Force the name of the molecule
        """
        if filename is None: raise RuntimeError("You should provide an input pickle filename")
        pickle_file=open(filename,'rb')
        raw=pickle.load(pickle_file, encoding='latin1')
        pickle_file.close()
        
        if mol is not None: 
            self.mol=mol
        else:
            self.mol=raw['name']
            if self.mol=='H2OP': self.mol='H2O'

        self.pgrid=raw['p']
        self.logpgrid=np.log10(self.pgrid)
        self.tgrid=raw['t']
        self.wns=raw['bin_centers']
        self.wnedges=raw['bin_edges']
        if 'weights' in raw.keys():
            self.weights=raw['weights']
        else:
            raise RuntimeError('No weights keyword. This file is probably a cross section file.')
        self.ggrid=raw['samples']
        self.gedges=np.insert(np.cumsum(self.weights),0,0.)
        self.kdata=raw['kcoeff']
    
        if 'p_unit' in raw.keys():
            self.p_unit=raw['p_unit']
        else:
            self.p_unit='bar'

        if 'kdata_unit' in raw.keys():
            self.kdata_unit=raw['kdata_unit']
        else:
            self.kdata_unit='cm^2/molec'

        if 'wn_unit' in raw.keys(): self.wn_unit=raw['wn_unit']

        self.Np,self.Nt,self.Nw,self.Ng=self.kdata.shape

    def write_pickle(self, filename):
        """Saves data in a pickle format

        Parameters
        ----------
            filename: str
                Name of the file to be created and saved
        """
        fullfilename=filename
        if not filename.lower().endswith('.pickle'): fullfilename=filename+'.pickle'
        pickle_file=open(fullfilename,'wb')
        dictout={'name':self.mol,
                 'p':self.pgrid,
                 'p_unit':self.p_unit,
                 't':self.tgrid,
                 'bin_centers':self.wns,
                 'bin_edges':self.wnedges,
                 'wn_unit':self.wn_unit,
                 'weights':self.weights,
                 'samples':self.ggrid,
                 'kcoeff':self.kdata,
                 'kdata_unit':self.kdata_unit}
        #print(dictout)
        pickle.dump(dictout,pickle_file,protocol=-1)
        pickle_file.close()

    def read_hdf5(self, filename=None, mol=None):
        """Initializes k coeff table and supporting data from an Exomol hdf5 file

        Parameters
        ----------
            file : str
                Name of the input hdf5 file
        """
        if (filename is None or not filename.lower().endswith(('.hdf5', '.h5'))):
            raise RuntimeError("You should provide an input hdf5 file")
        f = h5py.File(filename, 'r')
        if 'mol_name' in f:
            self.mol=f['mol_name'][()]
        elif 'mol_name' in f.attrs:
            self.mol=f.attrs['mol_name']
        else:
            self.mol=os.path.basename(filename).split(self._settings._delimiter)[0]
        if isinstance(self.mol, np.ndarray): self.mol=self.mol[0]
        if mol is not None:
            self.mol=mol
        self.wns=f['bin_centers'][...]
        self.wnedges=f['bin_edges'][...]
        if 'units' in f['bin_edges'].attrs:
            self.wn_unit=f['bin_edges'].attrs['units']
        else:
            if 'units' in f['bin_centers'].attrs:
                self.wn_unit=f['bin_centers'].attrs['units']
        self.kdata=f['kcoeff'][...]
        self.kdata_unit=f['kcoeff'].attrs['units']
        self.tgrid=f['t'][...]
        self.pgrid=f['p'][...]
        self.logpgrid=np.log10(self.pgrid)
        self.p_unit=f['p'].attrs['units']
        if 'weights' in f.keys():
            self.weights=f['weights'][...]
        else:
            raise RuntimeError('No weights keyword. This file is probably a cross section file.')
        self.ggrid=f['samples'][...]
        self.gedges=np.insert(np.cumsum(self.weights),0,0.)
        self.logk=False
        f.close()  
        self.Np,self.Nt,self.Nw,self.Ng=self.kdata.shape

    def write_hdf5(self, filename):
        """Saves data in a hdf5 format

        Parameters
        ----------
            filename: str
                Name of the file to be created and saved
        """
        fullfilename=filename
        if not filename.lower().endswith(('.hdf5', '.h5')):
            fullfilename=filename+'.hdf5'
        compression="gzip"
        f = h5py.File(fullfilename, 'w')
        f.create_dataset("mol_name", data=self.mol)
        f.create_dataset("p", data=self.pgrid,compression=compression)
        f["p"].attrs["units"] = self.p_unit
        f.create_dataset("t", data=self.tgrid,compression=compression)
        f.create_dataset("kcoeff", data=self.kdata,compression=compression)
        f["kcoeff"].attrs["units"] = self.kdata_unit
        f.create_dataset("samples", data=self.ggrid,compression=compression)
        f.create_dataset("weights", data=self.weights,compression=compression)
        f.create_dataset("bin_edges", data=self.wnedges,compression=compression)
        f.create_dataset("bin_centers", data=self.wns,compression=compression)
        f["bin_edges"].attrs["units"] = self.wn_unit
        f.close()    

    def read_LMDZ(self, path=None, res=None, band=None, mol=None):
        """Initializes k coeff table and supporting data from a .dat file in a gcm friendly format.
        
        Units are assumed to be cm^2 for kdata and mbar for pressure. 

        Parameters
        ----------
            path: str
                Name of the directory with the various input files
            res: str
                "IRxVI" where IR and VI are the numbers of bands
                in the infrared and visible of the k table to load.
            band: str
                "IR" or "VI" to specify which band to load.
            mol: str
                Name of the molecule to be saved in the Ktable object. 
        """        
        if (path is None) or (res is None): \
            raise TypeError("You should provide an input directory name and a resolution")

        self.filename=path
        if mol is not None:
            self.mol=mol
        else:
            self.mol=os.path.basename(self.filename).split(self._settings._delimiter)[0]

        self.weights=np.loadtxt(os.path.join(path,'g.dat'),skiprows=1)[:-1]
            # we remove the last point that is always zero.
            # in the gcm this last point is intended to take care of continuum
        self.Ng=self.weights.size
        self.gedges=np.insert(np.cumsum(self.weights),0,0.)
        self.ggrid=(self.gedges[1:]+self.gedges[:-1])*0.5

        self.p_unit='mbar'
        self.logpgrid=np.loadtxt(os.path.join(path,'p.dat'),skiprows=1)*1.
        self.Np=self.logpgrid.size
        self.pgrid=10**self.logpgrid

        self.tgrid=np.loadtxt(os.path.join(path,'T.dat'),skiprows=1)
        self.Nt=self.tgrid.size

        if band is None:
            raw=np.loadtxt(os.path.join(path,res,'narrowbands.in'), skiprows=1, unpack=True)
        else:
            raw=np.loadtxt(os.path.join(path,res,'narrowbands_'+band+'.in'), \
                skiprows=1, unpack=True)
        self.wnedges=np.append(raw[0],raw[1,-1])
        self.wns=(self.wnedges[1:]+self.wnedges[:-1])*0.5
        self.Nw=self.wns.size
        
        self.kdata_unit='cm^2/molecule'
        if band is None:
            file_to_load=os.path.join(path,res,'corrk_gcm.dat')
        else:
            file_to_load=os.path.join(path,res,'corrk_gcm_'+band+'.dat')        
        tmp=np.loadtxt(file_to_load).reshape((self.Nt,self.Np,self.Nw,self.Ng+1),order='F')
        self.kdata=tmp[:,:,:,:-1].transpose((1,0,2,3)) 
        # also removing the last g point which is equal to 0.
        self.logk=False        
        return None

    def write_LMDZ(self, path, band='IR', fmt='%22.15e', write_only_metadata=False):
        """Saves data in a LMDZ friendly format.
        
        The gcm requires p in mbar and kdata in cm^2/molec.
        The conversion is done automatically.

        Parameters
        ----------
            path: str
                Name of the directory to be created and saved,
                the one that will contain all the necessary files
            band: str
                The band you are computing: 'IR' or 'VI'
            fmt: str
                Fortran format for the corrk file. 
        """
        try:
            os.mkdir(path)
        except FileExistsError:
            print('Directory was already there '+path)
        file = open(os.path.join(path,'p.dat'), "w")
        file.write(str(self.Np)+'\n')
        lp_to_write=self.logpgrid+np.log10(u.Unit(self.p_unit).to(u.Unit('mbar')))
        for lp in lp_to_write:
            file.write(str(lp)+'\n')
        file.close()

        file = open(os.path.join(path,'T.dat'), "w")
        file.write(str(self.Nt)+'\n')
        for t in self.tgrid:
            file.write(str(t)+'\n')
        file.close()

        file = open(os.path.join(path,'g.dat'), "w")
        file.write(str(self.Ng+1)+'\n')
        for g in self.weights:
            file.write(str(g)+'\n')
        file.write(str(0.)+'\n')
        file.close()

        dirname=os.path.join(path,band+str(self.Nw))
        try:
            os.mkdir(dirname)
        except FileExistsError:
            print('Directory was already there: '+dirname)

        file = open(os.path.join(dirname,'narrowbands_'+band+'.in'), "w")
        file.write(str(self.Nw)+'\n')
        for iw in range(self.Nw):
            file.write(str(self.wnedges[iw])+' '+str(self.wnedges[iw+1])+'\n')
        file.close()

        if not write_only_metadata:
            #file = open(dirname+'/corrk_gcm_IR.in', "w")
            data_to_write=self.kdata.transpose((1,0,2,3)).flatten(order='F')
            data_to_write=data_to_write*u.Unit(rm_molec(self.kdata_unit)).to(u.Unit('cm^2'))
            data_to_write=np.append(data_to_write, \
                np.zeros(self.Np*self.Nt*self.Nw)).reshape((1,self.Np*self.Nt*self.Nw*(self.Ng+1)))
            np.savetxt(os.path.join(dirname,'corrk_gcm_'+band+'.dat'),data_to_write,fmt=fmt)

    def read_nemesis(self, filename=None, mol=None):
        """Initializes k coeff table and supporting data from a Nemesis binary file (.kta)

        Parameters
        ----------
            file: str
                Name of the input Nemesis binary file.
            mol: str, optional
                Name of the molecule.
        """
        if (filename is None or not filename.lower().endswith('.kta')):
            raise RuntimeError("You should provide an input nemesis (.kta) file")
        
        if mol is not None:
            self.mol=mol
        else:
            self.mol=os.path.basename(self.filename).split(self._settings._delimiter)[0]

        self.Np, self.Nt, self.Nw, self.Ng, \
            self.pgrid, self.tgrid, self.wns, \
            self.ggrid, self.weights, self.kdata = read_nemesis_binary(filename)

        self.logpgrid=np.log10(self.pgrid)
        self.wnedges=np.concatenate(  \
            ([self.wns[0]],(self.wns[:-1]+self.wns[1:])*0.5,[self.wns[-1]]))
        self.gedges=np.insert(np.cumsum(self.weights),0,0.)
        self.p_unit='bar'
        self.kdata_unit='cm^2/molecule'

    def write_nemesis(self, filename):
        """Saves data in a nemesis format.

        base on a routine provided by K. Chubb.

        Parameters
        ----------
            filename: str
                Name of the file to be created and saved
        """
        int_format = np.int32
        float_format = np.float32
        with open(filename,'wb') as o:
            o.write(int_format(11+2*self.Ng+2+self.Np+self.Nt+self.Nw).tostring())
            #int for irec0 (5221)
            o.write(int_format(self.Nw).tostring()) #int
            o.write(float_format(self.wls[-1]).tostring()) #float for VMIN (0.30015)
            o.write(float_format(-1.).tostring()) #float
            o.write(float_format(0.).tostring()) #float for FWHM
            o.write(int_format(self.Np).tostring())
            o.write(int_format(self.Nt).tostring())
            o.write(int_format(self.Ng).tostring())
            o.write(int_format(nemesis_hitran_id_numbers[self.mol]).tostring()) #IDGAS1 FROM HITRAN
            o.write(int_format(0).tostring())
            o.write(float_format(self.ggrid).tostring())
            o.write(float_format(self.weights).tostring())
            o.write(float_format(0.).tostring()) #float for 0
            o.write(float_format(0.).tostring()) #float for 0
            o.write(float_format(self.pgrid).tostring())
            o.write(float_format(self.tgrid).tostring())
            o.write(float_format(self.wls[::-1]).tostring())
            o.write(float_format(self.kdata[:,:,::-1,:].transpose(2,0,1,3)*1.e20).tostring())

    def xtable_to_ktable(self, xtable=None, wnedges=None, weights=None, ggrid=None,
        quad='legendre', order=20, mid_dw=True, write=0):
        """Fills the :class:`~exo_k.ktable.Ktable` object with a k-coeff table computed
        from a :class:`~exo_k.xtable.Xtable` object.

        The p and kcorr units are inherited from the :class:`~exo_k.xtable.Xtable` object.

        Parameters
        ----------
            xtable: :class:`~exo_k.xtable.Xtable`
                input Xtable object instance
            wnedges: Array
                edges of the wavenumber bins to be used to compute the corrk

        Other Parameters
        ----------------
            weights: array, optional
                If weights are provided, they are used instead of the legendre quadrature. 
            quad: string, optional
                Type of quadrature used. Default is 'legendre'
            order: Integer, optional
                Order of the Gauss legendre quadrature used. Default is 20.
            mid_dw: boolean, optional

                * If True, the Xsec values in the high resolution xtable data are assumed to
                  cover a spectral interval that is centered around
                  the corresponding wavenumber value.
                  The first and last Xsec values are discarded. 
                * If False, each interval runs from the wavenumber value to the next one.
                  The last Xsec value is dicarded.
        """        
        if xtable is None: raise TypeError("You should provide an input Xtable object")
        if wnedges is None: raise TypeError("You should provide an input wavenumber array")

        self.copy_attr(xtable,cp_kdata=False)
        self.wnedges=np.array(wnedges)
        if self.wnedges.size<2: raise TypeError('wnedges should at least have two values')
        self.wns=0.5*(self.wnedges[1:]+self.wnedges[:-1])

        if weights is not None:
            self.weights=np.array(weights)
            self.gedges=np.concatenate(([0],np.cumsum(self.weights)))
            if ggrid is not None: 
                self.ggrid=np.array(ggrid)
            else:
                self.ggrid=(self.gedges[1:]+self.gedges[:-1])*0.5
        else:
            if quad=='legendre':
                self.weights,self.ggrid,self.gedges=gauss_legendre(order)
            else:
                raise NotImplementedError("Type of quadrature (quad keyword) not known.")
        self.Ng=self.weights.size
        if write >6 : print(self.weights,self.ggrid,self.gedges)
        self.logk=False
        self.Nw=self.wns.size
        self.kdata=np.zeros(self.shape)
        if mid_dw:
            dwn_hr=(xtable.wns[2:]-xtable.wns[:-2])*0.5
            wn_hr=xtable.wns[1:-1]
        else:
            dwn_hr=(xtable.wns[1:]-xtable.wns[:-1])
            wn_hr=xtable.wns[:-1]
        for iP in range(self.Np):
          for iT in range(self.Nt):
            if mid_dw:
                k_hr=xtable.kdata[iP,iT,1:-1]
            else:
                k_hr=xtable.kdata[iP,iT,:-1]
            #print(self.ggrid)
            #print(wn_hr[[0,-1]])
            #print(self.wnedges[[0,-1]])
            self.kdata[iP,iT]=spectrum_to_kdist(k_hr,wn_hr,dwn_hr,self.wnedges,self.ggrid)

    def hires_to_ktable(self, path=None, filename_grid=None,
        logpgrid=None, tgrid=None, wnedges=None,
        quad='legendre', order=20, weights=None, ggrid=None,
        mid_dw=True, write=0, mol=None,
        grid_p_unit='Pa', p_unit='unspecified',
        kdata_unit='unspecified', file_kdata_unit='unspecified',
        **kwargs):
        """Computes a k coeff table from :class:`~exo_k.util.hires_spectrum.Hires_spectrum`
        objects.

        .. warning::
            By default, log pressures are specified in Pa in logpgrid!!! If you want
            to use another unit, do not forget to specify it with the grid_p_unit keyword.

        Parameters
        ----------
            path : String
                directory with the input files
            filename_grid : array of str with shape (logpgrid.size,tgrid.size)
                Names of the input high-res spectra. If None, the files are assumed to
                follow Kspectrum/LMDZ convention, i.e.
                be of the type 'k001', 'k002', etc. 
                See :func:`~exo_k.util.filenames.create_fname_grid_Kspectrum_LMDZ`
                for possible additional keyword arguments. 
            logpgrid: array
                Grid in log(pressure) of the input. Default unit is Pa, but can be changed
                with the `grid_p_unit` keyword.
            tgrid: array
                Grid in temperature of the input.
            wnedges : array
                Edges of the wavenumber bins to be used to compute the corr-k

        Other Parameters
        ----------------
            weights: array, optional
                If weights are provided, they are used instead of the legendre quadrature. 
            quad : string, optional
                Type of quadrature used. Default is 'legendre'
            order : Integer, optional
                Order of the Gauss legendre quadrature used. Default is 20.
            mid_dw: boolean, optional
                * If True, the Xsec values in the high resolution xsec data are assumed to
                  cover a spectral interval that is centered around
                  the corresponding wavenumber value.
                  The first and last Xsec values are discarded. 
                * If False, each interval runs from the wavenumber value to the next one.
                  The last Xsec value is dicarded.
            mol: string, optional
                Give a name to the molecule. Useful when used later in a Kdatabase
                to track molecules.
            p_unit: str, optional
                Pressure unit to convert to.
            grid_p_unit : str, optional
                Unit of the specified `logpgrid`.
            kdata_unit : str, optional
                Kdata unit to convert to.
            file_kdata_unit : str, optional
                Kdata unit in input files.

        See :func:`~exo_k.util.filenames.create_fname_grid_Kspectrum_LMDZ`
        or :class:`~exo_k.util.hires_spectrum.Hires_spectrum`
        for a list of additional arguments that can be provided to those funtion
        through `**kwargs`.
        """        
        if path is None: raise TypeError("You should provide an input hires_spectrum directory")
        if wnedges is None: raise TypeError("You should provide an input wavenumber array")

        self.filename=path
        if mol is not None:
            self.mol=mol
        else:
            self.mol=os.path.basename(self.filename).split(self._settings._delimiter)[0]

        if weights is not None:
            self.weights=np.array(weights)
            self.gedges=np.concatenate(([0],np.cumsum(self.weights)))
            if ggrid is not None: 
                self.ggrid=np.array(ggrid)
            else:
                self.ggrid=(self.gedges[1:]+self.gedges[:-1])*0.5
        else:
            if quad=='legendre':
                self.weights,self.ggrid,self.gedges=gauss_legendre(order)
            else:
                raise NotImplementedError("Type of quadrature (quad keyword) not known.")
        self.Ng=self.weights.size

        conversion_factor=u.Unit(grid_p_unit).to(u.Unit('Pa'))
        self.logpgrid=np.array(logpgrid, dtype=float)+np.log10(conversion_factor)
        self.pgrid=10**self.logpgrid #in Pa
        self.p_unit='Pa'
        self.Np=self.logpgrid.size
        if write >= 3 : print(self.Np,self.pgrid)

        self.tgrid=np.array(tgrid)
        self.Nt=self.tgrid.size
        if write >= 3 : print(self.Nt,self.tgrid)

        self.wnedges=np.array(wnedges)
        if self.wnedges.size<2: raise TypeError('wnedges should at least have two values')
        self.wns=(self.wnedges[1:]+self.wnedges[:-1])*0.5
        self.Nw=self.wns.size
        
        self.kdata=np.zeros(self.shape)
        if filename_grid is None:
            filename_grid=create_fname_grid_Kspectrum_LMDZ(self.Np, self.Nt,
                **select_kwargs(kwargs,['suffix','nb_digit']))
        else:
            filename_grid=np.array(filename_grid)
        
        for iP in range(self.Np):
          for iT in range(self.Nt):
            filename=filename_grid[iP,iT]
            fname=os.path.join(path,filename)
            if write >= 3 : print(fname)

            spec_hr=Hires_spectrum(fname, file_kdata_unit=file_kdata_unit,
                **select_kwargs(kwargs,['skiprows','wn_column','mult_factor',
                    'kdata_column','data_type']))
            # for later conversion, the real kdata_unit is in spec_hr.kdata_unit
            self.kdata_unit=spec_hr.kdata_unit
            was_xsec=(spec_hr.data_type=='xsec')
            wn_hr=spec_hr.wns
            k_hr=spec_hr.kdata
            if mid_dw:
                dwn_hr=(wn_hr[2:]-wn_hr[:-2])*0.5
                wn_hr=wn_hr[1:-1]
                k_hr=k_hr[1:-1]
            else:
                dwn_hr=(wn_hr[1:]-wn_hr[:-1])
                wn_hr=wn_hr[:-1]
                k_hr=k_hr[:-1]
            self.kdata[iP,iT]=spectrum_to_kdist(k_hr,wn_hr,dwn_hr,self.wnedges,self.ggrid)
            if not was_xsec:
                self.kdata[iP,iT]=self.kdata[iP,iT]*KBOLTZ*self.tgrid[iT]/self.pgrid[iP]
        if not was_xsec:
            self.kdata=self.kdata*u.Unit(self.kdata_unit).to(u.Unit('m^-1'))
            self.kdata_unit='m^2/molecule'
            # Accounts for the conversion of the abs_coeff to m^-1, so we know that
            #  self.kdata is now in m^2/molecule. Can now convert to the desired unit.
        if self._settings._convert_to_mks and kdata_unit is 'unspecified':
            kdata_unit='m^2/molecule'
        self.convert_kdata_unit(kdata_unit=kdata_unit)
        # converts from self.kdata_unit wich is either:
        #   - 'm^2/molecule' data_type was 'abs_coeff'
        #   - The units after Hires_spectrum (which have already taken
        #        file_kdata_unit into account)
        #  to the desired unit.
        self.convert_p_unit(p_unit=p_unit)


    def copy(self, cp_kdata=True, ktab5d=False):
        """Creates a new instance of :class:`Ktable` object and (deep) copies data into it

        Parameters
        ----------
            cp_kdata: bool, optional
                If false, the kdata table is not copied and
                only the structure and metadata are. 
            ktab5d: bool, optional
                If true, creates a Ktable5d object with the same structure. 
                Data are not copied.

        Returns
        -------
            :class:`Ktable` or :class:`Ktable5d`
                A new :class:`Ktable` or :class:`Ktable5d` 
                instance with the same structure as self.
        """
        if ktab5d:
            res=Ktable5d()
            cp_kdata=False
        else:    
            res=Ktable()
        res.copy_attr(self,cp_kdata=cp_kdata)
        res.weights = np.copy(self.weights)
        res.ggrid   = np.copy(self.ggrid)
        res.gedges  = np.copy(self.gedges)
        return res

    def gindex(self, g):
        """Finds the index corresponding to the given g
        """
        return min(np.searchsorted(self.ggrid,g),self.Ng-1)

    def spectrum_to_plot(self, p=1.e-5, t=200., x=1., g=None):
        """provide the spectrum for a given point to be plotted

        Parameters
        ----------
            p : float
                Pressure (Ktable pressure unit)
            t : float
                Temperature(K)
            g: float
                Gauss point
            x: float
                Mixing ratio of the species
        """
        if g is None: raise RuntimeError('A gauss point should be provided with the g= keyword.')
        gindex=self.gindex(g)
        return self.interpolate_kdata(log10(p), t, x)[0,:,gindex]

    def plot_distrib(self, ax, p=1.e-5, t=200., wl=1., x=1., xscale=None, yscale='log', **kwarg):
        """Plot the distribution for a given point

        Parameters
        ----------
            p : float
                Pressure (Ktable pressure unit)
            t : float
                Temperature(K)
            wl: float
                Wavelength (micron)
        """
        wlindex=self.wlindex(wl)
        toplot=self.interpolate_kdata(log10(p),t)[0,wlindex]*x
        if xscale is not None: 
            ax.set_xscale(xscale)
            ax.plot(1.-self.ggrid,toplot,**kwarg)
            ax.set_xlabel('1-g')
        else:
            ax.plot(self.ggrid,toplot,**kwarg)
            ax.set_xlabel('g')
        ax.set_ylabel('Cross section ('+self.kdata_unit+')')
        ax.grid(True)
        if yscale is not None: ax.set_yscale(yscale)

    def __repr__(self):
        """Method to output header
        """
        out1=super().__repr__()
        output=out1+"""
        weights      : {wg}
        data oredered following p, t, wn, g
        shape        : {shape}
        """.format(wg=self.weights, shape=self.shape)
        return output

    def RandOverlap(self, other, x_self, x_other, write=0, use_rebin=False):
        """Method to randomly mix the opacities of 2 species (self and other).

        Parameters
        ----------
            other: :class:`Ktable`
                A :class:`Ktable` object to be mixed with. Dimensions should be the same as self.
            x_self: float or array
                Volume mixing ratio of the first species
            x_other: float or array
                Volume mixing ratio of the species to be mixed with.

        If one of these is None, the kcoeffs of the species in question
        are considered to be already normalized with respect to the mixing ratio.

        Returns
        -------
            array
                array of k-coefficients for the mix. 
        """
        if not np.array_equal(self.shape,other.shape):
            raise TypeError("""in RandOverlap: kdata tables do not have the same dimensions.
                I'll stop now!""")
        Ng=self.Ng
        weights=self.weights
        try:
            kdatas=self.vmr_normalize(x_self)
        except TypeError:
            raise TypeError('Gave bad mixing ratio format to vmr_normalize')
        try:
            kdatao=other.vmr_normalize(x_other)
        except TypeError:
            raise TypeError('Gave bad mixing ratio format to vmr_normalize')
        kdataconv=np.zeros((self.Np,self.Nt,self.Nw,Ng**2))
        weightsconv=np.zeros(Ng**2)
        newkdata=np.zeros(self.shape)
#        newkdata2=np.zeros((self.shape[0]*self.shape[1]*self.shape[2],Ng))

        for jj in range(Ng):
            for ii in range(Ng):
                weightsconv[ii*Ng+jj]=weights[jj]*weights[ii]

        if write >= 3 : start=time.time()
        kdata_conv_loop(kdatas,kdatao,kdataconv,self.shape)
        if write >= 3 : end=time.time();print("kdata conv: ",end - start)        

        if write >= 3 : start=time.time()
        for iP in range(self.Np):
          for iT in range(self.Nt):
            for iW in range(self.Nw):
              tmp=kdataconv[iP,iT,iW,:]
              indices=np.argsort(tmp)
              kdatasort=tmp[indices]
              weightssort=weightsconv[indices]
              newggrid=np.cumsum(weightssort)
              if use_rebin:
                  newggrid=np.concatenate(([0.],newggrid))
                  kdatasort=np.concatenate(([kdatasort[0]],kdatasort))
                  newkdata[iP,iT,iW,:]=rebin(kdatasort,newggrid,self.gedges)
              else:
                  newkdata[iP,iT,iW,:]=np.interp(self.ggrid,newggrid,kdatasort)
        if write >= 3 : end=time.time();print("kdata rebin with loop: ",end - start)        

        return newkdata
        
    def bin_down(self, wnedges=None, weights=None, ggrid=None,
        remove_zeros=False, num=300, use_rebin=False, write=0):
        """Method to bin down a kcoeff table to a new grid of wavenumbers

        Parameters
        ----------
            wnedges: array
                Edges of the new bins of wavenumbers (cm-1)
                onto which the kcoeff should be binned down.
                if you want Nwnew bin in the end, wnedges.size must be Nwnew+1
                wnedges[0] should be greater than self.wnedges[0] (JL20 not sure anymore)
                wnedges[-1] should be lower than self.wnedges[-1]
            weights: array, optional
                Desired weights for the resulting Ktable.
            ggrid: array, optional
                Desired g-points for the resulting Ktable.
                Must be consistent with provided weights.
                If not given, they are taken at the midpoints of the array
                given by the cumulative sum of the weights
        """
        old_ggrid=self.ggrid
        if weights is not None:
            self.weights=np.array(weights)
            self.gedges=np.concatenate(([0],np.cumsum(self.weights)))
            if ggrid is not None: 
                self.ggrid=np.array(ggrid)
            else:
                self.ggrid=(self.gedges[1:]+self.gedges[:-1])*0.5
            self.Ng=self.ggrid.size
        wnedges=np.array(wnedges)
        if wnedges.size<2: raise TypeError('wnedges should at least have two values')
        wngrid_filter = np.where((wnedges <= self.wnedges[-1]) & (wnedges >= self.wnedges[0]))[0]
        if not is_sorted(wnedges):
            raise RuntimeError('wnedges should be sorted.')
        indicestosum,wn_weigths=rebin_ind_weights(self.wnedges, wnedges[wngrid_filter])
        if write> 10 :print(indicestosum);print(wn_weigths)

        self.kdata=bin_down_corrk_numba((self.Np,self.Nt,wnedges.size-1,self.Ng), \
            self.kdata, old_ggrid, self.ggrid, self.gedges, indicestosum, wngrid_filter, \
            wn_weigths, num, use_rebin)

        self.wnedges=wnedges
        self.wns=(wnedges[1:]+wnedges[:-1])*0.5
        self.Nw=self.wns.size
        if remove_zeros : self.remove_zeros(deltalog_min_value=10.)


########### OBSOLETE METHODS #################

    def remap_logPT_grid(self, logp_array=None, t_array=None):
        """Obsolete method. Replaced by remap_logPT in the parent class: Data_table.
        """
        res=Ktable()
        res.mol     =self.mol
        res.logpgrid=logp_array
        res.pgrid   =10**res.logpgrid
        res.tgrid   =t_array
        res.wns     =self.wns
        res.wnedges =self.wnedges
        res.weights =self.weights
        res.ggrid   =self.ggrid
        res.gedges  =self.gedges
        res.Np      =logp_array.shape
        res.Nt      =t_array.shape
        res.Nw      =self.Nw
        res.Ng      =self.Ng
        res.logk    =self.logk
        res.p_unit  =self.p_unit
        res.kdata_unit=self.kdata_unit
        #res.kdata=np.zeros((res.Np,res.Nt,res.Nw ,res.Ng))
        tind,tweight=interp_ind_weights(t_array,self.tgrid)
        lpind,lpweight=interp_ind_weights(logp_array,self.logpgrid)
        #print(tind,tweight)
        #print(lpind,lpweight)
        lpindextended=lpind.reshape((lpind.size,1))
        tw=tweight.reshape((1,tweight.size,1,1))    
        pw=lpweight.reshape((lpweight.size,1,1,1))
        # trick to broadcast over Nw and Ng a few lines below
        kc_p1t1=self.kdata[lpindextended,tind]
        kc_p0t1=self.kdata[lpindextended-1,tind]
        kc_p1t0=self.kdata[lpindextended,tind-1]
        kc_p0t0=self.kdata[lpindextended-1,tind-1]
        # kdata_tmp=kc_p1t1*pw*tw+kc_p0t1*(1.-pw)*tw+kc_p1t0*pw*(1.-tw)+kc_p0t0*(1.-pw)*(1.-tw)
        kdata_tmp=  np.log10(kc_p1t1)*pw*tw          \
                        +np.log10(kc_p0t1)*(1.-pw)*tw \
                        +np.log10(kc_p1t0)*pw*(1.-tw) \
                        +np.log10(kc_p0t0)*(1.-pw)*(1.-tw)
        res.kdata=10**kdata_tmp
        #kdata_tmp=  (kc_p1t1)*pw*tw          \
        #                +(kc_p0t1)*(1.-pw)*tw \
        #                +(kc_p1t0)*pw*(1.-tw) \
        #                +(kc_p0t0)*(1.-pw)*(1.-tw)
        #res.kdata=kdata_tmp
        return res
