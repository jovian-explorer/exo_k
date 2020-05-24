# -*- coding: utf-8 -*-
"""
@author: jeremy leconte
"""
import os.path
from math import log10
import pickle
import h5py
import numpy as np
from .data_table import Data_table
from .util.interp import rebin_ind_weights
from .util.cst import KBOLTZ


class Xtable(Data_table):
    """A class that handles tables of cross sections. Based on Data_table.
    """

    def __init__(self, *filename_filters, filename=None,
        p_unit='unspecified', old_p_unit='unspecified',
        kdata_unit='unspecified', old_kdata_unit='unspecified',
        remove_zeros=False, search_path=None, mol=None):
        """Initializes cross section table and
        supporting data from a file based on its extension.

        Parameters
        ----------
            filename: str
                Relative or absolute path to the input file.
            filename_filters: sequence of string
                As many strings as necessary to uniquely define
                a file in the Settings()._search_path
                The Settings()._search_path will be searched for a file
                with all the filename_filters in the name.
                The filename_filters can contain *.

        If there is no filename or filename_filters provided,
        just creates an empty object to be filled later

        Other Parameters
        ----------------
            search_path: str, optional
                If search_path is provided,
                it locally overrides the global _search_path settings
                and only files in search_path are returned.

        For unit options, see Ktable init function.
        """
        super().__init__()
        if filename is not None:
            self.filename=filename
        elif len(filename_filters)>=1:
            self.filename=self._settings.list_files(*filename_filters,
                                 only_one=True, search_path=search_path)[0]
        if self.filename is not None:
            if self.filename.lower().endswith('pickle'):
                self.read_pickle(filename=self.filename, mol=mol)
            elif self.filename.lower().endswith(('.hdf5', '.h5')):
                self.read_hdf5(filename=self.filename, mol=mol)
            elif self.filename.lower().endswith(('.dat')):
                self.read_exo_transmit(filename=self.filename, mol=mol)
            else:
                raise NotImplementedError("""Requested format not recognized.
                Currently recognized formats are Exomol .pickle, .hdf5, and Exo_Transmit .dat.""")

        if self.kdata is not None:
            if self._settings._convert_to_mks:
                if p_unit is 'unspecified': p_unit='Pa'
                if kdata_unit is 'unspecified': kdata_unit='m^2/molecule'
            self.convert_p_unit(p_unit=p_unit,old_p_unit=old_p_unit)
            self.convert_kdata_unit(kdata_unit=kdata_unit,old_kdata_unit=old_kdata_unit)
            if remove_zeros : self.remove_zeros(deltalog_min_value=10.)

    @property
    def shape(self):
        """Returns the shape of self.kdata
        """
        return np.array([self.Np,self.Nt,self.Nw])

    def read_pickle(self, filename=None, mol=None):
        """Initializes xsec table and supporting data from an Exomol pickle file

        Parameters
        ----------
            filename : str
                Relative or absolute name of the input pickle file
            mol : str, optional
                Overrides the name of the molecule to be put in the :class:`Xtable` object.
        """
        if filename is None: raise TypeError("You should provide an input pickle filename")
        pickle_file=open(filename,'rb')
        raw=pickle.load(pickle_file, encoding='latin1')
        pickle_file.close()
        
        if mol is not None: 
            self.mol=mol
        else:
            self.mol=raw['name']
            if self.mol=='H2OP': self.mol='H2O'

        self.tgrid=raw['t']
        self.wns=raw['wno']
        self.wnedges=np.concatenate( \
            ([self.wns[0]],(self.wns[:-1]+self.wns[1:])*0.5,[self.wns[-1]]))
        self.logk=False

        #deals with the p grid and units
        if 'p_unit' in raw.keys():
            self.p_unit=raw['p_unit']
        else:
            self.p_unit='unspecified'
        self.pgrid=raw['p']
        self.logpgrid=np.log10(self.pgrid)

        #deals with the k grid and units
        if 'kdata_unit' in raw.keys():
            self.kdata_unit=raw['kdata_unit']
        else:
            self.kdata_unit='unspecified'
        self.kdata=raw['xsecarr']

        self.Np,self.Nt,self.Nw=self.kdata.shape

    def write_pickle(self, filename):
        """Saves data in a pickle format

        Parameters
        ----------
            filename: str
                Relative or absolute name of the file to be created and saved
        """
        fullfilename=filename
        if not filename.lower().endswith('.pickle'): fullfilename=filename+'.pickle'
        pickle_file=open(fullfilename,'wb')
        dictout={'name':self.mol,
                 'p':self.pgrid,
                 'p_unit':self.p_unit,
                 't':self.tgrid,
                 'wno':self.wns,
                 'xsecarr':self.kdata,
                 'kdata_unit':self.kdata_unit}
        #print(dictout)
        pickle.dump(dictout,pickle_file,protocol=-1)
        pickle_file.close()

    def read_hdf5(self, filename=None, mol=None):
        """Initializes k coeff table and supporting data from an Exomol hdf5 file

        Parameters
        ----------
            filename : str
                Name of the input hdf5 file
            mol : str, optional
                Overrides the name of the molecule to be put in the :class:`Xtable` object.
        """
        if (filename is None or not filename.lower().endswith(('.hdf5', '.h5'))):
            raise RuntimeError("You should provide an input hdf5 file")
        f = h5py.File(filename, 'r')
        if 'mol_name' in f.attrs:
            self.mol=f.attrs['mol_name']
        else:
            self.mol=os.path.basename(filename).split(self._settings._delimiter)[0]
        if mol is not None:
            self.mol=mol
        self.wns=f['bin_edges'][...]
        self.wnedges=np.concatenate(  \
            ([self.wns[0]],(self.wns[:-1]+self.wns[1:])*0.5,[self.wns[-1]]))
        self.kdata=f['xsecarr'][...]
        self.kdata_unit=f['xsecarr'].attrs['units']
        self.tgrid=f['t'][...]
        self.pgrid=f['p'][...]
        self.logpgrid=np.log10(self.pgrid)
        self.p_unit=f['p'].attrs['units']
        self.logk=False
        f.close()  
        self.Np,self.Nt,self.Nw=self.kdata.shape

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
        f.attrs["mol_name"] = self.mol
        f.create_dataset("p", data=self.pgrid,compression=compression)
        f["p"].attrs["units"] = self.p_unit
        f.create_dataset("t", data=self.tgrid,compression=compression)
        f.create_dataset("xsecarr", data=self.kdata,compression=compression)
        f["xsecarr"].attrs["units"] = self.kdata_unit
        f.create_dataset("bin_edges", data=self.wns,compression=compression)
        f.close()

    def read_exo_transmit(self, filename, mol=None):
        """Creates an xsec object from an exo_transmit like spectra.
        See https://github.com/elizakempton/Exo_Transmit or Kempton et al. (2016) for details.
        Pressures are expected to be in Pa and cross sections in m^2/molecule

        Parameters
        ----------
            filename : str
                Name of the input file.
            mol: str
                Overrides the name of the molecule to be put in the :class:`Xtable` object.
        """
        tmp_wlgrid=[]
        tmp_kdata=[]
        with open(filename, 'r') as file:
            tmp = file.readline().split()
            self.tgrid=np.array([float(ii) for ii in tmp])
            tmp = file.readline().split()
            self.pgrid=np.array([float(ii) for ii in tmp])
            self.logpgrid=np.log10(self.pgrid)
            self.Np=self.pgrid.size
            self.Nt=self.tgrid.size
            while True:
                line=file.readline()
                if line is None or line=='': break
                tmp_wlgrid.append(float(line.split()[0]))
                tmp_kdata.append([])
                for _ in range(self.Np):
                    tmp_kdata[-1].append(file.readline().split()[1:])
        tmp_wlgrid=np.array(tmp_wlgrid)
        self.Nw=tmp_wlgrid.size
        self.kdata=np.zeros((self.Np,self.Nt,self.Nw))
        for iP in range(self.Np):
            for iT in range(self.Nt):
                for iW in range(self.Nw):
                    self.kdata[iP,iT,iW]=tmp_kdata[self.Nw-iW-1][iP][iT]
        self.wns=0.01/tmp_wlgrid[::-1]
        self.wnedges=np.concatenate( \
            ([self.wns[0]],(self.wns[:-1]+self.wns[1:])*0.5,[self.wns[-1]]))
        self.p_unit='Pa' 
        self.kdata_unit='m^2/molecule'
        if mol is not None:
            self.mol=mol
        else:
            self.mol=os.path.basename(filename).split(self._settings._delimiter)[0]
  

    def hires_to_xsec(self, path=None, filename_grid=None, logpgrid=None, tgrid=None,
        write=0, mol=None, kdata_unit='unspecified', old_kdata_unit='unspecified',
        k_to_xsec=True):
        """Computes a k coeff table from high resolution cross sections
        in the usual k-spectrum format.

        Parameters
        ----------
            path : String
                directory with the input files
            filename_grid : Numpy Array of strings with shape (logpgrid.size,tgrid.size)
                Names of the input high-res spectra.
            logpgrid: Array
                Grid in log(pressure/Pa) of the input
            tgrid: Array
                Grid intemperature of the input
            mol: str, optional
                give a name to the molecule. Useful when used later in a Kdatabase
                to track molecules.
            k_to_xsec : boolean, optional
                If true, performs a conversion from absorption coefficient (m^-1) to xsec.
        """        
        first=True

        if path is None: raise TypeError("You should provide an input hires_spectrum directory")
        filename_grid=np.array(filename_grid)
        self.filename=path
        if mol is not None:
            self.mol=mol
        else:
            self.mol=os.path.basename(self.filename).split(self._settings._delimiter)[0]

        self.p_unit='Pa'
        self.logpgrid=np.array(logpgrid)
        self.Np=self.logpgrid.size
        self.pgrid=10**self.logpgrid
        if write >= 3 : print(self.Np,self.pgrid)

        self.tgrid=np.array(tgrid)
        self.Nt=self.tgrid.size
        if write >= 3 : print(self.Nt,self.tgrid)
        
        for iP in range(self.Np):
          for iT in range(self.Nt):
            fname=os.path.join(path,filename_grid[iP,iT])
            if fname.lower().endswith('.dat'):
                wn_hr,k_hr=np.loadtxt(fname,skiprows=0,unpack=True)  
            elif fname.lower().endswith(('.hdf5', '.h5')):
                f = h5py.File(fname, 'r')
                wn_hr=f['wns'][...]
                k_hr=f['k'][...]
                f.close()
            else:
                raise NotImplementedError('Input file format not recognized.')
            if k_to_xsec: 
                k_hr=k_hr*KBOLTZ*self.tgrid[iT]/self.pgrid[iP]
            if first:
                self.wns=wn_hr[1:-1]  #to be consistent with kcorr
                self.Nw=self.wns.size
                self.kdata=np.zeros((self.Np,self.Nt,self.Nw))
                self.kdata[iP,iT]=k_hr[1:-1]
                first=False
            else:
                self.kdata[iP,iT]=np.interp(self.wns,wn_hr[1:-1],k_hr[1:-1])
        self.kdata_unit='m^2' #default unit assumed for the input file
        if self._settings._convert_to_mks and kdata_unit is 'unspecified': kdata_unit='m^2/molecule'
        self.convert_kdata_unit(kdata_unit=kdata_unit,old_kdata_unit=old_kdata_unit)

    def bin_down(self, wnedges=None, remove_zeros=False, write=0):
        """Method to bin down a xsec table to a new grid of wavenumbers (in place)

        Parameters
        ----------
            wnedges : array
                Edges of the new bins of wavenumbers (cm-1) onto which the xsec
                should be binned down.
                if you want Nwnew bin in the end, wngrid.size must be Nwnew+1
                wnedges[0] should be greater than self.wnedges[0]
                wnedges[-1] should be lower than self.wnedges[-1]
            remove_zeros: bool, optional
                If True, remove zeros in kdata. 
        """
        wnedges=np.array(wnedges)
        if wnedges.size<2: raise TypeError('wnedges should at least have two values')
        wngrid_filter = np.where((wnedges <= self.wnedges[-1]) & (wnedges >= self.wnedges[0]))[0]
        if write>=3:
            print(self.wnedges);print(wnedges);print(wngrid_filter);print(wnedges[wngrid_filter])

        indicestosum,weights=rebin_ind_weights(self.wnedges, wnedges[wngrid_filter])
        if write>=3: print(indicestosum[0]);print(weights[0])
        Nnew=wnedges.size-1
        Ntmp=wnedges[wngrid_filter].size-1
        newxsec=np.zeros((self.Np,self.Nt,Nnew))
        for iP in range(self.Np):
            for iT in range(self.Nt):
                tmp=self.kdata[iP,iT,:]
                newxsec[iP,iT,wngrid_filter[0:-1]]= \
                    [np.dot(tmp[indicestosum[ii]-1:indicestosum[ii+1]],weights[ii]) \
                        for ii in range(Ntmp)]        
        self.kdata=newxsec
        self.wnedges=wnedges
        self.wns=(wnedges[1:]+wnedges[:-1])*0.5
        self.Nw=self.wns.size
        if remove_zeros : self.remove_zeros(deltalog_min_value=10.)

    def sample(self, wngrid, remove_zeros=False, log_interp=None):
        """Method to re sample a xsec table to a new grid of wavenumbers (in place)

        Parameters
        ----------
            wngrid : array
                Location of the new wavenumbers points (cm-1)
        """
        wngrid=np.array(wngrid)
        wngrid_filter = np.where((wngrid <= self.wnedges[-1]) & (wngrid >= self.wnedges[0]))[0]
        Nnew=wngrid.size
        newxsec=np.zeros((self.Np,self.Nt,Nnew))
        if log_interp is None: log_interp=self._settings._log_interp
        if log_interp:
            for iP in range(self.Np):
                for iT in range(self.Nt):
                    tmp=np.log(self.kdata[iP,iT,:])
                    newxsec[iP,iT,wngrid_filter]=np.interp(wngrid[wngrid_filter],self.wns,tmp)
            self.kdata=np.exp(newxsec)
        else:
            for iP in range(self.Np):
                for iT in range(self.Nt):
                    tmp=self.kdata[iP,iT,:]
                    newxsec[iP,iT,wngrid_filter]=np.interp(wngrid[wngrid_filter],self.wns,tmp)
            self.kdata=newxsec
        self.wns=wngrid
        self.wnedges=np.concatenate( \
            ([self.wns[0]],0.5*(self.wns[1:]+self.wns[:-1]),[self.wns[-1]]))
        self.Nw=Nnew
        if remove_zeros : self.remove_zeros(deltalog_min_value=10.)

    def sample_cp(self, wngrid, **kwargs):
        """Creates a copy of the instance before resampling it.

        Parameters
        ----------
            See sample method for details. 

        Returns
        -------
            :class:`Xtable` object
                the re-sampled :class:`Xtable`
        """
        res=self.copy()
        res.sample(wngrid, **kwargs)
        return res


    def spectrum_to_plot(self, p=1.e-5, t=200., x=1., g=None):
        """provide the spectrum for a given point to be plotted

        Parameters
        ----------
            p : float
                Pressure (Ktable pressure unit)
            t : float
                Temperature(K)
            x: float
                Volume mixing ratio of the species
            g: is unused but here to be consistent with the method in data_table
        """
        return self.interpolate_kdata(log10(p),t)[0]*x


    def copy(self, cp_kdata=True):
        """Creates a new instance of :class:`Xtable` object and (deep) copies data into it

        Parameters
        ----------
            cp_kdata: bool, optional
                If false, the kdata table is not copied and
                only the structure and metadata are. 

        Returns
        -------
            :class:`Xtable`
                A new :class:`Xtable` instance with the same structure as self.
        """
        res=Xtable()
        res.copy_attr(self, cp_kdata=cp_kdata)
        return res

    def __repr__(self):
        """Method to output header
        """
        out1=super().__repr__()
        output=out1+"""
        data oredered following p, t, wl
        shape        : {shape}
        wl (microns) : {wl}
        xsec unit    : {kdata_unit}
        """.format(shape=self.shape, wl=self.wls, kdata_unit=self.kdata_unit)
        return output
