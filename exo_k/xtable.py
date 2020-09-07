# -*- coding: utf-8 -*-
"""
@author: jeremy leconte
"""
import os.path
from math import log10
import pickle
import h5py
import numpy as np
import astropy.units as u
from .data_table import Data_table
from .util.interp import rebin_ind_weights
from .util.cst import KBOLTZ
from .hires_spectrum import Hires_spectrum
from .util.filenames import create_fname_grid_Kspectrum_LMDZ, select_kwargs

class Xtable(Data_table):
    """A class that handles tables of cross sections.
    """

    def __init__(self, *filename_filters, filename=None,
        p_unit='unspecified', file_p_unit='unspecified',
        kdata_unit='unspecified', file_kdata_unit='unspecified',
        remove_zeros=False, search_path=None, mol=None):
        """Initializes cross section table and
        supporting data from a file based on its extension.

        Parameters
        ----------
            filename: str
                Relative or absolute path to the input file.
            filename_filters: sequence of string
                As many strings as necessary to uniquely define
                a file in the global search path defined in
                :class:`~exo_k.settings.Settings`.
                This path will be searched for a file
                with all the filename_filters in the name.
                The filename_filters can contain '*'.

        If there is no filename or filename_filters provided,
        just creates an empty object to be filled later

        See :class:`~exo_k.ktable.Ktable` __init__ mehthod for documentation on
        `p_unit`, `file_p_unit`, `kdata_unit`, `file_kdata_unit`, `remove_zeros`,
        `search_path`, and `mol` keywords.
        """
        super().__init__()
        if filename is not None:
            self.filename=filename
        elif filename_filters or mol is not None:
            self.filename=self._settings.list_files(*filename_filters, molecule = mol,
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

        super().finalize_init(p_unit=p_unit, file_p_unit=file_p_unit,
            kdata_unit=kdata_unit, file_kdata_unit=file_kdata_unit,
            remove_zeros=remove_zeros)

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
            self.p_unit='bar'
        self.pgrid=raw['p']
        self.logpgrid=np.log10(self.pgrid)

        #deals with the k grid and units
        if 'kdata_unit' in raw.keys():
            self.kdata_unit=raw['kdata_unit']
        else:
            self.kdata_unit='cm^2/molec'
        self.kdata=raw['xsecarr']

        if 'wn_unit' in raw.keys(): self.wn_unit=raw['wn_unit']

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
                 'wn_unit':self.wn_unit,
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
        if 'mol_name' in f:
            self.mol=f['mol_name'][()]
        elif 'mol_name' in f.attrs:
            self.mol=f.attrs['mol_name']
        else:
            self.mol=os.path.basename(filename).split(self._settings._delimiter)[0]
        if isinstance(self.mol, np.ndarray): self.mol=self.mol[0]
        if mol is not None:
            self.mol=mol
        if 'bin_edges' in f:
            self.wns=f['bin_edges'][...]
            if 'units' in f['bin_edges'].attrs:
                self.wn_unit=f['bin_edges'].attrs['units']
        else:
            self.wns=f['bin_centers'][...]
            if 'units' in f['bin_centers'].attrs:
                self.wn_unit=f['bin_centers'].attrs['units']
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
        f.create_dataset("mol_name", data=self.mol)
        f.create_dataset("p", data=self.pgrid,compression=compression)
        f["p"].attrs["units"] = self.p_unit
        f.create_dataset("t", data=self.tgrid,compression=compression)
        f.create_dataset("xsecarr", data=self.kdata,compression=compression)
        f["xsecarr"].attrs["units"] = self.kdata_unit
        f.create_dataset("bin_edges", data=self.wns,compression=compression)
        f["bin_edges"].attrs["units"] = self.wn_unit
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

    def hires_to_xtable(self, path=None, filename_grid=None,
        logpgrid=None, tgrid=None,
        write=0, mol=None,
        grid_p_unit='Pa', p_unit='unspecified',
        kdata_unit='unspecified', file_kdata_unit='unspecified',
        **kwargs):
        """Loads an `Xtable` from high-resolution spectra.

        .. warning::
            By default, log pressures are specified in Pa in logpgrid!!! If you want
            to use another unit, do not forget to specify it with the grid_p_unit keyword.

        see :func:`exo_k.ktable.Ktable.hires_to_ktable` method for details
        on the arguments and options.
        """        
        first=True

        if path is None: raise TypeError("You should provide an input hires_spectrum directory")
        self.filename=path
        if mol is not None:
            self.mol=mol
        else:
            self.mol=os.path.basename(self.filename).split(self._settings._delimiter)[0]

        conversion_factor=u.Unit(grid_p_unit).to(u.Unit('Pa'))
        self.logpgrid=np.array(logpgrid, dtype=float)+np.log10(conversion_factor)
        self.pgrid=10**self.logpgrid #in Pa
        self.p_unit='Pa'
        self.Np=self.logpgrid.size
        if write >= 3 : print(self.Np,self.pgrid)

        self.tgrid=np.array(tgrid)
        self.Nt=self.tgrid.size
        if write >= 3 : print(self.Nt,self.tgrid)
        
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
            if not was_xsec:
                k_hr=k_hr*KBOLTZ*self.tgrid[iT]/self.pgrid[iP]
            if first:
                self.wns=wn_hr[1:-1]  #to be consistent with kcorr
                self.Nw=self.wns.size
                self.kdata=np.zeros((self.Np,self.Nt,self.Nw))
                self.kdata[iP,iT]=k_hr[1:-1]
                first=False
            else:
                self.kdata[iP,iT]=np.interp(self.wns,wn_hr[1:-1],k_hr[1:-1])
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
        data oredered following p, t, wn
        shape        : {shape}
        """.format(shape=self.shape)
        return output
