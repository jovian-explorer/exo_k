# -*- coding: utf-8 -*-
"""
@author: jeremy leconte
"""
from math import log10
import os
import h5py
import numpy as np
import astropy.units as u
from scipy.interpolate import RegularGridInterpolator
from .data_table import Data_table
from .util.interp import rm_molec, rebin_ind_weights, \
        gauss_legendre, spectrum_to_kdist, is_sorted, bin_down_corrk_numba
from .util.cst import KBOLTZ
from .hires_spectrum import Hires_spectrum
from .util.filenames import create_fname_grid_Kspectrum_LMDZ, select_kwargs


class Ktable5d(Data_table):
    """A class that handles tables of k-coefficients with a variable gas.
    Based on the Data_table class that handles basic operations common to Ktable and Xtable.

    This class is specifically designed to deal with
    LMDZ type ktable where there is a variable gas.
    """
    
    def __init__(self, *filename_filters, filename=None, path=None,
        p_unit='unspecified', file_p_unit='unspecified',
        kdata_unit='unspecified', file_kdata_unit='unspecified',
        remove_zeros=False, search_path=None, mol=None,
        **kwargs):
        """Initializes k coeff table with variable gas and
        supporting data from various sources (see below by order of precedence)

        Parameters
        ----------
            filename: str (optional)
                Relative or absolute name of the file to be loaded. 
            path: str
                If none of the above is specifed,
                path can point to a directory with a LMDZ type k coeff table.
                In this case, see read_LMDZ for the keywords to specify.

        If there is no input, just creates an empty object to be filled later

        See :class:`~exo_k.ktable.Ktable` __init__ mehthod for documentation on
        `p_unit`, `file_p_unit`, `kdata_unit`, `file_kdata_unit`, `remove_zeros`,
        `search_path`, and `mol` arguments.

        """
        super().__init__()

        if filename is not None:
            self.filename=filename
        elif filename_filters or (mol is not None and path is None): 
        # a none empty sequence returns a True in a conditional statement
            self.filename=self._settings.list_files(*filename_filters, molecule = mol,
                only_one=True, search_path=search_path)[0]
        if self.filename is not None:
            if self.filename.lower().endswith(('.hdf5', '.h5')):
                self.read_hdf5(filename=self.filename, mol=mol)
            else:
                raise NotImplementedError( \
                    'Requested format not recognized. Should end with .pickle, .hdf5, or .h5')
        elif path is not None:
            self.read_LMDZ(path=path, mol=mol, **kwargs)
        else:                  #if there is no input file, just create an empty object 
            self.wnedges=None
            self.weights=None
            self.ggrid=None
            self.gedges=None
            self.xgrid=None
            self._finterp_kdata=None

        super().finalize_init(p_unit=p_unit, file_p_unit=file_p_unit,
            kdata_unit=kdata_unit, file_kdata_unit=file_kdata_unit,
            remove_zeros=remove_zeros)

        if self.kdata is not None:
            self.setup_interpolation()

    @property
    def shape(self):
        """Returns the shape of self.kdata
        """
        return np.array([self.Np,self.Nt,self.Nx,self.Nw,self.Ng])

    def read_hdf5(self, filename=None, mol=None):
        """Initializes k coeff table and supporting data from an Exomol hdf5 file

        Parameters
        ----------
            filename : str
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
        if mol is not None: self.mol=mol
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
        self.xgrid=f['x'][...]
        if 'weights' in f.keys():
            self.weights=f['weights'][...]
        else:
            raise RuntimeError('No weights keyword. This file is probably a cross section file.')
        self.ggrid=f['samples'][...]
        self.gedges=np.insert(np.cumsum(self.weights),0,0.)
        self.logk=False
        f.close()  
        self.Np,self.Nt,self.Nx,self.Nw,self.Ng=self.kdata.shape

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
        f.create_dataset("x", data=self.xgrid,compression=compression)
        f.create_dataset("kcoeff", data=self.kdata,compression=compression)
        f["kcoeff"].attrs["units"] = self.kdata_unit
        f.create_dataset("samples", data=self.ggrid,compression=compression)
        f.create_dataset("weights", data=self.weights,compression=compression)
        f.create_dataset("bin_edges", data=self.wnedges,compression=compression)
        f.create_dataset("bin_centers", data=self.wns,compression=compression)
        f["bin_edges"].attrs["units"] = self.wn_unit
        f.close()    

    def read_LMDZ(self, path=None, res=None, band=None, mol=None):
        """Initializes k coeff table and supporting data from a .dat file
        in a gcm friendly format.

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
        """        
        if (path is None) or (res is None): \
            raise TypeError("You should provide an input directory name and a resolution")

        self.filename=path

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

        _, self.mol, self.Nx, self.xgrid = read_Qdat(os.path.join(path,'Q.dat'))
        if mol is not None: self.mol=mol

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
        tmp=np.loadtxt(file_to_load) \
            .reshape((self.Nt,self.Np,self.Nx,self.Nw,self.Ng+1),order='F')
        self.kdata=tmp[:,:,:,:,:-1].transpose((1,0,2,3,4))  
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
            write_only_metadata: bool, optional
                If `True`, only supporting files are written (T.dat, p.dat, etc.)
        """
        try:
            os.mkdir(path)
        except FileExistsError:
            print('Directory was already there: '+path)
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
            data_to_write=self.kdata.transpose((1,0,2,3,4)).flatten(order='F')
            data_to_write=data_to_write*u.Unit(rm_molec(self.kdata_unit)).to(u.Unit('cm^2'))
            data_to_write=np.append(data_to_write, \
                np.zeros(self.Np*self.Nt*self.Nx*self.Nw)) \
                    .reshape((1,self.Np*self.Nt*self.Nx*self.Nw*(self.Ng+1)))
            np.savetxt(os.path.join(dirname,'corrk_gcm_'+band+'.dat'),data_to_write,fmt=fmt)

    def hires_to_ktable(self, path=None, filename_grid=None,
        logpgrid=None, tgrid=None, xgrid=None, wnedges=None,
        quad='legendre', order=20, weights=None, ggrid=None,
        mid_dw=True, write=0, mol=None,
        grid_p_unit='Pa', p_unit='unspecified',
        kdata_unit='unspecified', file_kdata_unit='unspecified',
        **kwargs):
        """Computes a k coeff table from :class:`~exo_k.util.hires_spectrum.Hires_spectrum`
        objects.

        see :func:`exo_k.ktable.Ktable.hires_to_ktable` method for details
        on the arguments and options.

        Other Parameters
        ----------------
            xgrid: array
                Input grid in vmr of the variable gas. Needed for a Ktable5d.

        .. warning::
            By default, log pressures are specified in Pa in logpgrid!!! If you want
            to use another unit, do not forget to specify it with the grid_p_unit keyword.

        """        
        if path is None: raise TypeError("You should provide an input hires_spectrum directory")
        if wnedges is None: raise TypeError("You should provide an input wavenumber array")

        self.filename=path
        if mol is not None:
            self.mol=mol
        else:
            self.mol=os.path.basename(self.filename).split(self._settings._delimiter)[0]

        if weights is not None:
            self.weights=weights
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

        self.xgrid=np.array(xgrid)
        self.Nx=self.xgrid.size
        if write >= 3 : print(self.Nx,self.xgrid)

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
            for iX in range(self.Nx):
                filename=filename_grid[iP,iT,iX]
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
                self.kdata[iP,iT,iX]=spectrum_to_kdist(k_hr,wn_hr,dwn_hr,self.wnedges,self.ggrid)
                if not was_xsec:
                    self.kdata[iP,iT,iX]=self.kdata[iP,iT,iX]*KBOLTZ*self.tgrid[iT]/self.pgrid[iP]
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

    def setup_interpolation(self, log_interp=None):
        """Creates interpolating functions to be called later on.
        """
        if log_interp is None: log_interp=self._settings._log_interp
        self._local_log_interp=log_interp
        if self._local_log_interp:
            self._finterp_kdata=RegularGridInterpolator( \
                (self.logpgrid,self.tgrid,np.log(self.xgrid)), np.log(self.kdata), \
                bounds_error=True)
        else:
            self._finterp_kdata=RegularGridInterpolator( \
                (self.logpgrid,self.tgrid,np.log(self.xgrid)), self.kdata, \
                bounds_error=True)

    def set_kdata(self, new_kdata):
        """Changes kdata. this is preferred to directly accessing kdata because one
        could forget to run setup_interpolation().

        Parameters
        ----------
            new_kdata: array
                New array of kdata.
        """
        self.kdata=new_kdata
        self.setup_interpolation()

    def interpolate_kdata(self, logp_array=None, t_array=None, x_array= None,
            log_interp=None, wngrid_limit=None):
        """interpolate_kdata interpolates the kdata at on a given temperature and
        log pressure profile. 

        Parameters
        ----------
            logp_array: Array
                log 10 pressure array to interpolate to
            t_array: Array, same size a logp_array
                Temperature array to interpolate to
            x_array: Array, same size a logp_array
                vmr of variable gas array to interpolate to
            If floats are given, they are interpreted as arrays of size 1.
            wngrid_limit: list or array, optional
                if an array is given, interpolates only within this array
            log_interp: bool, dummy
                Dummy variable to be consistent with interpolate_kdata in data_table.
                Whether the interpolation is linear in kdata or in log(kdata) is actually
                controlled by self._settings._log_interp but only when the ktable is loaded.
                If you change that after the loading, you should rerun setup_interpolation().

        """
        coord_to_interp=np.array([logp_array,t_array,np.log(x_array)]).transpose()
        tmp_res=self._finterp_kdata(coord_to_interp)
        if wngrid_limit is None:
            wngrid_filter = slice(None)
        else:
            wngrid_limit=np.array(wngrid_limit)
            wngrid_filter = np.where((self.wnedges > wngrid_limit.min()) & (
                self.wnedges <= wngrid_limit.max()))[0][:-1]
        if self._local_log_interp:
            return np.exp(tmp_res[:,wngrid_filter])
        else:
            return tmp_res[:,wngrid_filter]


    def remap_logPT(self, logp_array=None, t_array=None, x_array= None):
        """remap_logPT re-interpolates the kdata on a new temprature and log pressure grid. 

        Parameters
        ----------
            logp_array: Array
                log 10 pressure array to interpolate to
            t_array: Array
                temperature array to interpolate to
            x_array: Array
                vmr of variable gas array to interpolate to

        Whether the interpolation is linear in kdata or in log(kdata)
        is controlled by self._settings._log_interp but only when the ktable is loaded.
        If you change that after the loading, you should rerun setup_interpolation().
        """
        coord=np.array(np.meshgrid(logp_array, t_array, np.log(x_array))).transpose((2,1,3,0))
        if self._local_log_interp:
            tmp_res=np.exp(self._finterp_kdata(coord))
        else:
            tmp_res=self._finterp_kdata(coord)
        self.logpgrid= logp_array
        self.pgrid   = 10**self.logpgrid
        self.tgrid   = t_array
        self.xgrid   = x_array
        self.Np      = logp_array.size
        self.Nt      = t_array.size
        self.Nx      = x_array.size
        self.set_kdata(tmp_res)
        self.setup_interpolation()

    def copy(self,cp_kdata=True):
        """Creates a new instance of :class:`Ktable5d` object and (deep) copies data into it

        Parameters
        ----------
            cp_kdata: bool, optional
                If false, the kdata table is not copied and
                only the structure and metadata are. 

        Returns
        -------
            :class:`Ktable`
                A new :class:`Ktable5d` instance with the same structure as self.
        """
        res=Ktable5d()
        res.copy_attr(self,cp_kdata=cp_kdata)
        res._local_log_interp=self._local_log_interp
        res.xgrid   = np.copy(self.xgrid)
        res.weights = np.copy(self.weights)
        res.ggrid   = np.copy(self.ggrid)
        res.gedges  = np.copy(self.gedges)
        if res.kdata is not None: res.setup_interpolation()
        return res

    def gindex(self, g):
        """Finds the index corresponding to the given g
        """
        return min(np.searchsorted(self.ggrid,g),self.Ng-1)

    def xindex(self, x):
        """Finds the index corresponding to the given x
        """
        return min(np.searchsorted(self.xgrid,x),self.Nx-1)

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
        return self.interpolate_kdata(log10(p),t,x)[0,:,gindex]

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
        toplot=self.interpolate_kdata(log10(p),t,x)[0,wlindex]
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
        x      (vmr) : {xgrid}
        data oredered following p, t, x, wn, g
        shape        : {shape}
        """.format(wg=self.weights, xgrid=self.xgrid, shape=self.shape)
        return output

    def combine_with(self, other, x_self=None, x_other=None, **kwargs):
        """Method to create a new :class:`Ktable5d` where the kdata of 'self' are
        randomly mixed with 'other' (that must be a :class:`Ktable`).

        The main purpose is to add the opacity of a trace species to the background gas
        of the :class:`Ktable5d` instance. 

        .. warning::
            Because:
            
              * the opacity from the background and variable gases cannot be
                isolated,
              * The values of the array for the vmr of the variable gas (self.xgrid)
                are not modified (diluted),

            the treatment here is valid only if `x_other` << 1.
            
            For this reason, `x_self` should be either left to None, or 1-`x_other` depending
            exactly what you want to do. But if you put `x_self`<<1, you are on your own.

        Parameters
        ----------
            other: :class:`~exo_k.ktable.Ktable`
                A :class:`Ktable` object to be mixed with. Dimensions should be the same as self
                (except for xgrid).
            x_self: float only, optional
                Volume mixing ratio of self.
            x_other: float or array, optional
                Volume mixing ratio of the species to be mixed with (other).

        If either x_self or x_other are set to `None` (default),
        the cross section of the species in question
        are considered to be already normalized with respect to the mixing ratio.

        Returns
        -------
            :class:`Ktable5d`
                A new Ktable5d with the opacity of the new species added. 
        """
        if not ((self.Np == other.Np) and (self.Nt == other.Nt) and (self.Nw == other.Nw) \
            and (self.Ng == other.Ng)):
            raise TypeError("""in combine_with: kdata tables do not have the same dimensions.
                I'll stop now!""")
        if other.Nx is not None:
            raise TypeError("""in combine_with: cannot combine 2 Ktable5d""")
        if (self.p_unit!=other.p_unit) or \
            (rm_molec(self.kdata_unit)!=rm_molec(other.kdata_unit)):
            raise RuntimeError("""in combine_with: tables do not have the same units.
                I'll stop now!""")
        res=self.copy(cp_kdata=True)
        tmp=other.copy()
        if x_other is None:
            x_other=1.
        for iX in range(self.Nx):
            tmp.kdata=self.kdata[:,:,iX,:,:]
            res.kdata[:,:,iX,:,:]= \
                tmp.RandOverlap(other, x_self, x_other*(1.-self.xgrid[iX]), **kwargs)
        res.setup_interpolation()
        return res

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

        newkdata=np.zeros((self.Np,self.Nt,self.Nx,wnedges.size-1,self.Ng), dtype=np.float64)
        for iX in range(self.Nx):
            newkdata[:,:,iX]=bin_down_corrk_numba((self.Np,self.Nt,wnedges.size-1,self.Ng), \
                self.kdata[:,:,iX], old_ggrid, self.ggrid, self.gedges, indicestosum, \
                wngrid_filter, wn_weigths, num, use_rebin)
        self.kdata=newkdata
        self.wnedges=wnedges
        self.wns=(wnedges[1:]+wnedges[:-1])*0.5
        self.Nw=self.wns.size
        if remove_zeros : self.remove_zeros(deltalog_min_value=10.)
        self.setup_interpolation()
    
    def clip_spectral_range(self, wn_range=None, wl_range=None):
        """Limits the data to the provided spectral range:

           * Wavenumber in cm^-1 if using wn_range argument
           * Wavelength in micron if using wl_range
        """
        super().clip_spectral_range(wn_range=wn_range, wl_range=wl_range)
        if self.kdata is not None:
            self.setup_interpolation()


def read_Qdat(filename):
    """Reads Q.dat files LMDZ style and extract the vmr grid.

    Parameters
    ----------
        filename: str
            Path to file to read.
    Returns
    -------
        background_mol_names: list
            list of names of molecules in background gas
        var_mol: str
            Name of variable molecule
        Nx: int
            Size of xgrid
        xgrid: array
            grid of vmr for variable gas
    """
    file = open(filename, "r")
    Nmol=int(file.readline())
    background_mol_names=[]
    for ii in range(Nmol-1):
        #print(file.readline().split()[0])
        background_mol_names.append(file.readline().split()[0])
    var_mol=file.readline().split()[0]
    Nx=int(file.readline())
    xgrid=np.zeros(Nx)
    for ii in range(Nx):
        xgrid[ii]=float(file.readline())
    return background_mol_names,var_mol,Nx,xgrid
