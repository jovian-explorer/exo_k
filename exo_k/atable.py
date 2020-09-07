# -*- coding: utf-8 -*-
"""
@author: jeremy leconte

A class to handle continuum absorption (CIA)
"""

import os.path
import h5py
import numpy as np
from .settings import Settings
from .util.interp import linear_interpolation, interp_ind_weights
from .util.spectral_object import Spectral_object
from .util.cst import PI


class Atable(Spectral_object):
    """A class to handle aerosol optical properties in table form.
    """

    def __init__(self, *filename_filters, filename=None, aerosol_name=None, search_path=None,
            mks=False, remove_zeros=False):
        """Initialization for Cia_tables.

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

        Other Parameters
        ----------------
            search_path: str, optional
                If search_path is provided,
                it locally overrides the global _search_path
                in :class:`~exo_k.settings.Settings`
                and only files in search_path are returned.
        """
        self._init_empty()
        self._settings=Settings()
        if filename is not None:
            self.filename=filename
        elif filename_filters :
            # a none empty sequence returns a True in a conditional statement
            self.filename=self._settings.list_files(*filename_filters,
                only_one=True, search_path=search_path)[0]

        if self.filename is not None:
            if self.filename.lower().endswith(('h5','hdf5')):
                self.read_hdf5(self.filename, aerosol_name=aerosol_name)
            elif self.filename.lower().endswith('.dat'):
                self.read_LMDZ(self.filename, aerosol_name=aerosol_name)
            else:
                raise RuntimeError('Aerosol optical property file extension not known.')
        if self.ext_coeff is not None:
            if self._settings._convert_to_mks or mks: self.convert_to_mks()
            if remove_zeros: self.remove_zeros()

    def _init_empty(self):
        """Initializes attributes to none.
        """
        self.aerosol_name=None
        self.wns=None
        self.wnedges=None
        self.wn_unit='cm^-1'
        self.ext_coeff=None
        self.single_scat_alb=None
        self.asymmetry_factor=None
        self.r_eff_grid=None
        self.r_eff_unit=None
        self.Nw=None
        self.Nr=None
        self.filename=None


    def read_LMDZ(self, filename, aerosol_name=None):
        """Reads LMDZ like optical propertt files.

        Parameters
        ----------
            filename: str
                Name of the file to be read.
        """
        with open(filename, 'r') as file:
            _=file.readline()
            self.Nw=int(file.readline())
            _=file.readline()
            self.Nr=int(file.readline())
            _=file.readline()
            self.wns=self._read_array(file, self.Nw)
            self.wns=.01/self.wns # conversion from m to cm^-1
            _=file.readline()
            self.r_eff_grid=self._read_array(file, self.Nr, revert=False)
            # read ext coeff
            _=file.readline()
            self.ext_coeff=self._read_arrays(file, self.Nw, self.Nr)
            # read albedo
            _=file.readline()
            self.single_scat_alb=self._read_arrays(file, self.Nw, self.Nr)
            # read asymmetry factor
            _=file.readline()
            self.asymmetry_factor=self._read_arrays(file, self.Nw, self.Nr)
        self.r_eff_unit='m'
        if aerosol_name:
            self.aerosol_name=aerosol_name
        else:
            self.aerosol_name=os.path.basename(filename).split('.')[0]




    def _read_array(self, file, Nvalue, Nline=None, revert=True):
        """Reads an array in a optical property LMDZ file. 
        Assumes that the arrays are arranged 5 values per line. 

        Parameters
        ----------
            file: file stream
                File to be read.
            Nvalue: int
                Number of values to be read. 
            revert: boolean (optional)
                Whether or not to revert the array (because most arrays are
                ordered by decreasing wavenumber).

        Returns
        -------
            Array
                A numpy array with the values.
        """
        
        if Nline is None:
            Nline=Nvalue//5
            if Nvalue%5 != 0:
                Nline+=1
        new_array=[]
        for _ in range(Nline):
            line=[float(i) for i in file.readline().split()]
            new_array+=line
        if revert :
            new_array=np.array(new_array)[::-1]
        else:
            new_array=np.array(new_array)
        return new_array

    def _read_arrays(self, file, Nvalue, Narray):
        """Reads an array in a optical property LMDZ file. 
        Assumes that the arrays are arranged 5 values per line. 

        Parameters
        ----------
            file: file stream
                File to be read.
            Nvalue: int
                Number of values to be read in each array. 
            Narray: int
                Number of arrays to be read. 

        Returns
        -------
            Array
                A numpy array with the values.
        """
        
        Nline=Nvalue//5
        if Nvalue%5 != 0:
            Nline+=1
        new_array=np.zeros((Narray,Nvalue))
        for ii in range(Narray):
            file.readline()
            new_array[ii]=self._read_array(file, Nvalue, Nline=Nline)
        return new_array

    def read_hdf5(self, filename, aerosol_name=None):
        """Reads hdf5 cia files and load temperature, wavenumber, and absorption coefficient grid.

        Parameters
        ----------
            filename: str
                Name of the file to be read.
        """
        f = h5py.File(filename, 'r')
        if aerosol_name:
            self.aerosol_name=aerosol_name
        else:
            self.aerosol_name=f['aerosol_name'][()]
        self.wns=f['bin_centers'][...]
        if 'units' in f['bin_centers'].attrs:
            self.wn_unit=f['bin_centers'].attrs['units']
        self.ext_coeff=f['ext_coeff'][...]
        self.single_scat_alb=f['single_scat_alb'][...]
        self.asymmetry_factor=f['asymmetry_factor'][...]
        self.r_eff_grid=f['r_eff'][...]
        self.r_eff_unit=f['r_eff'].attrs['units']
        f.close()  
        self.wnedges=np.concatenate(([self.wns[0]],0.5*(self.wns[1:]+self.wns[:-1]),[self.wns[-1]]))
        self.Nr=self.r_eff_grid.size
        self.Nw=self.wns.size
          
    def write_hdf5(self, filename):
        """Writes hdf5 cia files.

        Parameters
        ----------
            filename: str
                Name of the file to be written.
        """
        if not filename.lower().endswith(('.hdf5', '.h5')):
            filename=filename+'.hdf5'
        f = h5py.File(filename, 'w')
        compression="gzip"
        f.create_dataset("aerosol_name", data=self.aerosol_name)
        f.create_dataset("bin_centers", data=self.wns, compression=compression)
        f.create_dataset("r_eff", data=self.r_eff_grid, compression=compression)
        f.create_dataset("ext_coeff", data=self.ext_coeff, compression=compression)
        f.create_dataset("single_scat_alb", data=self.single_scat_alb, compression=compression)
        f.create_dataset("asymmetry_factor", data=self.asymmetry_factor, compression=compression)
        f["bin_centers"].attrs["units"] = self.wn_unit
        f["r_eff"].attrs["units"] = self.r_eff_unit
        f.close()    

    def sample(self, wngrid, remove_zeros=False, use_grid_filter=False,
            sample_all_vars=False, **kwargs):
        """Method to re sample a Atable to a new grid of wavenumbers (in place)

        Parameters
        ----------
            wngrid : array
                new wavenumber grid (cm-1)
            use_grid_filter: boolean, optional
                If true, the table is sampled only within the boundaries
                of its current wavenumber grid. The coefficients are set to zero elswere
                (except if remove_zeros is set to True).
                If false, the values at the boundaries are used when sampling outside the grid.
            sample_all_vars: boolean, optional
                Whether to sample the single_scattering albedo and asymmetry_factor as well. 
        """
        wngrid=np.array(wngrid)
        Nnew=wngrid.size
        if use_grid_filter:
            wngrid_filter = np.where((wngrid <= self.wns[-1]) & (wngrid >= self.wns[0]))[0]
        else:
            wngrid_filter = np.ones(Nnew,dtype=bool)
        new_values=np.zeros((self.Nr, Nnew))
        for iR in range(self.Nr):
            tmp=self.ext_coeff[iR,:]
            new_values[iR,wngrid_filter]=np.interp(wngrid[wngrid_filter],self.wns,tmp)
        self.ext_coeff=new_values
        if sample_all_vars :
            new_values=np.zeros((self.Nr, Nnew))
            new_values2=np.zeros((self.Nr, Nnew))
            for iR in range(self.Nr):
                tmp=self.single_scat_alb[iR,:]
                new_values[iR,wngrid_filter]=np.interp(wngrid[wngrid_filter],self.wns,tmp)
                tmp=self.asymmetry_factor[iR,:]
                new_values2[iR,wngrid_filter]=np.interp(wngrid[wngrid_filter],self.wns,tmp)
            self.single_scat_alb=new_values
            self.asymmetry_factor=new_values2
        self.wns=wngrid
        self.wnedges=np.concatenate(([self.wns[0]],0.5*(self.wns[1:]+self.wns[:-1]),[self.wns[-1]]))
        self.Nw=Nnew
        if remove_zeros : self.remove_zeros(**kwargs)

    def sample_cp(self, wngrid, **kwargs):
        """Creates a copy of the object before resampling it.

        Parameters
        ----------
            See sample method for details. 

        Returns
        -------
            :class:`Atable` object
                the re-sampled :class:`Atable`
        """
        res=self.copy()
        res.sample(wngrid, **kwargs)
        return res

    def interpolate_optical_properties(self, r_array=None, var_type=0, 
            log_interp=None, wngrid_limit=None):
        """interpolate_cia interpolates the kdata at on a given temperature profile. 

        Parameters
        ----------
            r_array: float or array
                Effective radius array to interpolate to.
                If a float is given, it is interpreted as an array of size 1.
            data_id: int
                type of data to interpolate:
                  * 0 is extinction coefficient
                  * 1 is single scattering albedo
                  * 2 is asymmetry factor
            wngrid_limit: array, optional
                If an array is given, interpolates only within this array.
            log_interp: bool, optional
                Whether the interpolation is linear in kdata or in log(kdata).
        """
        if hasattr(r_array, "__len__"):
            r_array=np.array(r_array)
        else:
            r_array=np.array([r_array])
        rind,rweight=interp_ind_weights(r_array,self.r_eff_grid)
        if wngrid_limit is None:
            wngrid_filter = slice(None)
            Nw=self.Nw
        else:
            wngrid_limit=np.array(wngrid_limit)
            wngrid_filter = np.where((self.wnedges > wngrid_limit.min()) & (
                self.wnedges <= wngrid_limit.max()))[0][:-1]
            Nw=wngrid_filter.size
        res=np.zeros((rind.size,Nw))
        if var_type == 0:
            data_to_interp=self.ext_coeff
        elif var_type ==1:
            data_to_interp=self.single_scat_alb
        else:
            data_to_interp=self.asymmetry_factor
        if log_interp is None: log_interp=self._settings._log_interp
        if log_interp:
            for ii in range(rind.size):
                kc_t1=np.log(data_to_interp[rind[ii]][wngrid_filter].ravel())
                kc_t0=np.log(data_to_interp[rind[ii]-1][wngrid_filter].ravel())
                res[ii]=linear_interpolation(kc_t0, kc_t1, rweight[ii])
            return np.exp(res)
        else:
            for ii in range(rind.size):
                kc_t1=data_to_interp[rind[ii]][wngrid_filter].ravel()
                kc_t0=data_to_interp[rind[ii]-1][wngrid_filter].ravel()
                res[ii]=linear_interpolation(kc_t0, kc_t1, rweight[ii])
            return res

    def cross_section(self, r_array, wngrid_limit=None, log_interp=None):
        """Computes the cross section due to the aerosol in area per particles.

        Parameters
        ----------
            r_array: float or array
                Effective radius array to interpolate to.
                If a float is given, it is interpreted as an array of size 1.

        Other Parameters
        ----------------
            wngrid_limit: array, optional
                If an array is given, interpolates only within this array.
            log_interp: bool, optional
                Whether the interpolation is linear in kdata or in log(kdata).
        """
        if hasattr(r_array, "__len__"):
            r_array=np.array(r_array)
        else:
            r_array=np.array([r_array])
        area=PI*r_array**2
        return (area*self.interpolate_optical_properties(r_array=r_array,
            wngrid_limit=wngrid_limit, log_interp=log_interp).transpose()).transpose()

    def absorption_coefficient(self, r_array, n_density, wngrid_limit=None, log_interp=None):
        """Computes the cross section due to the aerosol in area per particles.

        Parameters
        ----------
            r_array: float or 1d array
                Effective radius array to interpolate to.
                If a float is given, it is interpreted as an array of size 1.
            n_density: float or 1d array (same dim as r_array)
                Number density of aerosol

        Other Parameters
        ----------------
            wngrid_limit: array, optional
                If an array is given, interpolates only within this array.
            log_interp: bool, optional
                Whether the interpolation is linear in kdata or in log(kdata).
        """
        if hasattr(r_array, "__len__"):
            r_array=np.array(r_array)
        else:
            r_array=np.array([r_array])
        factor=np.array(n_density)*PI*r_array**2
        return (factor*self.interpolate_optical_properties(r_array=r_array,
            wngrid_limit=wngrid_limit, log_interp=log_interp).transpose()).transpose()

    def plot_spectrum(self, ax, r=1.e-6, x_axis='wls', xscale=None, yscale=None, 
            var_type=0, **kwarg):
        """Plot the spectrum for a given point

        Parameters
        ----------
            ax : :class:`pyplot.Axes`
                A pyplot axes instance where to put the plot.
            r: float
                Effective radius (m)
            x_axis: str, optional
                If 'wls', x axis is wavelength. Wavenumber otherwise.
            x/yscale: str, optional
                If 'log' log axes are used.
        """
        toplot=self.interpolate_optical_properties(r_array=r, var_type=var_type)[0]
        if var_type==0:
            ax.set_ylabel('Ext. coeff')
        elif var_type==1:
            ax.set_ylabel('Single Scat. Albedo')
        else:
            ax.set_ylabel('Ass. factor')
        if x_axis == 'wls':
            ax.plot(self.wls,toplot,**kwarg)
            ax.set_xlabel('Wavelength (micron)')
        else:
            ax.plot(self.wns,toplot,**kwarg)
            ax.set_xlabel('Wavenumber (cm$^{-1}$)')
        ax.grid(True)
        if xscale is not None: ax.set_xscale(xscale)
        if yscale is not None: ax.set_yscale(yscale)

    def convert_to_mks(self):
        """Converts units to MKS
        """
        pass

    def rindex(self, r):
        """Finds the index corresponding to the given radius r
        (units must be the same as the ktable)
        """
        return min(np.searchsorted(self.r_eff_grid,r),self.Nr-1)

    def remove_zeros(self, deltalog_min_value=0.):
        """Finds zeros in the ext_coeff and set them to (10.^-deltalog_min_value)
        times the minimum positive value in the table.
        This is to be able to work in logspace. 
        """
        mask = np.zeros(self.ext_coeff.shape,dtype=bool)
        mask[np.nonzero(self.ext_coeff)] = True
        minvalue=np.amin(self.ext_coeff[mask])
        self.ext_coeff[~mask]=minvalue/(10.**deltalog_min_value)
        
    def __repr__(self):
        """Method to output header
        """
        output="""
        file          : {file}
        aerosol name  : {name}
        reff grid (mu): {r}
        wavenumber    : {w}
        ext coeff     : {ext_coeff}
        albedo        : {alb}
        asymmetry     : {ass}
        """.format(file=self.filename,name=self.aerosol_name, r=self.r_eff_grid,\
            w=self.wns,ext_coeff=self.ext_coeff, \
            alb=self.single_scat_alb, ass=self.asymmetry_factor)
        return output

    def __getitem__(self, key):
        """Overrides getitem.
        """
        return self.ext_coeff[key]

    def copy(self):
        """Creates a new instance of :class:`Atable` object and (deep) copies data into it
        """
        res=Atable()
        res.aerosol_name = self.aerosol_name
        res.wns = np.copy(self.wns)
        res.wn_unit = self.wn_unit
        res.wnedges = np.copy(self.wnedges)
        res.r_eff_grid = np.copy(self.r_eff_grid)
        res.r_eff_unit = self.r_eff_unit
        res.ext_coeff = np.copy(self.ext_coeff)
        res.single_scat_alb = np.copy(self.single_scat_alb)
        res.asymmetry_factor = np.copy(self.asymmetry_factor)
        res.Nr = self.Nr
        res.Nw = self.Nw
        res.filename = self.filename
        return res

def combine_tables(aerosol_name, *atables):
    """Combine several tables representing the same aerosol type for different wavelength range.
    """
    r_eff_grid=atables[0].r_eff_grid
    order=np.argsort([atable.wns[0] for atable in atables])
    ordered_atables=np.array(atables)[order]
    #print(ordered_atables)
    res=ordered_atables[0].copy()
    r_eff_grid=res.r_eff_grid
    wns_max=res.wns[-1]
    for atable in ordered_atables[1:]:
        if np.any(atable.r_eff_grid != r_eff_grid):
            raise RuntimeError('Not the same r_eff grid for all tables')
        print(atable.wns[0], wns_max)
        if atable.wns[0]<wns_max:
            raise RuntimeError('Overlapping wavenumber range')
        else:
            wns_max=atable.wns[-1]
    res.wns=np.concatenate([atable.wns[:-1] for atable in ordered_atables])
    res.Nw=res.wns.size
    res.ext_coeff=np.concatenate([atable.ext_coeff[:,:-1] for atable in ordered_atables], axis=1)
    res.single_scat_alb=np.concatenate(
        [atable.single_scat_alb[:,:-1] for atable in ordered_atables], axis=1)
    res.asymmetry_factor=np.concatenate(
        [atable.asymmetry_factor[:,:-1] for atable in ordered_atables], axis=1)
    res.aerosol_name=aerosol_name
    return res
