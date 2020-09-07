# -*- coding: utf-8 -*-
"""
@author: jeremy leconte

A class to handle continuum absorption (CIA)
"""

import os.path
import h5py
import numpy as np
from exo_k.util.filenames import EndOfFile
from .util.interp import linear_interpolation, interp_ind_weights, unit_convert
from .util.cst import KBOLTZ
from .settings import Settings


class Cia_table(object):
    """A class to handle CIA opacity data tables.
    """

    def __init__(self, *filename_filters, filename=None, molecule_pair=None, search_path=None,
            mks=False, remove_zeros=False, old_cia_unit='cm^5'):
        """Initialization for Cia_tables.

        Parameters
        ----------
            filename: str, optional
                Relative or absolute name of the file to be loaded. 
            filename_filters: sequence of string
                As many strings as necessary to uniquely define
                a file in the global search path defined with
                :func:`exo_k.settings.Settings.set_cia_search_path`.
                This path will be searched for a file
                with all the filename_filters in the name.
                The filename_filters can contain '*'.
            molecule_pair: list of size 2, optional
                The molecule pair we want to consider, 
                specified as an array with two strings (like ['H2','H2'] or ['N2','H2O']).
                The order of the molecules in the pair is irrelevant. 

        Other Parameters
        ----------------
            old_cia_unit : str, optional
                String to specify the current cia unit if it is unspecified or if 
                you have reasons to believe it is wrong (e.g. you just read a file where
                you know that the cia grid and the cia unit do not correspond).
                Available units are: 'cm^5', 'cm^2' that stand for cm^2/amagat,
                and 'cm^-1' that stand for cm^-1/amagat^2.
            remove_zeros: boolean, optional
                If True, the zeros in the kdata table are replaced by
                    a value 10 orders of magnitude smaller than the smallest positive value
            search_path: str, optional
                If search_path is provided,
                it locally overrides the global search path
                defined with :func:`exo_k.settings.Settings.set_cia_search_path`
                and only files in `search_path` are returned.
        """
        self._init_empty()
        self._settings=Settings()
        if filename is not None:
            self.filename=filename
        elif filename_filters or molecule_pair is not None:
            # a none empty sequence returns a True in a conditional statement
            self.filename=self._settings.list_cia_files(*filename_filters,
                molecule_pair=molecule_pair, only_one=True, search_path=search_path)[0]

        if self.filename is not None:
            if self.filename.lower().endswith(('h5','hdf5')):
                self.read_hdf5(self.filename)
            elif self.filename.lower().endswith('cia'):
                self.read_hitran_cia(self.filename, old_cia_unit=old_cia_unit)
            elif self.filename.lower().endswith('.dat'):
                self.read_CKD_cia(self.filename, old_cia_unit=old_cia_unit)
            else:
                raise RuntimeError('Cia file extension not known.')
        if self.abs_coeff is not None:
            if self._settings._convert_to_mks or mks: self.convert_to_mks()
            if remove_zeros : self.remove_zeros()

    def _init_empty(self):
        """Initializes attributes to none.
        """
        self.mol1=None
        self.mol2=None
        self.wns=None
        self.wnedges=None
        self.tgrid=None
        self.abs_coeff=None
        self.abs_coeff_unit='unspecified'
        self.wn_unit='cm^-1'
        self.Nt=None
        self.Nw=None
        self.filename=None


    def read_hitran_cia(self, filename, old_cia_unit='cm^5'):
        """Reads hitran cia files and load temperature, wavenumber, and absorption coefficient grid.

        Parameters
        ----------
            filename: str
                Name of the file to be read.
            old_cia_unit: str, optional
                Units found in the file.
        """
        tmp_tgrid=[]
        tmp_abs_coeff=[]
        First=True
        with open(filename, 'r') as file:
            while True:
                try:
                    Nw, Temp = self._read_header(file)
                except EndOfFile:
                    break
                tmp_tgrid.append(Temp)
                tmp_abs_coeff.append([])
                tmp_wns=[]
                if First: 
                    for _ in range(Nw): #as iterator not used, can be replaced by _
                        line=file.readline()
                        tmp=line.split()
                        tmp_wns.append(float(tmp[0]))
                        tmp_abs_coeff[-1].append(float(tmp[1]))
                    self.wns=np.array(tmp_wns)
                else:
                    for _ in range(Nw):
                        line=file.readline()
                        tmp=line.split()
                        tmp_abs_coeff[-1].append(float(tmp[1]))
        self.wnedges=np.concatenate(([self.wns[0]],0.5*(self.wns[1:]+self.wns[:-1]),[self.wns[-1]]))
        self.tgrid=np.array(tmp_tgrid)
        self.abs_coeff=np.array(tmp_abs_coeff)
        self.abs_coeff_unit=old_cia_unit
        self.Nt=self.tgrid.size
        self.Nw=self.wns.size


    def _read_header(self, file):
        """Reads the header lines in a Hitran CIA file.

        Parameters
        ----------
            file: file stream
                file to be read.
        """
        line=file.readline()
        if line is None or line=='':
            raise EndOfFile
        tmp = line.split()
        self.mol1,self.mol2 = tmp[0].split('-')
        Nw = int(tmp[3])
        Temp = float(tmp[4])

        return Nw, Temp

    def read_hdf5(self, filename):
        """Reads hdf5 cia files and load temperature, wavenumber, and absorption coefficient grid.

        Parameters
        ----------
            filename: str
                Name of the file to be read.
        """
        f = h5py.File(filename, 'r')
        self.wns=f['bin_centers'][...]
        self.abs_coeff=f['abs_coeff'][...]
        self.abs_coeff_unit=f['abs_coeff'].attrs['units']
        self.tgrid=f['t'][...]
        if 'cia_pair' in f:
            self.mol1,self.mol2=f['cia_pair'][()].split('-')
        elif 'cia_pair' in f.attrs:
            self.mol1,self.mol2=f.attrs['cia_pair'].split('-')
        f.close()  
        self.wnedges=np.concatenate(([self.wns[0]],0.5*(self.wns[1:]+self.wns[:-1]),[self.wns[-1]]))
        self.Nt=self.tgrid.size
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
        f.create_dataset("bin_centers", data=self.wns,compression="gzip")
        f.create_dataset("t", data=self.tgrid,compression="gzip")
        f.create_dataset("abs_coeff", data=self.abs_coeff,compression="gzip")
        f["abs_coeff"].attrs["units"] = self.abs_coeff_unit
        f.create_dataset("cia_pair", data=self.mol1+'-'+self.mol2)
        f.close()    

    def sample(self, wngrid, remove_zeros=False, use_grid_filter=False, **kwargs):
        """Method to re sample a cia table to a new grid of wavenumbers (in place)

        Parameters
        ----------
            wngrid : array
                new wavenumber grid (cm-1)
            use_grid_filter: boolean, optional
                If true, the table is sampled only within the boundaries
                of its current wavenumber grid. The coefficients are set to zero elswere
                (except if remove_zeros is set to True).
                If false, the values at the boundaries are used when sampling outside the grid.
        """
        wngrid=np.array(wngrid)
        #min_val=np.amin(self.abs_coeff)
        Nnew=wngrid.size
        #wngrid_filter = np.where((wngrid <= self.wnedges[-1]) & (wngrid >= self.wnedges[0]))[0]
        if use_grid_filter:
            wngrid_filter = np.where((wngrid <= self.wns[-1]) & (wngrid >= self.wns[0]))[0]
        else:
            wngrid_filter = np.ones(Nnew,dtype=bool)
        new_abs_coeff=np.zeros((self.Nt,Nnew))
        for iT in range(self.Nt):
            tmp=self.abs_coeff[iT,:]
            new_abs_coeff[iT,wngrid_filter]=np.interp(wngrid[wngrid_filter],self.wns,tmp)
        self.abs_coeff=new_abs_coeff
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
            :class:`Cia_table` object
                the re-sampled :class:`Cia_table`
        """
        res=self.copy()
        res.sample(wngrid, **kwargs)
        return res

    def interpolate_cia(self, t_array=None, log_interp=None, wngrid_limit=None):
        """interpolate_cia interpolates the kdata at on a given temperature profile. 

        Parameters
        ----------
            t_array: float or array
                Temperature array to interpolate to.
                If a float is given, it is interpreted as an array of size 1.
            wngrid_limit: array, optional
                If an array is given, interpolates only within this array.
            log_interp: bool, optional
                Whether the interpolation is linear in kdata or in log(kdata).
        """
        if hasattr(t_array, "__len__"):
            t_array=np.array(t_array)
        else:
            t_array=np.array([t_array])
        tind,tweight=interp_ind_weights(t_array,self.tgrid)
        if wngrid_limit is None:
            wngrid_filter = slice(None)
            Nw=self.Nw
        else:
            wngrid_limit=np.array(wngrid_limit)
            wngrid_filter = np.where((self.wnedges > wngrid_limit.min()) & (
                self.wnedges <= wngrid_limit.max()))[0][:-1]
            Nw=wngrid_filter.size
        res=np.zeros((tind.size,Nw))
        if log_interp is None: log_interp=self._settings._log_interp
        if log_interp:
            for ii in range(tind.size):
                kc_t1=np.log(self.abs_coeff[tind[ii]][wngrid_filter].ravel())
                kc_t0=np.log(self.abs_coeff[tind[ii]-1][wngrid_filter].ravel())
                res[ii]=linear_interpolation(kc_t0, kc_t1, tweight[ii])
            return np.exp(res)
        else:
            for ii in range(tind.size):
                #kc_t1=self.abs_coeff[tind[ii]]
                kc_t1=self.abs_coeff[tind[ii]][wngrid_filter].ravel()
                kc_t0=self.abs_coeff[tind[ii]-1][wngrid_filter].ravel()
                res[ii]=linear_interpolation(kc_t0, kc_t1, tweight[ii])
            return res

    def equivalent_xsec(self, logP, T, x_mol2, wngrid_limit=None):
        """Computes the cross section due to CIA in area per molecule of type 1.
        """
        n_density=10**logP/(KBOLTZ*T)
        return self.interpolate_cia(t_array=T,wngrid_limit=wngrid_limit)*n_density*x_mol2

    def effective_cross_section(self, logP, T, x_mol1, x_mol2, wngrid_limit=None):
        """Computes the total cross section for a molecule pair
        (in m^2 per total number of molecules; assumes data in MKS).

        Parameters
        ----------
            logP: float or array
                Log10 of the pressure (Pa).
            T: float or array
                Temperature (K).
            x_mol1/2: float or array
                Volume mixing ratio of the 1st and 2nd molecule of the pair.
            wngrid_limit: array, optional  
                If an array is given, interpolates only within this array. 

        Returns
        -------
            float or array
                total cross section for the molecule pair in m^2 per total number of molecules.
        """
        x_x_n_density=10**logP/(KBOLTZ*T)*x_mol1*x_mol2
        #return self.interpolate_cia( \
        # T,wngrid_limit=wngrid_limit)*n_density[:,None]*n_density[:,None]*x_mol1*x_mol2
        tmp=self.interpolate_cia(t_array=T, wngrid_limit=wngrid_limit)
        return (x_x_n_density*tmp.transpose()).transpose() 
        # trick for the broadcasting to work whether x_x_n_density is a float or an array

    def plot_spectrum(self, ax, t=200., x_axis='wls', xscale=None, yscale=None, **kwarg):
        """Plot the spectrum for a given point

        Parameters
        ----------
            ax : :class:`pyplot.Axes`
                A pyplot axes instance where to put the plot.
            t: float
                temperature(K)
            x_axis: str, optional
                If 'wls', x axis is wavelength. Wavenumber otherwise.
            x/yscale: str, optional
                If 'log' log axes are used.
        """
        toplot=self.interpolate_cia(t)[0]
        if x_axis == 'wls':
            ax.plot(self.wls,toplot,**kwarg)
            ax.set_xlabel('Wavelength (micron)')
        else:
            ax.plot(self.wns,toplot,**kwarg)
            ax.set_xlabel('Wavenumber (cm$^{-1}$)')
        ax.set_ylabel('Abs. coeff')
        ax.grid(True)
        if xscale is not None: ax.set_xscale(xscale)
        if yscale is not None: ax.set_yscale(yscale)

    def convert_abs_coeff_unit(self,abs_coeff_unit='unspecified',old_abs_coeff_unit='unspecified'):
        """Converts abs_coeff to a new unit (in place)

        Parameters
        ----------
            abs_coeff_unit: str
                String to identify the units to convert to.
                Accepts 'cm^5', 'm^5'. or any length^5 unit recognized by the 
                astropy.units library. If ='unspecified', no conversion is done.
            old_abs_coeff_unit : str, optional
                String to specify the current kdata unit if it is unspecified or if 
                you have reasons to believe it is wrong (e.g. you just read a file where
                you know that the kdata grid and the kdata unit do not correspond)
        """
        if abs_coeff_unit==old_abs_coeff_unit: return
        tmp_k_u_in=old_abs_coeff_unit
        tmp_k_u_out=abs_coeff_unit
        tmp_k_u_file=self.abs_coeff_unit
        self.abs_coeff_unit,conversion_factor=unit_convert( \
            'abs_coeff_unit',unit_file=tmp_k_u_file,unit_in=tmp_k_u_in,unit_out=tmp_k_u_out)
        self.abs_coeff=self.abs_coeff*conversion_factor

    def convert_to_mks(self):
        """Converts units to MKS
        """
        if self.abs_coeff_unit=='cm^5':
            print('Conversion from cm^5 to m^5')
            self.convert_abs_coeff_unit(abs_coeff_unit='m^5')
        elif self.abs_coeff_unit in ['cm^2','m^2']:
            print('Conversion from '+self.abs_coeff_unit+'/amagat to m^5')
            self.convert_abs_coeff_unit(abs_coeff_unit='m^2')
            self.abs_coeff=self.abs_coeff*(KBOLTZ*273.15/101325.0)
            #conversion from m^2/amagat to m^5
        elif self.abs_coeff_unit in ['cm^-1','m^-1']:
            print('Conversion from '+self.abs_coeff_unit+'/amagat^2 to m^5')
            self.convert_abs_coeff_unit(abs_coeff_unit='m^-1')
            self.abs_coeff=self.abs_coeff*(KBOLTZ*273.15/101325.0)**2
            #conversion from m^-1/amagat^2 to m^5
        else:
            pass
        self.abs_coeff_unit='m^5'
        return

    @property
    def wls(self):
        """Returns the wavelength array for the bin centers
        """
        return 10000./self.wns

    def remove_zeros(self, deltalog_min_value=0.):
        """Finds zeros in the abs_coeff and set them to (10.^-deltalog_min_value)
        times the minimum positive value in the table.
        This is to be able to work in logspace. 
        """
        mask = np.zeros(self.abs_coeff.shape,dtype=bool)
        mask[np.nonzero(self.abs_coeff)] = True
        minvalue=np.amin(self.abs_coeff[mask])
        self.abs_coeff[~mask]=minvalue/(10.**deltalog_min_value)
        
    def __repr__(self):
        """Method to output header
        """
        output="""
        file          : {file}
        molecule pair : {mol1} - {mol2}
        t grid   (K)  : {t}
        abs coeff     : {abs_coeff}
        """.format(file=self.filename,mol1=self.mol1, mol2=self.mol2, \
            t=self.tgrid,abs_coeff=self.abs_coeff)
        return output

    def __getitem__(self,key):
        """Overrides getitem.
        """
        return self.abs_coeff[key]

    def copy(self):
        """Creates a new instance of :class:`CIA_table` object and (deep) copies data into it
        """
        res=Cia_table()
        res.mol1 = self.mol1
        res.mol2 = self.mol2
        res.wns = np.copy(self.wns)
        res.wnedges = np.copy(self.wnedges)
        res.tgrid = np.copy(self.tgrid)
        res.abs_coeff = np.copy(self.abs_coeff)
        res.abs_coeff_unit = self.abs_coeff_unit
        res.Nt = self.Nt
        res.Nw = self.Nw
        res.filename = self.filename
        res._settings=Settings()
        return res

    def read_CKD_cia(self, filename, old_cia_unit='cm^2'):
        """Reads hitran cia files and load temperature, wavenumber, and absorption coefficient grid.

        Parameters
        ----------
            filename: str
                Name of the file to be read.
        """
        self.tgrid=np.array([200., 250., 300., 350., 400.,
            450., 500., 550., 600., 650., 700.])
        self.Nt=self.tgrid.size
        self.abs_coeff=np.loadtxt(filename, skiprows=1, unpack=True)
        self.mol1='H2O'
        if 'SELF' in os.path.basename(filename):
            self.mol2='H2O'
        else:
            self.mol2='N2'
        nu_name=os.path.join(os.path.dirname(filename),'H2O_CONT_NU.dat')
        self.wns=np.loadtxt(nu_name, skiprows=1, unpack=True)
        self.wnedges=np.concatenate(([self.wns[0]],0.5*(self.wns[1:]+self.wns[:-1]),[self.wns[-1]]))
        self.abs_coeff_unit=old_cia_unit
        self.Nw=self.wns.size

#    #  amagatS=(273.15/temp)*(presS/101325.0)
#    #  amagatF=(273.15/temp)*(presF/101325.0)
#
#    #  abcoef = abcoefS*amagatS + abcoefF*amagatF ! Eq. (15) in Clough (1989)
#    #  abcoef = abcoef*(presS/(presF+presS))      ! take H2O mixing ratio into account
#    #                                                ! abs coeffs are given per molecule of H2O
#
#    #  Nmolec = (presS+presF)/(kB*temp)           ! assume ideal gas
#!      print*,'Total number of molecules per m^3 is',Nmolec
#
#      abcoef = (abcoefS*x_H2O + abcoefF*(1-x_other))*(k_b*273.15/101325.0)*x_H2O
#                     *Nmolec*2/(100.0**2)          ! convert to m^-1

    def effective_cross_section2(self, logP, T, x_mol1, x_mol2, wngrid_limit=None):
        """Obsolete.
        
        Computes the total cross section for a molecule pair
        (in m^2 per total number of molecules; assumes data in MKS).
        """
        x_x_n_density=10**logP/(KBOLTZ*T)*x_mol1*x_mol2
        #return self.interpolate_cia( \
        # T,wngrid_limit=wngrid_limit)*n_density[:,None]*n_density[:,None]*x_mol1*x_mol2
        tmp=self.interpolate_cia(t_array=T, wngrid_limit=wngrid_limit)
        return x_x_n_density[:,None]*tmp

