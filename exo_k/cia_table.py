# -*- coding: utf-8 -*-
"""
@author: jeremy leconte

A class to handle continuum absorption (CIA)
"""

import h5py
import numpy as np
from .util.filenames import EndOfFile
from .util.interp import linear_interpolation,interp_ind_weights
from .util.cst import KBOLTZ
from .settings import Settings


class Cia_table(object):
    """A class to handle CIA opacity data tables.
    """

    def __init__(self,filename=None, mks=False, remove_zeros=False):
        self.filename=filename
        self.mol1=None
        self.mol2=None
        self.wns=None
        self.wnedges=None
        self.tgrid=None
        self.abs_coeff=None
        self.abs_coeff_unit='unspecified'
        self.Nt=None
        self.Nw=None
        self._settings=Settings()
        if filename is not None:
            if filename.lower().endswith(('h5','hdf5')):
                self.read_hdf5(filename)
            elif filename.lower().endswith('cia'):
                self.read_hitran_cia(filename)
            else:
                raise RuntimeError('Cia file extension not known.')
        if self.abs_coeff is not None:
            if self._settings._convert_to_mks or mks: self.convert_to_mks()
            if remove_zeros : self.remove_zeros()



    def read_hitran_cia(self,filename):
        """Reads hitran cia files and load temperature, wavenumber, and absorption coefficient grid.
        Parameters:
            filename: str
                Name of the file to be read.
        """
        tmp_tgrid=[]
        tmp_abs_coeff=[]
        First=True
        with open(filename, 'r') as file:
            while True:
                try:
                    Nw, Temp = self.read_header(file)
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
        self.abs_coeff_unit='cm^5'
        self.Nt=self.tgrid.size
        self.Nw=self.wns.size


    def read_header(self,file):
        """Reads the header lines in a Hitran CIA file.
        """
        line=file.readline()
        if line is None or line=='':
            raise EndOfFile
        tmp = line.split()
        self.mol1,self.mol2 = tmp[0].split('-')
        Nw = int(tmp[3])
        Temp = float(tmp[4])

        return Nw, Temp


    def read_hdf5(self,filename):
        """Reads hdf5 cia files and load temperature, wavenumber, and absorption coefficient grid.
        Parameters:
            filename: str
                Name of the file to be read.
        """
        f = h5py.File(filename, 'r')
        self.wns=f['bin_centers'][...]
        self.abs_coeff=f['abs_coeff'][...]
        self.abs_coeff_unit=f['abs_coeff'].attrs['units']
        self.tgrid=f['t'][...]
        self.mol1,self.mol2=f.attrs['cia_pair'].split('-')
        f.close()  
        self.wnedges=np.concatenate(([self.wns[0]],0.5*(self.wns[1:]+self.wns[:-1]),[self.wns[-1]]))
        self.Nt=self.tgrid.size
        self.Nw=self.wns.size
          
    def write_hdf5(self,filename):
        """Writes hdf5 cia files.
        Parameters:
            filename: str
                Name of the file to be written.
        """
        f = h5py.File(filename, 'w')
        f.create_dataset("bin_centers", data=self.wns,compression="gzip")
        f.create_dataset("t", data=self.tgrid,compression="gzip")
        f.create_dataset("abs_coeff", data=self.abs_coeff,compression="gzip")
        f["abs_coeff"].attrs["units"] = self.abs_coeff_unit
        f.attrs["cia_pair"] = self.mol1+'-'+self.mol2
        #f.create_dataset("cia_pair", data=self.mol1+'-'+self.mol2)
        f.close()    

    def sample(self,wngrid,remove_zeros=False,use_grid_filter=False,**kwargs):
        """Method to re sample a cia table to a new grid of wavenumbers
        Parameters:
            wngrid : new wavenumber grid (cm-1)
        Output    :
            Nothing returned. Directly changes the xsec self instance attributes.
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


    def interpolate_cia(self,t_array=None,log_interp=None,wngrid_limit=None):
        """interpolate_cia interpolates the kdata at on a given temperature profile. 
        Parameters:
            t_array: float or Array
                Temperature array to interpolate to.
                If a float is given, it is interpreted as an array of size 1.
        Options:
            wngrid_limit: if an array is given, interpolates only within this array
            log_interp: whether the interpolation is linear in kdata or in log(kdata)
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
                tmp=linear_interpolation(kc_t0, kc_t1, tweight[ii])
                res[ii]=np.reshape(tmp,(Nw,-1)).squeeze()
            return np.exp(res)
        else:
            for ii in range(tind.size):
                #kc_t1=self.abs_coeff[tind[ii]]
                kc_t1=self.abs_coeff[tind[ii]][wngrid_filter].ravel()
                kc_t0=self.abs_coeff[tind[ii]-1][wngrid_filter].ravel()
                tmp=linear_interpolation(kc_t0, kc_t1, tweight[ii])
                res[ii]=np.reshape(tmp,(Nw,-1)).squeeze()
            return res

    def equivalent_xsec(self, logP, T, x_mol2, wngrid_limit=None):
        """Computes the cross section due to CIA in area per molecule of type 1.
        """
        n_density=10**logP/(KBOLTZ*T)
        return self.interpolate_cia(t_array=T,wngrid_limit=wngrid_limit)*n_density*x_mol2

    def effective_cross_section(self, logP, T, x_mol1, x_mol2, wngrid_limit=None):
        """Computes the total cross section for a molecule pair
        (in m^2 per total number of molecules; assumes data in MKS).
        """
        x_x_n_density=10**logP/(KBOLTZ*T)*x_mol1*x_mol2
        #return self.interpolate_cia( \
        # T,wngrid_limit=wngrid_limit)*n_density[:,None]*n_density[:,None]*x_mol1*x_mol2
        tmp=self.interpolate_cia(t_array=T,wngrid_limit=wngrid_limit)
        return tmp*x_x_n_density[:,None]

    def plot_spectrum(self, ax, t=200., x_axis='wls', xscale=None, yscale=None, **kwarg):
        """Plot the spectrum for a given point
        Parameters:
            t: temperature(K)
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
        Parameters:
            abs_coeff_unit: str
                String to identify the units to convert to.
                Accepts 'cm^5', 'm^5'. or any length^5 unit recognized by the 
                astropy.units library. If ='unspecified', no conversion is done.
        Option:
            old_abs_coeff_unit : str
                String to specify the current kdata unit if it is unspecified or if 
                you have reasons to believe it is wrong (e.g. you just read a file where
                you know that the kdata grid and the kdata unit do not correspond)
        """
        from .util.interp import unit_convert
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
        self.convert_abs_coeff_unit(abs_coeff_unit='m^5')

    @property
    def wls(self):
        """Returns the wavelength array for the bin centers
        """
        return 10000./self.wns

    def remove_zeros(self,deltalog_min_value=0.):
        """Finds zeros in the abs_coeff and set them to (10.**-deltalog_min_value)
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
        molecule pair : {mol}
        t grid   (K)  : {t}
        abs coeff     : {abs_coeff}
        """.format(file=self.filename,mol=self.mol1+'-'+self.mol2, \
            t=self.tgrid,abs_coeff=self.abs_coeff)
        return output

    def __getitem__(self,key):
        """Overrides getitem.
        """
        return self.abs_coeff[key]
