# -*- coding: utf-8 -*-
"""
@author: jeremy leconte
"""
import numpy as np
from .util.interp import rm_molec,unit_convert,interp_ind_weights,bilinear_interpolation
from .settings import Settings

class Data_table(object):
    """An abstract class that will serve as a basis for Ktable and Xtable.
    This class includes all the interpolation and remapping methods.
    """

    __array_priority__=1000
    # this allows __rmul__ to take precedence over numpy array broadcasting in multiplications
    # See https://stackoverflow.com/questions/40694380/forcing-multiplication-to-use-rmul-instead-of-numpy-array-mul-or-byp/44634634#44634634

    def __init__(self):
        """Initializes all attributes to `None`"""
        self.filename=None
        self.mol=None
        self.pgrid=None
        self.logpgrid=None
        self.tgrid=None
        self.wns=None
        self.wnedges=None
        self.kdata=None
        self.logk=None
        self.Np=None
        self.Nt=None
        self.Nw=None
        self.Ng=None   # If the table is for xsec, Ng will stay at None.
                       # This will allow us to differentiate xsec from corrk
                       # when needed in the interpolation routines.
                       # Especially to reshape the output.
        self.Nx=None   # If we are dealing with a Qtable with a variable gas, Nx is the size of the
                       # grid along the dimension of the Volume mixing ratio of the variable gas.
                       # A None value means that we have either a regular Ktable or a Xtable.
        self.p_unit='unspecified'
        self.kdata_unit='unspecified'
        self.wn_unit='cm^-1'
        self._settings=Settings()

    def finalize_init(self, p_unit='unspecified', file_p_unit='unspecified',
        kdata_unit='unspecified', file_kdata_unit='unspecified',
        remove_zeros=False):
        """Common code at the end of the initialization of
        inheriting classes put here to avoid duplicates
        """
        if self.kdata is not None:
            if self._settings._convert_to_mks:
                if p_unit is 'unspecified': p_unit='Pa'
                if kdata_unit is 'unspecified': kdata_unit='m^2/molecule'
            self.convert_p_unit(p_unit=p_unit,file_p_unit=file_p_unit)
            self.convert_kdata_unit(kdata_unit=kdata_unit,file_kdata_unit=file_kdata_unit)
            if remove_zeros : self.remove_zeros(deltalog_min_value=10.)

    def copy_attr(self, other, cp_kdata=False):
        """Copy attributes from other

        Parameters
        ----------
            other: :class:`Data_table`
                :class:`Data_table` object that will be copied
            cp_kdata: bool, optional
                If `False`, only metadata are copied
        """
        self.filename=other.filename
        self.mol     =other.mol
        self.pgrid   =np.copy(other.pgrid)
        self.logpgrid=np.copy(other.logpgrid)
        self.tgrid   =np.copy(other.tgrid)
        self.wns     =np.copy(other.wns)
        self.wnedges =np.copy(other.wnedges)
        self.Np      =other.Np
        self.Nt      =other.Nt
        self.Nw      =other.Nw
        self.Ng      =other.Ng
        self.Nx      =other.Nx
        self.logk    =other.logk
        self.p_unit  =other.p_unit
        self.kdata_unit= other.kdata_unit
        if cp_kdata:
          self.kdata=np.copy(other.kdata)
        else:
          self.kdata=None

    def remove_zeros(self,deltalog_min_value=10.):
        """Finds zeros in the kdata and set them to (10.**-deltalog_min_value)
        times the minimum positive value in the table.

        This is to be able to work in logspace. 

        Parameters
        ----------
            deltalog_min_value: float        
        """
        mask = np.zeros(self.kdata.shape,dtype=bool)
        mask[np.nonzero(self.kdata)] = True
        minvalue=np.amin(self.kdata[mask])
        self.kdata[~mask]=minvalue/(10.**deltalog_min_value)

    def convert_p_unit(self, p_unit='unspecified', file_p_unit='unspecified'):
        """Converts pressure to a new unit (in place)

        Parameters
        ----------
            p_unit: str
                String identifying the pressure units to convert to (e.g. 'bar', 'Pa', 'mbar', 
                or any pressure unit recognized by the astropy.units library).
                If ='unspecified', no conversion is done.
            file_p_unit : str, optional
                String to specify the current pressure unit if it is unspecified or if 
                you have reasons to believe it is wrong (e.g. you just read a file where
                you know that the pressure grid and the pressure unit do not correspond)
        """
        #if p_unit==file_p_unit and self.p_unit != 'unspecified': return
        current_p_unit=self.p_unit
        self.p_unit,conversion_factor=unit_convert( \
            'p_unit',unit_file=current_p_unit,unit_in=file_p_unit,unit_out=p_unit)
        self.pgrid=self.pgrid*conversion_factor
        self.logpgrid=np.log10(self.pgrid)

    def convert_kdata_unit(self, kdata_unit='unspecified', file_kdata_unit='unspecified'):
        """Converts kdata to a new unit (in place)

        Parameters
        ----------
            kdata_unit: str
                String to identify the units to convert to. Accepts 'cm^2', 'm^2'
                or any surface unit recognized by the astropy.units library.
                If ='unspecified', no conversion is done.
                In general, kdata should be kept in 'per number' or 'per volume'
                units (as opposed to 'per mass' units) as composition will
                always be assumed to be a number or volume mixing ratio.
                Opacities per unit mass are not supported yet.
                Note that you do not need to specify the '/molec' or
                '/molecule' in the unit.
            file_kdata_unit : str
                String to specify the current kdata unit if it is unspecified or if 
                you have reasons to believe it is wrong (e.g. you just read a file where
                you know that the kdata grid and the kdata unit do not correspond)
        """
        #if kdata_unit==file_kdata_unit and self.kdata_unit != 'unspecified': return
        tmp_k_u_in=rm_molec(file_kdata_unit)
        tmp_k_u_out=rm_molec(kdata_unit)
        tmp_k_u_file=rm_molec(self.kdata_unit)
        self.kdata_unit,conversion_factor=unit_convert(  \
            'kdata_unit',unit_file=tmp_k_u_file,unit_in=tmp_k_u_in,unit_out=tmp_k_u_out)
        self.kdata_unit=self.kdata_unit+'/molecule'
        self.kdata=self.kdata*conversion_factor
    
    def convert_to_mks(self):
        """Converts units to MKS
        """
        self.convert_kdata_unit(kdata_unit='m^2')
        self.convert_p_unit(p_unit='Pa')

    def interpolate_kdata(self, logp_array=None, t_array=None, x_array=1.,
            log_interp=None, wngrid_limit=None):
        """interpolate_kdata interpolates the kdata at on a given temperature and
        log pressure profile. If a volume mixing ratio profile (`x_array`) is given,
        the cross section computed for the species is multiplied by `x_array` to account for the
        'dilution' of the opacity.

        Parameters
        ----------
            logp_array: array
                log 10 pressure array to interpolate to
            t_array: array, same size a logp_array
                Temperature array to interpolate to
            x_array: None
                Volume mixing ratio array used to renormalize the cross section.

        If floats are given, they are interpreted as arrays of size 1.

        Other Parameters
        ----------------
            wngrid_limit: list or array
                if an array of to values is given, interpolates only within this array
            log_interp: bool, optional
                Whether the interpolation is linear in kdata or in log(kdata)
        """
        if hasattr(logp_array, "__len__"):
            logp_array=np.array(logp_array, dtype=float)
        else:
            logp_array=np.array([logp_array], dtype=float)
        if hasattr(t_array, "__len__"):
            t_array=np.array(t_array)
        else:
            t_array=np.array([t_array])
        if hasattr(x_array, "__len__"):
            x_array=np.array(x_array)
        else:
            x_array=np.array([x_array])
        tind,tweight=interp_ind_weights(t_array,self.tgrid)
        lpind,lpweight=interp_ind_weights(logp_array,self.logpgrid)
        if wngrid_limit is None:
            wngrid_filter = slice(None)
            Nw=self.Nw
        else:
            wngrid_limit=np.array(wngrid_limit)
            wngrid_filter = np.where((self.wnedges > wngrid_limit.min()) & (
                self.wnedges <= wngrid_limit.max()))[0][:-1]
            Nw=wngrid_filter.size
        if self.Ng is None:
            res=np.zeros((tind.size,Nw))
        else:
            res=np.zeros((tind.size,Nw,self.Ng))
        if log_interp is None: log_interp=self._settings._log_interp
        if log_interp:
            for ii in range(tind.size):
                kc_p1t1=np.log(self.kdata[lpind[ii],tind[ii]][wngrid_filter].ravel())
                kc_p0t1=np.log(self.kdata[lpind[ii]-1,tind[ii]][wngrid_filter].ravel())
                kc_p1t0=np.log(self.kdata[lpind[ii],tind[ii]-1][wngrid_filter].ravel())
                kc_p0t0=np.log(self.kdata[lpind[ii]-1,tind[ii]-1][wngrid_filter].ravel())
                res[ii]=np.reshape(bilinear_interpolation(kc_p0t0, kc_p1t0, 
                    kc_p0t1, kc_p1t1, lpweight[ii], tweight[ii]),(Nw,-1)).squeeze()
            return (x_array*np.exp(res).transpose()).transpose()
            # trick for the broadcasting to work whatever the shape of x_array
        else:
            for ii in range(tind.size):
                kc_p1t1=self.kdata[lpind[ii],tind[ii]][wngrid_filter].ravel()
                kc_p0t1=self.kdata[lpind[ii]-1,tind[ii]][wngrid_filter].ravel()
                kc_p1t0=self.kdata[lpind[ii],tind[ii]-1][wngrid_filter].ravel()
                kc_p0t0=self.kdata[lpind[ii]-1,tind[ii]-1][wngrid_filter].ravel()
                res[ii]=np.reshape(bilinear_interpolation(kc_p0t0, kc_p1t0, 
                    kc_p0t1, kc_p1t1, lpweight[ii], tweight[ii]),(Nw,-1)).squeeze()
            return (x_array*res.transpose()).transpose()
            # trick for the broadcasting to work whatever the shape of x_array

    def remap_logPT(self, logp_array=None, t_array=None, x_array=None):
        """remap_logPT re-interpolates the kdata on a new temprature and log pressure grid. 

        Parameters
        ----------
            logp_array: Array
                log 10 pressure array to interpolate to
            t_array: Array
                temperature array to interpolate to
            x_array: dummy argument to be consistent with interpolate_kdata in Ktable5d

        Whether the interpolation is linear in kdata or in log10(kdata)
        is controlled by self._settings._log_interp
        """
        if x_array is not None: print('be careful, providing an x_array is usually for Ktable5d')
        t_array=np.array(t_array)
        logp_array=np.array(logp_array, dtype=float)
        tind,tweight=interp_ind_weights(t_array,self.tgrid)
        lpind,lpweight=interp_ind_weights(logp_array,self.logpgrid)
        lpindextended=lpind[:,None]
        if self.Ng is None:
            tw=tweight[None,:,None]    # trick to broadcast over Nw and Ng a few lines below
            pw=lpweight[:,None,None]
        else:
            tw=tweight[None,:,None,None]    # trick to broadcast over Nw and Ng a few lines below
            pw=lpweight[:,None,None,None]
        #tw=tweight.reshape((1,tweight.size,1,1))  
        ## trick to broadcast over Nw and Ng a few lines below
        #pw=lpweight.reshape((lpweight.size,1,1,1))
        kc_p1t1=self.kdata[lpindextended,tind]
        kc_p0t1=self.kdata[lpindextended-1,tind]
        kc_p1t0=self.kdata[lpindextended,tind-1]
        kc_p0t0=self.kdata[lpindextended-1,tind-1]
        if self._settings._log_interp is True:
            kdata_tmp=  np.log10(kc_p1t1)*pw*tw          \
                            +np.log10(kc_p0t1)*(1.-pw)*tw \
                            +np.log10(kc_p1t0)*pw*(1.-tw) \
                            +np.log10(kc_p0t0)*(1.-pw)*(1.-tw)
            self.kdata=10**kdata_tmp
        else:
            kdata_tmp=  (kc_p1t1)*pw*tw          \
                        +(kc_p0t1)*(1.-pw)*tw \
                        +(kc_p1t0)*pw*(1.-tw) \
                        +(kc_p0t0)*(1.-pw)*(1.-tw)
            self.kdata=kdata_tmp
        self.logpgrid=logp_array
        self.pgrid   =10**self.logpgrid
        self.tgrid   =t_array
        self.Np      =logp_array.size
        self.Nt      =t_array.size

    def pindex(self, p):
        """Finds the index corresponding to the given pressure p
        (units must be the same as the ktable)
        """
        return min(np.searchsorted(self.pgrid,p),self.Np-1)

    def tindex(self, t):
        """Finds the index corresponding to the given temperature t (in K)
        """
        return min(np.searchsorted(self.tgrid,t),self.Nt-1)

    def wlindex(self, wl):
        """Finds the index corresponding to the given wavelength (in microns)
        """
        return min(np.searchsorted(self.wns,10000./wl),self.Nw-1)-1
        #return min(np.searchsorted(self.wnedges,10000./wl),self.Nw-1)-1

    def __repr__(self):
        """Method to output header
        """
        output="""
        file         : {file}
        molecule     : {mol}
        p grid       : {p}
        p unit       : {p_unit}
        t grid   (K) : {t}
        wn grid      : {wn}
        wn unit      : {wnu}
        kdata unit   : {kdata_unit}""".format(file=self.filename,mol=self.mol,
            p=self.pgrid, p_unit=self.p_unit, t=self.tgrid,
            wn=self.wns, wnu=self.wn_unit,
            kdata_unit=self.kdata_unit)
        return output

    def plot_spectrum(self, ax, p=1.e-5, t=200., x=1., g=None,
            x_axis='wls', xscale=None, yscale=None, **kwarg):
        """Plot the spectrum for a given point

        Parameters
        ----------
            ax : :class:`pyplot.Axes`
                A pyplot axes instance where to put the plot.
            p : float
                Pressure (Ktable pressure unit)
            t : float
                Temperature (K)
            g: float
                Gauss point
            x: float
                Mixing ratio of the species

        Other Parameters
        ----------------
            x_axis: str, optional
                If 'wls', x axis is wavelength. Wavenumber otherwise.
            x/yscale: str, optional
                If 'log' log axes are used.
        """
        toplot=self.spectrum_to_plot(p=p, t=t, x=x, g=g)
        if x_axis == 'wls':
            ax.plot(self.wls,toplot,**kwarg)
            ax.set_xlabel('Wavelength (micron)')
        else:
            ax.plot(self.wns,toplot,**kwarg)
            ax.set_xlabel('Wavenumber (cm$^{-1}$)')
        ax.set_ylabel('Cross section ('+self.kdata_unit+')')
        ax.grid(True)
        if xscale is not None: ax.set_xscale(xscale)
        if yscale is not None: ax.set_yscale(yscale)

    def spectrum_to_plot(self, p=1.e-5, t=200., x=1., g=None):
        """Dummy function to be defined in inheriting classes
        """
        raise NotImplementedError()

    def vmr_normalize(self, x_self):
        """Rescales kdata to account for the fact that the gas is not a pure species

        Parameters
        ----------
            x_self: float or array of shape (`self.Np,self.Nt`)
                The volume mixing ratio of the species.
        
        Returns
        -------
            array
                The vmr normalized kdata (x_self*self.kdata)
        """
        if x_self is None : return self.kdata
        if isinstance(x_self, (float,int)): return x_self*self.kdata
        if not isinstance(x_self,np.ndarray):
            print("""in vmr_normalize:
            x_self should be a float or a numpy array: I'll probably stop now!""")
            raise TypeError('bad mixing ratio type')
        if np.array_equal(x_self.shape,self.kdata.shape[0:2]):
            if self.Ng is None:
                return x_self[:,:,None]*self.kdata
            else:                
                return x_self[:,:,None,None]*self.kdata
        else:
            print("""in vmr_normalize:
            x_self shape should be (pgrid.size,tgrid.size): I'll stop now!""")
            raise TypeError('bad mixing ratio type')

    def __rmul__(self, vmr):
        """Defines the right "*" operator with a vmr
        """
        res=self.copy()
        res.kdata=self.vmr_normalize(vmr)
        return res

    __mul__ = __rmul__

    def combine_with(self, other, x_self=None, x_other=None, **kwargs):
        """Method to create a new `Data_table` where the kdata of 'self' are:

            * randomly mixed with 'other' for a `Ktable`.
            * added to 'other' for an `Xtable`

        Parameters
        ----------
            other : same class as self
                A :class:`Data_table` object to be mixed with.
                Dimensions should be the same as self.
            x_self : float or array, optional
                Volume mixing ratio of self.
            x_other : float or array, optional
                Volume mixing ratio of the species to be mixed with (other).

        If either x_self or x_other are set to `None` (default),
        the cross section of the species in question
        are considered to be already normalized with respect to the mixing ratio.

        Returns
        -------
            :class:`Data_table`
                A new table for the mix
        """
        if other.Nx is not None:
            # if other is a Ktable5d, use the method for this class instead. 
            return other.combine_with(self, x_self=x_other, x_other=x_self, **kwargs)
        if not np.array_equal(self.shape,other.shape):
            raise TypeError("""in combine_with: kdata tables do not have the same dimensions.
                I'll stop now!""")
        if (self.p_unit!=other.p_unit) or \
            (rm_molec(self.kdata_unit)!=rm_molec(other.kdata_unit)):
            raise RuntimeError("""in combine_with: tables do not have the same units.
                I'll stop now!""")
        res=self.copy(cp_kdata=False)
        if self.Ng is None:
            res.kdata=self.vmr_normalize(x_self) + other.vmr_normalize(x_other)
        else:
            res.kdata=self.RandOverlap(other, x_self, x_other, **kwargs)
        return res

    def __add__(self, other):
        """Defines the "+" operator with another Data_table

        .. warning::
            __radd__ is not implmented because we want to use the __add__ (or combine_with)
            method of the left object.
        """
        return self.combine_with(other)

    @property
    def wls(self):
        """Returns the wavelength array for the bin centers
        """
        if self.wns is not None: return 10000./self.wns

    @property
    def wledges(self):
        """Returns the wavelength array for the bin edges
        """
        if self.wnedges is not None: return 10000./self.wnedges

    def __getitem__(self, key):
        """To access the data without typing self.kdata[]

        Parameters
        ----------
            key: can be slices, like for a numpy array.
        """
        return self.kdata[key]

    def clip_spectral_range(self, wn_range=None, wl_range=None):
        """Limits the data to the provided spectral range:

           * Wavenumber in cm^-1 if using wn_range argument
           * Wavelength in micron if using wl_range
        """
        if (wn_range is None) and (wl_range is None): return
        if wl_range is not None:
            if wn_range is not None:
                raise RuntimeError('Should provide either wn_range or wl_range, not both!')
            _wn_range=np.sort(10000./np.array(wl_range))
        else:
            _wn_range=np.sort(np.array(wn_range))
        iw_min, iw_max=np.searchsorted(self.wnedges, _wn_range, side='left')
        iw_max-=1
        self.wnedges=self.wnedges[iw_min:iw_max+1]
        self.wns=self.wns[iw_min:iw_max]
        self.Nw=self.wns.size
        if self.Nx is None:
            self.kdata=self.kdata[:,:,iw_min:iw_max]
        else:
            self.kdata=self.kdata[:,:,:,iw_min:iw_max]

    def bin_down_cp(self, wnedges=None, **kwargs):
        """Creates a copy of the instance and bins it down using the methods in 
        Ktable or Xtable.

        See :func:`exo_k.ktable.Ktable.bin_down` or :func:`exo_k.xtable.Xtable.bin_down`
        for details on parameters.

        Returns
        -------
            :class:`~exo_k.ktable.Ktable` or :class:`~exo_k.xtable.Xtable` object
                the binned down table
        """
        res=self.copy()
        res.bin_down(wnedges=wnedges, **kwargs)
        return res

    def toLogK(self):
        """Changes kdata to log 10.
        """
        if not self.logk:
            self.logk=True
            self.kdata=np.log10(self.kdata)
        return
            
    def toLinK(self):
        """Changes kdata back from log to linear scale.
        """
        if self.logk:
            self.logk=False
            self.kdata=np.power(10.,self.kdata)
        return

############# OLD METHODS ################

    def interpolate_kdata2(self, logp_array=None, t_array=None, x_array=None,
            log_interp=None,wngrid_limit=None):
        """Deprecated, should NOT be used.

        interpolate_kdata interpolates the kdata at on a given temperature and
        log pressure profile. 

        Parameters
        ----------
            logp_array: array
                log 10 pressure array to interpolate to
            t_array: array, same size a logp_array
                Temperature array to interpolate to
        If floats are given, they are interpreted as arrays of size 1.

        Parameters
        ----------
            x_array: None
                Dummy argument to be consistent with interpolate_kdata in Ktable5d
            wngrid_limit: list or array
                if an array of to values is given, interpolates only within this array
            log_interp: bool, optional
                Whether the interpolation is linear in kdata or in log(kdata)
        """
        if x_array is not None: print('be careful, providing an x_array is usually for Ktable5d')
        if hasattr(logp_array, "__len__"):
            logp_array=np.array(logp_array, dtype=float)
        else:
            logp_array=np.array([logp_array], dtype=float)
        if hasattr(t_array, "__len__"):
            t_array=np.array(t_array)
        else:
            t_array=np.array([t_array])
        tind,tweight=interp_ind_weights(t_array,self.tgrid)
        lpind,lpweight=interp_ind_weights(logp_array,self.logpgrid)
        if wngrid_limit is None:
            wngrid_filter = slice(None)
            Nw=self.Nw
        else:
            wngrid_limit=np.array(wngrid_limit)
            wngrid_filter = np.where((self.wnedges > wngrid_limit.min()) & (
                self.wnedges <= wngrid_limit.max()))[0][:-1]
            Nw=wngrid_filter.size
        if self.Ng is None:
            res=np.zeros((tind.size,Nw))
        else:
            res=np.zeros((tind.size,Nw,self.Ng))
        if log_interp is None: log_interp=self._settings._log_interp
        if log_interp:
            for ii in range(tind.size):
                kc_p1t1=np.log(self.kdata[lpind[ii],tind[ii]][wngrid_filter].ravel())
                kc_p0t1=np.log(self.kdata[lpind[ii]-1,tind[ii]][wngrid_filter].ravel())
                kc_p1t0=np.log(self.kdata[lpind[ii],tind[ii]-1][wngrid_filter].ravel())
                kc_p0t0=np.log(self.kdata[lpind[ii]-1,tind[ii]-1][wngrid_filter].ravel())
                res[ii]=np.reshape(bilinear_interpolation(kc_p0t0, kc_p1t0, 
                    kc_p0t1, kc_p1t1, lpweight[ii], tweight[ii]),(Nw,-1)).squeeze()
            return np.exp(res)
        else:
            for ii in range(tind.size):
                kc_p1t1=self.kdata[lpind[ii],tind[ii]][wngrid_filter].ravel()
                kc_p0t1=self.kdata[lpind[ii]-1,tind[ii]][wngrid_filter].ravel()
                kc_p1t0=self.kdata[lpind[ii],tind[ii]-1][wngrid_filter].ravel()
                kc_p0t0=self.kdata[lpind[ii]-1,tind[ii]-1][wngrid_filter].ravel()
                res[ii]=np.reshape(bilinear_interpolation(kc_p0t0, kc_p1t0, 
                    kc_p0t1, kc_p1t1, lpweight[ii], tweight[ii]),(Nw,-1)).squeeze()
            return res
