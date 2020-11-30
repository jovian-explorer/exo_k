# -*- coding: utf-8 -*-
"""
@author: jeremy leconte
"""
import os
from math import log10
import time
import numpy as np
import astropy.units as u
from .ktable_io import Ktable_io
from .ktable5d import Ktable5d
from .util.interp import rebin_ind_weights, rebin, is_sorted, \
        gauss_legendre, split_gauss_legendre, spectrum_to_kdist, \
        kdata_conv_loop, bin_down_corrk_numba, g_sample_4d
from .hires_spectrum import Hires_spectrum
from .util.filenames import create_fname_grid_Kspectrum_LMDZ, select_kwargs
from .util.cst import KBOLTZ


class Ktable(Ktable_io):
    """A class that handles 4D tables of k-coefficients.
    
    Based on the :class:`~exo_k.data_table.Data_table` class
    that handles basic operations common to Xtable.

    Based on the :class:`~exo_k.ktable_io.Ktable_io` class
    that incorporates all io routines (read/write_xxx for all the xxx supported formats).
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
                only_one=True, search_path=search_path, path_type='ktable')[0]

        if self.filename is not None:
            if self.filename.lower().endswith('pickle'):
                self.read_pickle(filename=self.filename)
            elif self.filename.lower().endswith(('.hdf5', '.h5')):
                self.read_hdf5(filename=self.filename, mol=mol)
            elif self.filename.lower().endswith('.kta'):
                self.read_nemesis(filename=self.filename, mol=mol)
            elif self.filename.lower().endswith('ktable.exorem.txt'):
                self.read_exorem(filename=self.filename, mol=mol)
            elif self.filename.lower().endswith('.fits'):
                self.read_arcis(filename=self.filename, mol=mol)
            else:
                raise NotImplementedError("""Requested format not recognized.
            Should end with .pickle, .hdf5, .h5, .fits, ktable.exorem.txt, or .kta""")
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

    def xtable_to_ktable(self, xtable=None, wnedges=None, weights=None, ggrid=None,
        quad='legendre', order=20, g_split=0.9, mid_dw=True, write=0, remove_zeros=False):
        """Fills the :class:`~exo_k.ktable.Ktable` object with a k-coeff table computed
        from a :class:`~exo_k.xtable.Xtable` object (inplace).

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
                Type of quadrature used. Default is 'legendre'.
                Also available: 'split-legendre' which uses half quadrature
                points between 0. and `g_split` and half between `g_split` and 1.
            order: Integer, optional
                Order of the Gauss legendre quadrature used. Default is 20.
            g_split: float, optional
                Used only if quad='split-legendre'. See above.
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
            self.sampling_method='custom'
        else:
            if quad=='legendre':
                self.weights,self.ggrid,self.gedges=gauss_legendre(order)
            elif quad=='split-legendre':
                self.weights,self.ggrid,self.gedges=split_gauss_legendre(order, g_split)
            else:
                raise NotImplementedError("Type of quadrature (quad keyword) not known.")
            self.sampling_method=quad
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
        if remove_zeros : self.remove_zeros(deltalog_min_value=10.)

    def hires_to_ktable(self, path=None, filename_grid=None,
        logpgrid=None, tgrid=None, wnedges=None,
        quad='legendre', order=20, g_split=0.9, weights=None, ggrid=None,
        mid_dw=True, write=0, mol=None,
        grid_p_unit='Pa', p_unit='unspecified',
        kdata_unit='unspecified', file_kdata_unit='unspecified',
        remove_zeros=False,
        **kwargs):
        """Computes a k coeff table from :class:`~exo_k.util.hires_spectrum.Hires_spectrum`
        objects (inplace).

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
            quad: string, optional
                Type of quadrature used. Default is 'legendre'.
                Also available: 'split-legendre' which uses half quadrature
                points between 0. and `g_split` and half between `g_split` and 1.
            order: Integer, optional
                Order of the Gauss legendre quadrature used. Default is 20.
            g_split: float, optional
                Used only if quad='split-legendre'. See above.
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
            self.sampling_method='custom'
        else:
            if quad=='legendre':
                self.weights,self.ggrid,self.gedges=gauss_legendre(order)
            elif quad=='split-legendre':
                self.weights,self.ggrid,self.gedges=split_gauss_legendre(order, g_split)
            else:
                raise NotImplementedError("Type of quadrature (quad keyword) not known.")
            self.sampling_method=quad
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
        if remove_zeros : self.remove_zeros(deltalog_min_value=10.)



    def copy(self, cp_kdata=True, ktab5d=False):
        """Creates a new instance of :class:`Ktable` object and (deep) copies data into it.

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
        """Method to bin down a kcoeff table to a new grid of wavenumbers (inplace).

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

    def remap_g(self, ggrid=None, weights=None):
        """Method to resample a kcoeff table to a new g grid (inplace).

        Parameters
        ----------
            ggrid: array
                New grid of abcissas for quadrature
            weights: array
                New grid of weights
        """
        weights=np.array(weights)
        ggrid=np.array(ggrid)
        self.Ng=weights.size
        newkdata=np.zeros((self.Np, self.Nt, self.Nw, self.Ng))
        g_sample_4d(ggrid, newkdata, self.ggrid, self.kdata)
        self.kdata=newkdata
        self.ggrid=ggrid
        self.weights=weights
        self.gedges=np.concatenate(([0],np.cumsum(self.weights)))

