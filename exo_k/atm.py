# -*- coding: utf-8 -*-
"""
@author: jeremy leconte

This module contain classes to handle atmosphes
and their radiative properties.

This alows us to compute the transmission and emission spectra
of those atmospheres.
"""
import numpy as np
import numba
from scipy.special import expn
import astropy.units as u
from numba.typed import List
from .gas_mix import Gas_mix
from .util.cst import N_A, PI, RGP, KBOLTZ, RSOL, RJUP, SIG_SB
#from .util.radiation import Bnu_integral_num, Bnu, rad_prop_corrk, rad_prop_xsec,\
#    Bnu_integral_array, path_integral_corrk, path_integral_xsec, dBnudT_array
from .util.radiation import path_integral_corrk, path_integral_xsec
from .util.radiation2 import Bnu_integral_num, Bnu, rad_prop_corrk, rad_prop_xsec,\
    Bnu_integral_array, dBnudT_array
from .two_stream import two_stream_toon as toon
from .two_stream import two_stream_lmdz as lmdz
from .two_stream import two_stream_crisp as crisp
from .util.interp import gauss_legendre
from .util.spectrum import Spectrum


class Atm_profile(object):
    """A class defining an atmospheric PT profile with some global data
    (gravity, etc.)

    The derived class :class:`~exo_k.atm.Atm` handles
    radiative transfer calculations.

    """
    
    def __init__(self, composition={}, psurf=None, ptop=None, logplay=None, tlay=None,
            Tsurf=None, Tstrat=None, grav=None, Rp=None, Mgas=None, rcp=0.28, Nlay=20,
            logplev=None, 
            ## old parameters that should be None. THey are left here to catch
            ## exceptions and warn the user that their use is obsolete
            Nlev=None, tlev=None,
            **kwargs):
        """Initializes atmospheric profiles

        Parameters
        ----------
            composition: dict
                Keys are molecule names and values the vmr.
                Vmr can be arrays of size Nlev-1 (i.e. the number of layers).
            grav: float
                Planet surface gravity (gravity constant with altitude for now).
            Rp: float or Astropy.unit quantity
                Planet radius. If float, Jupiter radii are assumed.
            rcp: float
                Adiabatic lapse rate for the gas (R/cp)
            Mgas: float, optional
                Molar mass of the gas (kg/mol). If given, overrides the molar mass computed
                from composition.
        
        There are two ways to define the profile.
        You can define:

        * Nlay: int
          Number of layers
        * psurf, Tsurf: float
          Surface pressure (Pa) and temperature 
        * ptop: float
          Pressure at the top of the model (Pa) 
        * Tstrat: float
          Stratospheric temperature        

        This way you will have an adiabatic atmosphere with Tsurf at the ground that
        becomes isothermal wherever T=<Tstrat.
        You can also specify:

        * logplay or play: array
        * tlay: array (same size)
          These will become the pressures (Pa; the log10 if you give
          logplay) and temperatures of the layers.
          This will be used to define the surface and top pressures.
          Nlay becomes the size of the arrays. 

        .. warning::
            Layers are counted from the top down (increasing pressure order).
            All methods follow the same convention.
        """
        if (Nlev is not None) or (tlev is not None):
            print("""
                since version 1.1.0, Nlev, tlev, and plev have been renamed
                Nlay, tlay, and logplay for consistency with other codes.
                Just change the name of the variables in the method call
                and you should be just fine!
                """)
            raise RuntimeError('Unknown keyword argument in __init__')
        self.gas_mix=Gas_mix(composition)
        self.rcp=rcp
        self.logplev=None
        if logplay is None:
            self.Nlay=Nlay
            self.Nlev=Nlay+1
            self.logplay=np.linspace(np.log10(ptop),np.log10(psurf),num=self.Nlay)
            self.compute_pressure_levels()
            self.set_adiab_profile(Tsurf=Tsurf, Tstrat=Tstrat)
        else:
            self.set_logPT_profile(logplay, tlay, logplev=logplev)
        self.set_Rp(Rp)        
        self.set_grav(grav)
        self.set_Mgas(Mgas)

    def set_logPT_profile(self, logplay, tlay, logplev=None):
        """Set the logP-T profile of the atmosphere with a new one

        Parameters
        ----------
            logplay: array
                Log pressure (in Pa) of the layer
            tlay: array (same size)
                temperature of the layers.
        
        Other Parameters
        ----------------
            logplev: array (size Nlay+1)
                If provided, allows the user to choose the location
                of the level surfaces separating the layers.
        """
        self.logplay=np.array(logplay, dtype=float)
        self.Nlay=self.logplay.size
        self.Nlev=self.Nlay+1
        if logplev is not None:
            if logplev.size == self.Nlev:
                self.logplev=np.array(logplev, dtype=float)
            else:
                raise RuntimeError('logplev does not have the size Nlay+1')
        self.compute_pressure_levels()
        self.set_T_profile(tlay)

    def set_T_profile(self, tlay):
        """Reset the temperature profile without changing the pressure levels
        """
        tlay=np.array(tlay, dtype=float)
        if tlay.shape != self.logplay.shape:
            raise RuntimeError('tlay and logplay should have the same size.')
        self.tlay=tlay
        self.t_rad=(self.tlay[:-1]+self.tlay[1:])*0.5
#        self.gas_mix.set_logPT(logp_array=self.logplay, t_array=self.tlay)
        self.gas_mix.set_logPT(logp_array=self.logp_rad, t_array=self.t_rad)

    def compute_pressure_levels(self):
        """Computes various pressure related quantities
        """
        if self.logplay[0] >= self.logplay[-1]:
            print("""
            Atmospheres are modelled from the top down.
            All arrays should be ordered accordingly
            (first values correspond to top of atmosphere)""")
            raise RuntimeError('Pressure grid is in decreasing order!')
        self.play=10**self.logplay
        if self.logplev is None:
        ## case where the levels are halfway between layer centers
        #    self.plev=np.zeros(self.Nlev)
        #    self.plev[1:-1]=(self.play[:-1]+self.play[1:])*0.5
        #    # we choose to use mid point so that there is equal mass in the bottom half
        #    # of any top layer and the top half of the layer below. 
        #    self.plev[0]=self.play[0]
        #    self.plev[-1]=self.play[-1]
        #    ## WARNING: Top and bottom pressure levels are set equal to the
        #    #  pressure in the top and bottom layers. If you change that,
        #    #  some assumptions here and there in the code may break down!!!
        #    self.logplev=np.log10(self.plev)

            self.logplev=np.zeros(self.Nlev)
            self.logplev[1:-1]=(self.logplay[:-1]+self.logplay[1:])*0.5
            self.logp_rad=self.logplev[1:-1]
            self.logplev[0]=self.logplay[0]
            self.logplev[-1]=self.logplay[-1]
            self.plev=np.power(10.,self.logplev)
        else:
            self.plev=10**self.logplev
        self.psurf=self.plev[-1]
        self.dp_lay=np.diff(self.plev) ### probably redundant with dmass
        self.exner=(self.play/self.psurf)**self.rcp

    def set_adiab_profile(self, Tsurf=None, Tstrat=None):
        """Initializes the logP-T atmospheric profile with an adiabat with index R/cp=rcp

        Parameters
        ----------
            Tsurf: float
                Surface temperature.
            Tstrat: float, optional
                Temperature of the stratosphere. If None is given,
                an isothermal atmosphere with T=Tsurf is returned.
        """
        if Tstrat is None: Tstrat=Tsurf
        self.tlay=Tsurf*self.exner
        self.tlay=np.where(self.tlay<Tstrat,Tstrat,self.tlay)
        self.t_rad=(self.tlay[:-1]+self.tlay[1:])*0.5
#        self.gas_mix.set_logPT(logp_array=self.logplay, t_array=self.tlay)
        self.gas_mix.set_logPT(logp_array=self.logp_rad, t_array=self.t_rad)

    def set_grav(self, grav=None):
        """Sets the surface gravity of the planet

        Parameters
        ----------
            grav: float
                surface gravity (m/s^2)
        """
        if grav is None: raise RuntimeError('A planet needs a gravity!')
        self.grav=grav
    
    def set_gas(self, composition_dict):
        """Sets the composition of the atmosphere

        Parameters
        ----------
            composition_dict: dictionary
                Keys are molecule names, and values are volume mixing ratios.
                A 'background' value means that the gas will be used to fill up to vmr=1
                If they do not add up to 1 and there is no background gas_mix,
                the rest of the gas_mix is considered transparent.
        """
        self.gas_mix.set_composition(composition_dict)
        self.set_Mgas()

    def set_Mgas(self, Mgas=None):
        """Sets the mean molar mass of the atmosphere.

        Parameters
        ----------
            Mgas: float or array of size Nlay
                Mean molar mass (kg/mol).
                If None is give, the mmm is computed from the composition.
        """
        if Mgas is not None:
            self.Mgas=Mgas
        else:
            self.Mgas=self.gas_mix.molar_mass()
        self.Mgas=self.Mgas*np.ones(self.Nlay-1, dtype=np.float)

    def set_rcp(self,rcp):
        """Sets the adiabatic index of the atmosphere

        Parameters
        ----------
            rcp: float
                R/c_p
        """
        self.rcp=rcp

    def set_Rp(self, Rp):
        """Sets the radius of the planet

        Parameters
        ----------
            Rp: float
                radius of the planet (m)
        """
        if Rp is None:
            self.Rp = None
            return
        if isinstance(Rp,u.quantity.Quantity):
            self.Rp=Rp.to(u.m).value
        else:
            self.Rp=Rp

    def set_Rstar(self, Rstar):
        """Sets the radius of the star

        Parameters
        ----------
            Rstar: float
                radius of the star (m)
        """
        if Rstar is None:
            self.Rstar = None
            return
        if isinstance(Rstar,u.quantity.Quantity):
            self.Rstar=Rstar.to(u.m).value
        else:
            self.Rstar=Rstar

    def compute_density(self):
        """Computes the number density (m^-3) profile of the atmosphere
        """
        self.density=self.play/(KBOLTZ*self.tlay)

    def compute_layer_col_density(self):
        """Computes the column number density (molecules/m^-2) per
        radiative layer of the atmosphere.

        There are Nlay-1 radiative layers as they go from the midle of a layer to the next.

        dcol_density_radlay_up (size Nlay-1) is the column density in the upper half of each layer
        dcol_density_radlay_dw (size Nlay-1) is the column density in the lower half of each layer

        """
        factor=N_A/(self.grav * self.Mgas)
        #self.dcol_density_radlay_up=(self.plev[1:-1]-self.play[:-1])*factor[:-1]
        #self.dcol_density_radlay_dw=(self.play[1:]-self.plev[1:-1])*factor[1:]
        self.dcol_density_rad = np.diff(self.play)*factor[:]

        if self.Rp is not None: #includes the altitude effect if radius is known
            self.compute_altitudes()
            self.dcol_density_radlay_up*=(1.+self.zlay[:-1]/self.Rp)**2
            self.dcol_density_radlay_dw*=(1.+self.zlay[1:]/self.Rp)**2

    def compute_altitudes(self):
        """Compute altitudes of the level surfaces (zlev) and mid layers (zlay).
        """
        H=RGP*self.tlay/(self.grav*self.Mgas)
        dlnP=np.diff(self.logplev)*np.log(10.)
        self.zlev=np.zeros_like(self.logplev)
        if self.Rp is None:
            self.dz=H*dlnP
            self.zlev[:-1]=np.cumsum(self.dz[::-1])[::-1]
            self.zlay=0.5*(self.zlev[1:]+self.zlev[:-1])
            ## assumes layer centers at the middle of the two levels
            ## which is not completely consistent with play, but should be
            ## a minor error.
        else:
            for i in range(H.size)[::-1]:
                z1=self.zlev[i+1]
                H1=H[i]
                dlnp=dlnP[i]
                self.zlev[i]=z1+( (H1 * (self.Rp + z1)**2 * dlnp) \
                    / (self.Rp**2 + H1 * self.Rp * dlnp + H1 * z1 * dlnp) )
        self.zlay=0.5*(self.zlev[1:]+self.zlev[:-1])
        ## assumes layer centers at the middle of the two levels
        ## which is not completely consistent with play, but should be
        ## a minor error.
        
    def compute_area(self):
        """Computes the area of the annulus covered by each layer in a transit setup. 
        """
        self.area=PI*(self.Rp+self.zlev[:-1])**2
        self.area[:-1]-=self.area[1:]
        self.area[-1]-=PI*self.Rp**2

    def compute_tangent_path(self):
        """Computes a triangular array of the tangent path length (in m) spent in each layer.
        self.tangent_path[ilay][jlay] is the length that the ray that is tangent to the ilay layer 
        spends in the jlay>=ilay layer (accounting for a factor of 2 due to symmetry)
        """
        if self.Rp is None: raise RuntimeError('Planetary radius should be set')
        self.compute_altitudes()
        self.tangent_path=List()
        # List() is a new numba.typed list to comply with new numba evolution after v0.50
        for ilay in range(self.Nlay):
            z0square=(self.Rp+self.zlay[ilay])**2
            dl=np.sqrt((self.Rp+self.zlev[:ilay+1])**2-z0square)
            dl[:-1]-=dl[1:]
            self.tangent_path.append(2.*dl)

    def __repr__(self):
        """Method to output header
        """
        output="""
    gravity (m/s^2) : {grav}
    Planet Radius(m): {rad}
    Ptop (Pa)       : {ptop}
    Psurf (Pa)      : {psurf}
    Tsurf (K)       : {tsurf}
    composition     :
        {comp}""".format(grav=self.grav, rad=self.Rp, comp=self.gas_mix,
            ptop=self.plev[0], psurf=self.psurf, tsurf=self.tlay[-1])
        return output



class Atm(Atm_profile):
    """Class based on Atm_profile that handles radiative trasnfer calculations.

    Radiative data are accessed through the :any:`gas_mix.Gas_mix` class.
    """

    def __init__(self, k_database=None, cia_database=None,
        wn_range=None, wl_range=None, **kwargs):
        """Initialization method that calls Atm_Profile().__init__() and links
        to Kdatabase and other radiative data. 
        """
        super().__init__(**kwargs)
        self.set_k_database(k_database)
        self.set_cia_database(cia_database)
        self.set_spectral_range(wn_range=wn_range, wl_range=wl_range)
        self.flux_net_nu=None
        self.kernel=None
        self.tlay_kernel=self.tlay

    def set_k_database(self, k_database=None):
        """Change the radiative database used by the
        :class:`Gas_mix` object handling opacities inside
        :class:`Atm`.

        See :any:`gas_mix.Gas_mix.set_k_database` for details.

        Parameters
        ----------
            k_database: :class:`Kdatabase` object
                New Kdatabase to use.
        """
        self.gas_mix.set_k_database(k_database=k_database)
        self.kdatabase=self.gas_mix.kdatabase
        self.Ng=self.gas_mix.Ng
        # to know whether we are dealing with corr-k or not and access some attributes. 

    def set_cia_database(self, cia_database=None):
        """Change the CIA database used by the
        :class:`Gas_mix` object handling opacities inside
        :class:`Atm`.

        See :any:`gas_mix.Gas_mix.set_cia_database` for details.

        Parameters
        ----------
            cia_database: :class:`CIAdatabase` object
                New CIAdatabase to use.
        """
        self.gas_mix.set_cia_database(cia_database=cia_database)

    def set_spectral_range(self, wn_range=None, wl_range=None):
        """Sets the spectral range in which computations will be done by specifying
        either the wavenumber (in cm^-1) or the wavelength (in micron) range.

        See :any:`gas_mix.Gas_mix.set_spectral_range` for details.
        """
        self.gas_mix.set_spectral_range(wn_range=wn_range, wl_range=wl_range)

    def spectral_integration(self, spectral_var):
        """Spectrally integrate an array, taking care of whether
        we are dealing with corr-k or xsec data.

        Parameters
        ----------
            spectral_var: array
                array to integrate

        Returns
        -------
            var: array
                array integrated over wavenumber (and g-space if relevant)
        """
        if self.Ng is None:
            var=np.sum(spectral_var*self.dwnedges,axis=-1)
        else:
            var=np.sum(np.sum(spectral_var*self.weights,axis=-1)*self.dwnedges,axis=-1)
        return var

    def opacity(self, rayleigh = False, **kwargs):
        """Computes the opacity of each of the layers.

        See :any:`gas_mix.Gas_mix.cross_section` for details.
        """
        self.kdata = self.gas_mix.cross_section(rayleigh=rayleigh, **kwargs)
        if rayleigh: self.kdata_scat=self.gas_mix.kdata_scat
        self.Nw=self.gas_mix.Nw
        self.wns=self.gas_mix.wns
        self.wnedges=self.gas_mix.wnedges
        self.dwnedges=self.gas_mix.dwnedges

    def source_function(self, integral=True, source=True):
        """Compute the blackbody source function (Pi*Bnu) for each layer of the atmosphere.

        Parameters
        ----------
            integral: boolean, optional
                * If true, the black body is integrated within each wavenumber bin.
                * If not, only the central value is used.
                  False is faster and should be ok for small bins,
                  but True is the correct version. 
            source: boolean, optional
                If False, the source function is put to 0 (for solar absorption calculations)
        """
        if source:
            if integral:
                #JL2020 What is below is slightly faster for very large resolutions
                #piBatm=np.empty((self.Nlay,self.Nw))
                #for ii in range(self.Nlay):
                #    piBatm[ii]=PI*Bnu_integral_num(self.wnedges,self.tlay[ii])/dw
                #JL2020 What is below is much faster for moderate to low resolutions
                piBatm=PI*Bnu_integral_array(self.wnedges,self.tlay,self.Nw,self.Nlay) \
                    /self.dwnedges
            else:
                piBatm=PI*Bnu(self.wns[None,:],self.tlay[:,None])
        else:
            piBatm=np.zeros((self.Nlay,self.Nw))
        return piBatm

    def setup_emission_caculation(self, mu_eff=0.5, rayleigh=False, integral=True,
            source=True, **kwargs):
        """Computes all necessary quantities for emission calculations
        (opacity, source, etc.)
        """
        self.opacity(rayleigh=rayleigh, **kwargs)
        self.piBatm = self.source_function(integral=integral, source=source)
        self.compute_layer_col_density()
        if self.Ng is None:
#            self.tau, self.dtau=rad_prop_xsec(self.dcol_density_radlay_up,
#                self.dcol_density_radlay_dw, self.kdata, mu_eff)
            self.tau, self.dtau=rad_prop_xsec(self.dcol_density_rad,
                self.kdata, mu_eff)
        else:
#            self.tau, self.dtau=rad_prop_corrk(self.dcol_density_radlay_up,
#                self.dcol_density_radlay_dw, self.kdata, mu_eff)
            self.tau, self.dtau=rad_prop_corrk(self.dcol_density_rad,
                self.kdata, mu_eff)
            self.weights=self.kdatabase.weights
            #print('trad:',self.t_rad.shape)
            #print('kdata:',self.kdata.shape)
            #print('tau:',self.tau.shape)
            #print('dtau:',self.dtau.shape)
            #print('piBatm:',self.piBatm.shape)

    def emission_spectrum(self, integral=True, mu0=0.5, mu_quad_order=None, **kwargs):
        """Returns the emission flux at the top of the atmosphere (in W/m^2/cm^-1)

        Parameters
        ----------
            integral: boolean, optional
                * If true, the black body is integrated within each wavenumber bin.
                * If not, only the central value is used.
                  False is faster and should be ok for small bins,
                  but True is the correct version. 

        Other Parameters
        ----------------
            mu0: float
                Cosine of the quadrature angle use to compute output flux
            mu_quad_order: int
                If an integer is given, the emission intensity is computed
                for a number of angles and integrated following a gauss legendre
                quadrature rule of order `mu_quad_order`.

        Returns
        -------
            Spectrum object 
                A spectrum with the Spectral flux at the top of the atmosphere (in W/m^2/cm^-1)
        """
        if mu_quad_order is not None:
            # if we want quadrature, use the more general method.
            return self.emission_spectrum_quad(integral=integral,
                mu_quad_order=mu_quad_order, **kwargs)

        try:
            self.setup_emission_caculation(mu_eff=mu0, rayleigh=False, integral=integral, **kwargs)
        except TypeError:
            raise RuntimeError("""
            Cannot use rayleigh option with emission_spectrum.
            If you meant to include scattering, you should use emission_spectrum_2stream.
            """)
        # self.tau and self.dtau include the 1/mu0 factor.
        expdtau=np.exp(-self.dtau)
        expdtauminone=np.where(self.dtau<1.e-13,-self.dtau,expdtau-1.)
        # careful: due to numerical limitations, 
        # the limited development of Exp(-dtau)-1 needs to be used for small values of dtau
        exptau=np.exp(-self.tau)
        if self.Ng is None:
            timesBatmTop=(1.+expdtauminone/self.dtau)*exptau[:-1]
            timesBatmBottom=(-expdtau-expdtauminone/self.dtau)*exptau[:-1]
            timesBatmBottom[-1]+=exptau[-1]
        else:
            timesBatmTop=np.sum((1.+expdtauminone/self.dtau)*exptau[:-1]*self.weights,axis=-1)
            timesBatmBottom=np.sum((-expdtau-expdtauminone/self.dtau)*exptau[:-1] \
                *self.weights,axis=-1)
            timesBatmBottom[-1]+=np.sum(exptau[-1]*self.weights,axis=-1)
        IpTop=np.sum(self.piBatm[:-1]*timesBatmTop+self.piBatm[1:]*timesBatmBottom,axis=0)

        return Spectrum(IpTop,self.wns,self.wnedges)

    def emission_spectrum_quad(self, integral=True, mu_quad_order=3, **kwargs):
        """Returns the emission flux at the top of the atmosphere (in W/m^2/cm^-1)
        using gauss legendre qudrature of order `mu_quad_order`

        Parameters
        ----------
            integral: boolean, optional
                * If true, the black body is integrated within each wavenumber bin.
                * If not, only the central value is used.
                  False is faster and should be ok for small bins,
                  but True is the correct version. 
                  
        Returns
        -------
            Spectrum object 
                A spectrum with the Spectral flux at the top of the atmosphere (in W/m^2/cm^-1)
        """
        self.setup_emission_caculation(mu_eff=1., rayleigh=False, integral=integral, **kwargs)
        # angle effect dealt with later

        IpTop=np.zeros(self.kdata.shape[1])
        mu_w, mu_a, _ = gauss_legendre(mu_quad_order)
        mu_w = mu_w * mu_a * 2.# takes care of the mu factor in last integral => int(mu I d mu)
                               # Factor 2 takes care of the fact that the source function is pi*Batm
                               # but we want 2*Pi*Batm
        for ii, mu0 in enumerate(mu_a):
            tau=self.tau/mu0
            dtau=self.dtau/mu0
            expdtau=np.exp(-dtau)
            expdtauminone=np.where(dtau<1.e-13,-dtau,expdtau-1.)
            exptau=np.exp(-tau)
            if self.Ng is None:
                timesBatmTop=(-expdtau-expdtauminone/dtau)*exptau[:-1]
                timesBatmBottom=(1.+expdtauminone/dtau)*exptau[:-1]
                timesBatmBottom[-1]=timesBatmBottom[-1]+exptau[-1]
            else:
                timesBatmTop=np.sum((-expdtau-expdtauminone/dtau)*exptau[:-1] \
                    *self.weights,axis=-1)
                timesBatmBottom=np.sum((1.+expdtauminone/dtau)*exptau[:-1] \
                    *self.weights,axis=-1)
                timesBatmBottom[-1]=timesBatmBottom[-1]+np.sum(exptau[-1]*self.weights,axis=-1)
            IpTop+=np.sum(self.piBatm[:-1]*timesBatmTop+self.piBatm[1:]*timesBatmBottom,axis=0) \
                *mu_w[ii]

        return Spectrum(IpTop,self.wns,self.wnedges)

    def emission_spectrum_2stream(self, integral=True, mu0=0.5,
            method='toon', dtau_min=1.e-10, flux_at_level=False, rayleigh=False,
            flux_top_dw=None, source=True, compute_kernel=False, **kwargs):
        """Returns the emission flux at the top of the atmosphere (in W/m^2/cm^-1)

        Parameters
        ----------
            integral: boolean, optional
                * If true, the black body is integrated within each wavenumber bin.
                * If not, only the central value is used.
                  False is faster and should be ok for small bins,
                  but True is the correct version. 

        Other Parameters
        ----------------
            mu0: float
                Cosine of the quadrature angle use to compute output flux
            dtau_min: float
                If the optical depth in a layer is smaller than dtau_min,
                dtau_min is used in that layer instead. Important as too
                transparent layers can cause important numerical rounding errors.

        Returns
        -------
            Spectrum object 
                A spectrum with the Spectral flux at the top of the atmosphere (in W/m^2/cm^-1)
        """
        self.setup_emission_caculation(mu_eff=1., rayleigh=rayleigh, integral=integral,
            source=source, **kwargs)
        # mu_eff=1. because the mu effect is taken into account in solve_2stream_nu
                 # we must compute the vertical optical depth here.
        self.dtau=np.where(self.dtau<dtau_min,dtau_min,self.dtau)

        module_to_use=globals()[method]
        # globals()[method] converts the method string into a module name
        #  if the module has been loaded
        if self.Ng is None:
            solve_2stream_nu=module_to_use.solve_2stream_nu_xsec
        else:
            solve_2stream_nu=module_to_use.solve_2stream_nu_corrk

        if rayleigh:
            if self.Ng is None:
                self.single_scat_albedo = self.kdata_scat / self.kdata
            else:
                self.single_scat_albedo = self.kdata_scat[:,:,None] / self.kdata
        else:
            self.single_scat_albedo = np.zeros_like(self.dtau)
        self.single_scat_albedo=np.clip(self.single_scat_albedo,None,0.9999999999999)
        self.asym_param = np.zeros_like(self.dtau)

        if flux_top_dw is None:
            self.flux_top_dw_nu = np.zeros((self.Nw))
        else:
            Tstar=5700.
            self.flux_top_dw_nu = Bnu_integral_num(self.wnedges,Tstar)
            self.flux_top_dw_nu = self.flux_top_dw_nu * flux_top_dw \
                / (np.sum(self.flux_top_dw_nu)*self.dwnedges)
            #self.flux_top_dw_nu=Spectrum(flux_top_dw_nu,self.wns,self.wnedges)

        if method == 'crisp':
            self.flux_up_nu, self.flux_down_nu, self.flux_net_nu, kernel_nu = \
                solve_2stream_nu(self.piBatm, self.dtau,
                self.single_scat_albedo, self.asym_param, self.flux_top_dw_nu, mu0 = mu0)
            if compute_kernel:
                dB = PI * dBnudT_array(self.wns, self.tlay, self.Nw, self.Nlay)
                if self.Ng is None:
                    self.kernel=np.sum(kernel_nu*dB*self.dwnedges,axis=2)
                else:
                    self.kernel=np.sum(np.sum(kernel_nu*self.weights,axis=3) \
                        *dB*self.dwnedges,axis=2)
                self.kernel[:,:-1]=(self.kernel[:,:-1]-self.kernel[:,1:])
                self.kernel*=self.grav/self.dp_lay
        else:
            self.flux_up_nu, self.flux_down_nu, self.flux_net_nu = \
                solve_2stream_nu(self.piBatm, self.dtau, self.single_scat_albedo, self.asym_param,
                    self.flux_top_dw_nu, mu0 = mu0, flux_at_level=flux_at_level)
            if compute_kernel: self.compute_kernel(solve_2stream_nu, mu0=mu0, flux_at_level=flux_at_level,
                        per_unit_mass=True, integral=True, **kwargs)

        if self.Ng is None:
            return Spectrum(self.flux_up_nu[0],self.wns,self.wnedges)
        else:
            return Spectrum(np.sum(self.flux_up_nu[0]*self.weights,axis=1),self.wns,self.wnedges)

    def compute_kernel(self, solve_2stream_nu, epsilon=0.01, flux_at_level=False, mu0 = 0.5,
            per_unit_mass=True, integral=True, **kwargs):
        """Compute the Jacobian matrix d Heating[lay=i] / d T[lay=j]
        """
        net=self.spectral_integration(self.flux_net_nu)
        self.kernel=np.empty((self.Nlay,self.Nlay))
        #dB = PI * dBnudT_array(self.wns, self.tlay, self.Nw, self.Nlay)
        tlay=self.tlay
        dT = epsilon*tlay
        self.tlay = tlay + dT
        newpiBatm = self.source_function(integral=integral)

        #newpiBatm=PI*Bnu_integral_array(self.wnedges,self.tlay+dT,self.Nw,self.Nlay) \
        #            /self.dwnedges

        for ilay in range(self.Nlay):
            pibatm = np.copy(self.piBatm)
            #dT = epsilon*self.tlay[ilay]
            #pibatm[ilay] += dB[ilay]*dT
            pibatm[ilay] = newpiBatm[ilay]
            _, _, flux_net_tmp = \
                solve_2stream_nu(pibatm, self.dtau, self.single_scat_albedo, self.asym_param,
                    self.flux_top_dw_nu, mu0 = mu0, flux_at_level=flux_at_level)
            net_tmp = self.spectral_integration(flux_net_tmp)
            #self.kernel[ilay]=(net_tmp-net)/dT
            self.kernel[ilay]=(net-net_tmp)/dT[ilay]
        self.kernel[:,:-1]-=self.kernel[:,1:]
        self.tlay=tlay
        if per_unit_mass: self.kernel*=self.grav/self.dp_lay

    def flux_divergence(self, internal_flux=0., per_unit_mass = True):
        """Computes the divergence of the net flux in the layers
        (used to compute heating rates).

        :func:`emission_spectrum_2stream` needs to be ran first.

        Parameters
        ----------
            internal_flux: float
                flux coming from below in W/m^2
            per_unit_mass: bool
                If True, the heating rates are normalized by the
                mass of each layer (result in W/kg).

        Returns
        -------
            net: array
                Net fluxes at level surfaces
            H: array
                Heating rate in each layer (Difference of the net fluxes). Positive means heating.
                The last value is the net flux impinging on the surface + the internal flux.
        """
        if self.flux_net_nu is None:
            raise RuntimeError('should have ran emission_spectrum_2stream.')
        net=self.spectral_integration(self.flux_net_nu)
        H=-np.copy(net)
        H[:-1]-=H[1:]
        H[-1]+=internal_flux
        if per_unit_mass: H*=self.grav/self.dp_lay
        return net, H

    def bolometric_fluxes(self, internal_flux = 0., per_unit_mass = True):
        """Computes the bolometric fluxes at levels and the divergence of the net flux in the layers
        (used to compute heating rates).

        :func:`emission_spectrum_2stream` needs to be ran first.

        Parameters
        ----------
            internal_flux: float
                flux coming from below in W/m^2
            per_unit_mass: bool
                If True, the heating rates are normalized by the
                mass of each layer (result in W/kg).

        Returns
        -------
            up: array
                Upward fluxes at level surfaces
            dw: array
                Downward fluxes at level surfaces
            net: array
                Net fluxes at level surfaces
            H: array
                Heating rate in each layer (Difference of the net fluxes). Positive means heating.
                The last value is the net flux impinging on the surface + the internal flux.
        """
        net, H = self.flux_divergence(internal_flux=internal_flux, per_unit_mass = per_unit_mass)
        up=self.spectral_integration(self.flux_up_nu)
        dw=self.spectral_integration(self.flux_down_nu)
        return up, dw, net, H

    def exp_minus_tau(self):
        """Sums Exp(-tau) over gauss points
        """
        weights=self.kdatabase.weights
        return np.sum(np.exp(-self.tau[1:])*weights,axis=2)

    def exp_minus_tau_g(self, g_index):
        """Sums Exp(-tau) over gauss point
        """
        return np.exp(-self.tau[1:,:,g_index])

    def surf_bb(self, integral=True):
        """Computes the surface black body flux (in W/m^2/cm^-1)

        Parameters
        ----------
            integral: boolean, optional
                * If true, the black body is integrated within each wavenumber bin.
                * If not, only the central value is used.
                  False is faster and should be ok for small bins,
                  but True is the correct version. 
        Returns
        -------
            Spectrum object
                Spectral flux at the surface (in W/m^2/cm^-1)
        """
        if integral:
            piBatm=PI*Bnu_integral_num(self.wnedges,self.tlay[-1])/np.diff(self.wnedges)
        else:
            piBatm=PI*Bnu(self.wns[:],self.tlay[-1])
        return Spectrum(piBatm,self.wns,self.wnedges)

    def top_bb(self, integral=True):
        """Computes the top of atmosphere black body flux (in W/m^2/cm^-1)

        Parameters
        ----------
            integral: boolean, optional
                If true, the black body is integrated within each wavenumber bin.
                If not, only the central value is used.
                    False is faster and should be ok for small bins,
                    but True is the correct version. 
        Returns
        -------
            Spectrum object
                Spectral flux of a bb at the temperature at the top of atmosphere (in W/m^2/cm^-1)
        """
        if integral:
            piBatm=PI*Bnu_integral_num(self.wnedges,self.tlay[0])/np.diff(self.wnedges)
        else:
            piBatm=PI*Bnu(self.wns[:],self.tlay[0])
        return Spectrum(piBatm,self.wns,self.wnedges)

    def transmittance_profile(self, **kwargs):
        """Computes the transmittance profile of an atmosphere,
        i.e. Exp(-tau) for each layer of the model.
        Real work done in the numbafied function path_integral_corrk/xsec
        depending on the type of data.
        """
        self.opacity(**kwargs)
        self.compute_tangent_path()
        self.compute_density()
        if self.Ng is not None:
            self.weights=self.kdatabase.weights
            transmittance=path_integral_corrk( \
                self.Nlay,self.Nw,self.Ng,self.tangent_path,self.density,self.kdata,self.weights)
        else:
            transmittance=path_integral_xsec( \
                self.Nlay,self.Nw,self.tangent_path,self.density,self.kdata)
        return transmittance

    def transmission_spectrum(self, normalized=False, Rstar=None, **kwargs):
        r"""Computes the transmission spectrum of the atmosphere.
        In general (see options below), the code returns the transit depth:

        .. math::
            \delta_\nu=(\pi R_p^2+\alpha_\nu)/(\pi R_{star}^2),

        where

        .. math::
          \alpha_\nu=2 \pi \int_0^{z_{max}} (R_p+z)*(1-e^{-\tau_\nu(z)) d z.
        
        Parameters
        ----------
            Rstar: float, optional
                Radius of the host star. Does not need to be given here if
                as already been specified as an attribute of the self.Atm object.
                If specified, the result is the transit depth:

                .. math::
                  \delta_\nu=(\pi R_p^2+\alpha_\nu)/(\pi R_{star}^2).

            normalized: boolean, optional
                Used only if self.Rstar and Rstar are None:

                * If True,
                  the result is normalized to the planetary radius:

                  .. math::
                    \delta_\nu=1+\frac{\alpha_\nu}{\pi R_p^2}.
                * If False,
                                    
                  .. math::
                    \delta_\nu=\pi R_p^2+\alpha_\nu.

        Returns
        -------
            array
                The transit spectrum (see above for normalization options).
        """
        self.set_Rstar(Rstar)
        transmittance=self.transmittance_profile(**kwargs)
        self.compute_area()
        res=Spectrum((np.dot(self.area,(1.-transmittance))),self.wns,self.wnedges)
        if self.Rstar is not None:
            return (res+(PI*self.Rp**2))/(PI*self.Rstar**2)
        elif normalized:
            return res/(PI*self.Rp**2)+1
        else:
            return res+(PI*self.Rp**2)

    def __repr__(self):
        """Method to output header
        """
        output=super().__repr__()
        output+="""
    k_database      :
        {kdatab}
    cia_database    :
        {cdatab}""".format(kdatab=self.kdatabase, cdatab=self.gas_mix.cia_database)
        if self.gas_mix._wn_range is not None:
            output+='    wn range        : '+ self.gas_mix._wn_range +'\n'

        return output


@numba.jit(nopython=True,fastmath=True)
def convadj(timestep, Nlay, tlay, exner, dp_lay, verbose = False):
    r"""Computes the heating rates needed to adjust unstable regions 
    of a given atmosphere to a convectively neutral T profile on
    a given timestep.

    .. important::
        The *layers* here are not the same as the *layers* in the radiative transfer!

        In the radiative transfer, a layer is the volume between two pressure level boundaries
        (`plev`).

        Here, a layer is centered on a pressure level (`plev`) with its temperature `tlev`.
        Therefore, `dp_lay[0]=play[0]-plev[0]`, `dp_lay[i]=play[i]-play[i-1]`, and
        `dp_lay[-1]=psurf-play[-1]`.
    
    Parameters
    ----------
        timestep: float
            Duration of the adjustment in seconds.
        Nlay: int
            Number of atmospheric layers
        tlay: array
            Temperatures of the atmospheric layers
        exner: array
            Exner function computed at the layer centers ((p/psurf)**rcp)

            .. math::
              \Pi=(p / p_{s})^{R/c_p}

        dp_lay: array
            Pressure difference between the bottom and top of each layer

    Returns
    -------
        array
            Heating rate in each atmospheric layer (K/s).  
    """
    epsilon=-1.e-5
    theta_lay=tlay/exner
    new_theta_lay=np.copy(theta_lay)
    dsig=dp_lay
    sdsig=dsig*exner
    H_conv=np.zeros(Nlay)
    n_iter=0
    while True:
        conv=np.nonzero(new_theta_lay[:-1]-new_theta_lay[1:]<epsilon)[0]
        # find convective layers
        N_conv=conv.size
        if N_conv==0: # no more convective regions, can exit
            return H_conv
        i_conv=0
        i_top=conv[i_conv] #upper unstable layer
        while i_conv<N_conv-1: #search from the top of the 1st unstable layer for its bottom
            if conv[i_conv+1]==conv[i_conv]+1:
                i_conv+=1
                continue
            else:
                break
        i_bot=conv[i_conv]+1
        mass_conv=0.
        intexner=0.
        theta_mean=0.
        for ii in range(i_top,i_bot+1): # compute new mean potential temperature
            intexner+=sdsig[ii]
            mass_conv+=dsig[ii]
            theta_mean+=sdsig[ii] * (theta_lay[ii] - theta_mean) / intexner
        while i_top>0:#look for newly unstable layers above
            if theta_lay[i_top-1]<theta_mean:
                i_top-=1
                intexner+=sdsig[i_top]
                mass_conv+=dsig[i_top]
                theta_mean+=sdsig[i_top] * (theta_lay[i_top] - theta_mean) / intexner
            else:
                break
        while i_bot<Nlay-1: #look for newly unstable layers below
            if theta_lay[i_bot+1]>theta_mean:
                i_bot+=1
                intexner+=sdsig[i_bot]
                mass_conv+=dsig[i_bot]
                theta_mean+=sdsig[i_bot] * (theta_lay[i_bot] - theta_mean) / intexner
            else:
                break
        # compute heating and adjust before looking for a new potential unstable layer
        H_conv[i_top:i_bot+1]=(theta_mean-new_theta_lay[i_top:i_bot+1]) \
            *exner[i_top:i_bot+1]/timestep
        new_theta_lay[i_top:i_bot+1]=theta_mean
        n_iter+=1
        if n_iter>Nlay+1:
            if verbose : print('oops, went crazy in convadj')
            break
    return H_conv
