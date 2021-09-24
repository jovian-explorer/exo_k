# -*- coding: utf-8 -*-
"""
@author: jeremy leconte

This module contain classes to handle atmospheres
and their radiative properties.
This alows us to compute the transmission and emission spectra
of those atmospheres.

The nomenclature for layers, levels, etc, is as follows (example with 4 layers; Nlay=4):

::

    -------------------  Model top or first level (plev[0])
    - - - - - - - - - -  First atmopheric layer (play[0]=plev[0], tlay[0], xlay[0])

    -------------------  plev[1]

    - - - - - - - - - -  play[1], tlay[1], xlay[1]

    -------------------  plev[2]

    - - - - - - - - - - play[2], tlay[2], xlay[2]

    -------------------  plev[3]

    - - - - - - - - - -  bottom layer (play[Nlay-1]=psurf, tlay[3], xlay[3])
    ------------------- Surface (plev[Nlev-1=Nlay]=psurf)
    ///////////////////

.. image:: /images/atm_schema.png


Tempratures (`tlay`) and volume mixing ratios (`xlay`) are provided at the
*mid layer* point (`play`) (Note that this does not mean
that this point is the middle of the layer, in general, it is not).

If pressure levels are not specified by the user with `logplevel`,
they are at the mid point between
the pressure of the layers directly above and below. The pressure of the
top and bottom levels are counfounded with the pressure of the
top and bottom mid layer points.

For radiative calculations, the source function (temperature) needs to be known at
the boudaries of the *radiative* layers but the opacity needs to be known
inside the *radiative* layer.
For this reason, there are `Nlay-1` radiative layers and they are offset
with respect to atmospheric layers.
Opacities are computed at the center of those radiative layers, i.e. at the pressure levels.
The temperature is interpolated at these levels with an arithmetic average.
Volume mixing ratios are interpolated using a geometric average. 

"""
import numpy as np
import astropy.units as u
from numba.typed import List
from .gas_mix import Gas_mix
from .aerosols import Aerosols
from .util.cst import N_A, PI, RGP, KBOLTZ, SIG_SB
from .util.radiation import Bnu_integral_num, Bnu, rad_prop_corrk, rad_prop_xsec,\
    Bnu_integral_array, path_integral_corrk, path_integral_xsec
from .two_stream import two_stream_toon as toon
from .two_stream import two_stream_lmdz as lmdz
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
            logplev=None, aerosols={},
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
                Planet radius. If float, meters are assumed.
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
        self.gas_mix = None
        self.set_gas(composition, compute_Mgas=False)
        self.aerosols = None
        self.set_aerosols(aerosols)
        self.rcp = rcp
        self.logplev = None
        if logplay is None:
            self.Nlay = Nlay
            self.Nlev = Nlay+1
            self.logplay = np.linspace(np.log10(ptop),np.log10(psurf),num=self.Nlay)
            self.compute_pressure_levels()
            self.set_adiab_profile(Tsurf=Tsurf, Tstrat=Tstrat)
        else:
            self.set_logPT_profile(logplay, tlay, logplev=logplev)
        self.set_Rp(Rp)        
        self.set_grav(grav)
        self.set_Mgas(Mgas=Mgas)

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
        self.t_opac=(self.tlay[:-1]+self.tlay[1:])*0.5
        self.gas_mix.set_logPT(logp_array=self.logp_opac, t_array=self.t_opac)

    def compute_pressure_levels(self):
        """Computes various pressure related quantities
        """
        if self.logplay[0] >= self.logplay[-1]:
            print("""
            Atmospheres are modelled from the top down.
            All arrays should be ordered accordingly
            (first values correspond to top of atmosphere)""")
            raise RuntimeError('Pressure grid is in decreasing order!')
        self.play=np.power(10., self.logplay)
        if self.logplev is None:
        # case where the levels are halfway between layer centers
            self.plev=np.zeros(self.Nlev)
            self.plev[1:-1]=(self.play[:-1]+self.play[1:])*0.5
            # we choose to use mid point so that there is equal mass in the bottom half
            # of any top layer and the top half of the layer below. 
            self.plev[0]=self.play[0]
            self.plev[-1]=self.play[-1]
            ## WARNING: Top and bottom pressure levels are set equal to the
            #  pressure in the top and bottom layers. If you change that,
            #  some assumptions here and there in the code may break down!!!
            self.logplev=np.log10(self.plev)

        # case where the levels are halfway between layer centers in LOG10
            #self.logplev=np.zeros(self.Nlev)
            #self.logplev[1:-1]=(self.logplay[:-1]+self.logplay[1:])*0.5
            #self.logplev[0]=self.logplay[0]
            #self.logplev[-1]=self.logplay[-1]
            #self.plev=np.power(10.,self.logplev)
        else:
            self.plev=np.power(10., self.logplev)

        self.logp_opac=self.logplev[1:-1]
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
        self.t_opac=(self.tlay[:-1]+self.tlay[1:])*0.5
        self.gas_mix.set_logPT(logp_array=self.logp_opac, t_array=self.t_opac)

    def set_grav(self, grav=None):
        """Sets the surface gravity of the planet

        Parameters
        ----------
            grav: float
                surface gravity (m/s^2)
        """
        if grav is None: raise RuntimeError('A planet needs a gravity!')
        self.grav=grav
        self.dmass=self.dp_lay/self.grav
        self.inv_dmass=1./self.dmass
    
    def set_gas(self, composition_dict, Mgas=None, compute_Mgas=True):
        """Sets the composition of the atmosphere.

        The composition_dict gives the composition in the layers, but we will
        need the composition in the radiative layers, so the interpolation is
        done here. For the moment we do a geometrical average.

        .. important::
            For the first initialization, compute_Mgas must be False
            because we need Gas_mix to be initialized before we know the number of layers
            in the atmosphere, but the number of layers is needed by set_Mgas!
            So set_Mgas needs to be called at the end of the initialization.

        Parameters
        ----------
            composition_dict: dictionary
                Keys are molecule names, and values are volume mixing ratios.
                A 'background' value means that the gas will be used to fill up to vmr=1
                If they do not add up to 1 and there is no background gas_mix,
                the rest of the gas_mix is considered transparent.
            compute_Mgas: bool
                If False, the molar mass of the gas is not updated. 
        """
        for mol, vmr in composition_dict.items():
            if isinstance(vmr,(np.ndarray, list)):
                tmp_vmr=np.array(vmr)
                #geometrical average:
                composition_dict[mol]=np.sqrt(tmp_vmr[1:]*tmp_vmr[:-1])
        if self.gas_mix is None:
            self.gas_mix=Gas_mix(composition_dict)
        else:
            self.gas_mix.set_composition(composition_dict)
        if compute_Mgas: self.set_Mgas(Mgas=Mgas)

    def set_Mgas(self, Mgas=None):
        """Sets the mean molar mass of the atmosphere.

        Parameters
        ----------
            Mgas_rad: float or array of size Nlay-1
                Mean molar mass in the radiative layers (kg/mol).
                If None is given, the mmm is computed from the composition.
        """
        if Mgas is not None:
            self.Mgas_rad=Mgas
        else:
            self.Mgas_rad=self.gas_mix.molar_mass()
        if not isinstance(self.Mgas_rad, np.ndarray):
            self.Mgas_rad=self.Mgas_rad*np.ones(self.Nlay-1, dtype=np.float)

    def set_rcp(self,rcp):
        """Sets the adiabatic index of the atmosphere

        Parameters
        ----------
            rcp: float
                R/c_p
        """
        self.rcp=rcp

    def set_aerosols(self, aerosols):
        """Sets the aerosols dictionary

        performs the interlayer averaging so that we only have properties at
        the middle of radiative layers
        """
        for aer, [reff, densities] in aerosols.items():
            if isinstance(reff,(np.ndarray, list)):
                tmp_reff=np.array(reff)
                aerosols[aer][0]=0.5*(tmp_reff[1:]+tmp_reff[:-1])
            if isinstance(densities,(np.ndarray, list)):
                tmp_densities=np.array(densities)
                #geometrical average:
                aerosols[aer][1]=np.sqrt(tmp_densities[1:]*tmp_densities[:-1])
        if self.aerosols is None:
            self.aerosols = Aerosols(aerosols)
        else:
            self.aerosols.set_aer_reffs_densities(aer_reffs_densities=aerosols)


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
        in the radiative layers
        """
        self.density=np.power(10., self.logp_opac)/(KBOLTZ*self.t_opac)

    def compute_layer_col_density(self):
        """Computes the column number density (molecules/m^2) per
        radiative layer of the atmosphere.

        There are Nlay-1 radiative layers as they go from the midle of a layer to the next.
        """
        factor=N_A/(self.grav * self.Mgas_rad)
        self.dcol_density_rad = np.diff(self.play)*factor[:]

        if self.Rp is not None: #includes the altitude effect if radius is known
            self.compute_altitudes()
            self.dcol_density_rad*=(1.+self.zlev[1:-1]/self.Rp)**2

    def compute_altitudes(self):
        """Compute altitudes of the level surfaces (zlev) and mid layers (zlay).
        """
        Mgas = np.empty(self.Nlay, dtype=np.float)
        Mgas[1:-1] = 0.5*(self.Mgas_rad[:-1]+self.Mgas_rad[1:])
        Mgas[0] = self.Mgas_rad[0]
        Mgas[-1] = self.Mgas_rad[-1]
        H = RGP*self.tlay/(self.grav*Mgas)
        dlnP = np.diff(self.logplev)*np.log(10.)
        self.zlev = np.zeros_like(self.logplev)
        if self.Rp is None:
            self.dz = H*dlnP
            self.zlev[:-1] = np.cumsum(self.dz[::-1])[::-1]
        else:
            for i in range(H.size)[::-1]:
                z1 = self.zlev[i+1]
                H1 = H[i]
                dlnp = dlnP[i]
                self.zlev[i] = z1+( (H1 * (self.Rp + z1)**2 * dlnp) \
                    / (self.Rp**2 + H1 * self.Rp * dlnp + H1 * z1 * dlnp) )
        self.zlay = 0.5*(self.zlev[1:]+self.zlev[:-1])
        self.zlay[-1] = 0.
        self.zlay[0] = self.zlev[0]
        ## assumes layer centers at the middle of the two levels
        ## which is not completely consistent with play, but should be
        ## a minor error.
        
    def compute_area(self):
        """Computes the area of the annulus covered by each
        radiative layer (from a mid layer to the next) in a transit setup. 
        """
        self.area=PI*(self.Rp+self.zlay[:-1])**2
        self.area[:-1]-=self.area[1:]
        self.area[-1]-=PI*self.Rp**2

    def compute_tangent_path(self):
        """Computes a triangular array of the tangent path length (in m) spent in each
        radiative layer.
        
        self.tangent_path[ilay][jlay] is the length that the ray that is tangent to the ilay 
        radiative layer spends in the jlay>=ilay layer
        (accounting for a factor of 2 due to symmetry)
        """
        if self.Rp is None: raise RuntimeError('Planetary radius should be set')
        self.compute_altitudes()
        self.tangent_path=List()
        # List() is a new numba.typed list to comply with new numba evolution after v0.50
        for ilay in range(self.Nlay-1): #layers counted from the top
            z0square=(self.Rp+self.zlev[ilay+1])**2
            dl=np.sqrt((self.Rp+self.zlay[:ilay+1])**2-z0square)
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

    def __init__(self, k_database=None, cia_database=None, a_database=None,
        wn_range=None, wl_range=None, internal_flux=0., **kwargs):
        """Initialization method that calls Atm_Profile().__init__() and links
        to Kdatabase and other radiative data. 
        """
        super().__init__(**kwargs)
        self.set_k_database(k_database)
        self.set_cia_database(cia_database)
        self.set_a_database(a_database)
        self.set_spectral_range(wn_range=wn_range, wl_range=wl_range)
        self.set_internal_flux(internal_flux)
        self.flux_net_nu=None
        self.kernel=None

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
        if self.kdatabase is not None:
            self.Nw=self.kdatabase.Nw
            self.flux_top_dw_nu = np.zeros((self.Nw))

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

    def set_a_database(self, a_database=None):
        """Change the Aerosol database used by the
        :class:`Aerosols` object handling aerosol optical properties.

        See :any:`aerosols.Aerosols.set_a_database` for details.

        Parameters
        ----------
            a_database: :class:`Adatabase` object
                New Adatabase to use.
        """    
        self.aerosols.set_a_database(a_database=a_database)

    def set_spectral_range(self, wn_range=None, wl_range=None):
        """Sets the spectral range in which computations will be done by specifying
        either the wavenumber (in cm^-1) or the wavelength (in micron) range.

        See :any:`gas_mix.Gas_mix.set_spectral_range` for details.
        """
        self.gas_mix.set_spectral_range(wn_range=wn_range, wl_range=wl_range)

    def set_incoming_stellar_flux(self, flux=0., Tstar=5778., **kwargs):
        """Sets the stellar incoming flux integrated in each wavenumber
        channel.

        .. important::
            If your simulated range does not include the whole spectral range
            where the star emits, the flux seen by the model will be smaller
            than the input one. 

        Parameters
        ----------
            flux: float
                Bolometric Incoming flux (in W/m^2).
            Tstar: float
                Stellar temperature (in K) used to compute the spectral distribution
                of incoming flux using a blackbody.

        """
        self.flux_top_dw_nu = Bnu_integral_num(self.wnedges, Tstar)
        factor = flux * PI / (SIG_SB*Tstar**4 * self.dwnedges)
        self.flux_top_dw_nu = self.flux_top_dw_nu * factor
    
    def set_internal_flux(self, internal_flux):
        """Sets internal flux from the subsurface in W/m^2
        """
        self.internal_flux = internal_flux

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

    def opacity(self, rayleigh = False, compute_all_opt_prop = False, **kwargs):
        """Computes the opacity of each of the radiative layers (m^2/molecule).

        Parameters
        ----------
            rayleigh: bool
                If true, the rayleigh cross section is computed in
                self.kdata_scat and added to kdata(total extinction cross section)

        See :any:`gas_mix.Gas_mix.cross_section` for details.
        """
        self.kdata = self.gas_mix.cross_section(rayleigh=rayleigh, **kwargs)
        shape = self.kdata.shape
        if self.aerosols.adatabase is not None:
            [kdata_aer, k_scat_aer, g_aer] = self.aerosols.optical_properties(
                    compute_all_opt_prop=compute_all_opt_prop, **kwargs)
            if self.Ng is None:
                self.kdata += kdata_aer
            else:
                self.kdata += kdata_aer[:,:,None]
        if compute_all_opt_prop:
            if rayleigh:
                kdata_scat_tot = self.gas_mix.kdata_scat
            else:
                kdata_scat_tot = np.zeros(shape[0:2], dtype=np.float)
            if self.aerosols.adatabase is not None:
                kdata_scat_tot += k_scat_aer
                self.asym_param = np.where(kdata_scat_tot<=0., 0., k_scat_aer * g_aer / kdata_scat_tot)
                #the line below is for test. Should be removed
                #self.asym_param = np.zeros(shape[0:2], dtype=np.float)
                self.asym_param = np.ones(shape[0:2], dtype=np.float)
            else:
                self.asym_param = np.zeros(shape[0:2], dtype=np.float)

            if self.Ng is None:
                self.single_scat_albedo = kdata_scat_tot / self.kdata
            else:
                self.single_scat_albedo = kdata_scat_tot[:,:,None] / self.kdata
                self.asym_param = self.asym_param[:,:,None] * np.ones(self.Ng)
                #JL21 the line above could be removed as asym param does not seem to change
                # with g point, but then the dimensions should be change in 2stream routines.
            self.single_scat_albedo=np.core.umath.minimum(self.single_scat_albedo,0.9999999999)
            #self.single_scat_albedo=np.core.umath.maximum(self.single_scat_albedo,0.1)

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
            source=True, gas_vmr=None, **kwargs):
        """Computes all necessary quantities for emission calculations
        (opacity, source, etc.)
        """
        if gas_vmr is not None: self.set_gas(gas_vmr)
        self.opacity(rayleigh=rayleigh, **kwargs)
        self.piBatm = self.source_function(integral=integral, source=source)
        self.compute_layer_col_density()
        if self.Ng is None:
            self.tau, self.dtau=rad_prop_xsec(self.dcol_density_rad,
                self.kdata, mu_eff)
        else:
            self.tau, self.dtau=rad_prop_corrk(self.dcol_density_rad,
                self.kdata, mu_eff)
            self.weights=self.kdatabase.weights

    def emission_spectrum(self, integral=True, mu0=0.5, mu_quad_order=None,
            dtau_min=1.e-13, **kwargs):
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
            dtau_min: float
                If the optical depth in a layer is smaller than dtau_min,
                dtau_min is used in that layer instead. Important as too
                transparent layers can cause important numerical rounding errors.

        Returns
        -------
            Spectrum object 
                A spectrum with the Spectral flux at the top of the atmosphere (in W/m^2/cm^-1)
        """
        if mu_quad_order is not None:
            # if we want quadrature, use the more general method.
            return self.emission_spectrum_quad(integral=integral,
                mu_quad_order=mu_quad_order, dtau_min=dtau_min, **kwargs)

        try:
            self.setup_emission_caculation(mu_eff=mu0, rayleigh=False, integral=integral, **kwargs)
        except TypeError:
            raise RuntimeError("""
            Cannot use rayleigh option with emission_spectrum.
            If you meant to include scattering, you should use emission_spectrum_2stream.
            """)
        # self.tau and self.dtau include the 1/mu0 factor.
        expdtau=np.exp(-self.dtau)
        expdtauminone=np.where(self.dtau < dtau_min, -self.dtau, expdtau-1.)
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

    def emission_spectrum_quad(self, integral=True, mu_quad_order=3, dtau_min=1.e-13, **kwargs):
        """Returns the emission flux at the top of the atmosphere (in W/m^2/cm^-1)
        using gauss legendre qudrature of order `mu_quad_order`

        Parameters
        ----------
            integral: boolean, optional
                * If true, the black body is integrated within each wavenumber bin.
                * If not, only the central value is used.
                  False is faster and should be ok for small bins,
                  but True is the correct version. 
            dtau_min: float
                If the optical depth in a layer is smaller than dtau_min,
                dtau_min is used in that layer instead. Important as too
                transparent layers can cause important numerical rounding errors.
                  
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
            expdtauminone=np.where(dtau<dtau_min,-dtau,expdtau-1.)
            exptau=np.exp(-tau)
            if self.Ng is None:
                timesBatmTop=(1.+expdtauminone/dtau)*exptau[:-1]
                timesBatmBottom=(-expdtau-expdtauminone/dtau)*exptau[:-1]
                timesBatmBottom[-1]+=exptau[-1]
            else:
                timesBatmTop=np.sum((1.+expdtauminone/dtau)*exptau[:-1] \
                    *self.weights,axis=-1)
                timesBatmBottom=np.sum((-expdtau-expdtauminone/dtau)*exptau[:-1] \
                    *self.weights,axis=-1)
                timesBatmBottom[-1]+=np.sum(exptau[-1]*self.weights,axis=-1)
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
            source=source, flux_top_dw=flux_top_dw, compute_all_opt_prop=True, **kwargs)
        # mu_eff=1. because the mu effect is taken into account in solve_2stream_nu
                 # we must compute the vertical optical depth here.
        # JL21: shouldn't we remove flux_top_dw, it seems not to be used.
        self.single_scat_albedo=np.where(self.dtau<dtau_min,0.,self.single_scat_albedo)
        self.dtau=np.where(self.dtau<dtau_min,dtau_min,self.dtau)

        module_to_use=globals()[method]
        # globals()[method] converts the method string into a module name
        #  if the module has been loaded
        if self.Ng is None:
            solve_2stream_nu=module_to_use.solve_2stream_nu_xsec
        else:
            solve_2stream_nu=module_to_use.solve_2stream_nu_corrk

        if flux_top_dw is not None:
            self.set_incoming_stellar_flux(flux=flux_top_dw, **kwargs)

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
        tlay=self.tlay
        dT = epsilon*tlay
        self.tlay = tlay + dT
        newpiBatm = self.source_function(integral=integral)

        for ilay in range(self.Nlay):
            pibatm = np.copy(self.piBatm)
            pibatm[ilay] = newpiBatm[ilay]
            _, _, flux_net_tmp = \
                solve_2stream_nu(pibatm, self.dtau, self.single_scat_albedo, self.asym_param,
                    self.flux_top_dw_nu, mu0 = mu0, flux_at_level=flux_at_level)
            net_tmp = self.spectral_integration(flux_net_tmp)
            self.kernel[ilay]=(net-net_tmp)/dT[ilay]
        self.kernel[:,:-1]-=self.kernel[:,1:]
        self.tlay=tlay
        self.tlay_kernel=self.tlay
        if per_unit_mass: self.kernel*=self.inv_dmass

    def flux_divergence(self, per_unit_mass = True,
            compute_kernel=False, **kwargs):
        """Computes the divergence of the net flux in the layers
        (used to compute heating rates).

        :func:`emission_spectrum_2stream` needs to be ran first.

        Parameters
        ----------
            per_unit_mass: bool
                If True, the heating rates are normalized by the
                mass of each layer (result in W/kg).

        Returns
        -------
            H: array
                Heating rate in each layer (Difference of the net fluxes). Positive means heating.
                The last value is the net flux impinging on the surface + the internal flux.
            net: array
                Net fluxes at level surfaces
        """
        if self.flux_net_nu is None:
            raise RuntimeError('should have ran emission_spectrum_2stream.')
        net = self.spectral_integration(self.flux_net_nu)
        H = -np.copy(net)
        H[:-1] -= H[1:]
        H[-1] += self.internal_flux
        if per_unit_mass: H *= self.inv_dmass
        if compute_kernel: self.H_kernel = H
        return H, net

    def heating_rate(self, compute_kernel=False, dTmax_use_kernel=None, **kwargs):
        if (not compute_kernel) and (dTmax_use_kernel is not None):
            dT=self.tlay-self.tlay_kernel
            if np.amax(np.abs(dT)) < dTmax_use_kernel:
                try:
                    H = self.H_kernel + np.dot(dT,self.kernel)
                except:
                    raise RuntimeError("Kernel has not been precomputed")
                net = np.zeros_like(H)
                return H, net
        _ = self.emission_spectrum_2stream(flux_at_level=True, integral=True,
                compute_kernel=compute_kernel, **kwargs)
        H, net = self.flux_divergence(compute_kernel=compute_kernel, **kwargs)
        if compute_kernel:
            self.H_kernel = H
            self.tau_rad = 1./np.amax(np.abs(self.kernel.diagonal()))
        return H, net


    def bolometric_fluxes(self, per_unit_mass = True):
        """Computes the bolometric fluxes at levels and the divergence of the net flux in the layers
        (used to compute heating rates).

        :func:`emission_spectrum_2stream` needs to be ran first.

        Parameters
        ----------
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
        H, net = self.flux_divergence(per_unit_mass = per_unit_mass)
        up=self.spectral_integration(self.flux_up_nu)
        dw=self.spectral_integration(self.flux_down_nu)
        return up, dw, net, H

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
                self.Nlay-1,self.Nw,self.Ng,self.tangent_path,self.density,self.kdata,self.weights)
        else:
            transmittance=path_integral_xsec( \
                self.Nlay-1,self.Nw,self.tangent_path,self.density,self.kdata)
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
            Rstar: float or astropy.unit object, optional
                Radius of the host star. If a float is specified, meters are assumed.
                Does not need to be given here if
                it has already been specified as an attribute of the :class:`Atm` object.
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
        
    def exp_minus_tau(self):
        """Sums Exp(-tau) over gauss points
        """
        weights=self.kdatabase.weights
        return np.sum(np.exp(-self.tau[1:])*weights,axis=2)

    def exp_minus_tau_g(self, g_index):
        """Sums Exp(-tau) over gauss point
        """
        return np.exp(-self.tau[1:,:,g_index])

    def blackbody(self, layer_idx=-1, integral=True):
        """Computes the surface black body flux (in W/m^2/cm^-1) for the temperature
        of layer `layer_idx`.

        Parameters
        ----------
            layer_idx; int
                Index of layer used for the temperature.
            integral: boolean, optional
                * If true, the black body is integrated within each wavenumber bin.
                * If not, only the central value is used.
                  False is faster and should be ok for small bins,
                  but True is the correct version. 
        Returns
        -------
            Spectrum object
                Spectral flux in W/m^2/cm^-1
        """
        if integral:
            piBatm=PI*Bnu_integral_num(self.wnedges,self.tlay[layer_idx])/np.diff(self.wnedges)
        else:
            piBatm=PI*Bnu(self.wns[:],self.tlay[layer_idx])
        return Spectrum(piBatm,self.wns,self.wnedges)

    def surf_bb(self, **kwargs):
        """Computes the surface black body flux (in W/m^2/cm^-1).

        See :any:`blackbody` for options.

        Returns
        -------
            Spectrum object
                Spectral flux of a bb at the surface (in W/m^2/cm^-1)
        """
        return self.blackbody(layer_idx=-1, **kwargs)

    def top_bb(self, **kwargs):
        """Computes the top of atmosphere black body flux (in W/m^2/cm^-1).

        See :any:`blackbody` for options.

        Returns
        -------
            Spectrum object
                Spectral flux of a bb at the temperature at the top of atmosphere (in W/m^2/cm^-1)
        """
        return self.blackbody(layer_idx=0, **kwargs)

