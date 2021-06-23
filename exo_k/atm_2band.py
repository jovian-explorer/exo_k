# -*- coding: utf-8 -*-
"""
@author: jeremy leconte
"""
import numpy as np
from .atm import Atm
from .gas_mix import Gas_mix
from .util.cst import PI, SIG_SB
from .util.radiation import Bnu_integral_num, rad_prop_corrk, rad_prop_xsec
from .two_stream import two_stream_toon as toon
from .two_stream import two_stream_lmdz as lmdz
from .util.spectrum import Spectrum

class Atm_2band(Atm):
    """Class based on :class:`Atm` that handles radiative trasnfer calculations
    when we want to use different radiative data to treat
    stellar absorption/scattering and the atmospheric emission.

    Only the method that change are overloaded.
    """

    def __init__(self, k_database_stellar=None, cia_database_stellar=None, **kwargs):
        """Initialization method that calls Atm_Profile().__init__() and links
        to Kdatabase and other radiative data. 
        """
        super().__init__(**kwargs)
        self.set_k_database_stellar(k_database_stellar)
        self.set_cia_database_stellar(cia_database_stellar)

    def set_k_database_stellar(self, k_database=None):
        """Change the radiative database used by the
        :class:`Gas_mix` object handling opacities inside
        :class:`Atm`.

        See :any:`gas_mix.Gas_mix.set_k_database` for details.

        Parameters
        ----------
            k_database: :class:`Kdatabase` object
                New Kdatabase to use.
        """
        self.gas_mix_stellar.set_k_database(k_database=k_database)
        self.kdatabase_stellar=self.gas_mix_stellar.kdatabase
        self.Ng_stellar=self.gas_mix_stellar.Ng
        # to know whether we are dealing with corr-k or not and access some attributes. 
        if self.kdatabase_stellar is not None:
            self.Nw_stellar=self.kdatabase_stellar.Nw
            self.flux_top_dw_nu_stellar = np.zeros((self.Nw_stellar))

    def set_cia_database_stellar(self, cia_database=None):
        """Change the CIA database used by the
        :class:`Gas_mix` object handling opacities inside
        :class:`Atm`.

        See :any:`gas_mix.Gas_mix.set_cia_database` for details.

        Parameters
        ----------
            cia_database: :class:`CIAdatabase` object
                New CIAdatabase to use.
        """
        self.gas_mix_stellar.set_cia_database(cia_database=cia_database)

    def set_gas(self, composition_dict, Mgas=None, compute_Mgas=True):
        """Sets the composition of the atmosphere.
        """
        for mol, vmr in composition_dict.items():
            if isinstance(vmr,(np.ndarray, list)):
                tmp_vmr=np.array(vmr)
                #geometrical average:
                composition_dict[mol]=np.sqrt(tmp_vmr[1:]*tmp_vmr[:-1])
        if self.gas_mix is None:
            self.gas_mix=Gas_mix(composition_dict)
            self.gas_mix_stellar=Gas_mix(composition_dict)
        else:
            self.gas_mix.set_composition(composition_dict)
            self.gas_mix_stellar.set_composition(composition_dict)
        if compute_Mgas: self.set_Mgas(Mgas=Mgas)

    def set_T_profile(self, tlay):
        """Reset the temperature profile without changing the pressure levels
        """
        tlay=np.array(tlay, dtype=float)
        if tlay.shape != self.logplay.shape:
            raise RuntimeError('tlay and logplay should have the same size.')
        self.tlay=tlay
        self.t_opac=(self.tlay[:-1]+self.tlay[1:])*0.5
        self.gas_mix.set_logPT(logp_array=self.logp_opac, t_array=self.t_opac)
        self.gas_mix_stellar.set_logPT(logp_array=self.logp_opac, t_array=self.t_opac)

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
        self.gas_mix_stellar.set_logPT(logp_array=self.logp_opac, t_array=self.t_opac)

    def spectral_integration_stellar(self, spectral_var):
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
        if self.Ng_stellar is None:
            var=np.sum(spectral_var*self.dwnedges_stellar,axis=-1)
        else:
            var=np.sum(np.sum(spectral_var*self.weights_stellar,axis=-1)*self.dwnedges_stellar,axis=-1)
        return var

    def opacity_stellar(self, rayleigh = False, **kwargs):
        """Computes the opacity of each of the radiative layers (m^2/molecule).

        Parameters
        ----------
            rayleigh: bool
                If true, the rayleigh cross section is computed in
                self.kdata_scat and added to kdata(total extinction cross section)

        See :any:`gas_mix.Gas_mix.cross_section` for details.
        """
        self.kdata_stellar = self.gas_mix_stellar.cross_section(rayleigh=rayleigh, **kwargs)
        if rayleigh: self.kdata_scat_stellar=self.gas_mix_stellar.kdata_scat
        self.Nw_stellar=self.gas_mix_stellar.Nw
        self.wns_stellar=self.gas_mix_stellar.wns
        self.wnedges_stellar=self.gas_mix_stellar.wnedges
        self.dwnedges_stellar=self.gas_mix_stellar.dwnedges

    def source_function_stellar(self, **kwargs):
        """Dummy function to have zero source in stellar channel
        """
        piBatm=np.zeros((self.Nlay,self.Nw_stellar))
        return piBatm

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
        self.flux_top_dw_nu_stellar = Bnu_integral_num(self.wnedges_stellar, Tstar)
        factor = flux * PI / (SIG_SB*Tstar**4 * self.dwnedges_stellar)
        self.flux_top_dw_nu_stellar = self.flux_top_dw_nu_stellar * factor


    def setup_emission_caculation_stellar(self, mu_eff=0.5, rayleigh=False,
            **kwargs):
        """Computes all necessary quantities for emission calculations
        (opacity, source, etc.)
        """
        self.opacity_stellar(rayleigh=rayleigh, **kwargs)
        self.piBatm_stellar = self.source_function_stellar()
        self.compute_layer_col_density() #done twice
        if self.Ng_stellar is None:
            self.tau_stellar, self.dtau_stellar=rad_prop_xsec(self.dcol_density_rad,
                self.kdata_stellar, mu_eff)
        else:
            self.tau_stellar, self.dtau_stellar=rad_prop_corrk(self.dcol_density_rad,
                self.kdata_stellar, mu_eff)
            self.weights_stellar=self.kdatabase_stellar.weights

    def reflexion_spectrum_2stream(self, integral=True, mu0_stellar=0.5,
            method='toon', dtau_min=1.e-10, flux_at_level=False, rayleigh=False,
            flux_top_dw=None, **kwargs):
        """Returns the reflexion flux at the top of the atmosphere (in W/m^2/cm^-1)

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
        self.setup_emission_caculation_stellar(mu_eff=1., rayleigh=rayleigh, integral=integral,
            flux_top_dw=flux_top_dw, **kwargs)
        # mu_eff=1. because the mu effect is taken into account in solve_2stream_nu
                 # we must compute the vertical optical depth here.
        self.dtau_stellar=np.where(self.dtau_stellar<dtau_min,dtau_min,self.dtau_stellar)

        module_to_use=globals()[method]
        # globals()[method] converts the method string into a module name
        #  if the module has been loaded
        if self.Ng_stellar is None:
            solve_2stream_nu=module_to_use.solve_2stream_nu_xsec
        else:
            solve_2stream_nu=module_to_use.solve_2stream_nu_corrk

        if rayleigh:
            if self.Ng_stellar is None:
                self.single_scat_albedo_stellar = self.kdata_scat_stellar / self.kdata_stellar
            else:
                self.single_scat_albedo_stellar = self.kdata_scat_stellar[:,:,None] / self.kdata_stellar
        else:
            self.single_scat_albedo_stellar = np.zeros_like(self.dtau)
        self.single_scat_albedo_stellar=np.clip(self.single_scat_albedo_stellar,None,0.9999999999999)
        self.asym_param_stellar = np.zeros_like(self.dtau_stellar)

        if flux_top_dw is not None:
            self.set_incoming_stellar_flux(flux=flux_top_dw, **kwargs)

        self.flux_up_nu_stellar, self.flux_down_nu_stellar, self.flux_net_nu_stellar = \
            solve_2stream_nu(self.piBatm_stellar, self.dtau_stellar,
                self.single_scat_albedo_stellar, self.asym_param_stellar,
                self.flux_top_dw_nu_stellar, mu0 = mu0_stellar, flux_at_level=flux_at_level)

        if self.Ng_stellar is None:
            return Spectrum(self.flux_up_nu_stellar[0],self.wns_stellar,self.wnedges_stellar)
        else:
            return Spectrum(np.sum(self.flux_up_nu_stellar[0]*self.weights_stellar,axis=1),
                        self.wns_stellar,self.wnedges_stellar)

    def flux_divergence_stellar(self, per_unit_mass = True, **kwargs):
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
                The last value is the net flux impinging on the surface.
            net: array
                Net fluxes at level surfaces
        """
        if self.flux_net_nu_stellar is None:
            raise RuntimeError('should have ran reflexion_spectrum_2stream.')
        net=self.spectral_integration_stellar(self.flux_net_nu_stellar)
        H=-np.copy(net)
        #print(H)
        H[:-1]-=H[1:]
        if per_unit_mass: H*=self.inv_dmass
        return H, net

    def heating_rate(self, compute_kernel=False, dTmax_use_kernel=None,
            flux_top_dw=None, **kwargs):
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
                compute_kernel=compute_kernel, flux_top_dw=None, **kwargs)
        H, net = self.flux_divergence(compute_kernel=compute_kernel, **kwargs)
        _ = self.reflexion_spectrum_2stream(flux_at_level=True, integral=True,
                flux_top_dw=flux_top_dw, **kwargs)
        H_stellar, net_stellar = self.flux_divergence_stellar(**kwargs)
        #print(H, H_stellar)
        if compute_kernel:
            self.H_kernel = H + H_stellar
            self.tau_rad = 1./np.amax(np.abs(self.kernel.diagonal()))
        return H+H_stellar, net+net_stellar


    def bolometric_fluxes_stellar(self, per_unit_mass = True):
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
                The last value is the net flux impinging on the surface.
        """
        H, net = self.flux_divergence_stellar(per_unit_mass = per_unit_mass)
        up=self.spectral_integration_stellar(self.flux_up_nu_stellar)
        dw=self.spectral_integration_stellar(self.flux_down_nu_stellar)
        return up, dw, net, H

    def __repr__(self):
        """Method to output header
        """
        output=super().__repr__()
        output+="""
    k_database_ste  :
        {kdatab}
    cia_database_ste:
        {cdatab}""".format(kdatab=self.kdatabase, cdatab=self.gas_mix.cia_database)
        if self.gas_mix._wn_range is not None:
            output+='    wn range        : '+ self.gas_mix._wn_range +'\n'

        return output
        
