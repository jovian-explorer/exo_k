# -*- coding: utf-8 -*-
"""
@author: jeremy leconte
"""
import pickle
import copy
import numpy as np
import numba
import exo_k as xk
from exo_k.util.cst import DAY
from .settings import Settings
from .convection import dry_convective_adjustment, turbulent_diffusion, molecular_diffusion
from .condensation import Condensing_species, moist_adiabat, compute_condensation, Tsat_P,\
                Condensation_Thermodynamical_Parameters 

class Atm_evolution(object):
    """Model of atmospheric evolution.

    Uses exo_k.Atm class to compute radiative transfer
    """

    def __init__(self, bg_vmr=None, verbose=False, **kwargs):
        """Initializes atmospheric profiles.

        Most arguments are passed directly to exo_k.Atm class through **kwargs

        .. warning::
            Layers are counted from the top down (increasing pressure order).
            All methods follow the same convention.
        """
        self.settings = Settings()
        self.settings.set_parameters(**kwargs)
        self.header={'rad':0,'conv':1,'cond':2,'madj':3,'rain':4,'tot':5}

        self.bg_gas = xk.Gas_mix(bg_vmr)
        self.M_bg = self.bg_gas.molar_mass()
        self.M_bg = self.settings.get('M_bg', self.M_bg)
        self.cp = self.bg_gas.cp()
        self.cp = self.settings.get('cp', self.cp)
        self.rcp = xk.RGP/(self.M_bg*self.cp)
        self.rcp = self.settings.get('rcp', self.rcp)
        if verbose: print('cp, M_bg, rcp:', self.cp, self.M_bg, self.rcp)


        self.tracers=Tracers(self.settings, bg_vmr = self.bg_gas.composition,
            M_bg = self.M_bg, **self.settings.parameters)
        self.initialize_condensation(**self.settings.parameters)

        self.setup_radiative_model(gas_vmr = self.tracers.gas_vmr,
            **self.settings.parameters)
        self.Nlay = self.atm.Nlay
        self.tlay = self.atm.tlay
        if verbose: print(self.settings.parameters)
        self.evol_time = 0.

    def set_options(self, reset_rad_model=False, **kwargs):
        """This method is used to store the global options
        in the `Settings` object.

        Arguments are all passed through **kwargs.

        Sometimes, one needs to reset the radiative model to take into
        account some modifications (like the databases). when such options are changed,
        one should set `reset_rad_model=True`
        """
        if 'tlay' not in kwargs.keys():
            self.settings.set_parameters(tlay=self.tlay, **kwargs)
        else:
            self.settings.set_parameters(**kwargs)
        if 'Kzz' in kwargs.keys():
            self.tracers.Kzz = np.ones(self.Nlay)*kwargs['Kzz']
        if reset_rad_model: self.setup_radiative_model(gas_vmr = self.tracers.gas_vmr,
                **self.settings.parameters)

    def initialize_condensation(self, condensing_species={}, **kwargs):
        """This method initializes the condensation module by
        listing all the condensing vapors and linking them to their
        condensed form. 

        For each vapor-condensate pair, a :class:`Condensible_species` object is created
        with the thermodynamical data provided. 

        Here is an example of dictinoary to provide as input to include CH4 condensation
        ```
        condensing_species={'ch4':{'Latent_heat_vaporization': 5.25e5, 'cp_vap': 2.232e3, 'Mvap': 16.e-3,
            'T_ref': 90., 'Psat_ref': 0.11696e5}}
        ```
        """
        self.condensing_pairs = list()
        self.condensing_pairs_idx = list()
        self.condensing_species_idx = dict()
        self.condensing_species_params=list()
        self.condensing_species_thermo=list()
        idx=0
        for name in self.tracers.namelist:
            if 'type' in self.tracers.dico[name]:
                if self.tracers.dico[name]['type'] == 'vapor':
                    if 'condensed_form' not in self.tracers.dico[name]:
                        print("You should identify the 'condensed_form' of:", name)
                        raise RuntimeError()
                    elif self.tracers.dico[name]['condensed_form'] not in self.tracers.namelist:
                        print("The condensed form of a vapor should be a tracer.")
                        raise RuntimeError()
                    else:
                        if name in condensing_species.keys():
                            cond_name = self.tracers.dico[name]['condensed_form']
                            self.condensing_species_idx[name]=idx
                            self.condensing_pairs.append([name, cond_name])
                            self.condensing_pairs_idx.append(\
                                [self.tracers.idx[name], self.tracers.idx[cond_name]])
                            self.condensing_species_params.append(\
                                Condensing_species(**condensing_species[name]))
                            self.condensing_species_thermo.append(\
                                Condensation_Thermodynamical_Parameters(**condensing_species[name]))
                            idx+=1
        self.Ncond=idx

    def condensation(self, timestep, Htot):
        """This method computes the vapor and temperature tendencies do to
        large scale condensation in saturated layers.

        The tracer array in modified in place.

        Parameters
        ----------
            timestep: float
                physical timestep of the current step (in s/cp).
            Htot: array
                Total heating rate (in W/kg) of all physical processes
                already computed

        Return
        ------
            H_cond: array
                Heating rate due to large scale condensation (W/kg)
        """
        new_t = self.atm.tlay + timestep * Htot
        H_cond = np.zeros(self.Nlay)
        for i_cond in range(self.Ncond): #careful i_cond is a dumy loop index, idx_cond is position of species i_cond in tracers array.
            idx_vap, idx_cond = self.condensing_pairs_idx[i_cond]
            thermo_parameters = self.condensing_species_thermo[i_cond].th_params
            Lvap, qsat, dqsat_dt = compute_condensation(new_t, self.atm.play, self.tracers.Mgas, *thermo_parameters[1:])
            qvap=self.tracers.qarray[idx_vap]
            if self.settings['latent_heating']:
                H_cond += np.where(qsat <= qvap, Lvap * (qvap-qsat) \
                    / ((self.cp+Lvap*dqsat_dt)*timestep), 0.)
                dqvap= np.where(qsat <= qvap, qsat+dqsat_dt*H_cond*timestep-qvap, 0.)
            else:
                dqvap= np.where(qsat <= qvap, qsat-qvap, 0.)
            
            self.tracers.qarray[idx_vap] += dqvap
            self.tracers.qarray[idx_cond] -= dqvap
#            if self.settings['latent_heating']:
#                H_cond += - cond_species.Lvap(new_t) * dqvap / (timestep *self.cp)
        return H_cond

    def rainout(self, timestep, Htot):
        """This method computes rainout.

        For the moment, all condensates are taken down and reevaporated in the
        bottom layer.

        The tracer array in modified in place.

        Parameters
        ----------
            timestep: float
                physical timestep of the current step (in s/cp).
            Htot: array
                Total heating rate (in W/kg) of all physical processes
                already computed

        Return
        ------
            H_rain: array
                Heating rate due to re evaporation (W/kg)
        """
        new_t = self.atm.tlay + timestep * Htot
        H_rain=np.zeros(self.Nlay)
        for i_cond in range(self.Ncond):
            idx_vap, idx_cond = self.condensing_pairs_idx[i_cond]
            thermo_parameters = self.condensing_species_thermo[i_cond].th_params
            H_rain += compute_rainout(timestep, self.Nlay, new_t, self.atm.play,
                self.atm.dmass, self.cp, self.tracers.Mgas, self.tracers.qarray,
                idx_vap, idx_cond, thermo_parameters,
                self.settings['evap_coeff'], self.settings['qvap_deep'])

        return H_rain


    def setup_radiative_model(self, k_database=None, k_database_stellar=None,
            cia_database=None, cia_database_stellar=None, gas_vmr=None, **kwargs):
        """This method initializes the exo_k.Atm object that will be used
        to carry out radiative transfer calculations. 

        This is where the radiative data used are chosen and transmitted to the 
        radiative transfer module, along with many other
        parameters including the incoming stellar flux (`flux_top_dw`), the
        blackbody temperature of the star (`Tstar`), the 

        If a `k_database_stellar` is provided, then this is this database
        that will be used to treat the scattering and absorption of incoming radiation.
        In that case, `k_database` will be used to treat the emission of the atmosphere.
        The effective cos(zenith angle) for the incoming stellar radiation can then
        be specified independently with the `mu0_stellar` keyword.

        If no `k_database_stellar` is provided, `k_database` will be used to treat
        both the atmospheric emission and the stellar radiation. Running a model
        with `k_database_stellar=k_database` yields the same results at twice the cost.
        
        Parameters
        ----------
            k_database, k_database_stellar: `exo_k.Kdatabase` objects
                radiative database for the molecules in the atmospheres.
            cia_database, cia_database_stellar: `exo_k.CIA_database` object
                radiative database for the CIA ofmolecules in the atmospheres.
        """
        if k_database is None:
            raise RuntimeError('We need at least a k_database')
        if k_database_stellar is None:
            self.atm=xk.Atm(k_database=k_database, cia_database=cia_database, composition=gas_vmr, **kwargs)
        else:
            self.atm=xk.Atm_2band(k_database=k_database, cia_database=cia_database,
                k_database_stellar=k_database_stellar, cia_database_stellar=cia_database_stellar,
                composition=gas_vmr, **kwargs)
        H, net = self.atm.heating_rate(compute_kernel=True, **kwargs)

    def compute_average_fluxes(self):
        """Use the averaged heating rates to compute the various fluxes (W/m^2)
        at the level interfaces. These fluxes are positive when the energy flows
        upward.

        To be consistent with radiative fluxes, the first array value corresponds
        to the top of atmosphere and should be 0 in most cases. The last value corresponds
        to the flux between the deepest layer (considered to be the surface) and the layer just above.
        """
        self.Fnet = np.zeros((6, self.Nlay))
        self.Fnet[0] = self.Fnet_rad
        for ii in range(1,5):
            self.Fnet[ii]=np.concatenate([[0.],
                np.cumsum(self.H_ave[ii]*self.atm.dmass)[:-1]])
        self.Fnet[-1] = np.sum(self.Fnet, axis=0)


    def evolve(self, N_timestep=1, N_kernel=10000, timestep_factor=1., dT_max = 100., verbose = False, **kwargs):
        r"""The time variable used in the model is tau=t/cp. 
        The equation we are solving in each layer is thus

        .. math::
            c_p \frac{d T}{d t}= \frac{d T}{d tau} = \sum_i H_i

        For a given timestep `dtau`, the physical time elapsed in second can be computed using `dt=dtau*cp`

        To work, the heating rates (H) must be computed in W/kg. 

        This also means that if one needs the physical rate of change of another quantity (like dq/dt)
        from the delta q over a time step,
        one needs to do `dq/dt = delta q / (timestep * cp)`

        Parameters
        ----------
            N_timestep: int
                Number of timesteps to perform.
            N_kernel: int
                Maximal number of timesteps between two computations of the radiative kernel.
            timestep_factor: float
                Multiplicative factor applied to timestep computed automatically by
                the radiative module.
                timestep_factor > 1 can lead to unstabilities.
            dT_max: float
                Maximum temperature increment in a single timestep.
            
        """
        self.H_ave = np.zeros((6, self.Nlay))

        self.tlay_hist = np.zeros((N_timestep,self.Nlay))
        self.Fnet_top = np.zeros((N_timestep))
        self.timestep_hist = np.zeros((N_timestep))
        tau0 = self.evol_time
        self.N_last_ker = 0
        compute_kernel = False
        dTlay_max = 2. * self.settings['dTmax_use_kernel']
        self.tracers.update_gas_composition(update_vmr=True)
        for ii in range(N_timestep):
            if np.amax(np.abs(self.tlay-self.atm.tlay_kernel)) < self.settings['dTmax_use_kernel']:
                self.N_last_ker +=1
                if verbose: print(self.N_last_ker, self.N_last_ker%N_kernel)
                compute_kernel = (self.N_last_ker%N_kernel == 0)
            else:
                if dTlay_max < 0.5 * self.settings['dTmax_use_kernel']:
                    compute_kernel = True
                    self.N_last_ker = 0
                else:
                    compute_kernel = False
                    self.N_last_ker +=1
            if ii == N_timestep-1:
                if ii != 0:
                    compute_kernel=True
            self.H_tot=np.zeros(self.Nlay)
            if verbose: print('iter, compute_kernel:', ii, compute_kernel)
            if self.tracers.some_var_gases:
                gas_vmr_rad = self.tracers.gas_vmr
            else:
                gas_vmr_rad = None
            self.H_rad, self.Fnet_rad = self.atm.heating_rate(compute_kernel = compute_kernel,
                rayleigh = self.settings['rayleigh'], dTmax_use_kernel=self.settings['dTmax_use_kernel'],
                gas_vmr = gas_vmr_rad, **kwargs)
#            if verbose and compute_kernel: print('H_rad', self.H_rad)
            self.H_tot += self.H_rad
            self.timestep = timestep_factor * self.atm.tau_rad
            #if verbose: print('tau_rad, dt:', self.atm.tau_rad, self.timestep)
            self.evol_time += self.timestep
            if self.settings['convection']:
                self.H_conv = self.tracers.dry_convective_adjustment(self.timestep, self.H_tot, self.atm, verbose=verbose)
                self.H_tot += self.H_conv
            else:
                self.H_conv = np.zeros(self.Nlay)
            if self.settings['diffusion']:
                self.tracers.turbulent_diffusion(self.timestep, self.H_tot, self.atm, self.cp)
                self.tracers.update_gas_composition(update_vmr=False)
            if self.settings['molecular_diffusion']:
                self.H_diff = self.molecular_diffusion(self.timestep,
                    self.H_tot, self.atm, self.cp)
                self.H_tot += self.H_diff
            if self.settings['moist_convection']:
                self.H_madj = self.moist_convective_adjustment(self.timestep, self.H_tot,
                                    moist_inhibition=self.settings['moist_inhibition'], verbose=verbose)                
                self.H_tot += self.H_madj
            else:
                self.H_madj = np.zeros(self.Nlay)
            if self.settings['condensation']:
                self.H_cond = self.condensation(self.timestep, self.H_tot)
                self.H_tot += self.H_cond
            else:
                self.H_cond = np.zeros(self.Nlay)
            if self.settings['rain']:
                self.H_rain = self.rainout(self.timestep, self.H_tot)
                self.H_tot += self.H_rain
            else:
                self.H_rain = np.zeros(self.Nlay)
            if self.settings['radiative_acceleration']:
                self.acceleration_factor = self.radiative_acceleration()
                self.H_tot *= self.acceleration_factor
            dTlay= self.H_tot * self.timestep
            dTlay_max = np.amax(np.abs(dTlay))
            if dTlay_max > dT_max:
                print('dT > dTmax:',dTlay_max,' at k=',np.argmax(np.abs(dTlay))) 
            if verbose:
                print('heat rates (rad, dry conv), dTmax:', 
                    np.sum(self.H_rad*self.atm.dmass), np.sum(self.H_conv*self.atm.dmass),
                    dTlay_max)
            dTlay=np.clip(dTlay,-dT_max,dT_max)
            self.tlay = self.tlay + dTlay
            self.tlay_hist[ii] = self.tlay
            for jj, H in enumerate([self.H_rad, self.H_conv, self.H_cond, self.H_madj,
                  self.H_rain, self.H_tot]):
                self.H_ave[jj] += H * self.timestep
            self.Fnet_top[ii] = self.Fnet_rad[0]
            self.timestep_hist[ii] = self.timestep
            self.atm.set_T_profile(self.tlay)
            self.tracers.update_gas_composition(update_vmr=True)
        inv_delta_t = 1./(self.evol_time-tau0)
        self.H_ave *= inv_delta_t
        self.compute_average_fluxes()

    def equilibrate(self, Fnet_tolerance = None, N_iter_max = 10,
        N_timestep_ini = 100, N_timestep_max = 1000000, verbose = False, **kwargs):
        """Evolves an atmosphere until it is at equilibrium.
        
        Equilibrium is assumed to be reached when the net top of atmosphere
        flux remains within +-Fnet_tolerance of the internal flux
        for a whole evolution step.

        The number of timesteps per evolution step in multiplied by two 
        at each iteration, starting from N_timestep_ini, until the limit of
        N_timestep_max is reached.

        Parameters
        ----------
            Fnet_tolerance: float
                Tolerance on net flux in W/m^2 to identify convergence.
            N_iter_max: int
                Max number of successive calls to evolve
            N_timestep_ini: int
                Initial number of timesteps in a single evolution step
            N_timestep_max: int
                Max number of timesteps in a single evolution step
        """
        iter=1
        if Fnet_tolerance is None:
            raise RuntimeError('You should provide the maximum tolerance on the net flux: Fnet_tolerance (in W/m^2)') 
        N_timestep = N_timestep_ini
        while iter <= N_iter_max:
            time_init = self.evol_time
            self.evolve(N_timestep = N_timestep, **kwargs)
            net = self.Fnet_top - self.atm.internal_flux
            if verbose:
                print('iter: {iter}, N_timestep: {Nt}'.format(iter = iter, Nt = N_timestep))
                print('Fnet mean: {fme:.3g} W/m^2, (min:{fmi:.3g}, max:{fma:.3g})'.format( \
                    fme = net.mean(), fmi = net.min(), fma = net.max()))   
                print('timestep: {ts1:.3g} d | {ts2:.3g} s, total time: {ev_t:.3g} yr'.format( \
                    ts1 = self.timestep*self.cp/(DAY), ts2 = self.timestep*self.cp,
                    ev_t = (self.evol_time-time_init)*self.cp/(DAY*365.)))   
                #print('timestep:',self.timestep*self.cp/(DAY),'days, ',
                #    self.timestep*self.cp,'s, evol_time(yr):',(self.evol_time-time_init)*self.cp/(DAY*365.))
            if np.all(np.abs(net) < Fnet_tolerance): break
            N_timestep = min( N_timestep * 2, N_timestep_max)
            iter += 1
            
    def moist_convective_adjustment(self, timestep, Htot, moist_inhibition = True,
            verbose = False):
        """This method computes the vapor and temperature tendencies do to
        moist convectoin in saturated layers.

        The tracer array in modified in place.

        Parameters
        ----------
            timestep: float
                physical timestep of the current step (in s/cp).
            Htot: array
                Total heating rate (in W/kg) of all physical processes
                already computed

        Return
        ------
            H_madj: array
                Heating rate due to large scale condensation (W/kg)
        """
        new_t = self.atm.tlay + timestep * Htot
        H_madj = np.zeros(self.Nlay)
        for i_cond in range(self.Ncond): #careful i_cond is the index of the condensing pair
            # in the list of condensing species, idx_cond is the position of the
            # condensate linked to i_cond in the tracers array.
            idx_vap, idx_cond = self.condensing_pairs_idx[i_cond]
            thermo_parameters = self.condensing_species_thermo[i_cond].th_params
            H, qarray, new_t = moist_convective_adjustment(timestep, self.Nlay,
                new_t, self.atm.play, self.atm.dmass, self.cp, self.tracers.Mgas, self.tracers.qarray, idx_vap, idx_cond,
                thermo_parameters,
                moist_inhibition = moist_inhibition, verbose = verbose)
            #print('t after madj:', new_t)
            H_madj += H
            self.tracers.qarray = qarray
            if verbose: print(qarray[idx_cond])
        return H_madj

    def molecular_diffusion(self, timestep, Htot, atm, cp):
        """Mixes energy following a diffusion equation
        with a constant Dmol parameter (self.Dmol in m^2/s).

        Parameters
        ----------
            timestep: float
                physical timestep of the current step (in s/cp).
                (needs to be converted before it is sent to `turbulent diffusion)
            Htot: array
                Total heating rate (in W/kg) of all physical processes
                already computed
            atm: :class:`Atm` object
                The Atm object used in the radiative transfer which
                contains many state variables. 
        """
        new_t = atm.tlay + timestep * Htot
        H_diff = molecular_diffusion(timestep*cp, self.Nlay,
                    atm.play, atm.plev,
                    atm.dmass, new_t, self.tracers.Mgas,
                    atm.grav, self.tracers.Dmol)
        return H_diff

    def radiative_acceleration(self):
        """"Computes an acceleration factor to speed up convergence in
        layers where only radiation plays a role.
        """
        acceleration_factor = np.ones_like(self.H_tot)
        rad_layers =  np.isclose(self.H_tot, self.H_rad, atol=0.e0, rtol=1.e-10)
        if np.all(rad_layers):
            self.base_timescale = self.atm.tau_rad
        else:
            self.base_timescale = np.amax(self.atm.tau_rads[np.logical_not(rad_layers)])
        acceleration_factor[rad_layers] = np.core.umath.maximum(
            self.atm.tau_rads[rad_layers] \
            / self.base_timescale * self.settings['radiative_acceleration_reducer'],
            1.)
        return acceleration_factor

    def write_pickle(self, filename, data_reduction_level = 1):
        """Saves the instance in a pickle file

        Parameters
        ----------
            filename: str
                Path to pickle file
            data_reduction_level: int
                Level of data to delete.
                0: keep everything (results in big files).
                1: removes some arrays, should not affect subsequent evolution.
                2: removes the k and cia databases. The radiative model will need to be reset.
                This can be done with
                `set_options(k_database=, cia_database=, reset_rad_model=True)` after
                re-loading the `Atm_evolution` instance.
        """
        other = copy.deepcopy(self)
        if data_reduction_level >=1 :
            other.tlay_hist = None
            other.atm.asym_param = None
            other.atm.kdata = None
            other.atm.tau = None
            other.atm.dtau = None
            other.atm.flux_down_nu = None
            other.atm.flux_net_nu = None
            other.atm.flux_up_nu = None
            other.atm.piBatm = None
            other.atm.single_scat_albedo = None
            other.atm.gas_mix.kdata_scat = None
        if data_reduction_level >=2 :
            other.settings['k_database'] = None
            other.settings['cia_database'] = None
            other.atm.k_database = None
            other.atm.gas_mix.k_database = None
            other.atm.gas_mix.cia_database = None
            other.atm.kernel = None
            other.atm.tlay_kernel = None
            other.atm.H_kernel = None
        with open(filename, 'wb') as filehandler:
            pickle.dump(other, filehandler)


@numba.jit(nopython=True, fastmath=True, cache=True)
def moist_convective_adjustment(timestep, Nlay, tlay, play, dmass, cp, Mgas, q_array,
        i_vap, i_cond, thermo_parameters, 
        moist_inhibition = True, verbose = False):
    r"""Computes the heating rates needed to adjust unstable regions 
    of a given atmosphere to a moist adiabat.

    Parameters
    ----------
        timestep
        Nlay: float
            Number of layers
        tlay: array
            Layer temperatures
        play:array
            Pressure at layer centers
        dmass: array
            mass of layers in kg/m^2
        cp: float
            specific heat capacity at constant pressure
        q_array: array
            mass mixing ratio of tracers
        i_vap: int
            index of vapor tracer in qarray
        i_cond: int
            index of condensate tracer in qarray
        qsat: array
            Saturation mmr for each layer
        dqsat_dt: array
            d qsat / dT in each layer
        Lvap: array
            Latent heat of vaporization (can have different values
            in each layer if Lvap=f(T))
        dlnt_dlnp_moist: array
            threshold thermal gradient (d ln T / d ln P) for a moist atmosphere
            computed at the layer centers.
        q_crit: array
            Critical mass mixing ratio for the inhibition of moist convection
            (Eq. 17 of Leconte et al. 2017)

    Returns
    -------
        H_madj: array
            Heating rate in each atmospheric layer (W/kg). 
        new_q: array
            tracer mmr array after adjustment.
        new_t: array
            Temperature of layers after adjustment. 
    """
    if verbose: print('enter moist adj')
    H_madj=np.zeros(Nlay)
    new_q = q_array.copy()
    new_t = tlay.copy()
    dp = np.diff(play)

    dlnt_dlnp_moist, Lvap, psat, qsat, dqsat_dt, q_crit = \
        moist_adiabat(new_t, play, cp, Mgas, thermo_parameters[0],
            thermo_parameters[1], thermo_parameters[2], thermo_parameters[3], thermo_parameters[4],
            thermo_parameters[5], thermo_parameters[6], thermo_parameters[7], thermo_parameters[8],
            thermo_parameters[9])

    nabla_interlayer = tlay * dlnt_dlnp_moist /play
    nabla_interlayer = 0.5*(nabla_interlayer[:-1]+nabla_interlayer[1:])
    dTmoist_array = nabla_interlayer * dp
    dT_inter_lay = np.diff(tlay)
    qvap = new_q[i_vap]
    mvap_sursat_array = (qvap-qsat) * dmass
    if moist_inhibition:
        q_crit_criterion = qvap<q_crit # convection possible if True
    else:
        q_crit_criterion = qvap<2. #should always be true
    #print('dT:', dT_inter_lay)
    #print('dTmoist:', dTmoist_array)
    #dT_unstab = np.nonzero(dT_inter_lay>dTmoist_array)[0]
    #saturated = np.nonzero(mvap_sursat_array>0.)[0]
    conv = np.nonzero((dT_inter_lay>dTmoist_array)*(mvap_sursat_array[:-1]>0.) \
            *q_crit_criterion[:-1])[0]# find convective layers
    if verbose: 
        print(conv)
        print(np.nonzero(dT_inter_lay>dTmoist_array)[0])
        print(np.nonzero(mvap_sursat_array[:-1]>0.)[0])
        print(np.nonzero(q_crit_criterion)[0])
    N_conv=conv.size
    if N_conv==0: # no more convective regions, can exit
        return H_madj, new_q, new_t
    i_top=conv[0] #upper unstable layer
    T_top = new_t[i_top]
    mvap_sursat = mvap_sursat_array[i_top]
    dqsdm = dqsat_dt[i_top]*dmass[i_top]
    int_dqsdm = dqsdm
    C = cp*dmass[i_top] + Lvap[i_top]*dqsdm
    B = C * new_t[i_top] + Lvap[i_top] * mvap_sursat_array[i_top]
    dT_moist = 0.
    int_m_cond = mvap_sursat_array[i_top] + dqsdm*(new_t[i_top] - dT_moist)
    i_bot=i_top+1
    while i_bot<Nlay: #search for the bottom of the 1st unstable layer from its top
        tmp_sursat = mvap_sursat + mvap_sursat_array[i_bot]
        tmp_dT_moist = dT_moist + dTmoist_array[i_bot-1]
        dqsdm = dqsat_dt[i_bot] * dmass[i_bot]
        tmp_int_dqsdm = int_dqsdm + dqsdm
        tmp_int_m_cond = int_m_cond + mvap_sursat_array[i_bot] + dqsdm * (new_t[i_bot] - tmp_dT_moist)
        tmp = cp *dmass[i_bot] + Lvap[i_bot]* dqsdm
        tmp_C = C + tmp
        tmp_B = B + tmp * (new_t[i_bot]-tmp_dT_moist) + Lvap[i_bot] * mvap_sursat_array[i_bot]
        tmp_new_Ttop = tmp_B / tmp_C
        tmp_m_cond = tmp_int_m_cond - tmp_int_dqsdm * tmp_new_Ttop
        if tmp_sursat>0. and tmp_dT_moist<new_t[i_bot]-T_top and q_crit_criterion[i_bot] and tmp_m_cond>0.:
            dT_moist = tmp_dT_moist
            mvap_sursat = tmp_sursat
            int_dqsdm = tmp_int_dqsdm
            int_m_cond = tmp_int_m_cond
            C = tmp_C
            B = tmp_B
            m_cond = tmp_m_cond
            i_bot += 1
            continue
        else:
            i_bot -= 1
            break
    if verbose: print('it,ib=', i_top, i_bot)
    if i_top == i_bot: # need at least 2 layers to convect, so exit
        return H_madj, new_q, new_t
    new_Ttop = B / C
    if verbose: print(new_Ttop, m_cond, dT_moist)
    dT = new_Ttop - new_t[i_top]
    qvap[i_top] = qsat[i_top] + dqsat_dt[i_top] * dT
    #if verbose: print('top i, dT, qv, qs, dqs, qf', i_top, dT, q_array[i_vap, i_top], qsat[i_top], dqsat_dt[i_top], qvap[i_top])
    new_t[i_top] = new_Ttop
    H_madj[i_top] = dT / timestep
    for ii in range(i_top+1, i_bot+1):
        dT = new_t[ii-1] + dTmoist_array[ii-1] - new_t[ii]
        #print(ii, new_t[ii-1], dTmoist_array[ii-1],  new_t[ii], new_t[ii-1] + dTmoist_array[ii-1] - new_t[ii])
        qvap[ii] = qsat[ii] + dqsat_dt[ii] * dT
        new_t[ii] += dT
        # compute heating and adjust before looking for a new potential unstable layer
        H_madj[ii] = dT / timestep
        #if verbose: print('i, dT, qv, qs, dqs, qf', ii, dT, q_array[i_vap, ii], qsat[ii], dqsat_dt[ii], qvap[ii])
    # put ice
    m_cond_2 = np.sum((q_array[i_vap, i_top:i_bot+1]-qvap[i_top:i_bot+1])*dmass[i_top:i_bot+1])
    dTmoist_array[i_top-1]=0.
    if m_cond<0.:
        print('Negative condensates in moist adj, i:', i_top, i_bot)
        print(q_array[i_vap, i_top:i_bot+1], qvap[i_top:i_bot+1], q_array[i_vap, i_top:i_bot+1]-qvap[i_top:i_bot+1])
    m_conv = np.sum(dmass[i_top:i_bot+1])
    new_q[i_cond, i_top:i_bot+1] += m_cond / m_conv
    if verbose: 
        print('m_cond, m_conv, m_cond2', m_cond, m_conv, m_cond_2)
    return H_madj, new_q, new_t

@numba.jit(nopython=True, fastmath=True, cache=True)
def compute_rainout(timestep, Nlay, tlay, play, dmass, cp, Mgas, qarray,
        idx_vap, idx_cond, thermo_parameters, evap_coeff, qvap_deep,
        verbose = False):
    r"""Computes the heating rates needed to adjust unstable regions 
    of a given atmosphere to a moist adiabat.

    Parameters
    ----------
        timestep
        Nlay: float
            Number of layers
        tlay: array
            Layer temperatures
    """
    H_rain=np.zeros(Nlay)
    Lvap, qsat, dqsat_dt = compute_condensation(tlay, play, Mgas, 
            thermo_parameters[1], thermo_parameters[2], thermo_parameters[3], thermo_parameters[4],
            thermo_parameters[5], thermo_parameters[6], thermo_parameters[7], thermo_parameters[8],
            thermo_parameters[9])

    Tsat_p = Tsat_P(play, thermo_parameters[5], thermo_parameters[8], thermo_parameters[9])

    mass_cond = 0.
    for i_lay in range(Nlay):
        #  if evap_coeff =1, rain vaporisation in an undersaturated layer can fill the layer up to the (estimated) saturation
        #  if 0 < evap_coeff < 1, rain vaporisation in one layer is limited to a fraction of the amount that would saturate the layer
        #  This allows not to exceed saturation, to spread rain vaporization in more and denser layers
        
        dqvap = evap_coeff*(qsat[i_lay] - qarray[idx_vap,i_lay])/(1 + Lvap[i_lay]*dqsat_dt[i_lay]/cp)
        mass_cond += qarray[idx_cond,i_lay]*dmass[i_lay]
        qarray[idx_cond,i_lay] = 0.
        # dqvap is the change in vapor content to reach saturation and accounting for the temperature change.
        # dqvap < 0 implies condensation, meaning that there is a remaining excess of vapor after the previous 
        # condensation step. In such case we apply this new change in vapor content and temperature and
        # increase the amount of falling condensed species (mass_cond).
        # dqvap > 0 implies evaporation. Here there are two possibilities:
        #   - the amount of condensates is lower than dqvap. All condensates are vaporised in the layer
        #    and the tempertaure change is -Ldqv/cp where dqv is the actual change in vapor.
        #   - the amount of condensates is larger than dqvap. We then apply dqvap and the corresponding 
        #    change in temperature and transfer the remaining condensate in the falling rain reservoir.
        if dqvap < 0: # more condensation in the layer
            qarray[idx_vap][i_lay] += dqvap
            H_rain[i_lay] = -Lvap[i_lay]*dqvap/(cp*timestep)
            mass_cond -= dqvap*dmass[i_lay]               
        else: # evaporation of rain
            mass_dvap = dqvap*dmass[i_lay]     
            if (mass_dvap > mass_cond) or (tlay[i_lay] >= Tsat_p[i_lay]): # evaporate everything
                qarray[idx_vap,i_lay] += mass_cond/dmass[i_lay]
                H_rain[i_lay] = - Lvap[i_lay] * mass_cond / (dmass[i_lay]*cp*timestep)
                mass_cond = 0.
            else:
                qarray[idx_vap,i_lay] += dqvap
                H_rain[i_lay] = -Lvap[i_lay]*dqvap/(cp*timestep)
                mass_cond -= dqvap*dmass[i_lay]
    if mass_cond != 0.:
        qarray[idx_cond,-1]+= mass_cond/dmass[-1]
        # Issue : qarray[idx_cond,-1] sometimes can become large (when starting with a hot atmosphere dominated with H2O that cools to saturation)
    if qvap_deep>=0.:
        qarray[idx_vap,-1] = qvap_deep
    return H_rain



class Tracers(object):

    def __init__(self, settings, tracers={}, tracer_values={}, Kzz=0., Dmol=0.,
            bg_vmr=None, M_bg=None, Nlay=None, **kwargs):
        """Deals with tracers. 

        Fills out the tracers.qarray and creates tracer_names, a table of correspondence
        between tracer names and indices in tracers.qarray.
        """
        self.settings=settings
        self.Ntrac = len(tracers)
        if self.Ntrac == 0:
            self.Ntrac=1
            tracers={'inactive_tracer':{}}
        if Nlay is not None:
            self.Nlay = Nlay
        else:
            raise RuntimeError("We need to know Nlay to initialize Tracers")
        self.bg_vmr = bg_vmr.copy()
        self.gas_vmr = bg_vmr.copy()
        self.var_gas_idx = list()
        self.var_gas_names = list()
        self.gas_molar_masses = list()
        self.some_var_gases = False
        self.dico=tracers
        self.namelist = list(tracers.keys())
        self.idx = dict()
        self.qarray = np.empty((self.Ntrac, self.Nlay))
        for ii, name in enumerate(self.namelist):
            self.idx[name]=ii
            if name in tracer_values.keys():
                self.qarray[ii]=np.copy(tracer_values[name])
            elif 'init_value' in self.dico[name].keys():
                self.qarray[ii]=np.ones(self.Nlay)*self.dico[name]['init_value']
            else:
                self.qarray[ii]=np.zeros(self.Nlay)
            if 'type' in self.dico[name]:
                if self.dico[name]['type'] in ('gas', 'vapor'):
                    self.some_var_gases = True
                    self.var_gas_idx.append(ii)
                    self.var_gas_names.append(name)
                    self.gas_molar_masses.append(xk.Molar_mass().fetch(name))
        self.var_gas_idx = np.array(self.var_gas_idx)
        self.var_gas_names = np.array(self.var_gas_names)
        self.gas_molar_masses = np.array(self.gas_molar_masses)
        self.M_bg = M_bg
        self.update_gas_composition()
        self.Kzz = np.ones(self.Nlay)*Kzz
        self.Dmol = Dmol

    def update_gas_composition(self, update_vmr=True):
        """Performs mass to volume mixing ratio conversion,
        computes the new molar mass of the total gas
        and transmits new composition to the radiative model. 

        !!!Small discrepancies with case without tracers!!!
        Maybe a problem with the background gas.

        Parameters
        ----------
            atm: exo_k.Atm object
                If atm is provided, its composition will be
                updated with the new composition.
        """
        qvar = np.zeros(self.Nlay)
        ovMvar = np.zeros(self.Nlay)
        for ii, idx in enumerate(self.var_gas_idx):
            qvar += self.qarray[idx]
            ovMvar += self.qarray[idx]/self.gas_molar_masses[ii]
        ovMbg=1./self.M_bg
        self.Mgas = 1./((1.-qvar)*ovMbg + ovMvar)
        if update_vmr:
            var_vmr_tot = 0.
            for ii, idx in enumerate(self.var_gas_idx):
                vmr_tmp = np.core.umath.maximum(
                    self.Mgas * self.qarray[idx] / self.gas_molar_masses[ii], 0.)
                    #conversion to vmr
                self.gas_vmr[self.var_gas_names[ii]] = vmr_tmp
                var_vmr_tot += vmr_tmp
                #print(vmr_tmp, var_vmr_tot, isinstance(vmr_tmp, (np.ndarray)))
            var_vmr_tot = 1. - var_vmr_tot    
            #print(var_vmr_tot)
            for mol, vmr in self.bg_vmr.items():
                #print(mol, vmr, var_vmr_tot, isinstance(var_vmr_tot, (np.ndarray)))
                self.gas_vmr[mol] = vmr * var_vmr_tot

    def turbulent_diffusion(self, timestep, Htot, atm, cp):
        """Mixes tracers following a diffusion equation
        with a constant Kzz parameter (self.Kzz in m^2/s).

        Parameters
        ----------
            timestep: float
                physical timestep of the current step (in s/cp).
                (needs to be converted before it is sent to `turbulent diffusion`)
            Htot: array
                Total heating rate (in W/kg) of all physical processes
                already computed
            atm: :class:`Atm` object
                The Atm object used in the radiative transfer which
                contains many state variables. 
        """
        new_t = atm.tlay + timestep * Htot
        if self.settings['moist_inhibition']:
            Mgas_tmp = self.Mgas
        else:
            Mgas_tmp = np.mean(self.Mgas)
            Mgas_tmp = np.full_like(self.Mgas, Mgas_tmp)
        qarray = turbulent_diffusion(timestep*cp, self.Nlay,
                    atm.play, atm.plev,
                    atm.dmass, new_t/Mgas_tmp,
                    atm.grav, self.Kzz, self.qarray)
        #self.dm_trac = (qarray - self.qarray) * atm.dmass / (timestep*cp)
        self.qarray = qarray


    def dry_convective_adjustment(self, timestep, Htot, atm, verbose=False):
        """Computes convective adjustement. 

        Parameters
        ----------
            timestep: float
                physical timestep of the current step (in s/cp).
                (needs to be converted before it is sent to `turbulent diffusion)
            Htot: array
                Total heating rate (in W/kg) of all physical processes
                already computed
            atm: :class:`Atm` object
                The Atm object used in the radiative transfer which
                contains many state variables. 
        """
        new_t = atm.tlay + timestep * Htot
        if self.settings['moist_inhibition']:
            Mgas_tmp = self.Mgas
        else:
            Mgas_tmp = np.mean(self.Mgas)
            Mgas_tmp = np.full_like(self.Mgas, Mgas_tmp)
        H_conv, q_array = dry_convective_adjustment(timestep, self.Nlay, new_t,
                    atm.exner, atm.dmass, self.qarray, Mgas_tmp, verbose=verbose)
        if self.settings['convective_transport']:
            self.qarray=q_array
        return H_conv

    def __getitem__(self, key):
        return self.qarray[self.idx[key]]
