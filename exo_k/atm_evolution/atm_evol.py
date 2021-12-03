# -*- coding: utf-8 -*-
"""
@author: jeremy leconte
"""
import numpy as np
import numba
import exo_k as xk
from .settings import Settings
from .convection import dry_convective_adjustment, turbulent_diffusion, molecular_diffusion
from .condensation import Condensing_species, moist_adiabat, compute_condensation,\
                Condensation_Thermodynamical_Parameters 

class Atm_evolution(object):
    """Model of atmospheric evolution.

    Uses exo_k.Atm class to compute radiative transfer
    """

    def __init__(self, bg_vmr={}, verbose=False, **kwargs):
        """Initializes atmospheric profiles.

        Most arguments are passed directly to exo_k.Atm class through **kwargs

        .. warning::
            Layers are counted from the top down (increasing pressure order).
            All methods follow the same convention.
        """
        self.settings = Settings()
        self.settings.set_parameters(**kwargs)
        self.header={'rad':0,'conv':1,'cond':2,'madj':3,'rain':4,'tot':5}

        self.bg_vmr=bg_vmr
        self.bg_gas = xk.Gas_mix(self.bg_vmr)
        self.M_bg = self.bg_gas.molar_mass()
        self.M_bg = self.settings.pop('M_bg', self.M_bg)
        self.cp = self.settings['cp']
        self.rcp = xk.RGP/(self.M_bg*self.cp)
        self.rcp = self.settings.pop('rcp', self.rcp)
        if verbose: print('cp, M_bg, rcp:', self.cp, self.M_bg, self.rcp)


        self.tracers=Tracers(self.settings, bg_vmr=self.bg_vmr,
            M_bg=self.M_bg, **self.settings.parameters)
        self.initialize_condensation(**self.settings.parameters)

        self.setup_radiative_model(rcp = self.rcp, gas_vmr=self.tracers.gas_vmr,
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
        if reset_rad_model: self.setup_radiative_model(**self.settings.parameters)

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

    def compute_condensation(self, timestep, Htot):
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
            mass_cond = np.sum(self.tracers.qarray[idx_cond]*self.atm.dmass)
            dqvap = mass_cond / self.atm.dmass[-1]
            self.tracers.qarray[idx_vap,-1] += dqvap
            #print(dqvap*self.atm.dmass[-1]/(self.atm.grav*timestep))
            if self.settings['latent_heating']:
                H_rain[-1] += - self.condensing_species_params[i_cond].Lvap(new_t[-1]) \
                    * dqvap / (timestep *self.cp)
            self.tracers.qarray[idx_cond] = np.zeros(self.Nlay)
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


    def evolve(self, Niter=1, N_kernel=10000., alpha=1., dT_max = 50., verbose = False, **kwargs):
        r"""The time variable used in the model is tau=t/cp. 
        The equation we are solving in each layer is thus

        .. math::
            c_p \frac{d T}{d t}= \frac{d T}{d tau} = \sum_i H_i

        For a given timestep `dtau`, the physical time elapsed in second can be computed using `dt=dtau*cp`

        To work, the heating rates (H) must be computed in W/kg. 

        THis also means that if one needs the physical rate of change of another quantity (like dq/dt)
        from the delta q over a time step,
        one needs to do `dq/dt = delta q / (timestep * cp)`
        """
        self.H_hist = np.zeros((6, Niter, self.Nlay))
        self.H_ave = np.zeros((6, self.Nlay))

        self.tlay_hist = np.zeros((Niter,self.Nlay))
        self.Fnet_top = np.zeros((Niter))
        self.timestep_hist = np.zeros((Niter))
        tau0 = self.evol_time
        self.N_last_ker = 0
        compute_kernel = False
        dTlay_max = 2. * self.settings['dTmax_use_kernel']
        self.tracers.update_gas_composition(update_vmr=True)
        for ii in range(Niter):
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
            if ii == Niter-1:
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
            if self.settings['radiative_acceleration']:
                self.H_rad =  self.H_rad * self.atm.tau_rads / self.atm.tau_rad
            self.H_tot += self.H_rad
            self.timestep = alpha * self.atm.tau_rad
            #if verbose: print('tau_rad, dt:', self.atm.tau_rad, self.timestep)
            self.evol_time += self.timestep
            if self.settings['diffusion']:
                self.tracers.compute_turbulent_diffusion(self.timestep, self.H_tot, self.atm, self.cp)
                self.tracers.update_gas_composition(update_vmr=False)
            if self.settings['molecular_diffusion']:
                self.H_diff = self.compute_molecular_diffusion(self.timestep,
                    self.H_tot, self.atm, self.cp)
                self.H_tot += self.H_diff
            if self.settings['convection']:
                self.H_conv = self.tracers.dry_convective_adjustment(self.timestep, self.H_tot, self.atm, verbose=verbose)
                self.H_tot += self.H_conv
            else:
                self.H_conv = np.zeros(self.Nlay)
            if self.settings['moist_convection']:
                self.H_madj = self.moist_convective_adjustment(self.timestep, self.H_tot,
                                    moist_inhibition=self.settings['moist_inhibition'], verbose=verbose)                
                self.H_tot += self.H_madj
            else:
                self.H_madj = np.zeros(self.Nlay)
            if self.settings['condensation']:
                self.H_cond = self.compute_condensation(self.timestep, self.H_tot)
                self.H_tot += self.H_cond
            else:
                self.H_cond = np.zeros(self.Nlay)
            if self.settings['rain']:
                self.H_rain = self.rainout(self.timestep, self.H_tot)
                self.H_tot += self.H_rain
            else:
                self.H_rain = np.zeros(self.Nlay)
            dTlay= self.H_tot * self.timestep
            dTlay_max = np.amax(np.abs(dTlay))
            if verbose: print('heat rates (rad, dry conv), dTmax:', 
                np.sum(self.H_rad*self.atm.dmass), np.sum(self.H_conv*self.atm.dmass), dTlay_max)
            #if dTlay_max > dT_max:
            #    print("got a big T step at iteration:", ii)
            #    print(dTlay)
            #    break
            dTlay=np.clip(dTlay,-dT_max,dT_max)
            self.tlay = self.tlay + dTlay
            self.tlay_hist[ii] = self.tlay
            for jj, H in enumerate([self.H_rad, self.H_conv, self.H_cond, self.H_madj,
                  self.H_rain, self.H_tot]):
                self.H_hist[jj,ii] = H
                self.H_ave[jj] += H * self.timestep
            self.Fnet_top[ii] = self.Fnet_rad[0]
            self.timestep_hist[ii] = self.timestep
            self.atm.set_T_profile(self.tlay)
            self.tracers.update_gas_composition(update_vmr=True)
        inv_delta_t = 1./(self.evol_time-tau0)
        self.H_ave *= inv_delta_t
        self.compute_average_fluxes()

    def moist_convective_adjustment2(self, timestep, Htot, moist_inhibition = True,
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
        for i_cond in range(self.Ncond): #careful i_cond is a dumy loop index, idx_cond is position of species i_cond in tracers array.
            idx_vap, idx_cond = self.condensing_pairs_idx[i_cond]
            cond_species = self.condensing_species_params[i_cond]
            #print(new_t, self.Mgas)
            #print(cond_species.Psat(new_t))
            dlnt_dlnp_moist, Lvap, psat, qsat, dqsat_dt, q_crit = \
                cond_species.moist_adiabat(new_t, self.atm.play, self.cp, self.tracers.Mgas)
            #print(dlnt_dlnp_moist, qsat, dqsat_dt)
            H, qarray, new_t = moist_convective_adjustment2(timestep, self.Nlay,
                new_t, self.atm.play, self.atm.dmass, self.cp, self.tracers.qarray, idx_vap, idx_cond,
                qsat, dqsat_dt, Lvap, dlnt_dlnp_moist, q_crit,
                moist_inhibition = moist_inhibition, verbose = verbose)
            #print('t after madj:', new_t)
            H_madj += H
            self.tracers.qarray = qarray
            if verbose: print(qarray[idx_cond])
        return H_madj

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

    def compute_molecular_diffusion(self, timestep, Htot, atm, cp):
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
def moist_convective_adjustment2(timestep, Nlay, tlay, play, dmass, cp, q_array,
        i_vap, i_cond, qsat, dqsat_dt, Lvap, dlnt_dlnp_moist, q_crit, 
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
    new_q = q_array.copy()
    new_t = tlay.copy()
    dp = np.diff(play)
    nabla_interlayer = tlay * dlnt_dlnp_moist /play
    nabla_interlayer = 0.5*(nabla_interlayer[:-1]+nabla_interlayer[1:])
    dTmoist_array = nabla_interlayer * dp
    dT_inter_lay = np.diff(tlay)
    qvap = new_q[i_vap]
    mvap_sursat_array = (qvap-qsat) * dmass
    H_madj=np.zeros(Nlay)
    if moist_inhibition:
        q_crit_criterion = qvap<q_crit # convection possible if True
    else:
        q_crit_criterion = qvap<2. #should always be true
    #print('dT:', dT_inter_lay)
    #print('dTmoist:', dTmoist_array)
    #dT_unstab = np.nonzero(dT_inter_lay>dTmoist_array)[0]
    #saturated = np.nonzero(mvap_sursat_array>0.)[0]
    #conv = np.intersect1d(dT_unstab, saturated, assume_unique = True)# find convective layers
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
    i_conv=0
    i_top=conv[i_conv] #upper unstable layer
    T_top = new_t[i_top]
    mvap_sursat = mvap_sursat_array[i_top]
    C = (cp + Lvap[i_top]*dqsat_dt[i_top])*dmass[i_top]
    B = C * new_t[i_top] + Lvap[i_top] * mvap_sursat_array[i_top]
    dT_moist = 0.
    i_bot=i_top
    while i_bot<Nlay-1: #search for the bottom of the 1st unstable layer from its top
        tmp1 = mvap_sursat + mvap_sursat_array[i_bot+1]
        tmp2 = dT_moist + dTmoist_array[i_bot]
        if tmp1>0. and tmp2<new_t[i_bot+1]-T_top and q_crit_criterion[i_bot+1]:
            i_bot += 1
            dT_moist = tmp2
            mvap_sursat = tmp1
            tmp = (cp + Lvap[i_bot]*dqsat_dt[i_bot])*dmass[i_bot]
            C += tmp
            B += tmp * (new_t[i_bot]-dT_moist) + Lvap[i_bot] * mvap_sursat_array[i_bot] 
            continue
        else:
            break
    if verbose: print('it,ib=', i_top, i_bot)
    if i_top==i_bot: # need at least 2 layers to convect, so exit
        return H_madj, new_q, new_t
    new_Ttop=(B)/C
    dT=new_Ttop-new_t[i_top]
    qvap[i_top] = qsat[i_top] + dqsat_dt[i_top] * dT
    #if verbose: print('top i, dT, qv, qs, dqs, qf', i_top, dT, q_array[i_vap, i_top], qsat[i_top], dqsat_dt[i_top], qvap[i_top])
    new_t[i_top] = new_Ttop
    H_madj[i_top]=dT/timestep
    for ii in range(i_top+1, i_bot+1):
        dT = new_t[ii-1] + dTmoist_array[ii-1] - new_t[ii]
        #print(ii, new_t[ii-1], dTmoist_array[ii-1],  new_t[ii], new_t[ii-1] + dTmoist_array[ii-1] - new_t[ii])
        qvap[ii] = qsat[ii] + dqsat_dt[ii] * dT
        new_t[ii] += dT
        # compute heating and adjust before looking for a new potential unstable layer
        H_madj[ii]=dT/timestep
        #if verbose: print('i, dT, qv, qs, dqs, qf', ii, dT, q_array[i_vap, ii], qsat[ii], dqsat_dt[ii], qvap[ii])
    # put ice
    m_cond = np.sum((q_array[i_vap, i_top:i_bot+1]-qvap[i_top:i_bot+1])*dmass[i_top:i_bot+1])
    if m_cond<0.:
        print('Negative condensates in moist adj, i:', i_top, i_bot)
        print(q_array[i_vap, i_top:i_bot+1], qvap[i_top:i_bot+1], q_array[i_vap, i_top:i_bot+1]-qvap[i_top:i_bot+1])
    m_conv = np.sum(dmass[i_top:i_bot+1])
    new_q[i_cond, i_top:i_bot+1] += m_cond / m_conv
    if verbose: 
        print('m_cond, m_conv, new_q', m_cond, m_conv)
    return H_madj, new_q, new_t

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
        self.Kzz = Kzz
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
        self.Mgas = 1./(qvar*(ovMvar - ovMbg)+ ovMbg)
        if update_vmr:
            for ii, idx in enumerate(self.var_gas_idx):
                self.gas_vmr[self.var_gas_names[ii]] = np.core.umath.maximum(
                    self.Mgas * self.qarray[idx] / self.gas_molar_masses[ii], 0.)
                    #conversion to vmr

    def compute_turbulent_diffusion(self, timestep, Htot, atm, cp):
        """Mixes tracers following a diffusion equation
        with a constant Kzz parameter (self.Kzz in m^2/s).

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
        qarray = turbulent_diffusion(timestep*cp, self.Nlay,
                    atm.play, atm.plev,
                    atm.dmass, new_t/self.Mgas,
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
        H_conv, q_array = dry_convective_adjustment(timestep, self.Nlay, new_t,
                    atm.exner, atm.dmass, self.qarray, self.Mgas, verbose=verbose)
        if self.settings['convective_transport']:
            self.qarray=q_array
        return H_conv

    def __getitem__(self, key):
        return self.qarray[self.idx[key]]
