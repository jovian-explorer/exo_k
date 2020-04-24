# -*- coding: utf-8 -*-
"""
@author: jeremy leconte

Contains classes for the atmospheric profile and its radiative properties.
That's where forward models are computed. 
"""
import numpy as np
import astropy.units as u
from .chemistry import gas_mix
from .util.cst import N_A,PI,RGP,KBOLTZ,RSOL,RJUP,SIG_SB
from .util.interp import RandOverlap_2_kdata_prof
from .util.radiation import Bnu_integral_num,Bnu,rad_prop_corrk,rad_prop_xsec,Bnu_integral_array
from .rayleigh import Rayleigh


class Atm_profile(object):
    """A class defining an atmospheric PT profile with some global data
    (gravity, Molar mass, etc.)
    """
    
    def __init__(self,composition={},psurf=None,ptop=None,
            Tsurf=None,Tstrat=None,grav=None,Rp=None,Mgas=None,rcp=0.28,Nlev=20):
        """Initializes atmospheric profiles
        Parameters:
        Option:
        """
        self.gas_mix=gas_mix(composition)
        self.rcp=rcp
        self.psurf=psurf
        self.ptop=ptop
        self.Nlev=Nlev
        self.Nlay=Nlev-1
        self.logplev=np.linspace(np.log10(ptop),np.log10(psurf),num=self.Nlev)
        self.plev=10**self.logplev
        self.logplay=(self.logplev[:-1]+self.logplev[1:])*0.5
        self.play=10**self.logplay
        self.set_Mgas(Mgas)
        self.set_grav(grav)
        self.set_Rp(Rp)        
        self.set_adiab_profile(Tsurf=Tsurf,Tstrat=Tstrat,rcp=rcp)
        self.tlay=(self.tlev[:-1]+self.tlev[1:])*0.5

    def set_logPT_profile(self,log_plev,tlev):
        """Set the logP-T profile of the atmosphere with a new one
        Parameters:
            log_plev: numpy array
                Log pressure (in Pa) at the level surfaces
            tlev: numpy array (same size)
                temperature at the level surfaces.
        """
        self.logplev=log_plev
        self.plev=10**self.logplev
        self.psurf=self.plev[-1]
        self.tlev=tlev
        self.logplay=(self.logplev[:-1]+self.logplev[1:])*0.5
        self.play=10**self.logplay
        self.tlay=(self.tlev[:-1]+self.tlev[1:])*0.5
        self.dmass=(self.plev[1:]-self.plev[:-1])/self.grav
        self.Nlev=log_plev.size
        self.Nlay=self.Nlev-1

    def set_adiab_profile(self,Tsurf=None,Tstrat=None,rcp=0.28):
        """Initializes atmospheric the logP-T profile with an adiabat with index R/cp=rcp
        Parameters:
            Tsurf: float
                Surface temperature
            Tstrat: float
                Temperature of the stratosphere
            rcp: float
                R/c_p of the atmosphere
        """
        self.tlev=Tsurf*(self.plev/self.psurf)**rcp
        self.tlev=np.where(self.tlev<Tstrat,Tstrat,self.tlev)

    def set_grav(self,grav=None):
        """Sets the surface gravity of the planet
        Parameters:
            grav: float
                surface gravity (m/s^2)
        """
        if grav is None: raise RuntimeError('A planet needs a gravity!')
        self.grav=grav
        self.dmass=(self.plev[1:]-self.plev[:-1])/self.grav
    
    def set_gas(self,composition_dict):
        """Sets the composition of the atmosphere
        Parameters:
            composition_dict: dictionary
                Keys are molecule names, and values are volume mixing ratios.
                A 'background' value means that the gas will be used to fill up to vmr=1
                If they do not add up to 1 and there is no background gas_mix,
                the rest of the gas_mix is considered transparent.
        """
        self.gas_mix=gas_mix(composition_dict)
        self.set_Mgas()

    def set_Mgas(self,Mgas=None):
        """Sets the mean molecular weight of the atmosphere
        Parameters:
            Mgas: float
                mean molecular weight (kg/mol)
        """
        if Mgas is not None:
            self.Mgas=Mgas
        else:
            self.Mgas=self.gas_mix.molar_mass()

    def set_rcp(self,rcp):
        """Sets the adiabatic index of the atmosphere
        Parameters:
            rcp: float
                R/c_p
        """
        self.rcp=rcp

    def set_Rp(self,Rp):
        """Sets the radius of the planet
        Parameters:
            Rp: float
                radius of the planet (m)
        """
        if Rp is None:
            self.Rp = None
            return
        if isinstance(Rp,u.quantity.Quantity):
            self.Rp=Rp.to(u.m).value
        else:
            self.Rp=Rp*RJUP

    def set_Rstar(self,Rstar):
        """Sets the radius of the star
        Parameters:
            Rstar: float
                radius of the star (m)
        """
        if Rstar is None:
            self.Rstar = None
            return
        if isinstance(Rstar,u.quantity.Quantity):
            self.Rstar=Rstar.to(u.m).value
        else:
            self.Rstar=Rstar*RSOL


    def compute_density(self):
        """Computes the number density (m^-3) profile of the atmosphere
        """
        self.density=self.play/(KBOLTZ*self.tlay)

    def compute_altitudes(self):
        """Compute altitudes of the level surfaces (zlev) and mid layers (zlay).
        """
        H=RGP*self.tlay/(self.grav*self.Mgas)*np.log(10)
        #print(RGP*self.tlay/(self.grav*self.Mgas))
        self.dz=H*np.diff(self.logplev)
        self.zlay=np.cumsum(self.dz[::-1])
        self.zlev=np.concatenate(([0.],self.zlay))[::-1]
        self.zlay-=0.5*self.dz[::-1]
        self.zlay=self.zlay[::-1]
        
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
        self.tangent_path=[]
        for ilay in range(self.Nlay):
            z0square=(self.Rp+self.zlay[ilay])**2
            dl=np.sqrt((self.Rp+self.zlev[:ilay+1])**2-z0square)
            dl[:-1]-=dl[1:]
            self.tangent_path.append(2.*dl)

class RadAtm(Atm_profile):
    """Class based on Atm_profile that contains the link to the radiative data.
    """

    def __init__(self,kdatabase=None,CIAdatabase=None,wl_range=None,**kwargs):
        """Initialization method that calls Atm_Profile().__init__() and links
        to Kdatabase and other radiative data. 
        """
        super().__init__(**kwargs)
        self.kdatabase=kdatabase
        self.CIAdatabase=CIAdatabase
        if wl_range is not None:
            self.set_wl_range(wl_range)
        else:
            self.wn_range=None
        if self.kdatabase is None:
            self.Ng=None
        else:
            self.Ng=self.kdatabase.Ng

    def set_wl_range(self,wl_range):
        """Sets the wavelength range in which computations will be done.
        Parameters:
            wl_range: Array of size 2
                Minimum and maximum wavelength
        """
        self.set_wn_range(np.sort(10000./np.array(wl_range)))

    def set_wn_range(self,wn_range):
        """Sets the wavenumber range in which computations will be done.
        Parameters:
            wn_range: Array of size 2
                Minimum and maximum wavenumber
        """
        self.wn_range=np.sort(np.array(wn_range))
        self.iw_min,self.iw_max=np.searchsorted(self.kdatabase.wnedges,wn_range,side='left')
        self.iw_max-=1

    def compute_wn_range_indices(self):
        """Compute the starting and ending indices to be used for current wn_range
        """
        if self.wn_range is None:
            self.iw_min=0
            self.iw_max=self.kdatabase.Nw
        else:
            self.iw_min,self.iw_max=  \
                np.searchsorted(self.kdatabase.wnedges,self.wn_range,side='left')
            self.iw_max-=1

    @property
    def wls(self):
        """Returns the wavelength array for the bin centers
        """
        return 10000./self.wns

    @property
    def wledges(self):
        """Returns the wavelength array for the bin edges.
        """
        return 10000./self.wnedges

    def set_database(self,kdatabase):
        """Change the radiative database attached to the current instance of AtmRad
        Parameters:
            kdatabase: Kdatabase object
                New Kdatabase to use.
        """
        self.kdatabase=kdatabase
        self.Ng=self.kdatabase.Ng

    def set_CIAdatabase(self,CIAdatabase):
        """Change the CIA database attached to the current instance of AtmRad
        Parameters:
            CIAdatabase: CIAdatabase object
                New CIAdatabase to use.
        """
        self.CIAdatabase=CIAdatabase

    def opacity(self, wl_range=None, rayleigh=True, write=0, random_overlap=False, **kwargs):
        """Computes the opacity of each of the layers for the composition given
        for every wavelength (and possibly g point).
        For the moment the kcoeff are added to each other (maximum recovery assumption).
        Option:
            wl_range: array of two values
                Wavelength range to cover
        Output:
            self.kdata: array
                opacity of each layer for each wavenumber (and potentially g point).
        """
        if self.kdatabase is None: raise RuntimeError("""kdatabase not provided. 
        Use the kdatabase keyword during initialization or use the set_database method.""")
        if not self.kdatabase.consolidated_wn_grid: raise RuntimeError("""
           All tables in the database should have the same wavenumber grid to proceed.
           You should probably use bin_down().""")
        kdatabase=self.kdatabase
        first_mol=True
        molecs=self.gas_mix.keys()
        mol_to_be_done=list(set(molecs).intersection(kdatabase.molecules))
        if all(elem in kdatabase.molecules for elem in molecs):
            if write>3 : print("""I have all the molecules present in the atmosphere
              in ktables provided:""")
        else:
            if write>3 : print("""Some missing molecules in my database,
             I ll compute opacites with the available ones:""")
        if write>3 : print(mol_to_be_done)

        if wl_range is not None: self.set_wl_range(wl_range)
        self.compute_wn_range_indices()
        self.wnedges=kdatabase.wnedges[self.iw_min:self.iw_max+1]
        self.wns=kdatabase.wns[self.iw_min:self.iw_max]
        self.Nw=self.wns.size
        for mol in mol_to_be_done:
            tmp_kdata=kdatabase.ktables[mol].interpolate_kdata( \
                logp_array=self.logplay,t_array=self.tlay,wngrid_limit=self.wn_range,**kwargs)
            if first_mol:
                self.kdata=self.gas_mix[mol]*tmp_kdata
                first_mol=False
            else:
                if random_overlap and (self.Ng is not None):
                    self.kdata=RandOverlap_2_kdata_prof(self.Nlay,self.Nw,self.Ng, \
                        self.kdata,self.gas_mix[mol]*tmp_kdata,self.kdatabase.weights,
                        self.kdatabase.ggrid)
                else:
                    self.kdata+=self.gas_mix[mol]*tmp_kdata
        if rayleigh or (self.CIAdatabase is not None):
            cont_sig=np.zeros((self.Nlay,self.Nw))
            if self.CIAdatabase is not None:
                cont_sig+=self.CIAdatabase.cia_cross_section( \
                    self.logplay,self.tlay,self.gas_mix,wngrid_limit=self.wn_range)
            if rayleigh:
                cont_sig+=Rayleigh().sigma(self.wns,self.gas_mix)
            if self.Ng is None:
                self.kdata+=cont_sig
            else:
                self.kdata+=cont_sig[:,:,None]

    def emission_spectrum(self, integral=True, mu0=0.5, **kwargs):
        """Computes the emission flux at the top of the atmosphere (in W/m**2/cm**-1)
        Parameters:
        Options:
            integral: Boolean
                If true, the black body is integrated within each wavenumber bin.
                If not, only the central value is used.
                    False is faster and should be ok for small bins,
                    but True is the correct version. 
        Output:
            A Spectrum object with the Spectral flux at the top of the atmosphere (in W/m**2/cm**-1)
        """
        from .util.spectrum import Spectrum
        #if self.dtau is None:
        #    self.rad_prop(**kwargs)
        self.opacity(rayleigh=False, **kwargs)
        NaOvMuMg=N_A/(mu0*self.Mgas)
        if integral:
            #JL2020 What is below is slightly faster for very large resolutions
            #dw=np.diff(self.wnedges)
            #piBatm=np.empty((self.Nlev,self.Nw))
            #for ii in range(self.Nlev):
            #    piBatm[ii]=PI*Bnu_integral_num(self.wnedges,self.tlev[ii])/dw
            #JL2020 What is below is much faster for moderate to low resolutions
            piBatm=PI*Bnu_integral_array(self.wnedges,self.tlev,self.Nw,self.Nlev) \
                /np.diff(self.wnedges)
        else:
            piBatm=PI*Bnu(self.wns[None,:],self.tlev[:,None])
        self.SurfFlux=piBatm[-1]
        if self.Ng is None:
            self.tau,dtau=rad_prop_xsec(self.dmass,self.kdata,NaOvMuMg)
        else:
            self.tau,dtau=rad_prop_corrk(self.dmass,self.kdata,NaOvMuMg)
            weights=self.kdatabase.weights
        expdtau=np.exp(-dtau)
        expdtauminone=np.where(dtau<1.e-14,-dtau,expdtau-1.)
        # careful: due to numerical limitations, 
        # the limited development of Exp(-dtau)-1 needs to be used for small values of dtau
        exptau=np.exp(-self.tau)
        if self.Ng is None:
            timesBatmTop=(-expdtau-expdtauminone/dtau)*exptau[:-1]
            timesBatmBottom=(1.+expdtauminone/dtau)*exptau[:-1]
            timesBatmBottom[-1]=timesBatmBottom[-1]+exptau[-1]
        else:
            timesBatmTop=np.sum((-expdtau-expdtauminone/dtau)*exptau[:-1]*weights,axis=2)
            timesBatmBottom=np.sum((1.+expdtauminone/dtau)*exptau[:-1]*weights,axis=2)
            timesBatmBottom[-1]=timesBatmBottom[-1]+np.sum(exptau[-1]*weights,axis=1)
        IpTop=np.sum(piBatm[:-1]*timesBatmTop+piBatm[1:]*timesBatmBottom,axis=0)

        return Spectrum(IpTop,self.wns,self.wnedges)


    def exp_minus_tau(self):
        """Sums Exp(-tau) over gauss points
        """
        weights=self.kdatabase.weights
        return np.sum(np.exp(-self.tau[1:])*weights,axis=2)


    def surf_bb(self, integral=True):
        """Computes the surface black body flux (in W/m**2/cm**-1)
        Parameters:
        Options:
            integral: Boolean
                If true, the black body is integrated within each wavenumber bin.
                If not, only the central value is used.
                    False is faster and should be ok for small bins,
                    but True is the correct version. 
        Output:
            surf_bb: Spectrum object
                Spectral flux at the surface (in W/m**2/cm**-1)
        """
        from .util.spectrum import Spectrum
        if integral:
            piBatm=PI*Bnu_integral_num(self.wnedges,self.tlev[-1])/np.diff(self.wnedges)
        else:
            piBatm=PI*Bnu(self.wns[:],self.tlev[-1])
        return Spectrum(piBatm,self.wns,self.wnedges)

    def top_bb(self, integral=True):
        """Computes the top of atmosphere black body flux (in W/m**2/cm**-1)
        Parameters:
        Options:
            integral: Boolean
                If true, the black body is integrated within each wavenumber bin.
                If not, only the central value is used.
                    False is faster and should be ok for small bins,
                    but True is the correct version. 
        Output:
            top_bb: Spectrum object
                Spectral flux of a bb at the temperature at the top of atmosphere (in W/m**2/cm**-1)
        """
        from .util.spectrum import Spectrum
        if integral:
            piBatm=PI*Bnu_integral_num(self.wnedges,self.tlev[0])/np.diff(self.wnedges)
        else:
            piBatm=PI*Bnu(self.wns[:],self.tlev[0])
        return Spectrum(piBatm,self.wns,self.wnedges)

    def transmittance_profile(self,**kwargs):
        """Computes the transmittance profile of an atmosphere,
        i.e. Exp(-tau) for each layer of the model.
        Real work done in the numbafied function path_integral_corrk/xsec
        depending on the type of data.
        """
        from .util.radiation import path_integral_corrk,path_integral_xsec
        self.opacity(**kwargs)
        self.compute_altitudes()
        self.compute_tangent_path()
        self.compute_density()
        if self.Ng is not None:
            weights=self.kdatabase.weights
            transmittance=path_integral_corrk( \
                self.Nlay,self.Nw,self.Ng,self.tangent_path,self.density,self.kdata,weights)
        else:
            transmittance=path_integral_xsec( \
                self.Nlay,self.Nw,self.tangent_path,self.density,self.kdata)
        return transmittance

    def transmission_spectrum(self,normalized=False,Rstar=None,**kwargs):
        """Computes the transmission spectrum of the atmosphere normalized
        to the radius of the planet:
        delta_sigma=(R_planet+h_sigma)^2/R_planet^2
        """
        from .util.spectrum import Spectrum
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

    def heating_rate(self, Fin=1., Tstar=5570., szangle=60., cp=1000., **kwargs):
        """Computes the heating rate in the atmosphere
        """
        mu0=np.cos(szangle*PI/180.)
        self.opacity(rayleigh=False, **kwargs)
        NaOvMuMg=N_A/(mu0*self.Mgas)
        Fstar=Fin*PI*Bnu_integral_array(self.wnedges,[Tstar],self.Nw,1)/(SIG_SB*Tstar**4)
        print(np.sum(Fstar))
        if self.Ng is None:
            self.tau, _ =rad_prop_xsec(self.dmass,self.kdata,NaOvMuMg)
            #the second returned variable is ignored
        else:
            self.tau, _ =rad_prop_corrk(self.dmass,self.kdata,NaOvMuMg)
            weights=self.kdatabase.weights
        exptau=np.exp(-self.tau)
        if self.Ng is None:
            tmp_heat_rate=Fstar[0,:]*exptau
            tmp_heat_rate[:-1]-=tmp_heat_rate[1:]
            heat_rate=tmp_heat_rate[:-1]
        else:
            tmp_heat_rate=Fstar[0,:,None]*exptau
            tmp_heat_rate[:-1]-=tmp_heat_rate[1:]
            heat_rate=np.sum(tmp_heat_rate[:-1]*weights,axis=2)
        heat_rate=np.sum(heat_rate,axis=1)
        heat_rate=heat_rate/self.dmass/cp
        return heat_rate





