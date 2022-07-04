# -*- coding: utf-8 -*-
"""
@author: jeremy leconte

This module contain classes to handle atmospheric profiles.
Radiative properties are handled in atm.py which contains a daughter class.

The nomenclature for layers, levels, etc., can be found in atm.py.
"""
import numpy as np
import warnings
import astropy.units as u
from numba.typed import List
from .gas_mix import Gas_mix
from .aerosols import Aerosols
from .util.cst import N_A, PI, RGP, KBOLTZ


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
        self.grav = None
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
        self.compute_layer_masses()

    def update_pressure_profile(self, play = None, plev = None):
        """Updates pressure levels without changing temperatures.
        
        To be used Atm_evolution class.
        """
        play = np.array(play, dtype=float)
        if play.size != self.Nlay:
            raise RuntimeError("You cannot change the number of layers in update_pressure_profile")
        self.play = play
        self.plev = np.array(plev, dtype=float)
        self.logplay = np.log10(self.play)
        self.logplev = np.log10(self.plev)
        self.logp_opac = self.logplev[1:-1]
        self.psurf = self.plev[-1]
        self.dp_lay = np.diff(self.plev) ### probably redundant with dmass
        self.exner = (self.play/self.psurf)**self.rcp
        self.compute_layer_masses()

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
        self.compute_layer_masses()

    def compute_layer_masses(self):
        """compute_layer_masses
        """
        if self.grav is not None:
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
        vmr_midlevel_dict=dict()
        for mol, vmr in composition_dict.items():
            if isinstance(vmr,(np.ndarray, list)):
                tmp_vmr=np.array(vmr)
                # geometrical average:
                with np.errstate(invalid='raise'):
                    try:
                        vmr_midlevel_dict[mol] = np.sqrt(tmp_vmr[1:]*tmp_vmr[:-1])
                    except FloatingPointError:
                        print('inset_gas, mol, tmp_vmr=',mol,tmp_vmr)
                        vmr_midlevel_dict[mol] = tmp_vmr[1:]
            else:
                vmr_midlevel_dict[mol] = vmr
        if self.gas_mix is None:
            self.gas_mix = Gas_mix(vmr_midlevel_dict)
        else:
            self.gas_mix.set_composition(vmr_midlevel_dict)
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
            self.Mgas_rad=self.Mgas_rad*np.ones(self.Nlay-1, dtype=float)

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
        Mgas = np.empty(self.Nlay, dtype=float)
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

    def plot_T_profile(self, ax, invert_p = True, use_altitudes = False,
            xscale=None, yscale=None, **kwarg):
        """Plot the T P profile
        
        Parameters
        ----------
            ax : :class:`pyplot.Axes`
                A pyplot axes instance where to put the plot.
            x/yscale: str, optional
                If 'log' log axes are used.
        """
        if use_altitudes:
            self.compute_altitudes()
            ax.plot(self.tlay,self.zlay,**kwarg)
            ax.set_ylabel('Altitude (m)')
        else:
            ax.plot(self.tlay,self.play,**kwarg)
            if invert_p: ax.invert_yaxis()
            ax.set_ylabel('Pressure (Pa)')
        ax.set_xlabel('Temperature (K)')
        if xscale is not None: ax.set_xscale(xscale)
        if yscale is not None: ax.set_yscale(yscale)


    def write_soundings(self, dirname='.', fmt='%.10e', cp=None, qvap=None, p_dry=None):
        """Writes sounding files that can be used to initiate the mesoscale model
        """
        self.compute_altitudes()
        mgas = np.zeros(self.Nlay)
        mgas[:-1] = self.Mgas_rad
        mgas[-1] = self.Mgas_rad[-1]
        r = RGP/mgas
        rho_lay = self.play/(r*self.tlay)
        if cp is not None:
            cp_array = cp * np.ones(self.Nlay)
        else:
            cp_array =  r / self.rcp
        zeros=np.zeros(self.Nlay)
        teta = self.tlay / self.exner
        if qvap is not None:
            q=qvap
        else:
            q=zeros
        filename=dirname+'/input_sounding'
        np.savetxt(filename, np.transpose([self.zlay[::-1], teta[::-1], q[::-1], zeros ,zeros]),
            fmt=fmt, header=str(self.psurf/100.)+'  '+str(self.tlay[-1])+'    0.0', comments=' ')
        # the last dummy column are q_vap (not used as of 2021, u, and v). We put it to zero.

        filename=dirname+'/input_therm'
        if p_dry is None:
            np.savetxt(filename,
                np.transpose([r[::-1], cp_array[::-1], self.play[::-1], rho_lay[::-1], self.tlay[::-1]]), fmt=fmt)
        else:
            np.savetxt(filename,
                np.transpose([r[::-1], cp_array[::-1], p_dry[::-1], rho_lay[::-1], self.tlay[::-1]]), fmt=fmt)
