# -*- coding: utf-8 -*-
"""
A dictionary like class that captures all the global settings used for an atmospheric evolution. 

@author: jeremy leconte
"""

class Settings(object):

    def __init__(self):
        """Initializes all global parameters to default values

        *Note for developers*: The initial self.parameters dictionary (the one below)
        should not contain any item whose key is a keyword argument for Atm.__init__().
        It should also contain every other possible option that does not require
        re-initializing the radiative transfer model. This is necessary for
        Atm_evolution.set_options() to know when to reset the radiative model. 
        """
        self.parameters={'rayleigh': True,
                        'convection': False,
                        'convective_transport': True,
                        'diffusion': False,
                        'molecular_diffusion': False,
                        'condensation': False,
                        'rain': False,
                        'latent_heating': True,
                        'moist_convection': False,
                        'moist_inhibition': False,
                        'surface_reservoir': False,
                        'mass_redistribution': False,
                        'compute_mass_fluxes': True,
                        'dTmax_use_kernel': 10.,
                        'qvap_deep': -1.,
                        'evap_coeff': 1.,
                        'acceleration_mode': 0,
                        'radiative_acceleration_reducer': 1.,
                        'condensation_timestep_reducer': .8,
                        'convective_acceleration_mode': 0,
                        'qcond_surf_layer': 0.1,
                        }
            
        self._forbidden_changes = ['logplay', 'play']

        self._non_radiative_parameters = set(self.parameters.keys())

    def set_parameters(self, **kwargs):
        """Sets various global options
        """
        for key, val in kwargs.items():
            if val is not None:
                self.parameters[key]=val
        if 'logplay' in self.keys():
            self['Nlay']=self['logplay'].size

    def use_or_set(self, key, value):
        """Returns the value stored in the parameters if available.
        Stores the given value if not. 
        """
        if key in self.parameters.keys():
            return self.parameters[key]
        else:
            self.parameters[key] = value
            return value

    def __getitem__(self, param):
        return self.parameters[param]

    def __setitem__(self, key, item):
        self.parameters[key] = item

    def pop(self, key, default):
        return self.parameters.pop(key, default)

    def get(self, key, default):
        return self.parameters.get(key, default)

    def keys(self):
        return self.parameters.keys()

    def items(self):
        return self.parameters.items()

    def values(self):
        return self.parameters.values()

    def __repr__(self):
        """Method to output parameters
        """
        return self.parameters.__repr__()
