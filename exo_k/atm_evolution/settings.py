# -*- coding: utf-8 -*-
"""
A dictionary like class that captures all the global settings used for an atmospheric evolution. 

@author: jeremy leconte
"""

class Settings(object):

    def __init__(self):
        """Initializes all global parameters to default values
        """
        self.parameters={'rayleigh': True,
                        'internal_flux': 0.,
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
                        'dTmax_use_kernel': 10.,
                        'qvap_deep': -1.,
                        'evap_coeff': 1.,
                        'acceleration_mode': 0,
                        'radiative_acceleration_reducer': 1.,
                        'condensation_timestep_reducer': .8,
                        'convective_acceleration_mode': 0,
                        'qcond_surf_layer': 0.1,
                        }

    def set_parameters(self, **kwargs):
        """Sets various global options
        """
        for key, val in kwargs.items():
            if val is not None:
                self.parameters[key]=val
        if 'logplay' in self.keys():
            self['Nlay']=self['logplay'].size
    
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
