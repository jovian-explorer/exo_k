# -*- coding: utf-8 -*-
"""
@author: jeremy leconte
"""

class Settings(object):

    def __init__(self):
        """Initializes all gloabal parameters to default values
        """
        self.parameters={'rayleigh': True,
                        'internal_flux': 0.,
                        'convection': False,
                        'radiative_acceleration_reducer': 1.,
                        'convective_transport': True,
                        'diffusion': False,
                        'molecular_diffusion': False,
                        'condensation': False,
                        'rain': False,
                        'latent_heating': True,
                        'moist_convection': False,
                        'moist_inhibition': False,
                        'dTmax_use_kernel': 10.,
                        'qvap_deep': -1.,
                        'evap_coeff': 1.,
                        'acceleration_mode': 0,
                        'condensation_timestep_reducer': .8,
                        }

    def set_parameters(self, **kwargs):
        """Sets various global options
        """
        for key, val in kwargs.items():
            if val is not None:
                self.parameters[key]=val
        if 'logplay' in kwargs.keys():
            self['Nlay']=self['logplay'].size
    
    def __getitem__(self, param):
        return self.parameters[param]

    def __setitem__(self, key, item):
        self.parameters[key] = item

    def pop(self, key, default):
        return self.parameters.pop(key, default)

    def get(self, key, default):
        return self.parameters.get(key, default)

    def __repr__(self):
        """Method to output parameters
        """
        return self.parameters.__repr__()
