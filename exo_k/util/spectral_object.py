"""
@author: jeremy leconte

A class with some basic functions for all objects with a spectral dimension
"""

class Spectral_object(object):
    """A class with some basic functions for all objects with a spectral dimension
    """

    @property
    def wls(self):
        """Returns the wavelength array for the bin centers
        """
        return 10000./self.wns
