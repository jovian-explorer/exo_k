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
        if self.wns is not None:
            return 10000./self.wns
        else:
            return None

    @property
    def wledges(self):
        """Returns the wavelength array for the bin edges
        """
        if self.wnedges is not None:
            return 10000./self.wnedges
        else:
            return None

