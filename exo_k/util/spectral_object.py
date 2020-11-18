"""
@author: jeremy leconte

A class with some basic functions for all objects with a spectral dimension
"""
import numpy as np

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

    @property
    def wnrange(self):
        """Returns the limits of the wavenumber range.

        First tries with wnedges (for Ktables) and then wns (Xtables).
        """
        if self.wnedges is not None:
            return self.wnedges[[0,-1]]
        elif self.wns is not None:
            return self.wns[[0,-1]]
        else:
            return None

    @property
    def wlrange(self):
        """Returns the limits of the wavelength range
        """
        if (self.wnedges is None) and (self.wns is None):
            return None
        else:
            return np.sort(10000./self.wnrange)
