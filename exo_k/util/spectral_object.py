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
        """Returns the wavelength array for the bin centers (in micron)
        """
        if self.wns is not None:
            return 10000./self.wns
        else:
            return None

    @property
    def wledges(self):
        """Returns the wavelength array for the bin edges (in micron)
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

    def select_spectral_range(self, wn_range=None, wl_range=None):
        """Select spectral range, without restricting the data. Should use either wn_range OR wl_range, not both.

        Parameters
        ----------
            wn_range: array
                Wavenumber range in cm^-1.
            wl_range: array
                Wavelength range in micron.

        Returns
        -------
            tuple:
                iw_min, iw_max the boundary indices of the spectral range
        """
        if (wn_range is None) and (wl_range is None): return None, None
        if wl_range is not None:
            if wn_range is not None:
                raise RuntimeError('Should provide either wn_range or wl_range, not both!')
            _wn_range=np.sort(10000./np.array(wl_range))
        else:
            _wn_range=np.sort(np.array(wn_range))
        iw_min, iw_max=np.searchsorted(self.wns, _wn_range, side='left')
        self.wnedges=self.wnedges[iw_min:iw_max+1]
        self.wns=self.wns[iw_min:iw_max]
        self.Nw=self.wns.size
        return iw_min, iw_max
