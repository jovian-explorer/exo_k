# -*- coding: utf-8 -*-
"""
@author: jeremy leconte

A module to handle ouputs rebinning and plotting
"""
import numpy as np
import h5py
from exo_k.util.interp import rebin

class Spectrum(object):
    """A class defining a Spectrum object to plot and manipulate.
    """

    def __init__(self, value=None, wns=None, wnedges=None, filename=None,
            dataset='native_spectrum'):
        """Instanciate with a value, bin centers, and bin edges.
        Can also load a Taurex spectrum if filename is provided.
        """
        if filename is None:
            self.value=value
            self.wns=wns
            self.wnedges=wnedges
        else:
            self.load_taurex(filename,dataset)

    @property
    def wls(self):
        """Returns the wavelength array for the bin centers
        """
        return 10000./self.wns

    @property
    def wledges(self):
        """Returns the wavelength array for the bin edges
        """
        return 10000./self.wnedges
    
    def copy(self):
        """Deep copy of the spectrum.
        """
        return Spectrum(self.value,self.wns,self.wnedges)

    def plot_spectrum(self, ax, per_wavenumber=True, x_axis='wls',
            xscale=None, yscale=None, **kwarg):
        """Plot the spectrum
        
        Parameters
        ----------
            ax : :class:`pyplot.Axes`
                A pyplot axes instance where to put the plot.
            per_wavenumber: bool, optional
                Defines the units of spectral flux density.
                False converts to per wavelength units.
            x_axis: str, optional
                If 'wls', x axis is wavelength. Wavenumber otherwise.
            x/yscale: str, optional
                If 'log' log axes are used.
        """
        if per_wavenumber:
            toplot=self.value
        else:
            toplot=self.value/self.wls**2*1.e4
        if x_axis == 'wls':
            ax.plot(self.wls,toplot,**kwarg)
            ax.set_xlabel('Wavelength (micron)')
        else:
            ax.plot(self.wns,toplot,**kwarg)
            ax.set_xlabel('Wavenumber (cm$^{-1}$)')
        ax.set_ylabel('Flux')
        ax.grid(True)
        if xscale is not None: ax.set_xscale(xscale)
        if yscale is not None: ax.set_yscale(yscale)

    def bin_down(self, wnedges):
        """Bins down the spectrum to a new grid of wnedges by conserving area.
        
        Parameters
        ----------
            wnedges: array
                Wavenumbers of the bin edges to be used
        """
        wnedges=np.array(wnedges)
        self.value=rebin(self.value,self.wnedges,wnedges)
        self.wnedges=wnedges
        self.wns=0.5*(self.wnedges[:-1]+self.wnedges[1:])

    def bin_down_cp(self, wnedges):
        """Returns a new binned down spectrum to a grid of wnedges by conserving area.
        
        Parameters
        ----------
            wnedges: array
                Wavenumbers of the bin edges to be used

        Returns
        -------
            :class:`Spectrum`
                Binned down spectrum
        """
        res=self.copy()
        res.bin_down(wnedges)
        return res

    def __add__(self,other):
        """Defines addition
        """
        if (isinstance(other,float) or isinstance(other,int)):
            return Spectrum(self.value+other,self.wns,self.wnedges)
        elif (self.wns.size==other.wns.size) and np.array_equal(self.wns,other.wns):
            val=self.value+other.value
            return Spectrum(val,self.wns,self.wnedges)
        else:
            raise RuntimeError('The two spectra do not have the same sampling.')

    def __sub__(self,other):
        """Defines substraction
        """
        if (isinstance(other,float) or isinstance(other,int)):
            return Spectrum(self.value-other,self.wns,self.wnedges)
        elif (self.wns.size==other.wns.size) and np.array_equal(self.wns,other.wns):
            val=self.value-other.value
            return Spectrum(val,self.wns,self.wnedges)
        else:
            raise RuntimeError('The two spectra do not have the same sampling.')

    def __mul__(self,other):
        """Defines multiplication
        """
        if (isinstance(other,float) or isinstance(other,int)):
            return Spectrum(self.value*other,self.wns,self.wnedges)
        elif (self.wns.size==other.wns.size) and np.array_equal(self.wns,other.wns):
            val=self.value*other.value
            return Spectrum(val,self.wns,self.wnedges)
        else:
            raise RuntimeError('The two spectra do not have the same sampling.')

    def __truediv__(self,other):
        """Defines division
        """
        if (isinstance(other,float) or isinstance(other,int)):
            return Spectrum(self.value/other,self.wns,self.wnedges)
        elif (self.wns.size==other.wns.size) and np.array_equal(self.wns,other.wns):
            val=self.value/other.value
            return Spectrum(val,self.wns,self.wnedges)
        else:
            raise RuntimeError('The two spectra do not have the same sampling.')

    def std(self):
        """Defines standard deviation
        """
        return self.value.std()

    def abs(self):
        """Defines absolute value
        """
        return Spectrum(np.abs(self.value),self.wns,self.wnedges)

    def log10(self):
        """Defines Log 10
        """
        return Spectrum(np.log10(self.value),self.wns,self.wnedges)

    @property
    def total(self):
        """Defines the weighted sum over the spectrum
        """
        dw=np.diff(self.wnedges)
        return np.dot(self.value,dw)

    def load_taurex(self,filename,dataset='native_spectrum'):
        """Loads a taurex file

        Parameters
        ----------
            filename: str
                Full name (path) of the hdf5 file to load
            dataset: str
                Name of the hdf5 dataset to load
        """
        f = h5py.File(filename, 'r')
        self.wns=f['Output/Spectra/native_wngrid'][...]
        self.value=f['Output/Spectra/'+dataset][...]
        f.close()
        self.wnedges=np.concatenate(([self.wns[0]],(self.wns[:-1]+self.wns[1:])*0.5,[self.wns[-1]]))


    def __repr__(self):
        """Method to output header
        """
        output="""
        value        : {val}
        wl (microns) : {wl}
        """.format(val=self.value,wl=self.wls)
        return output

