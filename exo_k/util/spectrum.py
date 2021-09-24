# -*- coding: utf-8 -*-
"""
@author: jeremy leconte

A module to handle ouputs rebinning and plotting
"""
import numpy as np
import h5py
from exo_k.util.interp import rebin
from exo_k.util.spectral_object import Spectral_object

class Spectrum(Spectral_object):
    """A class defining a Spectrum object to plot and manipulate.
    """

    def __init__(self, value=None, wns=None, wnedges=None, filename=None,
            from_taurex=False, dataset='native_spectrum', **kwargs):
        """Instanciate with a value, bin centers, and bin edges.
        Can also load a Taurex spectrum if filename is provided.
        """
        self.wn_unit='cm^-1'
        self.value=value
        self.wns=wns
        self.wnedges=wnedges
        if filename is not None:
            if from_taurex:
                self.load_taurex(filename, dataset)
            elif filename.lower().endswith(('.hdf5', '.h5')):
                self.read_hdf5(filename)
            elif filename.lower().endswith(('.dat', '.txt')):
                self.read_ascii(filename, **kwargs)
        if (self.wnedges is None) and (self.wns is not None):
            self.wnedges=np.concatenate(([self.wns[0]],
                (self.wns[:-1]+self.wns[1:])*0.5,[self.wns[-1]]))

    
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
            raise RuntimeError('The two spectra do not have the same spectral sampling.')

    def __sub__(self,other):
        """Defines substraction
        """
        if (isinstance(other,float) or isinstance(other,int)):
            return Spectrum(self.value-other,self.wns,self.wnedges)
        elif (self.wns.size==other.wns.size) and np.array_equal(self.wns,other.wns):
            val=self.value-other.value
            return Spectrum(val,self.wns,self.wnedges)
        else:
            raise RuntimeError('The two spectra do not have the same spectral sampling.')

    def __mul__(self,other):
        """Defines multiplication
        """
        if (isinstance(other,float) or isinstance(other,int)):
            return Spectrum(self.value*other,self.wns,self.wnedges)
        elif (self.wns.size==other.wns.size) and np.array_equal(self.wns,other.wns):
            val=self.value*other.value
            return Spectrum(val,self.wns,self.wnedges)
        else:
            raise RuntimeError('The two spectra do not have the same spectral sampling.')

    def __truediv__(self,other):
        """Defines division
        """
        if (isinstance(other,float) or isinstance(other,int)):
            return Spectrum(self.value/other,self.wns,self.wnedges)
        elif (self.wns.size==other.wns.size) and np.array_equal(self.wns,other.wns):
            val=self.value/other.value
            return Spectrum(val,self.wns,self.wnedges)
        else:
            raise RuntimeError('The two spectra do not have the same spectral sampling.')

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

    def read_hdf5(self, filename=None):
        """Reads data in a hdf5 file

        Parameters
        ----------
            filename: str
                Name of the file to be created and saved
        """
        if (filename is None or not filename.lower().endswith(('.hdf5', '.h5'))):
            raise RuntimeError("You should provide an input hdf5 file")
        f = h5py.File(filename, 'r')
        self.wns=f['bin_centers'][...]
        self.wnedges=f['bin_edges'][...]
        if 'units' in f['bin_edges'].attrs:
            self.wn_unit=f['bin_edges'].attrs['units']
        else:
            if 'units' in f['bin_centers'].attrs:
                self.wn_unit=f['bin_centers'].attrs['units']
        self.value=f['spectrum'][...]
        f.close()  

    def write_hdf5(self, filename):
        """Saves data in a hdf5 format

        Parameters
        ----------
            filename: str
                Name of the file to be created and saved
        """
        fullfilename=filename
        if not filename.lower().endswith(('.hdf5', '.h5')):
            fullfilename=filename+'.h5'
        compression="gzip"
        f = h5py.File(fullfilename, 'w')
        f.create_dataset("spectrum", data=self.value, compression=compression)
        f.create_dataset("bin_edges", data=self.wnedges, compression=compression)
        f["bin_edges"].attrs["units"] = 'cm^-1'
        f.create_dataset("bin_centers", data=self.wns, compression=compression)
        f["bin_centers"].attrs["units"] = 'cm^-1'
        f.close()

    def read_ascii(self, filename, spec_axis='wns', skip_header=None):
        """Saves data in a ascii format

        Parameters
        ----------
            filename: str
                Name of the file to be read
            spec_axis: str
                Whether the secptral axis in the file is
                wavenumber in cm^-1 ('wns') or wavelength in microns ('wls')
            skip_header: int
                Number of lines to skip
        """
        raw=np.genfromtxt(filename, skip_header=skip_header,
            usecols=(0,1), names=('spec_axis','value'))
        if spec_axis=='wls':
            raw['spec_axis']=1.e4/raw['spec_axis'][::-1] #conversion to wavenumber in cm-1
            raw['value']=raw['value'][::-1]
        self.value=raw['value']
        self.wns=raw['spec_axis']

    def write_ascii(self, filename, fmt='%.18e', spec_axis='wns', header=None):
        """Saves data in a ascii format

        Parameters
        ----------
            filename: str
                Name of the file to be created and saved
        """
        fullfilename=filename
        if not filename.lower().endswith(('.dat', '.txt')):
            fullfilename=filename+'.dat'
        head=header
        if spec_axis=='wns':
            if head is None: head='wavenumber(cm^-1)     spectrum'
            np.savetxt(fullfilename, np.array([self.wns,self.value]).transpose(),
                fmt=fmt, header=head)
        else:    
            if head is None: head='wavelength(micron)    spectrum'
            np.savetxt(fullfilename, np.array([self.wls[::-1],self.value[::-1]]).transpose(),
                fmt=fmt, header=head)

    def load_taurex(self, filename,dataset='native_spectrum'):
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

