# -*- coding: utf-8 -*-
"""
@author: jeremy leconte
"""
import os
import array
import pickle
import h5py
from astropy.io import fits
import numpy as np
import astropy.units as u
from .data_table import Data_table
from .util.interp import rm_molec, gauss_legendre
from .util.filenames import _read_array, _read_exorem_k_array
from .util.cst import nemesis_hitran_id_numbers

class Ktable_io(Data_table):
    """A class to handle the input-output methods of the :class:`~exo_k.ktable.Ktable` class.
    """

    def read_hdf5(self, filename=None, mol=None):
        """Initializes k coeff table and supporting data from an hdf5 file
        (compatible with Exomol format)

        Parameters
        ----------
            file : str
                Name of the input hdf5 file
        """
        if (filename is None or not filename.lower().endswith(('.hdf5', '.h5'))):
            raise RuntimeError("You should provide an input hdf5 file")
        with h5py.File(filename, 'r') as f:
            if 'mol_name' in f:
                self.mol=f['mol_name'][()]
            elif 'mol_name' in f.attrs:
                self.mol=f.attrs['mol_name']
            else:
                if mol is not None:
                    self.mol=mol
                else:
                    self.mol=os.path.basename(filename).split(self._settings._delimiter)[0]
            if isinstance(self.mol, np.ndarray): self.mol=self.mol[0]
            if isinstance(self.mol, bytes): self.mol=self.mol.decode('UTF-8')
            if 'method' in f:
                self.sampling_method=f['method'][()][0]
                if isinstance(self.sampling_method, bytes):
                    self.sampling_method=self.sampling_method.decode('UTF-8')
            if 'DOI' in f:
                self.DOI=f['DOI'][()][0]
                if isinstance(self.DOI, bytes): self.DOI=self.DOI.decode('UTF-8')
            self.wns=f['bin_centers'][...]
            self.wnedges=f['bin_edges'][...]
            if 'units' in f['bin_edges'].attrs:
                self.wn_unit=f['bin_edges'].attrs['units']
            else:
                if 'units' in f['bin_centers'].attrs:
                    self.wn_unit=f['bin_centers'].attrs['units']
            self.kdata=f['kcoeff'][...]
            self.kdata_unit=f['kcoeff'].attrs['units']
            self.tgrid=f['t'][...]
            self.pgrid=f['p'][...]
            self.logpgrid=np.log10(self.pgrid)
            self.p_unit=f['p'].attrs['units']
            if 'weights' in f.keys():
                self.weights=f['weights'][...]
            else:
                raise RuntimeError('No weights keyword. This file is probably a cross section file.')
            self.ggrid=f['samples'][...]
            self.gedges=np.insert(np.cumsum(self.weights),0,0.)
            self.logk=False
        self.Np,self.Nt,self.Nw,self.Ng=self.kdata.shape

    def write_hdf5(self, filename, compression="gzip", compression_level=9,
        kdata_unit=None, p_unit=None, exomol_units=False):
        """Saves data in a hdf5 format

        Parameters
        ----------
            filename: str
                Name of the file to be created and saved
            exomol_units: bool (optional)
                If True, data are converted back to
                cm^2 and bar units before being written.
        """
        dt = h5py.string_dtype(encoding='utf-8')
        fullfilename=filename
        if not filename.lower().endswith(('.hdf5', '.h5')):
            fullfilename=filename+'.h5'
        with h5py.File(fullfilename, 'w') as f:
            f.create_dataset("temperature_grid_type", (1,), data='regular', dtype=dt)
            f.create_dataset("key_iso_ll", (1,), data=self.isotopolog_id)
            f.create_dataset("mol_mass", (1,), data=self.molar_mass*1000.)
            f["mol_mass"].attrs["units"] = 'AMU'
            if exomol_units:
                kdata_unit='cm^2/molecule'
                p_unit='bar'
            if kdata_unit is not None:
                conv_factor=u.Unit(rm_molec(self.kdata_unit)).to(u.Unit(rm_molec(kdata_unit)))
                data_to_write=self.kdata*conv_factor
                f.create_dataset("kcoeff", data=data_to_write,
                    compression=compression, compression_opts=compression_level)
                f["kcoeff"].attrs["units"] = kdata_unit
            else:
                f.create_dataset("kcoeff", data=self.kdata,
                    compression=compression, compression_opts=compression_level)
                f["kcoeff"].attrs["units"] = self.kdata_unit
            f.create_dataset("method", (1,), data=self.sampling_method, dtype=dt)
            f.create_dataset("samples", data=self.ggrid,
                compression=compression, compression_opts=compression_level)
            f.create_dataset("weights", data=self.weights,
                compression=compression, compression_opts=compression_level)
            f.create_dataset("ngauss", data=self.Ng)
            f.create_dataset("bin_centers", data=self.wns,
                compression=compression, compression_opts=compression_level)
            f.create_dataset("bin_edges", data=self.wnedges,
                compression=compression, compression_opts=compression_level)

            # where most of the data is actually written
            self.write_hdf5_common(f, compression=compression, compression_level=compression_level,
            p_unit=p_unit)

    def read_LMDZ(self, path=None, res=None, band=None, mol=None):
        """Initializes k coeff table and supporting data from a .dat file in a gcm friendly format.
        
        Units are assumed to be cm^2 for kdata and mbar for pressure. 

        Parameters
        ----------
            path: str
                Name of the directory with the various input files
            res: str
                "IRxVI" where IR and VI are the numbers of bands
                in the infrared and visible of the k table to load.
            band: str
                "IR" or "VI" to specify which band to load.
            mol: str
                Name of the molecule to be saved in the Ktable object. 
        """        
        if (path is None) or (res is None): \
            raise TypeError("You should provide an input directory name and a resolution")

        self.filename=path
        if mol is not None:
            self.mol=mol
        else:
            self.mol=os.path.basename(self.filename).split(self._settings._delimiter)[0]

        self.weights=np.loadtxt(os.path.join(path,'g.dat'),skiprows=1)[:-1]
            # we remove the last point that is always zero.
            # in the gcm this last point is intended to take care of continuum
        self.Ng=self.weights.size
        self.gedges=np.insert(np.cumsum(self.weights),0,0.)
        self.ggrid=(self.gedges[1:]+self.gedges[:-1])*0.5

        self.p_unit='mbar'
        self.logpgrid=np.loadtxt(os.path.join(path,'p.dat'),skiprows=1)*1.
        self.Np=self.logpgrid.size
        self.pgrid=10**self.logpgrid

        self.tgrid=np.loadtxt(os.path.join(path,'T.dat'),skiprows=1)
        self.Nt=self.tgrid.size

        if band is None:
            raw=np.loadtxt(os.path.join(path,res,'narrowbands.in'), skiprows=1, unpack=True)
        else:
            raw=np.loadtxt(os.path.join(path,res,'narrowbands_'+band+'.in'), \
                skiprows=1, unpack=True)
        self.wnedges=np.append(raw[0],raw[1,-1])
        self.wns=(self.wnedges[1:]+self.wnedges[:-1])*0.5
        self.Nw=self.wns.size
        
        self.kdata_unit='cm^2/molecule'
        if band is None:
            file_to_load=os.path.join(path,res,'corrk_gcm.dat')
        else:
            file_to_load=os.path.join(path,res,'corrk_gcm_'+band+'.dat')        
        tmp=np.loadtxt(file_to_load).reshape((self.Nt,self.Np,self.Nw,self.Ng+1),order='F')
        self.kdata=tmp[:,:,:,:-1].transpose((1,0,2,3)) 
        # also removing the last g point which is equal to 0.
        self.logk=False        
        return None

    def write_LMDZ(self, path, band='IR', fmt='%22.15e', write_only_metadata=False):
        """Saves data in a LMDZ friendly format.
        
        The gcm requires p in mbar and kdata in cm^2/molec.
        The conversion is done automatically.

        Parameters
        ----------
            path: str
                Name of the directory to be created and saved,
                the one that will contain all the necessary files
            band: str
                The band you are computing: 'IR' or 'VI'
            fmt: str
                Fortran format for the corrk file. 
        """
        try:
            os.mkdir(path)
        except FileExistsError:
            print('Directory was already there '+path)
        file = open(os.path.join(path,'p.dat'), "w")
        file.write(str(self.Np)+'\n')
        lp_to_write=self.logpgrid+np.log10(u.Unit(self.p_unit).to(u.Unit('mbar')))
        for lp in lp_to_write:
            file.write(str(lp)+'\n')
        file.close()

        file = open(os.path.join(path,'T.dat'), "w")
        file.write(str(self.Nt)+'\n')
        for t in self.tgrid:
            file.write(str(t)+'\n')
        file.close()

        file = open(os.path.join(path,'g.dat'), "w")
        file.write(str(self.Ng+1)+'\n')
        for g in self.weights:
            file.write(str(g)+'\n')
        file.write(str(0.)+'\n')
        file.close()

        dirname=os.path.join(path,band+str(self.Nw))
        try:
            os.mkdir(dirname)
        except FileExistsError:
            print('Directory was already there: '+dirname)

        file = open(os.path.join(dirname,'narrowbands_'+band+'.in'), "w")
        file.write(str(self.Nw)+'\n')
        for iw in range(self.Nw):
            file.write(str(self.wnedges[iw])+' '+str(self.wnedges[iw+1])+'\n')
        file.close()

        if not write_only_metadata:
            #file = open(dirname+'/corrk_gcm_IR.in', "w")
            data_to_write=self.kdata.transpose((1,0,2,3)).flatten(order='F')
            data_to_write=data_to_write*u.Unit(rm_molec(self.kdata_unit)).to(u.Unit('cm^2'))
            data_to_write=np.append(data_to_write, \
                np.zeros(self.Np*self.Nt*self.Nw)).reshape((1,self.Np*self.Nt*self.Nw*(self.Ng+1)))
            np.savetxt(os.path.join(dirname,'corrk_gcm_'+band+'.dat'),data_to_write,fmt=fmt)

    def read_nemesis(self, filename=None, mol=None):
        """Initializes k coeff table and supporting data from a Nemesis binary file (.kta)

        Parameters
        ----------
            file: str
                Name of the input Nemesis binary file.
            mol: str, optional
                Name of the molecule.
        """
        if (filename is None or not filename.lower().endswith('.kta')):
            raise RuntimeError("You should provide an input nemesis (.kta) file")
        
        self.Np, self.Nt, self.Nw, self.Ng, \
            self.pgrid, self.tgrid, self.wns, \
            self.ggrid, self.weights, self.kdata, \
            self.mol, self.isotopolog_id = read_nemesis_binary(filename)

        if self.mol is None:
            if mol is not None:
                self.mol=mol
            else:
                self.mol=os.path.basename(self.filename).split(self._settings._delimiter)[0]

        self.logpgrid=np.log10(self.pgrid)
        self.wnedges=np.concatenate(  \
            ([self.wns[0]],(self.wns[:-1]+self.wns[1:])*0.5,[self.wns[-1]]))
        self.gedges=np.insert(np.cumsum(self.weights),0,0.)
        self.p_unit='bar'
        self.kdata_unit='cm^2/molecule'

    def write_nemesis(self, filename):
        """Saves data in a nemesis format.

        Based on a routine provided by K. Chubb.

        Parameters
        ----------
            filename: str
                Name of the file to be created and saved
        """
        int_format = np.int32
        float_format = np.float32
        fullfilename=filename
        if not filename.lower().endswith(('.kta')):
            fullfilename=filename+'.kta'
        with open(fullfilename,'wb') as o:
            o.write(int_format(11+2*self.Ng+2+self.Np+self.Nt+self.Nw).tostring())
            #int for irec0 (5221)
            o.write(int_format(self.Nw).tostring()) #int
            o.write(float_format(self.wls[-1]).tostring()) #float for VMIN (0.30015)
            o.write(float_format(-1.).tostring()) #float
            o.write(float_format(0.).tostring()) #float for FWHM
            o.write(int_format(self.Np).tostring())
            o.write(int_format(self.Nt).tostring())
            o.write(int_format(self.Ng).tostring())
            try:
                mol_id=nemesis_hitran_id_numbers[self.mol]
            except:
                mol_id=0
            o.write(int_format(mol_id).tostring()) #IDGAS1 FROM HITRAN
            o.write(int_format(self.isotopolog_id).tostring())
            o.write(float_format(self.ggrid).tostring())
            o.write(float_format(self.weights).tostring())
            o.write(float_format(0.).tostring()) #float for 0
            o.write(float_format(0.).tostring()) #float for 0
            conv_factor=u.Unit(self.p_unit).to(u.Unit('bar'))
            o.write(float_format(self.pgrid).tostring())
            o.write(float_format(self.tgrid).tostring())
            o.write(float_format(self.wls[::-1]).tostring())
            data_to_write=self.kdata[:,:,::-1,:].transpose(2,0,1,3)
            conv_factor=u.Unit(rm_molec(self.kdata_unit)).to(u.Unit('cm^2'))*1.e20
            data_to_write=data_to_write*conv_factor
            o.write(float_format(data_to_write).tostring())
    
    def read_exorem(self, filename, mol=None):
        """Reads data in an ExoREM .dat format

        Parameters
        ----------
            filename: str
                Name of the input file.
            mol: str, optional
                Name of the molecule.
        """
        if self.mol is None:
            if mol is not None:
                self.mol=mol
            else:
                self.mol=os.path.basename(self.filename).split(self._settings._delimiter)[0]
        with open(filename, 'r') as file:
            tmp = file.readline().split()
            self.Np=int(tmp[0])
            self.Nt=int(tmp[1])
            self.Nw=int(tmp[2])
            self.Ng=int(tmp[-1])
            #print(self.Np,self.Nt,self.Nw,self.Ng)
            self.pgrid=_read_array(file, self.Np,  N_per_line=5, revert=True)
            self.logpgrid=np.log10(self.pgrid)
            self.p_unit='Pa'
            for _ in range(self.Np):
                self.tgrid=_read_array(file, self.Nt,  N_per_line=5)
            self.ggrid=_read_array(file, self.Ng,  N_per_line=5, revert=False)
            self.weights=_read_array(file, self.Ng,  N_per_line=5, revert=False)
            self.gedges=np.concatenate(([0],np.cumsum(self.weights)))

            self.kdata=np.zeros((self.Np, self.Nt, self.Nw, self.Ng))
            self.wns=np.zeros((self.Nw))
            #read_new_line=True
            tmp = file.readline().split()
            for iW in range(self.Nw):
                self.wns[iW]=float(tmp[0])
                #if iW=self.Nw-1: read_new_line=False
                for iP in range(self.Np):
                    for iT in range(self.Nt):
                        self.kdata[-iP-1,iT,iW]=_read_exorem_k_array(file, self.Ng)
                        #if read_new_line: tmp = file.readline().split()
                        tmp = file.readline().split()
            self.kdata_unit='cm^2/molec'
            self.wnedges=np.concatenate( \
                ([self.wns[0]],(self.wns[:-1]+self.wns[1:])*0.5,[self.wns[-1]]))
                
    def read_arcis(self, filename=None, mol=None):
        """Initializes k coeff table and supporting data from an ARCI fits file (.fits)

        Parameters
        ----------
            file: str
                Name of the input fits file.
            mol: str, optional
                Name of the molecule.
        """
        with fits.open(filename) as f:
            self.kdata=f[0].data[:,:,:,::-1] #reverse spectral axis
            self.kdata=np.transpose(self.kdata, axes=(0,1,3,2))
            self.kdata_unit='cm^2/molecule'
            self.Np, self.Nt, self.Nw, self.Ng = self.kdata.shape
            self.tgrid=f[1].data
            self.pgrid=f[2].data
            self.logpgrid=np.log10(self.pgrid)
            self.p_unit='bar'
            self.wns=1./f[3].data[::-1] #data given in wavelength in cm
            self.wn_unit='cm^-1'
            wns_max=1./f[0].header['L_MIN']
            wns_min=1./f[0].header['L_MAX']
            self.wnedges=np.concatenate( \
                ([wns_min],(self.wns[:-1]+self.wns[1:])*0.5,[wns_max]))
            self.weights,self.ggrid,self.gedges=gauss_legendre(self.Ng)

        if mol is not None:
            self.mol=mol
        else:
            self.mol=os.path.basename(filename).split(self._settings._delimiter)[0]

    def write_arcis(self, filename):
        """Saves data in an ARCIS fits format.

        Parameters
        ----------
            filename: str
                Name of the file to be created and saved
        """
        hdr = fits.Header()
        p_conv_fac=u.Unit(self.p_unit).to(u.Unit('bar'))
        hdr['TMIN'] = self.tgrid[0]
        hdr['TMAX'] = self.tgrid[-1]
        hdr['PMIN'] = self.pgrid[0]*p_conv_fac
        hdr['PMAX'] = self.pgrid[-1]*p_conv_fac
        hdr['L_MIN'] = 1./self.wns[-1]
        hdr['L_MAX'] = 1./self.wns[0]
        hdr['NT'] = self.Nt
        hdr['NP'] = self.Np
        hdr['NLAM'] = self.Nw
        hdr['NG'] = self.Ng
        data_to_write=self.kdata[:,:,::-1,:].transpose(0,1,3,2)
        conv_factor=u.Unit(rm_molec(self.kdata_unit)).to(u.Unit('cm^2'))
        hdu0=fits.PrimaryHDU(data_to_write*conv_factor, header=hdr)
        hdu1=fits.ImageHDU(self.tgrid)
        hdu2=fits.ImageHDU(self.pgrid*p_conv_fac)
        hdu3=fits.ImageHDU(1./self.wns[::-1])
        hdul = fits.HDUList([hdu0, hdu1, hdu2, hdu3])

        fullfilename=filename
        if not filename.lower().endswith(('.fits')):
            fullfilename=filename+'.fits'
        hdul.writeto(fullfilename, overwrite=True)

    def read_pickle(self, filename=None):
        """Initializes k coeff table and supporting data from an Exomol pickle file

        Parameters
        ----------
            filename: str
                Name of the input pickle file
            mol: str, optional
                Force the name of the molecule
        """
        if filename is None: raise RuntimeError("You should provide an input pickle filename")
        pickle_file=open(filename,'rb')
        raw=pickle.load(pickle_file, encoding='latin1')
        pickle_file.close()
        
        self.mol=raw['name']
        if self.mol=='H2OP': self.mol='H2O'

        self.pgrid=raw['p']
        self.logpgrid=np.log10(self.pgrid)
        self.tgrid=raw['t']
        self.wns=raw['bin_centers']
        self.wnedges=raw['bin_edges']
        if 'weights' in raw.keys():
            self.weights=raw['weights']
        else:
            raise RuntimeError('No weights keyword. This file is probably a cross section file.')
        self.ggrid=raw['samples']
        self.gedges=np.insert(np.cumsum(self.weights),0,0.)
        self.kdata=raw['kcoeff']
    
        if 'p_unit' in raw.keys():
            self.p_unit=raw['p_unit']
        else:
            self.p_unit='bar'

        if 'kdata_unit' in raw.keys():
            self.kdata_unit=raw['kdata_unit']
        else:
            self.kdata_unit='cm^2/molec'

        if 'wn_unit' in raw.keys(): self.wn_unit=raw['wn_unit']

        self.Np,self.Nt,self.Nw,self.Ng=self.kdata.shape

    def write_pickle(self, filename):
        """Saves data in a pickle format

        Parameters
        ----------
            filename: str
                Name of the file to be created and saved
        """
        fullfilename=filename
        if not filename.lower().endswith('.pickle'): fullfilename=filename+'.pickle'
        pickle_file=open(fullfilename,'wb')
        dictout={'name':self.mol,
                 'p':self.pgrid,
                 'p_unit':self.p_unit,
                 't':self.tgrid,
                 'bin_centers':self.wns,
                 'bin_edges':self.wnedges,
                 'wn_unit':self.wn_unit,
                 'weights':self.weights,
                 'samples':self.ggrid,
                 'kcoeff':self.kdata,
                 'kdata_unit':self.kdata_unit}
        #print(dictout)
        pickle.dump(dictout,pickle_file,protocol=-1)
        pickle_file.close()

    def spectrum_to_plot(self, p=1.e-5, t=200., x=1., g=None):
        """Dummy function to be defined in inheriting classes
        """
        raise NotImplementedError()

def read_nemesis_binary(filename):
    """reads a nemesis binary file.
    """
    f = open(filename, mode='rb')
    int_array = array.array('i')
    float_array = array.array('f')
    int_array.fromfile(f, 2)
    irec0, Nw=int_array[-2:]
    float_array.fromfile(f, 3)
    wl_min, dwl, FWHM = float_array[-3:]
    int_array.fromfile(f, 5)
    Np, Nt, Ng, mol_id, isotopolog_id = int_array[-5:]
    try:
        mol_name=list(nemesis_hitran_id_numbers.keys())[ \
            list(nemesis_hitran_id_numbers.values()).index(mol_id)]
    except:
        mol_name=None
    float_array.fromfile(f, Ng)
    ggrid = np.array(float_array[-Ng:])
    float_array.fromfile(f, Ng)
    weights = np.array(float_array[-Ng:])
    float_array.fromfile(f, 2)
    float_array.fromfile(f, Np)
    pgrid = np.array(float_array[-Np:])
    float_array.fromfile(f, Nt)
    tgrid = np.array(float_array[-Nt:])
    ntot=Nw*Np*Nt*Ng
    if dwl>=0.: #regular grid
        wls=wl_min+np.arange(Nw)*dwl
    else:
        float_array.fromfile(f, Nw)
        wls = np.array(float_array[-Nw:])
    wns=10000./wls[::-1]
    f.close()
    f = open(filename, mode='rb') # restart to start reading kdata at record number irec0
    kdata=array.array('f')
    kdata.fromfile(f, irec0-1+ntot)
    kdata=np.reshape(kdata[-ntot:],(Nw,Np,Nt,Ng))[::-1]*1.e-20
    kdata=kdata.transpose(1,2,0,3)
    f.close()
    return Np, Nt, Nw, Ng, pgrid, tgrid, wns, ggrid, weights, \
        kdata, mol_name, isotopolog_id
