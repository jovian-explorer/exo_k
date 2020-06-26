# -*- coding: utf-8 -*-
"""
@author: jeremy leconte

A class based on a singleton to store global options only once for every instance. 
"""
import os.path
from glob import glob
from fnmatch import fnmatch,fnmatchcase
from .util.singleton import Singleton

class Settings(Singleton):
    """A class based on a singleton to store global options only once for every instance.

    So all the following methods can be called using the following syntax:

    >>> exo_k.Settings().method_name(agrs)

    In gerneal, they will change internal global attributes
    that change the global behavior of some routines in the library. 
    """

    def init(self, *args, **kwds):
        self.reset_search_path()
        self._log_interp = True
        self._convert_to_mks = False
        self._delimiter = '_'
        self._case_sensitive = False

    def reset_search_path(self):
        """Set default search path.
        """
        self._search_path = [os.path.abspath('.')]

    def add_search_path(self, *search_paths):
        """Add path(s) to the list of paths that will be searched for
        correlated-k and x-sec files.

        Parameters
        ----------
            search_path : string or list of strings
                Search path(s) to look for opacities.

        Examples
        --------

            >>> exo_k.Settings().add_search_path('data/xsec','data/corrk','data/cia')
            >>> exo_k.Settings().search_path()
            ['/Users/jleconte/atmosphere/RadiativeTransfer/exo_k',
             '/Users/jleconte/atmosphere/RadiativeTransfer/exo_k/data/xsec',
             '/Users/jleconte/atmosphere/RadiativeTransfer/exo_k/data/corrk',
             '/Users/jleconte/atmosphere/RadiativeTransfer/exo_k/data/cia']
        """
        for path in search_paths:
            if not os.path.isdir(path):
                raise NotADirectoryError("""The search_path you provided
                    does not exist or is not a directory""")
            else:
                if os.path.abspath(path) not in self._search_path:
                    self._search_path.append(os.path.abspath(path))

    def set_search_path(self, *search_paths):
        """Sets the path(s) that will be searched for correlated-k and x-sec files .

        Parameters
        ----------
            search_path : string or list of strings
                Search path(s) to look for opacities.
        """
        if not os.path.isdir(search_paths[0]):
            raise NotADirectoryError("""The search_path you provided
                does not exist or is not a directory""")
        else:
            self._search_path=[os.path.abspath(search_paths[0])]
        if len(search_paths)>1: self.add_search_path(*search_paths[1:])

    def search_path(self):
        """Returns the current value of the global search path (_search_path)
        """
        return self._search_path

    def set_delimiter(self, newdelimiter):
        """Sets the delimiter string used to separate molecule names in filenames.

        Parameters
        ----------
            newdelimiter: string
                New delimiter to use. Default is '_'.

        Exmaple
        -------
            If I have a file named 'H2O.R10000_xsec.hdf5'
            that I want to load in a `Kdatabase`, the default
            settings will result in an error:

            >>> database=xk.Kdatabase(['H2O'],'R10000')
             No file was found with these filters: 
             ('H2O_', 'R1000') in the following directories:
             ['/home/falco/xsec/xsec_sampled_R10000_0.3-15'] 

            Using

            >>> xk.Settings().set_delimiter('.')
            >>> database=xk.Kdatabase(['H2O'],'R10000')

            finds the file.
        """
        self._delimiter = newdelimiter

    def set_log_interp(self, log_interp):
        """Sets the default interpolation mode for kdata. Default is Log. 

        Parameters
        ----------
            log_interp: boolean
                If True, log interpolation. Linear if False.
        """
        self._log_interp = log_interp

    def set_case_sensitive(self, case_sensitive):
        """Set whether name matching is case sensitive. Default is False.

        Parameters
        ----------
            case_sensitive: boolean
                If True, name matching is case sensitive.
        """
        self._case_sensitive = case_sensitive

    def set_mks(self, set_mks):
        """Forces conversion to mks system.
        
        Parameters
        ----------
            set_mks: boolean
                If True, all datasets are converted to mks upon loading.
        """
        self._convert_to_mks = set_mks

    def list_files(self, *str_filters, only_one = False, search_path = None):
        """A routine that provides a list of all filenames containing
        a set of string filters in the global _search_path or a local one.

        Parameters
        ----------
            *str_filters: str
                A set of strings that need to be contained in the name of the file
            only_one: boolean, optional
                If true, only one filename is returned (the first one).
                If false, a list is returned. Default is False.
            search_path: str, optional
                If search_path is provided, it locally overrides
                the global _search_path settings
                and only files in search_path are returned.

        Returns
        -------
            list of strings
                List of filenames corresponding to all the str_filters
        """
        local_search_path=self._search_path
        if search_path is not None:
            local_search_path=[search_path]

        filenames = [f for path in local_search_path for f in glob(os.path.join(path,'*'))]
        finalnames=filenames[:]
        for filename in filenames:
            if os.path.isdir(filename):
                finalnames.remove(filename)
                continue
            if self._case_sensitive:
                for filt in str_filters:
                    if not fnmatchcase(os.path.basename(filename),'*'+filt+'*'):
                        finalnames.remove(filename)
                        break
            else:
                for filt in str_filters:
                    if not fnmatch(os.path.basename(filename.lower()),'*'+filt.lower()+'*'):
                        finalnames.remove(filename)
                        break
        if len(finalnames)>1 and only_one: 
            print("""be careful: {filt}
            String filters not specific enough, several corresponding files have been found.
            We will use the first one:
            {file}
            Other files are:
            {other_files}""".format(filt=str_filters,file=finalnames[0], \
                other_files=finalnames[1:]))
            finalnames=[finalnames[0]]
        if not finalnames:
            print("""No file was found with these filters: 
            {filt}
            in the following directories:
            {path}
            """.format(filt=str_filters,path=local_search_path))
            raise RuntimeError()
        # an empty sequence yields False in a conditional statement ! 
        return finalnames

