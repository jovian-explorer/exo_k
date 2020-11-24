# -*- coding: utf-8 -*-
"""
@author: jeremy leconte

A class based on a singleton to store global options only once for every instance. 
"""
import os.path
import re
from glob import glob
from string import Template
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
        self.set_delimiters(r"_.\-")

        self._case_sensitive = False

    def reset_search_path(self, path_type = 'all', no_path = False):
        """Set default search path.

        Parameters
        ----------
            path_type: str
                What type of path to reset. Possibilities are
                'all' (default), 'kdata', 'cia', and 'aerosol'
            no_path: bool (optional)
                If False (default), the search path is reset to '.'.
                If True, empties the search path.
        """
        if no_path:
            local_path=[]
        else:
            local_path=[os.path.abspath('.')]
        if path_type in ('kdata', 'ktable', 'all'):
            self._ktable_search_path = local_path.copy()
        if path_type in ('kdata', 'xtable', 'all'):
            self._xtable_search_path = local_path.copy()
        if path_type in ('cia', 'all'):
            self._cia_search_path = local_path.copy()
        if path_type in ('aerosol', 'all'):
            self._aerosol_search_path = local_path.copy()

    def add_search_path(self, *search_paths, path_type = 'kdata'):
        """Add path(s) to the list of paths that will be searched for
        various files.

        Parameters
        ----------
            search_path : string or list of strings
                Search path(s) to look for opacities.
            path_type: str
                What type of path to change. Possibilities are:
                  * 'kdata' (default): global xsec and corr-k search path
                  * 'ktable' : only corr-k files
                  * 'xtable' : only cross section files
                  * 'cia': search path for CIA files
                  * 'aerosol': search path for Aerosol optical property files

        Examples
        --------

            >>> exo_k.Settings().add_search_path('data/xsec','data/corrk')
            >>> exo_k.Settings().search_path()
            ['/your/path/to/exo_k',
             '/your/path/to/exo_k/data/xsec',
             '/your/path/to/exo_k/data/corrk']
        """
        path_to_change=[]
        if path_type == 'cia':
            path_to_change.append(self._cia_search_path)
        elif path_type == 'aerosol':
            path_to_change.append(self._aerosol_search_path)
        else:
            if path_type in ('kdata', 'ktable'):
                path_to_change.append(self._ktable_search_path)
            if path_type in ('kdata', 'xtable'):
                path_to_change.append(self._xtable_search_path)

        if path_to_change == []:
            raise RuntimeError('I did not understand the path type provided.')
        for path in search_paths:
            if not os.path.isdir(path):
                raise NotADirectoryError("""The search_path you provided
                    does not exist or is not a directory""")
            else:
                for to_change in path_to_change:
                    if os.path.abspath(path) not in to_change:
                        to_change.append(os.path.abspath(path))

    def add_cia_search_path(self, *search_paths):
        """Add path(s) to the list of paths that will be searched for
        cia files.

        Parameters
        ----------
            search_path : string or list of strings
                Search path(s) to look for opacities.
        """
        self.add_search_path(*search_paths, path_type = 'cia')

    def add_aerosol_search_path(self, *search_paths):
        """Add path(s) to the list of paths that will be searched for
        aerosol files.

        Parameters
        ----------
            search_path : string or list of strings
                Search path(s) to look for opacities.
        """
        self.add_search_path(*search_paths, path_type = 'aerosol')

    def set_search_path(self, *search_paths, path_type = 'kdata'):
        """Like :func:~`exo_k.settings.Settings.add_search_path` except for
        the fact that the path is reset first.
        """
        self.reset_search_path(path_type = path_type, no_path = True)
        self.add_search_path(*search_paths, path_type = path_type)

    def set_cia_search_path(self, *search_paths):
        """Like :func:~`exo_k.settings.Settings.add_cia_search_path` except for
        the fact that the path is reset first.
        """
        self.set_search_path(*search_paths, path_type = 'cia')

    def set_aerosol_search_path(self, *search_paths):
        """Like :func:~`exo_k.settings.Settings.add_aerosol_search_path` except for
        the fact that the path is reset first.
        """
        self.set_search_path(*search_paths, path_type = 'aerosol')

    @property
    def search_path(self):
        """Returns the current value of the global search path (_search_path)
        """
        return {'ktable': self._ktable_search_path, 'xtable': self._xtable_search_path,
             'cia': self._cia_search_path, 'aerosol': self._aerosol_search_path}

    def cia_search_path(self):
        """Returns the current value of the cia search path (_cia_search_path)
        """
        raise RuntimeError('Deprecated after v0.0.5. Use xk.Settings().search_path["cia"] instead.')

    def aerosol_search_path(self):
        """Returns the current value of the aerosol search path (_aerosol_search_path)
        """
        raise RuntimeError(
            'Deprecated after v0.0.5. Use xk.Settings().search_path["aerosol"] instead.')

    def set_delimiter(self, newdelimiter):
        """Sets the delimiter string used to separate molecule names in filenames.

        Parameters
        ----------
            newdelimiter: string
                New delimiter to use. Default is '_'.

        Example
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

    def set_delimiters(self, newdelimiters):
        """Sets the delimiter string used to separate molecule names in filenames.
        If you want to include the '-' character, you must put a '\' before to avoid
        it being interpreted as a special character by the `re` module.

        Parameters
        ----------
            newdelimiters: string
                New delimiters to use. Default is '_.\-'.
        """
        self._delimiters = newdelimiters
        self._search_mol_template = Template('(?=^$mol['+self._delimiters+  \
                ']|['+self._delimiters+']$mol['+self._delimiters+'])')
        self._search_cia_template = Template('(?=^$mol1\-$mol2['+self._delimiters+  \
                ']|^$mol2\-$mol1['+self._delimiters+'])')


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

    def list_files(self, *str_filters, molecule = None, 
            only_one = False, search_path = None, path_type = 'kdata'):
        """A routine that provides a list of all filenames containing
        a set of string filters in one of the global _search_path or a local one.

        Whether the search is case sensitive is specified through the
        Settings.set_case_sensitive() method. 

        .. warning::
            The pattern matching with the `str_filters` is done using regular expressions.
            If you want to match special characters (like a dot in a filename), do not
            forget to put a backslash in front of it.

        Parameters
        ----------
            *str_filters: str
                A set of strings that need to be contained in the name of the file
            molecule: str
                The name of a molecule to be looked for in the filename. It must be followed
                by one of the characters in self._delimiters and either at the begining of the name
                or just after one of the characters in self._delimiters.
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
        if search_path is not None:
            local_search_path=[search_path]
        else:
            if path_type == 'ktable':
                local_search_path=self._ktable_search_path
            elif path_type == 'xtable':
                local_search_path=self._xtable_search_path
            elif path_type == 'aerosol':
                local_search_path=self._aerosol_search_path

        filenames = [f for path in local_search_path for f in glob(os.path.join(path,'*'))]
        finalnames=filenames[:]
        if molecule is not None:
            if self._case_sensitive:
                template=self._search_mol_template.substitute(mol=molecule)
            else:
                template=self._search_mol_template.substitute(mol=molecule.lower())
        else:
            template='^'
        for filt in str_filters:
            if self._case_sensitive:
                template+='(?=.*'+filt+')'
            else:
                template+='(?=.*'+filt.lower()+')'
        #print(template)
        for filename in filenames:
            if os.path.isdir(filename):
                finalnames.remove(filename)
                continue
            if self._case_sensitive:
                basename=os.path.basename(filename)
            else:
                basename=os.path.basename(filename).lower()
            if re.search(template, basename) is None:
                finalnames.remove(filename)
                continue
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
            raise NoFileFoundError("""No file found to match this regular expression
            (based on your filters): 
            {filt}
            in the following directories:
            {path}
            """.format(filt=template,path=local_search_path))
        # an empty sequence yields False in a conditional statement ! 
        return finalnames

    def list_cia_files(self, *str_filters, molecule_pair = None, 
            only_one = False, search_path = None):
        """A routine that provides a list of all filenames containing
        a set of string filters in the global _search_path or a local one.

        Whether the search is case sensitive is specified through the
        Settings.set_case_sensitive() method. 

        .. warning::
            The pattern matching with the `str_filters` is done using regular expressions.
            If you want to match special characters (like a dot in a filename), do not
            forget to put a backslash in front of it.

        Parameters
        ----------
            *str_filters: str
                A set of strings that need to be contained in the name of the file
            molecule_pair: list of 2 str
                The name of the 2 molecules in the cia pair.
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
        if search_path is not None:
            local_search_path=[search_path]
        else:
            local_search_path=self._cia_search_path

        filenames = [f for path in local_search_path for f in glob(os.path.join(path,'*'))]
        finalnames=filenames[:]
        if molecule_pair is not None:
            if self._case_sensitive:
                template=self._search_cia_template.substitute(mol1=molecule_pair[0],
                    mol2=molecule_pair[1])
            else:
                template=self._search_cia_template.substitute(mol1=molecule_pair[0].lower(),
                    mol2=molecule_pair[1].lower())
        else:
            template='^'
        for filt in str_filters:
            if self._case_sensitive:
                template+='(?=.*'+filt+')'
            else:
                template+='(?=.*'+filt.lower()+')'
        #print(template)
        for filename in filenames:
            if os.path.isdir(filename):
                finalnames.remove(filename)
                continue
            if self._case_sensitive:
                basename=os.path.basename(filename)
            else:
                basename=os.path.basename(filename).lower()
            if re.search(template, basename) is None:
                finalnames.remove(filename)
                continue
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
            raise NoFileFoundError("""No file found to match this regular expression
            (based on your filters): 
            {filt}
            in the following directories:
            {path}
            """.format(filt=template,path=local_search_path))
        # an empty sequence yields False in a conditional statement ! 
        return finalnames

class NoFileFoundError(Exception):
    """Error when no file is found
    """
