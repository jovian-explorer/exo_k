#!/usr/bin/env python
import setuptools
from setuptools import find_packages
from numpy.distutils.core import setup
from numpy.distutils.core import Extension
from numpy.distutils import log
import re, os
import glob

packages = find_packages(exclude=('tests', 'doc'))
provides = ['exo_k', ]

requires = []

install_requires = ['numpy',
                    'cython',
                    'configobj',
                    'scipy',
                    'numba',
                    'astropy',
                    'numexpr',
                    'numpy',
                    'nestle',
                    'h5py',
                    'tabulate', ]

console_scripts = ['']

entry_points = {'console_scripts': console_scripts, }

classifiers = [
    'Development Status :: 4 - Beta',
    'Environment :: Console',
    'Environment :: Win32 (MS Windows)',
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Operating System :: Microsoft :: Windows',
    'Operating System :: POSIX',
    'Operating System :: POSIX :: Linux',
    'Operating System :: Unix',
    'Operating System :: OS Independent',
    'Programming Language :: Python',
    'Topic :: Scientific/Engineering',
    'Topic :: Software Development :: Libraries',
]

# Handle versioning
version = '1.0.0-Beta'

with open("README.md", "r") as fh:
    long_description = fh.read()


setup(name='Exo_k',
      author='J. Leconte',
      author_email='jeremy.leconte@u-bordeaux.fr',
      license="BSD",
      version=version,
      description='Library to handle radiative opacities from various sources for atmospheric applications',
      classifiers=classifiers,
      packages=packages,
      long_description=long_description,
      url='https://forge.oasu.u-bordeaux.fr/jleconte/exo_k',
      long_description_content_type="text/markdown",
      keywords = ['opacities','cross sections','corr-k','spectra','atmosphere','atmospheric'],
      include_package_data=True,
      entry_points=entry_points,
      provides=provides,
      requires=requires,
      install_requires=install_requires,
      extras_require={
        'Plot':  ["matplotlib"], },
      )
