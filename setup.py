#!/usr/bin/env python
"""
@author: Aurelien Falco
"""
from distutils.core import setup
from setuptools import find_packages

packages = find_packages(exclude=('tests', 'doc'))
provides = ['exo_k', ]

requires = []

install_requires = ['numpy',
                    'scipy',
                    'numba',
                    'astropy',
                    'h5py',
                    ]

console_scripts = ['']

entry_points = {'console_scripts': console_scripts, }

classifiers = [
    'Development Status :: 4 - Beta',
    'Environment :: Console',
    'Intended Audience :: Developers',
    'Intended Audience :: Education',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    'Operating System :: MacOS',
    'Operating System :: POSIX',
    'Operating System :: POSIX :: Linux',
    'Operating System :: Unix',
    'Operating System :: OS Independent',
    'Programming Language :: Python :: 3',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Astronomy',
    'Topic :: Scientific/Engineering :: Atmospheric Science',
    'Topic :: Software Development :: Libraries',
]

# Handle versioning
version = '0.0.4'

with open("README.md", "r") as fh:
    long_description = fh.read()


setup(name='exo_k',
      author='Jeremy Leconte',
      author_email='jeremy.leconte@u-bordeaux.fr',
      license="GPLv3",
      version=version,
      description='Library to handle radiative opacities from various sources for atmospheric applications',
      classifiers=classifiers,
      packages=packages,
      long_description=long_description,
      url='https://forge.oasu.u-bordeaux.fr/jleconte/exo_k-public',
      long_description_content_type="text/markdown",
      keywords = ['opacities','cross sections','corr-k','spectra','atmosphere','atmospheric',
        'radiative transfer'],
      include_package_data=True,
      entry_points=entry_points,
      provides=provides,
      requires=requires,
      install_requires=install_requires,
      extras_require={
        'Plot':  ["matplotlib"], },
      )
