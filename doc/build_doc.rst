Building the documentation
==========================

If you are reading this, you probably do not need to build the documentation. 

But, as we made the choice not to store previous versions of the documentation online,
you may want to build a local version of the documentation
for the version of the code you are currently using.
The documentation relies on `sphinx`. The necessary packages can be installed using::

    pip install nbsphinx sphinx-autoapi sphinx_rtd_theme
    conda install sphinx # installs more (required) dependencies than pip

Then the documentation can be produced with::

    python setup.py doc

This will create the documentation that you can check out by opening `doc/html/index.html`.