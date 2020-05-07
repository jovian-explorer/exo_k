Installation
============

The exo_k package
-----------------

The exo_k library is hosted at https://forge.oasu.u-bordeaux.fr/jleconte/exo_k-public. You can download using::

    git clone https://forge.oasu.u-bordeaux.fr/jleconte/exo_k-public.git

Once you have downloaded it, you can move into the downloaded folder, and install exo_k with::

    pip install -e .

You should now be able to import it into your own python module::

    import exo_k

Building the documentation
---------------------------

The documentation relies on `sphinx`. The necessary packages can be installed using::

    pip install nbsphinx sphinx-autoapi sphinx_rtd_theme
    conda install sphinx # installs more (required) dependencies than pip

Then the documentation can be produced with::

    python setup.py doc

This will create the documentation under, which you can check out by opening the file `doc/html/index.html`.