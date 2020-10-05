Getting `exo_k`
===============

There are several ways to install `exo_k`.


Installation using `pip`
------------------------

If the `pip` Python Packages manager is installed on your machine,
getting `exo_k` is as simple as running::

    pip install exo_k


Installing from the source files
--------------------------------

The `exo_k` library is hosted at https://forge.oasu.u-bordeaux.fr/jleconte/exo_k-public.
You can download it as follows::

    git clone https://forge.oasu.u-bordeaux.fr/jleconte/exo_k-public.git

Once you have downloaded it, you can move into the downloaded folder, and install `exo_k` with::

    pip install -e .

or::

    python setup.py install


Using the module
----------------

However you installed it,
you should now be able to import `exo_k` into your own python module using::

    import exo_k as xk

See the Getting Started section for a tutorial.