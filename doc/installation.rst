Getting `exo_k`
===============

There are several ways to install `exo_k`. Because the code
relies on several other python libraries that will be installed
automatically as well, it is strongly recommended
that you use your own python environment (not the one used
by your system), especially if you are
working on a cluster.

If you do not have such a python environment setup yet, the
`Anaconda <https://www.anaconda.com/products/individual>`_ or
`Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ (lighter)
distributions are good options. With any of these, you can create your
own environment with a specific version of python and get access to
the `pip` Python Packages manager and can directly install the library
as shown below. (See the documentation of `conda` on how to use it)


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


As an alternative to `conda` you can also easily install `exo_k` from source with `poetry <https://python-poetry.org>`_.
This tool allow you to create a virtual environment, and manage the packages inside it.
This way, you can be sure that the dependencies installed use the exact version that have been tested.



Using the module
----------------

However you installed it,
you should now be able to import `exo_k` into your own python module using::

    import exo_k as xk

Whether you are using this command in a jupyter notebook or a script, 
you should however make sure that it is using the environment where you have installed the package.
See the Getting Started section for a complete tutorial.