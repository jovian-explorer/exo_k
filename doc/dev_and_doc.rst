=================
Development and documentation
=================


As we made the choice not to store previous versions of the documentation online,
you may want to build a local version of the documentation for the version of the code you are currently using.

To do that you can use the following resources to quickly be able to do it.

Quickly use the source
======================
- If you want to build the documentation, 
  only using a virtual environment and the following should be enough.
  The documentation will be available at `doc/_build/html/index.html`::

    pip install nbsphinx sphinx-autoapi sphinx_rtd_theme
    cd doc && make html

- If you want to test `exo_k`, only using a virtual environment
  and the following should be enough::

    pip install pytest
    pytest

Using Poetry
============

Whether you want to contribute to exo_k, build the documentation for your local version
or ensure that the tests passed, another way is to use
`poetry <https://python-poetry.org>`_ to setup your virtual environment.
All requirements are listed inside `pyproject.toml` and `poetry.lock`.
They ensure that everyone use the same dependencies.

After cloning the repository, you only need to execute the following
command to create your environment and install the dependencies::

    poetry install

After that you can spawn a shell with `poetry shell` or just run a command with `poetry run command`.

Documentation
-------------
The documentation can be produced with::

    cd doc && poetry run make

If you are editing it, you can locally serve it and refreshed on change with the following::

    cd doc && poetry run make livehtml

Tests
-----
You can easily run the tests with `pytest`::

    poetry run pytest

You can enable a more verbose output with the flag `-vv`, `pytest-clarity`
provide a better output in this verbose mode.

To run a specific test, just add the path to the test file as an argument.

Poetry recap
------------

You can read the documentation of poetry at this `address <https://python-poetry.org>`_ .

We use the latest version available to this day (`1.2.0.b2`) of `poetry`.
It's in preview, so if you installed it, you can use the following command to migrate to it::

    poetry self update --preview 1.2.0.b2

You can get the version of poetry with::

    poetry --version


You can add a package in the following manner::

    poetry add <important-pkg>
    poetry add -G dev <pkg-by-the-devs>

After you modify `pyproject.toml`, to manually add dependencies for example, lock it like that::

    poetry lock

To install the environment defined::

    poetry install

