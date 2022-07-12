import random
import uuid
from typing import Optional

import h5py
import numpy as np
from numpy.testing import assert_array_equal, assert_almost_equal
from pytest import fixture

import exo_k as xk


@fixture()
def synthetic_table(tmp_path, set_up_xk):
    n_t = random.randint(a=5, b=10)
    t = np.linspace(start=random.randint(a=100, b=1000),
                    stop=random.randint(a=2000, b=4000),
                    num=n_t)

    n_p = random.randint(a=5, b=10)
    p = 10. ** np.linspace(-5., 2., n_p)

    mol_name = random.choice(['CO2', 'H2O'])

    samples = [0.0034357, 0.01801404, 0.04388279, 0.08044151, 0.12683405, 0.18197316,
               0.2445665, 0.31314696, 0.38610707, 0.46173674, 0.53826326, 0.61389293,
               0.68685304, 0.7554335, 0.81802684, 0.87316595, 0.91955849, 0.95611721,
               0.98198596, 0.9965643]

    weights = [0.008807003569576637, 0.02030071490019311, 0.03133602416705472, 0.041638370788352336,
               0.05096505990862013, 0.05909726598075912, 0.06584431922458826, 0.07104805465919094,
               0.07458649323630183, 0.07637669356536289, 0.07637669356536289, 0.07458649323630183,
               0.07104805465919094, 0.06584431922458826, 0.05909726598075912, 0.05096505990862013,
               0.041638370788352336, 0.03133602416705472, 0.02030071490019311, 0.008807003569576637]

    n_bin = 10
    bin_offset = 10
    bin_centers = np.linspace(0.5, 9.5, n_bin) + bin_offset
    bin_edges = np.linspace(0, 10, n_bin + 1) + bin_offset

    kcoeff = np.ones((len(p), len(t), len(bin_centers), len(weights)))

    filename = tmp_path / f'{str(uuid.uuid4().hex)}.h5'

    with h5py.File(name=filename, mode='w') as f:
        f.create_dataset('mol_name', data=mol_name)

        f.create_dataset('t', data=t)
        d = f.create_dataset('p', data=p)
        d.attrs['units'] = 'bar'

        f.create_dataset('bin_centers', data=bin_centers)
        d = f.create_dataset('bin_edges', data=bin_edges)
        d.attrs['units'] = 'cm^-1'

        f.create_dataset('samples', data=samples)
        f.create_dataset('weights', data=weights)

        d = f.create_dataset('kcoeff', data=kcoeff)
        d.attrs['units'] = 'cm^2/molecule'

    table = xk.Ktable(filename=filename.as_posix())


    yield table


@fixture(scope='session', params=['CO2', 'H2O'])
def table_molecule(request) -> str:
    yield request.param


@fixture()
def table_from_molecule(table_molecule) -> xk.Ktable:
    yield xk.Ktable(mol=table_molecule)


@fixture()
def table_from_regex(table_molecule) -> xk.Ktable:
    yield xk.Ktable(table_molecule)


@fixture()
def table_from_filename(table_from_molecule) -> xk.Ktable:
    yield xk.Ktable(filename=table_from_molecule.filename)

@fixture(scope='session', params=[[8,15], [13,18], [10,20]])
def wn_range(request) -> str:
    synthetic_table
    yield request.param

@fixture()
def table_from_range(synthetic_table, wn_range) -> xk.Ktable:
    yield xk.Ktable(filename=synthetic_table.filename, wn_range=wn_range)

# noinspection PyPep8Naming
def assert_same_ktable(K: xk.Ktable, L: xk.Ktable, same_source: bool = True,
                       decimal: Optional[int] = None) -> None:
    if same_source:
        assert K.filename == L.filename

    assert_array_equal(K.shape, L.shape)

    if K.mol != 'table.fits' and L.mol != 'table.fits':
        assert K.mol == L.mol
        assert K.molar_mass == L.molar_mass

    assert K.wn_unit == L.wn_unit
    if decimal is None:
        assert_array_equal(K.wls, L.wls)
    else:
        assert_almost_equal(K.wls, L.wls, decimal=decimal)

    assert K.p_unit == L.p_unit
    if decimal is None:
        assert_array_equal(K.pgrid, L.pgrid)
    else:
        assert_almost_equal(K.pgrid, L.pgrid, decimal=decimal)

    if decimal is None:
        assert_array_equal(K.tgrid, L.tgrid)
    else:
        assert_almost_equal(K.tgrid, L.tgrid, decimal=decimal)

    assert K.kdata_unit == L.kdata_unit
    if decimal is None:
        assert_array_equal(K.kdata, L.kdata)
    else:
        assert_almost_equal(K.kdata, L.kdata, decimal=decimal)

    if decimal is None:
        assert_array_equal(K.weights, L.weights)
    else:
        assert_almost_equal(K.weights, L.weights, decimal=decimal)


class TestKTable:
    @staticmethod
    def disabled_test_read_from_molecule(table_from_molecule):
        assert table_from_molecule is not None

        assert_same_ktable(table_from_molecule, table_from_molecule)

    @staticmethod
    def disabled_test_read_from_regex(table_from_regex, table_from_molecule):
        assert table_from_regex is not None

        assert_same_ktable(table_from_regex, table_from_molecule)

    @staticmethod
    def disabled_test_read_from_filename(table_from_filename, table_from_molecule):
        assert table_from_filename is not None

        assert_same_ktable(table_from_filename, table_from_molecule)

    @staticmethod
    def test_read_synthetic(synthetic_table):
        assert table_from_filename is not None

    @staticmethod
    def test_read_wn_range(table_from_range, wn_range):
        table = table_from_range
        assert table is not None
        assert table.wnedges.min() >= wn_range[0]
        assert table.wnedges.max() <= wn_range[1]
        assert table.wns.min() >= wn_range[0]
        assert table.wns.max() <= wn_range[1]
