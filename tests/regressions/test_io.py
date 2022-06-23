import pytest

from ..units.test_io import *  # noqa




@fixture(params=['h5', 'hdf5'])
def rw_h5(request, tmp_path, synthetic_table) -> xk.Ktable:
    filename = tmp_path / f'table.{request.param}'

    synthetic_table.write_hdf5(filename=filename.as_posix())
    return xk.Ktable(filename=filename.as_posix())


@fixture()
def rw_pickle(tmp_path, synthetic_table) -> xk.Ktable:
    filename = tmp_path / 'table.pickle'

    synthetic_table.write_pickle(filename=filename.as_posix())
    return xk.Ktable(filename=filename.as_posix())


@fixture()
def rw_nemesis(tmp_path, synthetic_table) -> xk.Ktable:
    filename = tmp_path / 'table.kta'

    synthetic_table.write_nemesis(filename=filename.as_posix())
    return xk.Ktable(filename=filename.as_posix())


@fixture()
def rw_exorem(tmp_path, synthetic_table) -> xk.Ktable:
    filename = tmp_path / 'table.ktable.exorem.txt.h5'

    synthetic_table.write_hdf5(filename=filename.as_posix(), exomol_units=True)
    return xk.Ktable(filename=filename.as_posix())


@fixture()
def rw_arcis(tmp_path, synthetic_table) -> xk.Ktable:
    filename = tmp_path / 'table.fits'

    synthetic_table.write_arcis(filename=filename.as_posix())
    return xk.Ktable(filename=filename.as_posix())


class TestKTable:
    @staticmethod
    def test_read_write_h5(rw_h5, synthetic_table):
        assert_same_ktable(rw_h5, synthetic_table, same_source=False)

    @staticmethod
    def test_read_write_pickle(rw_pickle, synthetic_table):
        assert_same_ktable(rw_pickle, synthetic_table, same_source=False)

    @staticmethod
    @pytest.mark.xfail(reason='Issue with the pression (difference of 1e5 -> bar to Pa issue)')
    def test_read_write_nemesis(rw_nemesis, synthetic_table):
        assert_same_ktable(rw_nemesis, synthetic_table, same_source=False, decimal=2)

    @staticmethod
    def test_read_write_exorem(rw_exorem, synthetic_table):
        assert_same_ktable(rw_exorem, synthetic_table, same_source=False, decimal=8)

    @staticmethod
    def test_read_write_arcis(rw_arcis, synthetic_table):
        assert_same_ktable(rw_arcis, synthetic_table, same_source=False, decimal=8)
