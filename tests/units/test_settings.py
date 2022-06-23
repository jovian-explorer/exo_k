from pathlib import Path

import pytest
from pytest import fixture

import exo_k as xk
from exo_k.settings import NoFileFoundError


@fixture()
def project_path() -> Path:
    return Path('.').resolve()


class TestSettings:
    @staticmethod
    def test_dataset_existing(dataset_root_path):
        """
        This test ensures that the dataset rootpath exist, ie that the submodule is pulled or
         that the user manually added data."""

        assert dataset_root_path.exists()
        assert dataset_root_path.is_dir()

        assert (dataset_root_path / 'gcm_data').is_dir()
        assert (dataset_root_path / 'radiative_data').is_dir()

    @staticmethod
    def test_default_path(dataset_root_path, settings):
        """
        This test ensures that the settings fixture is correctly initialized and,
         that `xk.Settings()` is in the same state as the fixtures
        """

        data_corrk: Path = dataset_root_path / 'radiative_data' / 'corrk'
        data_xsec: Path = dataset_root_path / 'radiative_data' / 'xsec'
        data_aerosol = dataset_root_path / 'radiative_data' / 'aerosol'
        data_cia: Path = dataset_root_path / 'radiative_data' / 'cia'

        assert data_corrk.as_posix() in xk.Settings().search_path['ktable']
        assert data_xsec.as_posix() in xk.Settings().search_path['xtable']
        assert data_cia.as_posix() in xk.Settings().search_path['cia']
        assert data_aerosol.as_posix() in xk.Settings().search_path['aerosol']

        assert len(xk.Settings().search_path['ktable']) == 1
        assert len(xk.Settings().search_path['xtable']) == 1
        assert len(xk.Settings().search_path['cia']) == 1
        assert len(xk.Settings().search_path['aerosol']) == 1

        assert len(settings.search_path['ktable']) == 1
        assert len(settings.search_path['xtable']) == 1
        assert len(settings.search_path['cia']) == 1
        assert len(settings.search_path['aerosol']) == 1

        assert settings.search_path == xk.Settings().search_path

    @staticmethod
    def test_list_files(settings):
        """This test check that we find the correct numbers of files"""
        with pytest.raises(NoFileFoundError):
            assert len(settings.list_cia_files()) == 6
        with pytest.raises(NoFileFoundError):
            assert len(settings.list_cia_files(molecule_pair=['H2', 'h2'])) == 1

        with pytest.raises(NoFileFoundError):
            assert len(settings.list_files(path_type='ktable')) == 2
        with pytest.raises(NoFileFoundError):
            assert len(settings.list_files(path_type='ktable', molecule='H2O')) == 1

        with pytest.raises(NoFileFoundError):
            assert len(settings.list_files(path_type='xtable')) == 1

        with pytest.raises(NoFileFoundError):
            assert len(settings.list_files(path_type='aerosol')) == 1
