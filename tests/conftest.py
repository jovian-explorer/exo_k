from pathlib import Path

from pytest import fixture

import exo_k as xk


@fixture()
def dataset_root_path(tmp_path) -> Path:
    (tmp_path / 'radiative_data' / 'corrk').mkdir(parents=True, exist_ok=True)
    (tmp_path / 'radiative_data' / 'xsec').mkdir(parents=True, exist_ok=True)
    (tmp_path / 'radiative_data' / 'aerosol').mkdir(parents=True, exist_ok=True)
    (tmp_path / 'radiative_data' / 'cia').mkdir(parents=True, exist_ok=True)

    (tmp_path / 'gcm_data').mkdir(parents=True, exist_ok=True)

    yield tmp_path.resolve()


@fixture(scope='function', autouse=True)
def ensure_Settings_not_created():
    xk.Settings.reset_singleton()


@fixture()
def settings(ensure_Settings_not_created, dataset_root_path):
    s = xk.Settings()

    data_corrk: Path = dataset_root_path / 'radiative_data' / 'corrk'
    data_xsec: Path = dataset_root_path / 'radiative_data' / 'xsec'
    data_aerosol = dataset_root_path / 'radiative_data' / 'aerosol'
    data_cia: Path = dataset_root_path / 'radiative_data' / 'cia'

    s.set_search_path(data_corrk.as_posix(), path_type='ktable')
    s.set_search_path(data_xsec.as_posix(), path_type='xtable')
    s.set_search_path(data_cia.as_posix(), path_type='cia')
    s.set_search_path(data_aerosol.as_posix(), path_type='aerosol')

    yield s


@fixture(autouse=True)
def set_up_xk(settings, dataset_root_path):
    settings.set_mks(True)
    settings.set_log_interp(True)

    yield settings
