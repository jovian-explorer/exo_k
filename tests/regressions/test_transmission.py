import astropy.units as u
import numpy as np
import pytest

import exo_k as xk
from exo_k import Atm, Spectrum


@pytest.mark.skip(reason='Should be treated as regression')
class TestTransmission:
    def test_model(self, ProblemHotJupiter):
        w_range = [5, 10] * u.Unit('um')

        database = xk.Kdatabase(['H2O'], 'R300', '0.3-50mu', remove_zeros=True)

        database.clip_spectral_range(
            np.sort(w_range.to_value(u.Unit('cm-1'), equivalencies=u.spectral())).tolist())

        cia_data = xk.CIAdatabase(molecules=['H2'], mks=True)
        cia_data.sample(database.wns)

        n_lay: int = 25

        atm_ck: Atm = xk.Atm(psurf=ProblemHotJupiter.p_surf.to_value('Pa'),
                             ptop=ProblemHotJupiter.p_strat.to_value('Pa'),
                             Tsurf=ProblemHotJupiter.T_surf.to_value('K'),
                             Tstrat=ProblemHotJupiter.T_strat.to_value('K'),
                             grav=ProblemHotJupiter.g.to_value('m* s**-2'),
                             composition={'H2': 'background', 'H2O': 1.e-3},
                             Nlay=n_lay,
                             Rp=ProblemHotJupiter.Rp.to_value('m'),
                             k_database=database,
                             cia_database=cia_data)

        s: Spectrum = atm_ck.transmission_spectrum()  # noqa

        assert s is not None
