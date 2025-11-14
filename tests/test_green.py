import numpy as np
import json
from typing import Dict
from frida.integral import (
    calc_complete_elliptic_integral,
    green_brbz_filament,
    green_flux_filament,
)


def laod_and_prepare_ref_data(path: str) -> Dict[str, np.ndarray]:
    data = json.load(open(path, "r"))
    for k in data.keys():
        data[k] = np.array(data[k])
    return data


def test_elliptic_integrals():
    data = laod_and_prepare_ref_data("tests/fixtures/elliptic_integrals.json")

    k, e = calc_complete_elliptic_integral(data["m"])
    assert np.allclose(k, data["k"], atol=1e-15)
    assert np.allclose(e, data["e"], atol=1e-15)

    # check for corner cases m = 1 and m = 0
    k, e = calc_complete_elliptic_integral(np.array(1))
    assert np.isinf(k).all()
    assert (e == 1).all()


def test_green_flux():
    data = laod_and_prepare_ref_data("tests/fixtures/green.json")
    flux = green_flux_filament(data["source"], data["target"])
    assert np.allclose(flux, data["flux"], atol=1e-15)


def test_green_br_bz():
    data = laod_and_prepare_ref_data("tests/fixtures/green.json")
    br, bz = green_brbz_filament(data["source"], data["target"])
    assert np.allclose(br, data["br"], atol=1e-15)
    assert np.allclose(bz, data["bz"], atol=1e-15)
    ...
