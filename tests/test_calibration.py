import numpy as np

from radarwow.calibration import RELATIONS, self_consistency_z_bias


def _bins(bias_db, relation="small drops", seed=0):
    a, b, c = RELATIONS[relation]
    rng = np.random.default_rng(seed)
    zc = np.arange(20.5, 48.5)
    zdr = 0.3 + 0.05 * (zc - 20)  # ZDR rising with Z
    kdp = 10 ** ((zc - a - c * zdr) / b) * np.exp(rng.normal(0, 0.05, zc.size))  # true Z
    n = np.full(zc.size, 1000.0)
    return np.c_[zc + bias_db, n, kdp, zdr]  # measured Z = true + bias


def test_recovers_known_bias_and_relation():
    for bias in (-3.0, 0.0, 2.5):
        r = self_consistency_z_bias(_bins(bias))
        assert abs(r["bias"] - bias) < 0.3 and r["relation"] == "small drops" and r["accepted"]


def test_too_few_bins_rejected():
    r = self_consistency_z_bias(_bins(1.0)[:3])
    assert not r["accepted"] and np.isnan(r["bias"])
