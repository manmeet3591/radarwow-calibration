"""Calibration from the weather itself: Ryzhkov, Giangrande, Melnikov & Schuur 2005, J. Atmos. Oceanic Technol. 22,
1138-1155, doi:10.1175/JTECH1772.1 (as summarised in radar_wow/literature/reading_notes/04_ryzhkov2005_dualpol_calibration.md).

Z (self-consistency, S band). In rain Z, ZDR and KDP are linked: Z = a + b log10 KDP + c ZDR (dBZ, deg/km, dB). Rain
gates pooled over an area and time are binned by 1 dB of Z; per bin the mean KDP and ZDR are taken; the Z bias is
the offset d for which the sum over bins of measured KDP equals the sum of KDP predicted from (Z - d, ZDR). The paper
gives two S-band relations for Oklahoma drop-size regimes, large drops (LD: a, b, c = 44.0, 12.2, 2.32) and small
drops (SD: 46.0, 9.59, 1.68), to use whichever fits better; it accepts a bias only if the mean |log10 KDP - model|
is below 0.07 after correction. ZDR must be calibrated first (it enters the relation). Other bands need their own
relations, which the paper does not give.

ZDR (natural targets). Dry aggregated snow just above the melting layer has an intrinsic ZDR of at most ~0.25 dB at
low elevation, falling roughly as cos^2(elevation) (the paper: ZDR(60 deg) ~ 0.24 ZDR(0 deg)); the offset is the
median measured ZDR there minus 0.2 cos^2(elevation) dB (0.2 dB: radar-wow choice inside the paper's <= 0.25 dB).
Light rain (20-28 dBZ) is a weaker target (the paper: intrinsic ZDR 0-1.5 dB); radar-wow uses an expected median of
0.25 dB there only as a cross-check.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from radarwow.core.noise import signal_mask

RELATIONS = {"large drops": (44.0, 12.2, 2.32), "small drops": (46.0, 9.59, 1.68)}


def _elev(ds):
    return float(np.nanmedian(ds["elevation"].values))


def _height_km(ds):
    from radarwow.estimation.vpr import beam_height_km

    return np.broadcast_to(beam_height_km(ds["range"].values, _elev(ds)), ds["DBZH"].shape)


def zdr_offset_snow(sweeps, ml_top_km, min_elevation: float = 10.0, depth_km: float = 1.0, z_range=(15.0, 35.0),
                    rhohv_min: float = 0.97, intrinsic_db: float = 0.2):
    """ZDR offset (dB) from dry snow in the ``depth_km`` above the melting-layer top, on sweeps >= ``min_elevation``.

    ``sweeps``: iterable of sweep Datasets (several volumes may be pooled). Returns (offset, number of gates)."""
    diffs = []
    for ds in sweeps:
        el = _elev(ds)
        if el < min_elevation or "ZDR" not in ds or "RHOHV" not in ds:
            continue
        h = _height_km(ds)
        z, zdr, rho = ds["DBZH"].values, ds["ZDR"].values, ds["RHOHV"].values
        m = signal_mask(ds) & (h > ml_top_km + 0.2) & (h < ml_top_km + 0.2 + depth_km) & (z >= z_range[0]) & (z <= z_range[1]) & (rho > rhohv_min)
        m &= np.isfinite(zdr)
        diffs.append(zdr[m] - intrinsic_db * np.cos(np.deg2rad(el)) ** 2)
    d = np.concatenate(diffs) if diffs else np.array([])
    return (float(np.median(d)) if d.size >= 500 else np.nan), int(d.size)


def zdr_offset_light_rain(sweeps, ml_bottom_km, z_range=(20.0, 28.0), rhohv_min: float = 0.98, expected_db: float = 0.25):
    """Cross-check: median ZDR of light rain below the melting layer minus ``expected_db``. Returns (offset, n)."""
    vals = []
    for ds in sweeps:
        if "ZDR" not in ds or "RHOHV" not in ds:
            continue
        h = _height_km(ds)
        z = ds["DBZH"].values
        m = signal_mask(ds) & (h < ml_bottom_km) & (z >= z_range[0]) & (z <= z_range[1]) & (ds["RHOHV"].values > rhohv_min)
        vals.append(ds["ZDR"].values[m & np.isfinite(ds["ZDR"].values)])
    v = np.concatenate(vals) if vals else np.array([])
    return (float(np.median(v) - expected_db) if v.size >= 500 else np.nan), int(v.size)


def rain_bins(sweeps, ml_bottom_km, zdr_offset_db: float = 0.0, z_range=(20.0, 50.0), rhohv_min: float = 0.98, max_pia_db: float = 1.0):
    """Pool rain gates and return per 1-dB Z bin: (centre Z, n, mean KDP, mean ZDR). Sweeps need KDP (radarwow-kdp)
    and, if available, PIA (radarwow-attenuation) so attenuated paths can be left out."""
    zs, ks, ds_ = [], [], []
    for ds in sweeps:
        h = _height_km(ds)
        z = ds["DBZH"].values
        m = signal_mask(ds) & (h < ml_bottom_km) & (z >= z_range[0]) & (z < z_range[1]) & (ds["RHOHV"].values > rhohv_min)
        m &= np.isfinite(ds["KDP"].values) & np.isfinite(ds["ZDR"].values)
        if "PIA" in ds:
            m &= ds["PIA"].values <= max_pia_db
        zs.append(z[m])
        ks.append(ds["KDP"].values[m])
        ds_.append(ds["ZDR"].values[m] - zdr_offset_db)
    z, k, d = (np.concatenate(a) if a else np.array([]) for a in (zs, ks, ds_))
    edges = np.arange(z_range[0], z_range[1] + 1e-9, 1.0)
    out = []
    for lo in edges[:-1]:
        s = (z >= lo) & (z < lo + 1)
        if s.sum() >= 50:
            out.append((lo + 0.5, int(s.sum()), float(np.mean(k[s])), float(np.mean(d[s]))))
    return np.array(out).reshape(-1, 4)


def self_consistency_z_bias(bins, relation: str | None = None, max_misfit: float = 0.07):
    """Z bias (dB, measured minus true) from binned rain (`rain_bins`), Ryzhkov et al. 2005.

    Returns dict(bias, relation, misfit, accepted, n). With ``relation`` None both are tried and the one with the
    smaller misfit is kept. misfit = n-weighted mean |log10 KDP - log10 KDP_model| after correction, over bins with
    mean KDP > 0.05 deg/km."""
    if len(bins) < 5:
        return dict(bias=np.nan, relation=None, misfit=np.nan, accepted=False, n=int(bins[:, 1].sum()) if len(bins) else 0)
    zc, n, kdp, zdr = bins.T
    best = None
    for name in ([relation] if relation else RELATIONS):
        a, b, c = RELATIONS[name]
        model = lambda dz: 10 ** ((zc - dz - a - c * zdr) / b)  # noqa: E731
        f = lambda dz: np.sum(n * model(dz)) - np.sum(n * kdp)  # noqa: E731
        try:
            dz = brentq(f, -15, 15)
        except ValueError:
            continue
        ok = kdp > 0.05
        mis = float(np.sum(n[ok] * np.abs(np.log10(kdp[ok]) - np.log10(model(dz)[ok]))) / np.sum(n[ok])) if ok.any() else np.nan
        if best is None or mis < best["misfit"]:
            best = dict(bias=float(dz), relation=name, misfit=mis)
    if best is None:
        return dict(bias=np.nan, relation=None, misfit=np.nan, accepted=False, n=int(n.sum()))
    best["accepted"] = bool(np.isfinite(best["misfit"]) and best["misfit"] < max_misfit)
    best["n"] = int(n.sum())
    return best
