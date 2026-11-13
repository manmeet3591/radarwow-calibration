"""radar-wow stage 1: calibration from natural targets (Ryzhkov et al. 2005)."""

from .natural import RELATIONS, rain_bins, self_consistency_z_bias, zdr_offset_light_rain, zdr_offset_snow

PLANNED = ["Warren2018GPMMatch", "Crisologo2018"]
__all__ = ["RELATIONS", "rain_bins", "self_consistency_z_bias", "zdr_offset_light_rain", "zdr_offset_snow"]
