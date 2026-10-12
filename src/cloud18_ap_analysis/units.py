"""Shared concentration conversions; loaded measurements are never modified."""
from __future__ import annotations

import numpy as np
from .models import Channel, Dataset
from .interpolation import resample

GAS_CONSTANT = 8.31446261815324  # Pa m3 / (mol K)
AVOGADRO_CONSTANT = 6.02214076e23  # molecules / mol (exact)
STANDARD_PRESSURE_PA = 101325.0
_MIXING_FRACTIONS = {"ppt": 1e-12, "pptv": 1e-12, "ppb": 1e-9, "ppbv": 1e-9}


def validate_concentration_unit(unit):
    if unit not in {"native", "molecules/cm3"}:
        raise ValueError('concentration_unit must be "native" or "molecules/cm3"')


def concentration_label(unit, native="pptv"):
    validate_concentration_unit(unit)
    return r"molecules/cm$^3$" if unit == "molecules/cm3" else native


def concentration_yscale(unit, override=None):
    validate_concentration_unit(unit)
    return override if override is not None else ("log" if unit == "molecules/cm3" else "linear")


def _condition_values(condition, times, name):
    if condition is None:
        raise ValueError(f"{name} is required for molecules/cm3 conversion")
    if isinstance(condition, Channel):
        values = resample(condition, times).values
        if not np.any(np.isfinite(values)):
            raise ValueError(f"No {name} measurements overlap the concentration timestamps")
        return values
    value = float(condition)
    if not np.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def convert_channel_concentration(channel, source_unit, target_unit, *,
                                  temperature_c=None, pressure_pa=STANDARD_PRESSURE_PA):
    """Return a new channel, converting pptv/ppbv to molecules/cm3 or mixing units.

    Legacy ppt/ppb spellings are accepted. Temperature (Celsius) and pressure
    (Pa) can be constants or Channels, interpolated without extrapolation.
    Missing conditions produce NaN; all-missing coverage raises an error.
    Temperature/pressure are treated as known; their uncertainties are not
    included. Names remain source identifiers, including suffixes like _ppt.
    """
    if source_unit not in _MIXING_FRACTIONS:
        raise ValueError("source_unit must be pptv or ppbv (ppt/ppb also accepted)")
    if target_unit == "native":
        target_unit = source_unit
    if target_unit in _MIXING_FRACTIONS:
        factor = _MIXING_FRACTIONS[source_unit] / _MIXING_FRACTIONS[target_unit]
    elif target_unit == "molecules/cm3":
        temperature_k = _condition_values(temperature_c, channel.time, "temperature_c") + 273.15
        pressure = _condition_values(pressure_pa, channel.time, "pressure_pa")
        if np.any(np.asarray(temperature_k) <= 0) or np.any(np.asarray(pressure) <= 0):
            raise ValueError("Absolute temperature and pressure must be positive")
        factor = _MIXING_FRACTIONS[source_unit] * pressure / (GAS_CONSTANT * temperature_k) * 1e-6 * AVOGADRO_CONSTANT
    else:
        raise ValueError("target_unit must be native, pptv, ppbv, or molecules/cm3")
    return Channel(channel.name, channel.time.copy(), channel.values * factor)


def convert_dataset_concentrations(data, source_units, target_unit, *,
                                   temperature_c=None, pressure_pa=STANDARD_PRESSURE_PA):
    """Copy a Dataset and convert only channels named in source_units."""
    missing = set(source_units).difference(data.channels)
    if missing:
        raise KeyError("Missing concentration channels: " + ", ".join(sorted(missing)))
    channels = {}
    for name, ch in data.channels.items():
        channels[name] = (convert_channel_concentration(
            ch, source_units[name], target_unit,
            temperature_c=temperature_c, pressure_pa=pressure_pa,
        ) if name in source_units else Channel(ch.name, ch.time.copy(), ch.values.copy()))
    return Dataset(channels)
