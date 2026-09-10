from dataclasses import dataclass
from typing import Dict, List
import numpy as np


@dataclass
class Channel:
    """A named time series.

    ``time`` is stored as Unix time in milliseconds, matching the Julia project.
    """
    name: str
    time: np.ndarray
    values: np.ndarray

    def __post_init__(self):
        self.time = np.asarray(self.time, dtype=float)
        self.values = np.asarray(self.values, dtype=float)


@dataclass
class Dataset:
    """Collection of named channels."""
    channels: Dict[str, Channel]

    def __getitem__(self, key: str) -> Channel:
        return self.channels[key]

    def __getattr__(self, key: str):
        if key != "channels" and "channels" in self.__dict__ and key in self.channels:
            return self.channels[key]
        raise AttributeError(key)

    def __len__(self):
        return len(self.channels)

@dataclass
class SMPSData:
    name: str
    time: np.ndarray          # datetime objects or np.nan for missing
    diameters: np.ndarray     # bin midpoints [nm]
    dNdlogDp: np.ndarray      # shape: time × bins

    def __post_init__(self):
        self.time = np.asarray(self.time, dtype=object)
        self.diameters = np.asarray(self.diameters, dtype=float)
        self.dNdlogDp = np.asarray(self.dNdlogDp, dtype=float)

        if self.dNdlogDp.ndim != 2:
            raise ValueError("SMPSData.dNdlogDp must be a 2D matrix")

        expected = (len(self.time), len(self.diameters))
        if self.dNdlogDp.shape != expected:
            raise ValueError(
                "SMPSData.dNdlogDp must have shape "
                f"(len(time), len(diameters)) = {expected}, "
                f"got {self.dNdlogDp.shape}"
            )