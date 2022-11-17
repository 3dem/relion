"""Simple data models for GUI views."""
from pathlib import Path
from typing import Optional, Dict, List

import numpy as np
from pydantic import BaseModel


class TiltSeries(BaseModel):
    name: str
    tilt_image_files: List[Path]
    tomogram_file: Optional[Path] = None

    def __len__(self):
        return len(self.tilt_image_files)


TiltSeriesSet = Dict[str, TiltSeries]


class LazyTiltSeriesData(BaseModel):
    name: str
    data: np.ndarray
    last_loaded_index: Optional[int]
    n_images_loaded: int = 0

    class Config:
        arbitrary_types_allowed = True

    @property
    def n_images(self):
        return len(self.data)

    @property
    def completed(self):
        return self.n_images_loaded == self.n_images
