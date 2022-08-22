"""Simple data models for GUI views."""
from pathlib import Path
from typing import Optional, Dict, List

from pydantic import BaseModel


class TiltSeries(BaseModel):
    name: str
    tilt_image_files: List[Path]
    tomogram_file: Optional[Path] = None

    def __len__(self):
        return len(self.tilt_image_files)


TiltSeriesSet = Dict[str, TiltSeries]
