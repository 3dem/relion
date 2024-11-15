from pathlib import Path
from typing import List, Optional

from pydantic import BaseModel


class GuiTiltSeries(BaseModel):
    name: str
    tilt_image_files: List[Path]
    tomogram_file: Optional[Path] = None
    denoised_tomogram_file: Optional[Path] = None

    def __len__(self):
        return len(self.tilt_image_files)
