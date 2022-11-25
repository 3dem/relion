from typing import Optional

import numpy as np
from pydantic import BaseModel


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
