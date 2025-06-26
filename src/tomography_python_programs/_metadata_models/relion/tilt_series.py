import os

import pandas as pd
import starfile
from pydantic import BaseModel


class RlnTiltSeries(BaseModel):
    name: str
    data: pd.DataFrame

    class Config:
        arbitrary_types_allowed = True

    @classmethod
    def from_star_file(cls, filename: os.PathLike):
        star = starfile.read(filename, always_dict=True, parse_as_string=['rlnTomoName'])
        tilt_series_id = list(star.keys())[0]
        return cls(name=tilt_series_id, data=starfile.read(filename, parse_as_string=['rlnTomoName']))

    def write_star_file(self, filename: os.PathLike):
        starfile.write({self.name: self.data}, filename, overwrite=True)

    def __len__(self):
        return len(self.data)
