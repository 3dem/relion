import os
from pathlib import Path
from typing import List, Optional, Tuple, Iterable

import pandas as pd
import starfile
from pydantic import BaseModel, validator

from tomography_python_programs.gui.components.data_models import TiltSeries, TiltSeriesSet


class RelionTiltSeries(BaseModel):
    name: str
    data: pd.DataFrame

    class Config:
        arbitrary_types_allowed = True

    @classmethod
    def from_star_file(cls, filename: os.PathLike):
        star = starfile.read(filename, always_dict=True)
        tilt_series_id = list(star.keys())[0]
        return cls(name=tilt_series_id, data=starfile.read(filename))

    def write_star_file(self, filename: os.PathLike):
        starfile.write({self.name: self.data}, filename, overwrite=True)

    def __len__(self):
        return len(self.data)


class RelionTiltSeriesSet(BaseModel):
    global_data: pd.DataFrame
    tilt_series: List[RelionTiltSeries]

    class Config:
        arbitrary_types_allowed = True

    @validator('global_data')
    def _index_on_tilt_series_id(cls, value: pd.DataFrame) -> pd.DataFrame:
        return value.set_index('rlnTomoName', drop=False)

    @validator('global_data')
    def _drop_tilt_series_star_file_column(cls, value: pd.DataFrame) -> pd.DataFrame:
        return value.drop(columns=['rlnTomoTiltSeriesStarFile'])

    @property
    def _tilt_series_names(self) -> List[str]:
        return [tilt_series.name for tilt_series in self.tilt_series]

    @classmethod
    def from_star_file(cls, filename: os.PathLike,
                       tilt_series_id: Optional[str] = None):
        star = starfile.read(filename, always_dict=True)
        global_df = star['global'].set_index('rlnTomoName', drop=False)
        if tilt_series_id is not None:  # get single tilt-series
            global_df = global_df.loc[tilt_series_id]
            tilt_series_star_file = global_df['rlnTomoTiltSeriesStarFile']
            tilt_series = [RelionTiltSeries.from_star_file(tilt_series_star_file)]
        else:  # get all tilt-series
            tilt_series = [
                RelionTiltSeries.from_star_file(tilt_series_star_file)
                for tilt_series_star_file
                in global_df['rlnTomoTiltSeriesStarFile']
            ]
        return cls(global_data=global_df, tilt_series=tilt_series)

    def write_star_file(self, filename: os.PathLike):
        tilt_series_directory = Path(filename).parent / 'tilt_series'
        tilt_series_directory.mkdir(parents=True, exist_ok=True)
        tilt_series_star_files = [
            tilt_series_directory / f'{tilt_series.name}.star'
            for tilt_series
            in self.tilt_series
        ]
        global_df = self.global_data.copy()
        global_df['rlnTomoTiltSeriesStarFile'] = tilt_series_star_files
        for star_file_name, tilt_series in zip(tilt_series_star_files, self.tilt_series):
            if len(tilt_series) > 0:
                tilt_series.write_star_file(star_file_name)
            else:
                global_df = global_df.drop(tilt_series.name)
        starfile.write({'global': global_df}, filename, overwrite=True)

    def as_gui_model(self) -> TiltSeriesSet:
        tilt_series_set = {
            tilt_series.name: TiltSeries(
                name=tilt_series.name,
                tilt_image_files=tilt_series.data['rlnMicrographName'].to_list()
            )
            for tilt_series
            in self.tilt_series
        }
        return tilt_series_set

    def __getitem__(self, tilt_series_id: str) -> RelionTiltSeries:
        return self.tilt_series[self._tilt_series_names.index(tilt_series_id)]

    def __iter__(self) -> Iterable[Tuple[pd.Series, RelionTiltSeries]]:
        return ((self.global_data.loc[tilt_series.name], tilt_series) for
                tilt_series in self.tilt_series)
