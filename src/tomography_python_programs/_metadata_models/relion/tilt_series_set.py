import os
from pathlib import Path
from typing import List, Optional, Iterable, Tuple

import pandas as pd
import starfile
from pydantic import BaseModel, validator

from .tilt_series import RlnTiltSeries
from ..gui.tilt_series_set import GuiTiltSeriesSet as GuiTiltSeriesSet
from ..gui.tilt_series import GuiTiltSeries as GuiTiltSeries


class RlnTiltSeriesSet(BaseModel):
    global_data: pd.DataFrame
    tilt_series: List[RlnTiltSeries]

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
    def from_star_file(
            cls,
            filename: os.PathLike,
            tilt_series_id: Optional[str] = None
    ):
        star = starfile.read(filename, always_dict=True, parse_as_string=['rlnTomoName'])
        global_df = star['global'].set_index('rlnTomoName', drop=False)
        if tilt_series_id is not None:  # get single tilt-series
            global_df = global_df.loc[tilt_series_id]
            tilt_series_star_file = global_df['rlnTomoTiltSeriesStarFile']
            tilt_series = [RlnTiltSeries.from_star_file(tilt_series_star_file)]
        else:  # get all tilt-series
            tilt_series = [
                RlnTiltSeries.from_star_file(tilt_series_star_file)
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

    def as_gui_model(self) -> GuiTiltSeriesSet:
        tilt_series_set = {
            tilt_series.name: GuiTiltSeries(
                name=tilt_series.name,
                tilt_image_files=tilt_series.data['rlnMicrographName'].to_list()
            )
            for tilt_series
            in self.tilt_series
        }
        if 'rlnTomoReconstructedTomogram' in self.global_data:
            for gui_tilt_series, tomogram_file in zip(
                    self.tilt_series, self.global_data['rlnTomoReconstructedTomogram']
            ):
                tilt_series_set[gui_tilt_series.name].tomogram_file = tomogram_file
        if 'rlnTomoReconstructedTomogramDenoised' in self.global_data:
            for gui_tilt_series, tomogram_file in zip(
                    self.tilt_series, self.global_data['rlnTomoReconstructedTomogramDenoised']
            ):
                tilt_series_set[gui_tilt_series.name].denoised_tomogram_file = tomogram_file
        return tilt_series_set

    def __getitem__(self, tilt_series_id: str) -> RlnTiltSeries:
        return self.tilt_series[self._tilt_series_names.index(tilt_series_id)]

    def __iter__(self) -> Iterable[Tuple[pd.Series, RlnTiltSeries]]:
        return ((self.global_data.loc[tilt_series.name], tilt_series) for
                tilt_series in self.tilt_series)
