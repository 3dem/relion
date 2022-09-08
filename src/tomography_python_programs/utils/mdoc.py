import os
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .file import basename


def calculate_pre_exposure_dose(
        mdoc_df: pd.DataFrame, dose_per_tilt: Optional[float] = None
) -> np.ndarray:
    if "exposure_dose" in mdoc_df.columns and dose_per_tilt is None:
        pre_exposure_dose = np.cumsum(mdoc_df["exposure_dose"].to_numpy())
    elif dose_per_tilt is not None:
        pre_exposure_dose = dose_per_tilt * np.arange(len(mdoc_df))
    else:
        pre_exposure_dose = [0] * len(mdoc_df)
    return pre_exposure_dose


def basename_from_sub_frame_path(filename: os.PathLike) -> str:
    """Get a basename from an mdoc 'SubFramePath' entry.

    Example of a 'SubFramePath' entry:
        "D:\\DATA\\Flo\\HGK149_20151130\\frames\\TS_01_000_0.0.mrc"
    """
    return basename(Path(str(filename).split("\\")[-1]))


def construct_tomogram_id(mdoc_file: Path, prefix: str) -> str:
    """Construct a tomogram ID."""
    if prefix == '':
        return f'{basename(mdoc_file)}'
    else:
        return f'{prefix}_{basename(mdoc_file)}'
