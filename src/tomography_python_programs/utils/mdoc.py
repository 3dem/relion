import os
from os import PathLike
from pathlib import Path
from typing import Optional, List, Sequence

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


def match_filenames(
        source: Sequence[PathLike], to_match: Sequence[PathLike]
) -> List:
    """Match filenames which have common basenames."""
    cache = {basename(f): f for f in to_match}
    mdoc_basenames = get_tilt_image_basenames(source)
    return [cache.get(name, None) for name in mdoc_basenames]


def get_tilt_image_basenames(mdoc_tilt_image_files: Sequence[os.PathLike]) -> map:
    return map(lambda x: Path(str(x).split("\\")[-1]).stem, mdoc_tilt_image_files)


def construct_tomogram_id(mdoc_file: Path, prefix: str) -> str:
    """Construct a tomogram ID."""
    if prefix == '':
        return f'{basename(mdoc_file)}'
    else:
        return f'{prefix}_{basename(mdoc_file)}'
