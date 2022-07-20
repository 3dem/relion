from os import PathLike
from pathlib import Path
from typing import Optional, List

import numpy as np
import pandas as pd

from .file import basename


def add_pre_exposure_dose(
        mdoc_df: pd.DataFrame, dose_per_tilt: Optional[float] = None
) -> pd.DataFrame:
    if dose_per_tilt is not None:
        pre_exposure_dose = dose_per_tilt * np.arange(len(mdoc_df))
    else:  # all zeros if exposure dose values not present in mdoc
        pre_exposure_dose = [0] * len(mdoc_df)

    if "exposure_dose" not in mdoc_df.columns or dose_per_tilt is not None:
        mdoc_df["pre_exposure_dose"] = pre_exposure_dose
    else:
        mdoc_df["pre_exposure_dose"] = np.cumsum(mdoc_df["exposure_dose"].to_numpy())
    return mdoc_df


def add_tilt_image_files(
        mdoc_df: pd.DataFrame, tilt_image_files: List[PathLike]
) -> pd.DataFrame:
    basename_filename_map = {Path(f).stem: f for f in tilt_image_files}
    mdoc_tilt_image_basenames = get_tilt_image_basenames(mdoc_df["sub_frame_path"])
    mdoc_df["tilt_image_file"] = [
        basename_filename_map.get(basename, None)
        for basename in mdoc_tilt_image_basenames
    ]
    mdoc_df = mdoc_df[mdoc_df["tilt_image_file"] != None]
    return mdoc_df


def get_tilt_image_basenames(mdoc_tilt_image_files: pd.Series) -> pd.Series:
    tilt_image_basenames = mdoc_tilt_image_files.apply(
        lambda x: Path(str(x).split("\\")[-1]).stem
    )
    return tilt_image_basenames


def construct_tomogram_id(mdoc_file: Path, prefix: str) -> str:
    """Construct a unique tomogram ID."""
    if prefix == '':
        return f'{utils.file.basename(mdoc_file)}'
    else:
        return f'{prefix}_{basename(mdoc_file)}'
