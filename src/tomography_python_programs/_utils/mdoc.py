import os
import warnings
import rich
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .file import basename


def calculate_pre_exposure_dose(
    df: pd.DataFrame,
    dose_per_tilt_image: Optional[float] = None,
    dose_per_movie_frame: Optional[float] = None,
) -> np.ndarray:
    """Assumes that mdoc dataframe is already sorted by datetime."""
    if dose_per_tilt_image is not None and dose_per_movie_frame is not None:
        raise ValueError(
            'only one of dose_per_tilt_image and dose_per_movie_frame can be set.'
        )
    dose_override_provided = (
        dose_per_tilt_image is not None
        or dose_per_movie_frame is not None
    )
    if "ExposureDose" in df.columns and dose_override_provided is False:
        pre_exposure_dose = np.cumsum(df["ExposureDose"].to_numpy())
    elif dose_per_tilt_image is not None:
        pre_exposure_dose = dose_per_tilt_image * np.arange(len(df))
    elif dose_per_movie_frame is not None:
        cumulative_frame_number = np.cumsum(df["NumSubFrames"])
        post_exposure_dose = dose_per_movie_frame * cumulative_frame_number
        pre_exposure_dose = np.pad(post_exposure_dose, pad_width=(1, 0))[:-1]
    else:
        warnings.warn(
            'no dose information found in mdoc or provided, defaulting to zero.'
        )
        pre_exposure_dose = [0] * len(df)
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


def remove_missing_images(df: pd.DataFrame, mdoc_file: Path) -> pd.DataFrame:
    """Removes images from dataframe if image specified in mdoc file cannot be found, and warns user"""
    console = rich.console.Console(record=True)
    missing_images = df['SubFramePath'].loc[df['tilt_image_file'].isnull()]
    for missing_image in missing_images:
        missing_image = Path(str(missing_image).split("\\")[-1])
        console.log(
            f'WARNING: Image {missing_image} from {mdoc_file} not found in image files.')
    return df[df['tilt_image_file'].notna()]
