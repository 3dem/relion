from __future__ import annotations

from typing import TYPE_CHECKING, Tuple, Optional
from rich.console import Console

import pandas as pd
import starfile

console = Console(record=True)

if TYPE_CHECKING:
    import os
    from typing import Iterable

# (tilt_series_id, tilt_series_df, tilt_image_df)
TiltSeriesMetadata = Tuple[str, pd.DataFrame, pd.DataFrame]


def _iterate_all_tilt_series_metadata(
        tilt_series_star_file: os.PathLike
) -> Iterable[TiltSeriesMetadata]:
    """Yield all metadata from a tilt-series data STAR file."""
    star = starfile.read(tilt_series_star_file, always_dict=True)
    for _, tilt_series_df in star['global'].iterrows():
        tilt_series_id = tilt_series_df['rlnTomoName']
        tilt_image_df = starfile.read(
            tilt_series_df['rlnTomoTiltSeriesStarFile'], always_dict=True
        )[tilt_series_id]
        yield tilt_series_id, tilt_series_df, tilt_image_df


def _extract_single_tilt_series_metadata(
        tilt_series_star_file: os.PathLike,
        tilt_series_id: str
) -> TiltSeriesMetadata:
    """Get metadata for a specific tilt-series from a tilt-series data STAR file."""
    star = starfile.read(tilt_series_star_file, always_dict=True)
    # Check if Tomogram Name provided is actually found in the star file
    if not star['global']['rlnTomoName'].str.contains(tilt_series_id).any():
        e = 'Specified Tomogram provided in Tomogram-Name is not found in rlnTomoName in the given star file.'
        console.log(f'ERROR: {e}')
        raise RuntimeError(e)
    tilt_series_df = star['global'].set_index('rlnTomoName').loc[tilt_series_id, :]
    tilt_image_df = starfile.read(
        tilt_series_df['rlnTomoTiltSeriesStarFile'], always_dict=True
    )[tilt_series_id]
    return tilt_series_id, tilt_series_df, tilt_image_df


def iterate_tilt_series_metadata(
        tilt_series_star_file: os.PathLike,
        tilt_series_id: Optional[str] = None,
) -> Iterable[TiltSeriesMetadata]:
    """Yield metadata from a tilt-series data STAR file."""
    if tilt_series_id is None:  # align all tilt-series
        all_tilt_series_metadata = _iterate_all_tilt_series_metadata(
            tilt_series_star_file)
    else:  # do single tilt-series alignment
        all_tilt_series_metadata = [
            _extract_single_tilt_series_metadata(tilt_series_star_file, tilt_series_id)
        ]
    for tilt_series_metadata in all_tilt_series_metadata:
        yield tilt_series_metadata
