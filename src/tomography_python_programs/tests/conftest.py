from pathlib import Path

import pytest

from tomography_python_programs.metadata_model import RelionTiltSeriesSet


@pytest.fixture
def test_data_directory() -> Path:
    return Path(__file__).parent / 'test_data'


@pytest.fixture
def tilt_series_star_file(test_data_directory) -> Path:
    return test_data_directory / 'tilt_series_metadata' / 'tilt_series' / '220728_TS_01.star'


@pytest.fixture
def global_relion_metadata_file(test_data_directory) -> Path:
    return test_data_directory / 'tilt_series_metadata' / 'tilt_series.star'


@pytest.fixture
def relion_metadata(global_relion_metadata_file) -> RelionTiltSeriesSet:
    return RelionTiltSeriesSet.from_star_file(global_relion_metadata_file)