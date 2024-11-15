import os
from pathlib import Path

import pytest

from tomography_python_programs._metadata_models.relion.tilt_series_set import RlnTiltSeriesSet


@pytest.fixture
def test_project_directory() -> Path:
    return Path(__file__).parent / 'test_project_directory'


@pytest.fixture
def tilt_series_star_file(test_project_directory) -> Path:
    return test_project_directory / 'tilt_series_metadata' / 'tilt_series' / '220728_TS_01.star'


@pytest.fixture
def global_relion_metadata_file(test_project_directory) -> Path:
    return test_project_directory / 'tilt_series_metadata' / 'tilt_series.star'


@pytest.fixture
def relion_metadata(global_relion_metadata_file) -> RlnTiltSeriesSet:
    return RlnTiltSeriesSet.from_star_file(global_relion_metadata_file)


@pytest.fixture
def run_in_project_directory(monkeypatch, test_project_directory):
    monkeypatch.chdir(test_project_directory)
