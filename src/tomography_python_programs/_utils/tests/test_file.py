from pathlib import Path
import pytest
import numpy as np

from tomography_python_programs._utils.file import basename, match_filenames


@pytest.mark.parametrize(
    "input, expected",
    [
        (Path('test/test2/test3.mrc'), 'test3'),  # multiple directories
        (Path('test.mrc.mdoc'), 'test'),  # multiple suffixes
    ]
)
def test_basename(input: Path, expected: str):
    """Check that basename is working as expected."""
    output = basename(input)
    assert output == expected


def test_match_filenames():
    """Check that matching of filenames is working as expected."""
    source_filenames = ['test.mrc', 'test2_[0.00].mrc']
    filenames_to_match = ['test.mrc', 'test2_[0.00]', 'test3']
    output = match_filenames(source_filenames, filenames_to_match)
    expected = ['test.mrc', 'test2_[0.00]']
    np.testing.assert_array_equal(output, expected)
    assert len(output) == 2

