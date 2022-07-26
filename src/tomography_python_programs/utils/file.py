import os
from os import PathLike
from pathlib import Path
from typing import List, Sequence


def basename(path: os.PathLike) -> str:
    while Path(path).stem != str(path):
        path = Path(path).stem
    return str(path)


def write_empty_file(filename: os.PathLike) -> None:
    open(filename, mode='w').close()


def match_filenames(
        source: Sequence[PathLike], to_match: Sequence[PathLike]
) -> List:
    """Match filenames which have common basenames."""
    cache = {basename(f): f for f in to_match}
    source_file_basenames = [basename(f) for f in source]
    return [cache.get(name, None) for name in source_file_basenames]
