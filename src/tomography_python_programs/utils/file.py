import os
from pathlib import Path
from typing import List


def basename(path: Path) -> str:
    while path.stem != str(path):
        path = Path(path.stem)
    return path


def glob(pattern: str) -> List[Path]:
    return list(Path().glob(pattern))


def write_empty_file(filename: os.PathLike) -> None:
    open(filename, mode='w').close()
