from pathlib import Path
from typing import Optional, List

import starfile
import rich
import typer
import subprocess
import os

from ._utils import (
    create_denoising_directory_structure,
    parse_training_tomograms,
    generate_training_tomograms_star,
    find_tomogram_halves,
    generate_train_config_json,
    generate_train_data_config_json,
    save_json,
    save_tilt_series_stars,
    save_global_star,
)
from .constants import TRAIN_CONFIG_PREFIX, TRAIN_DATA_CONFIG_PREFIX, MODEL_NAME
from .._cli import cli
from ..._utils.relion import relion_pipeline_job

console = rich.console.Console(record=True)
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '1'  


@cli.command(name='cryoCARE:train')
@relion_pipeline_job
def cryoCARE_train(
        tilt_series_star_file: Path = typer.Option(...),
        output_directory: Path = typer.Option(...),
        training_tomograms: Optional[str] = typer.Option(None),
        number_training_subvolumes: int = typer.Option(1200),
        subvolume_sidelength: int = typer.Option(72),
        gpu: Optional[List[int]] = typer.Option(None)
):
    """Train a denoising model using cryoCARE (>=v0.2.0)
    
    Requires that two tomograms have been generated using the same sample.
    These can be generated via taking odd/even frames during motion correction
    (optimal) or by taking odd/even tilts during tomogram reconstruction.

    The location of tomograms should be specified in the global star file for all tilt series with the headers:
    `rlnTomoReconstructedTomogramHalf1` and `rlnTomoReconstructedTomogramHalf2`.


    Parameters
    ----------
    tilt_series_star_file: Path
        RELION tilt-series STAR file.
    output_directory: Path
        directory in which results will be stored.
    training_tomograms: Optional[str]
        tomograms which are to be used for denoise neural network training.
        User should enter tomogram names as rlnTomoName, separated by colons, e.g. TS_1:TS_2
    number_training_subvolumes: int
        Number of sub-volumes to be extracted per training tomogram
        Corresponds to num_slices in cryoCARE_extract_train_data.py.
        Default is 1200.
        Number of normalisation samples will be 10% of this value.
    subvolume_sidelength: int
        Dimensions (for XYZ, in pixels) of the subvolumes extracted for training
        Default is 72. This number should not be lower than 64.
        Corresponds to patch_shape in cryoCARE_extract_train_data.py
    gpu: Optional[List[int]]
        Specify which GPU to use.
    """
    if not tilt_series_star_file.exists():
        e = 'Could not find tilt series star file'
        console.log(f'ERROR: {e}')
        raise RuntimeError(e)

    global_star = starfile.read(tilt_series_star_file, always_dict=True)['global']

    if not 'rlnTomoReconstructedTomogramHalf1' in global_star.columns:
        e = 'Could not find rlnTomoReconstructedTomogramHalf1 in tilt series star file.'
        console.log(f'ERROR: {e}')
        raise RuntimeError(e)

    training_dir, tomogram_dir, tilt_series_dir = \
        create_denoising_directory_structure(
            output_directory=output_directory,
            training_job=True,
        )

    training_tomograms = parse_training_tomograms(training_tomograms)

    training_tomograms_star = generate_training_tomograms_star(
        global_star=global_star,
        training_tomograms=training_tomograms,
    )

    even_tomos, odd_tomos = find_tomogram_halves(training_tomograms_star)

    console.log('Beginning to train denoise model.')

    train_data_config_json = generate_train_data_config_json(
        even_tomos=even_tomos,
        odd_tomos=odd_tomos,
        training_dir=training_dir,
        number_training_subvolumes=number_training_subvolumes,
        subvolume_dimensions=subvolume_sidelength,
    )

    save_json(
        training_dir=training_dir,
        output_json=train_data_config_json,
        json_prefix=TRAIN_DATA_CONFIG_PREFIX,
    )

    cmd = f"cryoCARE_extract_train_data.py --conf {training_dir}/{TRAIN_DATA_CONFIG_PREFIX}.json"
    subprocess.run(cmd, shell=True, stderr=subprocess.STDOUT)

    train_config_json = generate_train_config_json(
        training_dir=training_dir,
        output_directory=output_directory,
        model_name=MODEL_NAME,
        gpu=gpu,
    )

    save_json(
        training_dir=training_dir,
        output_json=train_config_json,
        json_prefix=TRAIN_CONFIG_PREFIX,
    )

    cmd = f"cryoCARE_train.py --conf {training_dir}/{TRAIN_CONFIG_PREFIX}.json"
    subprocess.run(cmd, shell=True)

    save_tilt_series_stars(
        global_star=global_star,
        tilt_series_dir=tilt_series_dir,
    )

    save_global_star(
        global_star=global_star,
        output_directory=output_directory,
    )

    if Path(f'{output_directory}/{MODEL_NAME}.tar.gz').exists():
        console.log(f'Finished training denoise model.')
        console.log(f'Denoising model can be found in {output_directory}/{MODEL_NAME}.tar.gz')
    else:
        e = f'Could not find denoise model ({MODEL_NAME}.tar.gz) in {output_directory}. Training has likely failed.'
        raise RuntimeError(e)

    console.save_html(str(output_directory / 'log.html'), clear=False)
    console.save_text(str(output_directory / 'log.txt'), clear=False)
