from pathlib import Path
from typing import Optional, Tuple, List

import starfile
import rich
import typer
import subprocess

from ._utils import (
    create_denoising_directory_structure,
    find_tomogram_halves,
    generate_predict_json,
    save_json,
    rename_predicted_tomograms,
    save_tilt_series_stars,
    add_denoised_tomo_to_global_star,
    save_global_star,
)
from .constants import PREDICT_CONFIG_PREFIX, EVEN_SUFFIX
from .._cli import cli
from ..._utils.relion import relion_pipeline_job

console = rich.console.Console(record=True)


@cli.command(name='cryoCARE:predict')
@relion_pipeline_job
def cryoCARE_predict(
        tilt_series_star_file: Path = typer.Option(...),
        output_directory: Path = typer.Option(...),
        model_file: Path = typer.Option(...),
        n_tiles: Optional[Tuple[int, int, int]] = typer.Option((1, 1, 1)),
        tomogram_name: Optional[str] = typer.Option(None),
        gpu: Optional[List[int]] = typer.Option(None)

):
    """Denoise tomograms using cryoCARE (>=v0.2.0).
    
    Requires that two tomograms have been generated using the same sample.
    These can be generated via taking odd/even frames during motion correction
    (optimal) or by taking odd/even tilts during tomogram reconstruction.

    The location of these tomograms should be specified in the
    global star file for all tilt series with the headers:
    `rlnTomoReconstructedTomogramHalf1` and `rlnTomoReconstructedTomogramHalf2`.

    Parameters
    ----------
    tilt_series_star_file: Path
        RELION tilt-series STAR file.
    output_directory: Path
        directory in which results will be stored.
    model_file: Path
        user should provide the path to the model.tar.gz produced by a
        previous cryoCARE denoise job.
    n_tiles: Tuple[int, int, int]
        Initial number of tiles per dimension during prediction step.
        Should get increased if the tiles do not fit on the GPU.
        However, sometimes this parameter is a source of errors causing out of
        memory problems and tomogram dimensions to be wrong. We have found
        other useful values for this include 4,4,2 and 2,4,4 (default = 1,1,1).
        e.g. `--n-tiles 4 4 2`
    tomogram_name: Optional[str]
         Name of tomogram to generate. Use the name in rlnTomoName to specify tomogram.
    gpu: Optional[List[int]]
        specify which GPU(s) to use. To use multiple GPUs provide the flag multiple
        times with a different GPU after each.
        e.g. `--gpu 0 --gpu 1`
    """
    if not tilt_series_star_file.exists():
        e = 'Could not find tilt series star file'
        console.log(f'ERROR: {e}')
        raise RuntimeError(e)

    global_star = starfile.read(tilt_series_star_file, always_dict=True)['global']

    if not model_file.exists():
        e = f'Could not find denoise model'
        console.log(f'ERROR: {e}')
        raise RuntimeError(e)

    if not 'rlnTomoReconstructedTomogramHalf1' in global_star.columns:
        e = 'Could not find rlnTomoReconstructedTomogramHalf1 in tilt series star file.'
        console.log(f'ERROR: {e}')
        raise RuntimeError(e)
    training_dir, tomogram_dir, tilt_series_dir = \
        create_denoising_directory_structure(
            output_directory=output_directory,
            training_job=False,
        )
    even_tomos, odd_tomos = find_tomogram_halves(global_star, tomogram_name)
    predict_json = generate_predict_json(
        even_tomos=even_tomos,
        odd_tomos=odd_tomos,
        training_dir=training_dir,
        model_name=model_file,
        output_directory=output_directory,
        n_tiles=n_tiles,
        gpu=gpu,
    )
    save_json(
        training_dir=training_dir,
        output_json=predict_json,
        json_prefix=PREDICT_CONFIG_PREFIX,
    )

    console.log('Generating denoised tomograms')
    cmd = f"cryoCARE_predict.py --conf {training_dir}/{PREDICT_CONFIG_PREFIX}.json"
    subprocess.run(cmd, shell=True)
    rename_predicted_tomograms(
        even_tomos=even_tomos,
        tomogram_dir=tomogram_dir,
        even_suffix=EVEN_SUFFIX,
    )
    console.log('Denoised tomograms successfully generated, finalising metadata')
    save_tilt_series_stars(
        global_star=global_star,
        tilt_series_dir=tilt_series_dir,
    )
    add_denoised_tomo_to_global_star(
        global_star=global_star,
        tomogram_dir=tomogram_dir,
        output_directory=output_directory,
        tomo_name=tomogram_name
    )
    save_global_star(
        global_star=global_star,
        output_directory=output_directory,
    )
    console.save_html(str(output_directory / 'log.html'), clear=False)
    console.save_text(str(output_directory / 'log.txt'), clear=False)
