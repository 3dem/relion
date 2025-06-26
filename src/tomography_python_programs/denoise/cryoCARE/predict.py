from pathlib import Path
from typing import Optional, Tuple, List
from rich.console import Console

import starfile
import typer
import subprocess
import os

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

console = Console(record=True)
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '1'  

@cli.command(name='cryoCARE:predict')
@relion_pipeline_job
def cryoCARE_predict(
        tomogram_star_file: Path = typer.Option(...),
        output_directory: Path = typer.Option(...),
        model_file: Path = typer.Option(...),
        cryocare_path: Optional[Path] = typer.Option(""),
        n_tiles: Optional[Tuple[int, int, int]] = typer.Option((1, 1, 1)),
        tomogram_name: Optional[str] = typer.Option(None),
        gpu: Optional[List[int]] = typer.Option(None)

):
    """Denoise tomograms using cryoCARE (>=v0.2.1).
    
    Requires that two tomograms have been generated using the same sample.
    These can be generated via taking odd/even frames during motion correction
    (optimal) or by taking odd/even tilts during tomogram reconstruction.

    The location of these tomograms should be specified in the
    global star file for all tomograms with the headers:
    `rlnTomoReconstructedTomogramHalf1` and `rlnTomoReconstructedTomogramHalf2`.

    Parameters
    ----------
    tomogram_star_file: Path
        RELION tomogram STAR file from a RELION Reconstruct Tomogram job.
    output_directory: Path
        directory in which results will be stored.
    model_file: Path
        user should provide the path to the model.tar.gz produced by a
        previous cryoCARE denoise job.
    cryocare_path: Optional[Path]
        directory where the cryoCARE executables can be found. This can be
        left empty if the cryoCARE executables are already in the PATH
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
    if not tomogram_star_file.exists():
        e = 'Could not find tomogram star file'
        console.log(f'ERROR: {e}')
        raise RuntimeError(e)

    cryocare_path = str(cryocare_path)
    if cryocare_path != "":
        predict_executable = os.path.join(cryocare_path,  "cryoCARE_predict.py")
        if not os.path.isfile(predict_executable):
            e = 'Could not find cryoCARE_predict.py executable in the input path: ' + cryocare_path
            console.log(f'ERROR: {e}')
            raise RuntimeError(e)
    else:
        predict_executable = "cryoCARE_predict.py"

    global_star = starfile.read(tomogram_star_file, always_dict=True, parse_as_string=['rlnTomoName'])['global']

    if not model_file.exists():
        e = f'Could not find denoise model'
        console.log(f'ERROR: {e}')
        raise RuntimeError(e)

    if not 'rlnTomoReconstructedTomogramHalf1' in global_star.columns:
        e = 'Could not find rlnTomoReconstructedTomogramHalf1 in tomogram star file.'
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
    cmd = f"{predict_executable} --conf {training_dir}/{PREDICT_CONFIG_PREFIX}.json"
    subprocess.run(cmd, shell=True, stderr=subprocess.STDOUT)
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
