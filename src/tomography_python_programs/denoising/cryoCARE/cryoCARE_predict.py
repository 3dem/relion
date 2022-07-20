from pathlib import Path
from typing import Optional, Tuple

import pandas as pd
import starfile
import rich
import typer
import subprocess
from rich.progress import track

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
from ...utils.relion import relion_pipeline_job

console = rich.console.Console(record=True)

@cli.command(name='cryoCARE:predict')
@relion_pipeline_job
def cryoCARE_predict(
    tilt_series_star_file: Path = typer.Option(...),
    output_directory: Path = typer.Option(...),
    model_name: Path = typer.Option(...),
    n_tiles: Optional[Tuple[int,int,int]] = typer.Option((1,1,1))

):
    """Generates denoised tomograms using cryoCARE from a previously trained denoising model (.tar.gz)
    (Euan Pyle version, https://github.com/EuanPyle/cryoCARE_mpido) branched from Thorsten Wagner 
    version, https://github.com/thorstenwagner/cryoCARE_mpido)
    
    Requires that two tomograms have been generated using the same sample. These can be generated via taking odd/even 
    frames during Motion Correction (optimal) or by taking odd/even tilts during tomogram reconstruction.
    The location of these tomograms should be specified in the global star file for all tilt series with the headers: 
    
     
    rlnTomoReconstructedTomogramHalf1
    
  rlnTomoReconstructedTomogramHalf2

    Parameters
    ----------
    
    tilt_series_star_file: RELION tilt-series STAR file.
    
    output_directory: directory in which results will be stored.
    	
    model_name: user should provide the path to the model.tar.gz produced by a previous cryoCARE denoising job. 
	
    n-tiles (optional): Initial tiles per dimension during prediction step. Should get increased 
        if the tiles do not fit on the GPU. However, sometimes this parameter is a source of errors causing out of memory 
        problems and tomogram dimensions to be wrong. We have found other useful values for this include 4,4,2 and 2,4,4 
        (default = 1,1,1). Enter as: `--n-tiles 4 4 2`
        
    Returns
    -------
    Denoised tomograms.
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
            training_job=False,
        )
    
    even_tomos, odd_tomos = find_tomogram_halves(global_star)

    predict_json = generate_predict_json(
        even_tomos=even_tomos,
        odd_tomos=odd_tomos,
	training_dir=training_dir,
	model_name=model_name,
        output_directory=output_directory,
        n_tiles=n_tiles,
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
    )
    
    save_global_star(
        global_star=global_star,
        output_directory=output_directory,
    )    
    
    console.save_html(str(output_directory / 'log.html'), clear=False)
    console.save_text(str(output_directory / 'log.txt'), clear=False)
