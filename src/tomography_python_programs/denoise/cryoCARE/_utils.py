from pathlib import Path
from typing import Tuple, List, Dict, Optional

import pandas as pd
import starfile
import typer
import json
import shutil

def create_denoising_directory_structure(
        output_directory: Path,
        training_job: bool,
) -> Tuple[Path, Path, Path]:
    """
    Creates directory structure for denoise jobs. Doe not create tomogram directory if the job is for training a
    denoise model as no tomograms are generated in this step.
    """
    training_dir = output_directory / 'external' / 'training' 
    training_dir.mkdir(parents=True, exist_ok=True)
    tomogram_dir = output_directory / 'tomograms'
    if not training_job:
        tomogram_dir.mkdir(parents=True, exist_ok=True)
    tilt_series_dir = output_directory / 'tilt_series'
    tilt_series_dir.mkdir(parents=True, exist_ok=True)
    return training_dir, tomogram_dir, tilt_series_dir
    
def parse_training_tomograms(
        training_tomograms: str
) -> List:
    """
    Reads the string given to the CLI to ascertain which tomograms to train on. String should
    be a list of : separated tomograms (name from rlnTomoName)
    """
    training_tomograms = training_tomograms.strip().split(':')
    return training_tomograms

def generate_training_tomograms_star(
        global_star: pd.DataFrame,
        training_tomograms: List,
) -> pd.DataFrame:
    """
    Generates a pandas dataframe of the tomograms the user has selected for training in global star format
    """
    training_tomograms_idx = pd.DataFrame(global_star.rlnTomoName.tolist()).isin(training_tomograms).values
    if not any(training_tomograms_idx):
        e = f"Could not user specified training tomograms ({', '.join(str(x) for x in training_tomograms)}) in tomogram star file"
        raise RuntimeError(e)
    training_tomograms_star = global_star[training_tomograms_idx]
    return training_tomograms_star
    
def find_tomogram_halves(
        training_tomograms_star: pd.DataFrame,
        tomo_name: Optional[str] = None,
) -> Tuple[List, List]:
    """
    Returns lists (even and odd) of the location of the the tomograms the user wishes to train on.
    """
    if tomo_name is None:
        return training_tomograms_star['rlnTomoReconstructedTomogramHalf1'].values.tolist(), training_tomograms_star['rlnTomoReconstructedTomogramHalf2'].values.tolist() 
    else:
        training_tomograms_idx = training_tomograms_star.index[(training_tomograms_star['rlnTomoName'] == tomo_name)]
        if len(training_tomograms_idx) == 0:
            e = f"Could not user specified tomogram rlnTomoName {tomo_name} in tomogram star file"
            raise RuntimeError(e)
        training_tomograms_star = training_tomograms_star.loc[training_tomograms_idx]
        return training_tomograms_star['rlnTomoReconstructedTomogramHalf1'].values.tolist(), training_tomograms_star['rlnTomoReconstructedTomogramHalf2'].values.tolist()
	
        
def generate_train_data_config_json(
        even_tomos: List,
        odd_tomos: List,
        training_dir: Path,
        number_training_subvolumes: int,
        subvolume_dimensions: int,
) -> Dict:
    """
    Creates a Dict which can be saved as a json file for train_data_config.json file
    """
    number_normalisation_subvolumes = round(number_training_subvolumes * 0.1)
    train_data_config_json = json.loads(f'{{"even": {json.dumps(even_tomos)}, "odd": {json.dumps(odd_tomos)}, "patch_shape": [{subvolume_dimensions}, {subvolume_dimensions}, {subvolume_dimensions}], \
    "num_slices": {number_training_subvolumes}, "split": 0.9, "tilt_axis": "Y", "n_normalization_samples": {number_normalisation_subvolumes}, "path": "{training_dir}", "overwrite": "True"}}')
    return train_data_config_json

def generate_train_config_json(
        training_dir: Path,
        output_directory: Path,
	model_name: str,
        gpu: Optional[List[int]] = None,
) -> Dict:
    """
    Creates a Dict which can be saved as a json file for train_config.json file
    """
    if gpu is None:
        gpu_json =  f''
    else:
        gpu_json = f', "gpu_id": {gpu}'
    train_config_json = json.loads(f'{{"train_data": "{training_dir}", "epochs": 100, "steps_per_epoch": 200, "batch_size": 16, "unet_kern_size": 3, \
    "unet_n_depth": 3, "unet_n_first": 16, "learning_rate": 0.0004, "model_name": "{model_name}", "path": "{output_directory}", "overwrite": "True"{gpu_json}}}')
    return train_config_json

def generate_predict_json(
        even_tomos: List,
        odd_tomos: List,
	training_dir: Path,
	model_name: Path,
        output_directory: Path,
        n_tiles: Tuple[int,int,int],
        gpu: Optional[List[int]] = None,
) -> Dict:
    """
    Creates a Dict which can be saved as a json file for predict_config.json file
    """
    if gpu is None:
        gpu_json =  f''
    else:
        gpu_json = f', "gpu_id": {gpu}'
    predict_json = json.loads(f'{{"path": "{model_name}", "even": {json.dumps(even_tomos)}, \
    "odd": {json.dumps(odd_tomos)}, "n_tiles": {list(n_tiles)}, "output": "{output_directory / "tomograms"}", "overwrite": "True"{gpu_json}}}')
    return predict_json

def save_json(
        training_dir: Path,
        output_json: Dict,
	json_prefix: str,
):
    """
    Saves json file in output directory with desired file name (prefix).
    """
    with open(f'{training_dir}/{json_prefix}.json', 'w') as outfile:
        json.dump(output_json, outfile, indent=4) 
	
def save_tilt_series_stars(
        global_star: pd.DataFrame,
        tilt_series_dir: Path,
):
    """
    Saves tilt series star files in output directory.
    """
    for idx,row in global_star.iterrows():
        shutil.copyfile(f"{row['rlnTomoTiltSeriesStarFile']}", f'{tilt_series_dir}/{row["rlnTomoName"]}.star')
    global_star['rlnTomoTiltSeriesStarFile'] = global_star.apply(lambda x: f'{tilt_series_dir}/{x["rlnTomoName"]}.star', axis=1)
    
def add_denoised_tomo_to_global_star(
        global_star: pd.DataFrame,
        tomogram_dir: Path,
        output_directory: Path,
        tomo_name: Optional[str] = None
):
    """
    Adds location of the denoise tomogram to the global star file.
    """
    if tomo_name == None:
        global_star['rlnTomoReconstructedTomogramDenoised'] = global_star.apply(lambda x: f'{tomogram_dir}/rec_{x["rlnTomoName"]}.mrc', axis=1)
        return global_star
    else:
        if 'rlnTomoReconstructedTomogramDenoised' not in global_star.columns:
            global_star['rlnTomoReconstructedTomogramDenoised'] = 'N/A'  
        global_star.loc[global_star.rlnTomoName == tomo_name, 'rlnTomoReconstructedTomogramDenoised'] = f"{tomogram_dir}/rec_{global_star.loc[global_star.rlnTomoName == tomo_name, 'rlnTomoName'].values[0]}.mrc"
        return global_star
    
def save_global_star(
        global_star: pd.DataFrame,
        output_directory: Path,
):
    """
    Saves global star file (tomograms.star) in output directory.
    """
    starfile.write({'global': global_star}, f'{output_directory}/tomograms.star')
    
def rename_predicted_tomograms(
    even_tomos: List,
    tomogram_dir: Path,
    even_suffix: str,
):
    """
    Gives denoised tomograms as cryoCARE likes to name them after the even tomograms.
    """
    even_tomos = [Path(tomo) for tomo in even_tomos]
    even_tomos = [Path(f"{tomogram_dir}/{tomo.name}") for tomo in even_tomos]
    for tomo in even_tomos:
        if not tomo.exists():
            raise RuntimeError(f'{tomo} not found and it should have been generated. If OOM error, try increasing number of tiles (e.g. 8,8,8).') 
    [tomo.rename(Path(f"{tomogram_dir}/{tomo.stem.replace(even_suffix,'')}{tomo.suffix}")) for tomo in even_tomos]
