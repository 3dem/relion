import os
import argparse
import sys


def import_error(module_name):
    print(
        "PYTHON ERROR: The required python module '",
        module_name,
        "' was not found.",
        sep=""
    )
    exit(1)


if sys.version_info >= (3, 5):
    import importlib
    if not importlib.util.find_spec("torch"):
        import_error('torch')
    if not importlib.util.find_spec("numpy"):
        import_error('numpy')
else:
    import imp
    try:
        imp.find_module('torch')
    except ImportError:
        import_error('torch')
    try:
        imp.find_module('numpy')
    except ImportError:
        import_error('numpy')

import torch
import numpy as np

parser = argparse.ArgumentParser()

parser.add_argument('model_path', type=str)
parser.add_argument('project_dir', type=str)
parser.add_argument('--features', type=str, default=None)
parser.add_argument('--images', type=str, default="images.npy")
args = parser.parse_args()

project_dir = args.project_dir
model_fn = args.model_path
feature_fn = os.path.join(project_dir, "features.npy")
images_fn = os.path.join(project_dir, "images.npy")

model = torch.jit.load(model_fn)
features = np.load(feature_fn)
images = np.load(images_fn)

model = model.to("cpu")

features_tensor = torch.Tensor(features).reshape(1, features.shape[0])
images_tensor = torch.Tensor(images).reshape(1, 1, images.shape[0], images.shape[1])
score = model(images_tensor, features_tensor).detach().cpu().numpy()

print(score[0, 0])


