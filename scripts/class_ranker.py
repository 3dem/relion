import os
import argparse
import sys

try:
    import torch
except ImportError:
    print("PYTHON ERROR: The required python module 'torch' was not found.")
    exit(1)
try:
    import numpy
except ImportError:
    print("PYTHON ERROR: The required python module 'numpy' was not found.")
    exit(1)

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


