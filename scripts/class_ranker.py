import os
import argparse
import sys
from packaging import version

try:
    import torch
    if version.parse(torch.__version__) < version.parse("1.0.1"):
        print("PYTHON ERROR: Module 'torch' is too old (must be >= 1.0.1).")
        exit(1)
except ImportError:
    print("PYTHON ERROR: The required python module 'torch' was not found.")
    exit(1)

try:
    import numpy as np
    if version.parse(np.__version__) < version.parse("1.0.0"):
        print("PYTHON ERROR: Module 'numpy' is too old (must be >= 1.0.0).")
        exit(1)
except ImportError:
    print("PYTHON ERROR: The required python module 'numpy' was not found.")
    exit(1)

parser = argparse.ArgumentParser()

parser.add_argument('model_path', type=str)
parser.add_argument('project_dir', type=str)
args = parser.parse_args()

project_dir = args.project_dir
model_fn = args.model_path
feature_fn = os.path.join(project_dir, "features.npy")
images_fn = os.path.join(project_dir, "images.npy")

model = torch.jit.load(model_fn)
features = np.load(feature_fn)
images = np.load(images_fn)

count = features.shape[0]

model = model.to("cpu")

features_tensor = torch.Tensor(features)
images_tensor = torch.unsqueeze(torch.Tensor(images), 1)
score = model(images_tensor, features_tensor).detach().cpu().numpy()

for i in range(count):
    print(score[i, 0], end=" ")

