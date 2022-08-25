import os
import argparse
import sys

if (sys.version_info < (3, 0)):
    raise Exception('This script supports Python 3 or above.')
    print(end="")  # This script requires Python 3. A Syntax error here means you are running it in Python 2.

try:
    import torch
except ImportError:
    print("PYTHON ERROR: The required python module 'torch' was not found.")
    exit(1)

try:
    import numpy as np
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

features_tensor = torch.Tensor(features)
images_tensor = torch.unsqueeze(torch.Tensor(images), 1)
score = model(images_tensor, features_tensor).detach().cpu().numpy()

for i in range(count):
    print(score[i, 0], end=" ")

