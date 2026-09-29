#!/public/EM/anaconda3/envs/relion-5.0/bin/python

import argparse
import math
import sys, os
from pathlib import Path

import numpy as np
import pandas as pd
import starfile
import torch
import torch.nn as nn
import torch.nn.functional as F

import mrcfile
import matplotlib
matplotlib.use("Agg")  # non-GUI backend, no Qt/X11 needed
import matplotlib.pyplot as plt

from skimage.morphology import remove_small_objects, skeletonize
from skimage.measure import regionprops, label as label_connected_component
from skan import Skeleton, summarize
from scipy.spatial.distance import cdist
import cv2

from tqdm import tqdm
import time

###############################################################################
#                               U-Net Model
###############################################################################
class DownCarbon(nn.Module):
    """Downscaling with maxpool then double conv"""

    def __init__(self, in_channels, out_channels):
        super().__init__()
        self.maxpool_conv = nn.Sequential(
            nn.MaxPool2d(2),
            DoubleConv(in_channels, out_channels)
        )

    def forward(self, x):
        return self.maxpool_conv(x)


class UpCarbon(nn.Module):
    """Upscaling then double conv"""

    def __init__(self, in_channels, out_channels):
        super().__init__()

        # 1) spatial upsample
        # 2) reduce channels 
        self.up = nn.Sequential(
            nn.Upsample(scale_factor=2, mode='bilinear', align_corners=True),
            nn.Conv2d(in_channels, in_channels // 2, kernel_size=1),
        )
        # after cat: (in_channels//2 + in_channels//2) == in_channels
        self.conv = DoubleConv(in_channels, out_channels)

    def forward(self, x1, x2):
        x1 = self.up(x1)
        # input is CHW
        diffY = x2.size()[2] - x1.size()[2]
        diffX = x2.size()[3] - x1.size()[3]

        x1 = F.pad(x1, [diffX // 2, diffX - diffX // 2,
                        diffY // 2, diffY - diffY // 2])
        # if you have padding issues, see
        # https://github.com/HaiyongJiang/U-Net-Pytorch-Unstructured-Buggy/commit/0e854509c2cea854e247a9c615f175f76fbb2e3a
        # https://github.com/xiaopeng-liao/Pytorch-UNet/commit/8ebac70e633bac59fc22bb5195e513d5832fb3bd
        x = torch.cat([x2, x1], dim=1)
        return self.conv(x)


class DoubleConv(nn.Module):
    """(convolution => [InstanceNorm] => GELU) * 2"""
    def __init__(self, in_channels, out_channels, mid_channels=None):
        super().__init__()
        if not mid_channels:
            mid_channels = out_channels
        self.double_conv = nn.Sequential(
            nn.Conv2d(in_channels, mid_channels, kernel_size=3, padding=1, bias=False),
            nn.InstanceNorm2d(mid_channels),
            nn.GELU(),
            nn.Conv2d(mid_channels, out_channels, kernel_size=3, padding=1, bias=False),
            nn.InstanceNorm2d(out_channels),
            nn.GELU()
        )

    def forward(self, x):
        return self.double_conv(x)


class Down(nn.Module):
    def __init__(self, in_channels, out_channels):
        super().__init__()
        self.maxpool_conv = nn.Sequential(
            nn.MaxPool2d(2),
            DoubleConv(in_channels, out_channels)
        )

    def forward(self, x):
        return self.maxpool_conv(x)


class Up(nn.Module):
    def __init__(self, in_channels, out_channels, bilinear=True):
        super().__init__()
        if bilinear:
            self.up = nn.Upsample(scale_factor=2, mode='bilinear', align_corners=True)
            self.conv = DoubleConv(in_channels, out_channels, in_channels // 2)
        else:
            self.up = nn.ConvTranspose2d(in_channels, in_channels // 2, kernel_size=2, stride=2)
            self.conv = DoubleConv(in_channels, out_channels)

    def forward(self, x1, x2):
        x1 = self.up(x1)
        diffY = x2.size()[2] - x1.size()[2]
        diffX = x2.size()[3] - x1.size()[3]

        x1 = F.pad(
            x1, 
            [diffX // 2, diffX - diffX // 2, diffY // 2, diffY - diffY // 2]
        )
        x = torch.cat([x2, x1], dim=1)
        return self.conv(x)


class OutConv(nn.Module):
    def __init__(self, in_channels, out_channels):
        super().__init__()
        self.conv = nn.Conv2d(in_channels, out_channels, kernel_size=1)

    def forward(self, x):
        return self.conv(x)


class UNet(nn.Module):
    def __init__(self, n_channels, n_classes, embedding_dim=32, bilinear=False):
        super().__init__()
        self.n_channels = n_channels
        self.n_classes = n_classes
        self.bilinear = bilinear

        self.base_channel_sz = 128

        self.inc = DoubleConv(n_channels, self.base_channel_sz)
        self.down1 = Down(self.base_channel_sz, 2 * self.base_channel_sz)
        self.down2 = Down(2 * self.base_channel_sz, 4 * self.base_channel_sz)
        self.down3 = Down(4 * self.base_channel_sz, 8 * self.base_channel_sz)
        factor = 2 if bilinear else 1
        self.down4 = Down(8 * self.base_channel_sz, 16 * self.base_channel_sz // factor)

        self.up1 = Up(16 * self.base_channel_sz, 8 * self.base_channel_sz // factor, bilinear)
        self.up2 = Up(8 * self.base_channel_sz, 4 * self.base_channel_sz // factor, bilinear)
        self.up3 = Up(4 * self.base_channel_sz, 2 * self.base_channel_sz // factor, bilinear)
        self.up4 = Up(2 * self.base_channel_sz, self.base_channel_sz, bilinear)

        self.outc = OutConv(self.base_channel_sz, n_classes)
        self.embedding_head = nn.Sequential(
            nn.Conv2d(self.base_channel_sz, embedding_dim, kernel_size=1),
            nn.Sigmoid()
        )

    def forward(self, x):
        x1 = self.inc(x)
        x2 = self.down1(x1)
        x3 = self.down2(x2)
        x4 = self.down3(x3)
        x5 = self.down4(x4)

        x = self.up1(x5, x4)
        x = self.up2(x, x3)
        x = self.up3(x, x2)
        x = self.up4(x, x1)

        seg_logits = self.outc(x)
        embedding = self.embedding_head(x)
        return seg_logits, embedding

    @staticmethod
    def load_from_checkpoint(model_path: str, n_channels: int, n_classes: int):
        model = UNet(n_channels=n_channels, n_classes=n_classes)
        lightning_dict = torch.load(model_path, map_location="cpu")
        state_dict = lightning_dict['state_dict']
        try:
            model.load_state_dict(state_dict, strict=False)
        except:
            model.load_state_dict(state_dict)
        return model


class UNetCarbon(nn.Module):
    def __init__(self, n_channels, n_classes, dropout=0.25):
        super(UNetCarbon, self).__init__()
        self.n_channels = n_channels
        self.n_classes = n_classes
        self.base_channel_sz = 64

        self.inc = (DoubleConv(n_channels, self.base_channel_sz))
        self.down1 = (DownCarbon(self.base_channel_sz, 2 * self.base_channel_sz))
        self.down2 = (DownCarbon(2 * self.base_channel_sz, 4 * self.base_channel_sz))
        self.down3 = (DownCarbon(4 * self.base_channel_sz, 8 * self.base_channel_sz))
        self.down4 = (DownCarbon(8 * self.base_channel_sz, 16 * self.base_channel_sz ))
        self.up1 = (UpCarbon(16 * self.base_channel_sz, 8 * self.base_channel_sz))
        self.up2 = (UpCarbon(8 * self.base_channel_sz, 4 * self.base_channel_sz))
        self.up3 = (UpCarbon(4 * self.base_channel_sz, 2 * self.base_channel_sz))
        self.up4 = (UpCarbon(2 * self.base_channel_sz, self.base_channel_sz))
        self.outc = (OutConv(self.base_channel_sz, n_classes))


    def forward(self, x):
        x1 = self.inc(x)
        x2 = self.down1(x1)
        x3 = self.down2(x2)
        x4 = self.down3(x3)
        x5 = self.down4(x4)

        x = self.up1(x5, x4)
        x = self.up2(x, x3)
        x = self.up3(x, x2)
        x = self.up4(x, x1)
        logits = self.outc(x)

        return logits

    @staticmethod
    def load_from_checkpoint(model_path: str, n_channels: int, n_classes: int):
        model = UNetCarbon(n_channels=n_channels, n_classes=n_classes)
        lightning_dict = torch.load(model_path, map_location="cpu")
        state_dict = lightning_dict['state_dict']
        try:
            model.load_state_dict(state_dict, strict=False)
        except:
            model.load_state_dict(state_dict)
        return model

   

###############################################################################
#                        MRC + instance segmentation
###############################################################################
def load_mrc(mrc_fn: str, multiply_global_origin: bool = True):
    with mrcfile.open(mrc_fn, "r") as mrc_file_handle:
        voxel_size = float(mrc_file_handle.voxel_size.x)
        if voxel_size <= 0:
            raise RuntimeError(f"MRC file {mrc_fn} does not have a valid header.")

        c = mrc_file_handle.header["mapc"]
        r = mrc_file_handle.header["mapr"]
        s = mrc_file_handle.header["maps"]

        global_origin = mrc_file_handle.header["origin"]
        global_origin = np.array([global_origin.x, global_origin.y, global_origin.z])
        global_origin[0] += mrc_file_handle.header["nxstart"]
        global_origin[1] += mrc_file_handle.header["nystart"]
        global_origin[2] += mrc_file_handle.header["nzstart"]

        if multiply_global_origin:
            global_origin *= mrc_file_handle.voxel_size.x

        raw_data = mrc_file_handle.data
        if c == 1 and r == 2 and s == 3:
            grid = raw_data
        elif c == 3 and r == 2 and s == 1:
            grid = np.moveaxis(raw_data, [0, 1, 2], [2, 1, 0])
        elif c == 2 and r == 1 and s == 3:
            grid = np.moveaxis(raw_data, [1, 2, 0], [2, 1, 0])
        else:
            raise RuntimeError("MRC file axis arrangement not supported!")
    return grid, voxel_size, global_origin


###############################################################################
#        Sample coordinates from the final pruned label, including IDs
###############################################################################
def sample_coordinates_with_ids(label_image, step=50):

    labels = np.unique(label_image)
    labels = labels[labels != 0]
    labels.sort()  # ensure ascending
    results = []

    filament_counter = 1
    for lbl in labels:
        bin_mask = (label_image == lbl)
        sk = Skeleton(bin_mask)
        if sk.n_paths == 0:
            continue

        # For each path in this label, sample coordinates
        for p_idx in range(sk.n_paths):
            coords = sk.path_coordinates(p_idx).astype(int)
            if len(coords) == 0:
                continue
            for i in range(0, len(coords), step):
                r, c = coords[i]
                results.append([c, r, filament_counter])
            # ensure final coordinate is included
            if (len(coords) - 1) % step != 0:
                r, c = coords[-1]
                results.append([c, r, filament_counter])

        filament_counter += 1

    return np.array(results, dtype=int)


###############################################################################
#                            Visualization Helpers
###############################################################################
import matplotlib.colors as mcolors

def color_nonrandom(label_image):
    H, W = label_image.shape
    colorized = np.full((H, W, 3), 127, dtype=np.uint8)

    labels = np.unique(label_image)
    labels = labels[labels != 0]
    i = 0
    for lbl in labels:
        if i%6==0:
            r= 255
            g= 0
            b =0
        elif i%6==1:
            r= 0
            g= 255
            b =0
        elif i%6==2:
            r= 0
            g= 0
            b =255
        elif i%6==3:
            r= 255
            g= 255
            b =0
        elif i%6==4:
            r= 0
            g= 255
            b =255
        elif i%6==5:
            r= 255
            g= 0
            b =255
        mask = (label_image == lbl)
        colorized[mask] = (r, g, b)
        i=i+1
    return colorized

def color_foreground(grey_image, threshold):
    grey_image_uint8 = (grey_image * 255).astype(np.uint8)
    colorized = np.stack([grey_image_uint8] * 3, axis=-1)
    mask = (grey_image > threshold)
    colorized[mask] = (255, 0, 0)
    return colorized


###############################################################################
#  Sjors go at pruning the raw skeleton image
###############################################################################

def trim_end(branch_coords, trim_length=10, trim_start=False):
    """Trims the start or end of a branch by `trim_length` pixels."""
    if len(branch_coords) > trim_length:
        return branch_coords[trim_length:] if trim_start else branch_coords[:-trim_length]
    return np.array([])  # Remove short branches completely

def angular_diff(a, b, period=180):
    """Minimal circular difference between angles a and b on [0, period)."""
    raw = abs(a - b)
    return min(raw, period - raw)

def segment_orientation(p0, p1):
    """
    Compute the orientation (in degrees) and fold it into [0, 180).
    p0, p1 are (row, col) integer tuples or arrays.
    """
    # Note: row = y, col = x
    dy = p1[0] - p0[0]
    dx = p1[1] - p0[1]
    # conventions of y-axis is reversed in RELION, so use -dy!
    ang = np.degrees(np.arctan2(-dy, dx))  # gives [-180, 180)

    # fold into [0, 180):
    if ang < 0:
        ang += 180
    elif ang >= 180:
        ang -= 180

    return ang

def split_branch_on_psi_mismatch(branch_coords, psi_image, threshold=20):
    
    segments = []
    current_seg = [tuple(branch_coords[0])]
    prev_pt = branch_coords[0]
    prev_ang = psi_image[prev_pt[0], prev_pt[1]]
    
    # step through each next point
    for i in range(1, len(branch_coords)):
        pt = branch_coords[i]
        psi_ang = psi_image[pt[0], pt[1]]

        # split if jump in neighbour psi image values > threshold
        do_split = False
        if angular_diff(prev_ang, psi_ang) > threshold:
            do_split = True
        # compute seg_ang between i-3 and i
        elif i >= 3:
            p0 = branch_coords[i - 3]
            seg_ang = segment_orientation(p0, pt)
            # split if mismatch > threshold
            if angular_diff(seg_ang, psi_ang) > threshold:
                do_split = True

        if do_split:
            segments.append(np.array(current_seg))
            current_seg = []
        
        current_seg.append(tuple(pt))
        prev_ang = psi_ang
        
    # append last segment
    if current_seg:
        segments.append(np.array(current_seg))

    return segments   

def prune_skeleton(skeleton, psi_image, min_length=100, min_distance=25, step=50, psi_threshold=45):
    """Prune skeleton branches that too short."""
    # Most of the compute goes into the Skeleton call; rest of function is cheap
    skel = Skeleton(skeleton)
    df = summarize(skel, separator='-')

    trim_length = int(min_distance / 2)
    labeled_skeleton = np.zeros_like(skeleton, dtype=np.int32)
    accepted_branches = []  # Stores (branch coordinates, start, end)
    sampled_coords = []

    # Pass 0: split branches with jumps in psi-angle values
    all_splits = []
    for idx, row in df.iterrows():
        coords = skel.path_coordinates(idx)
        pieces = split_branch_on_psi_mismatch(coords, psi_image, threshold=psi_threshold)
        all_splits.extend(pieces)

    # Pass 1: Collect branches and compute pairwise distances
    branch_data = []  # Store branch info for pass 2
    endpoints = []  # Store (start, end) tuples
    for branch_coords in all_splits:

        # Compute Euclidean length of the branch
        length = np.sum(np.sqrt(np.sum(np.diff(branch_coords, axis=0) ** 2, axis=1)))
        if length < min_length:
            continue  # Skip short branches

        start, end = tuple(branch_coords[0]), tuple(branch_coords[-1])
        branch_data.append((branch_coords, start, end))
        endpoints.append(start)
        endpoints.append(end)

    if not branch_data:
        return labeled_skeleton, np.array(sampled_coords, dtype=int) # No valid branches

    # Compute all pairwise distances between endpoints
    endpoints = np.array(endpoints)  # Convert to numpy array
    dists = cdist(endpoints, endpoints)

    # Mark close endpoint pairs for trimming
    close_pairs = np.where((dists < min_distance) )  # Exclude self-distances
    trim_flags = {}  # Dictionary to store which endpoints need trimming

    for i in range(len(close_pairs[0])):
        idx1, idx2 = close_pairs[0][i], close_pairs[1][i]
        branch_idx1, endpoint_type1 = divmod(idx1, 2)  # Get branch index & endpoint type
        branch_idx2, endpoint_type2 = divmod(idx2, 2)

        if branch_idx1 == branch_idx2:
            continue  # Ignore self-pairs

        # Mark both endpoints for trimming
        trim_flags.setdefault((branch_idx1, endpoint_type1), True)
        trim_flags.setdefault((branch_idx2, endpoint_type2), True)

    # **Pass 2: Trim branches based on identified close endpoints**
    final_branches = []

    for i, (branch_coords, start, end) in enumerate(branch_data):
        # Trim start or end if flagged
        if (i, 0) in trim_flags:
            branch_coords = trim_end(branch_coords, trim_length, trim_start=True)
        if (i, 1) in trim_flags:
            branch_coords = trim_end(branch_coords, trim_length, trim_start=False)

        # Ensure branch still exists after trimming
        if len(branch_coords) > 0:
            final_branches.append(branch_coords)

    # Assign unique labels to final pruned branches
    for path_id, branch in enumerate(final_branches, start=1):
        for x, y in branch:
            labeled_skeleton[int(x), int(y)] = path_id
            
        for i in range(0, len(branch), step):
            r, c = branch[i]
            sampled_coords.append([c, r, path_id])
            # ensure final coordinate is included
        if (len(branch) - 1) % step != 0:
            r, c = branch[-1]
            sampled_coords.append([c, r, path_id])

    return labeled_skeleton, np.array(sampled_coords, dtype=int)

def resize_and_norm(image, mysize, angpix, do_norm=True):
    ys, xs = image.shape
    if (xs > ys):
        scale = ys/mysize
    else:
        scale = xs/mysize

    xsp = int(xs/scale)
    ysp = int(ys/scale)
    resized_angpix = angpix *scale
    resized = cv2.resize(image.astype(np.float32), (xsp, ysp), interpolation=cv2.INTER_AREA)
    if do_norm:
        max = resized.max()
        min = resized.min()
        if (max > min):
            resized = (resized - min) / (max - min)

    return resized, resized_angpix


def install_model(
        name: str,
        verbose: bool = False,
):
    model_list = {
        "amytracer-v1.0": [
            "https://zenodo.org/records/17949642/files/amytracer-v1.0.ckpt.gz",
            "d83777816d899d64aef593c97888e88f56d49e2b652cdc460de554f59d6f94ed"
        ],
        "amytracer-v2.0": [
            "https://zenodo.org/records/17949642/files/amytracer-v2.0.ckpt.gz",
            "57f0cd566ca1779641691f7daf4a6147ccdb4f509f01c7472861e32bb7dbe14c"
        ],
        "carbonpicker-v1.0": [
            "https://zenodo.org/records/17949642/files/carbonpicker-v1.0.ckpt.gz",
            "e25702cc84850339b46ac5b7d8e19ebdd55ad420c15bda8c9f9a584f59e7fd6b"
        ]
    }

    if name in model_list.keys():
        dest_dir = os.path.join(torch.hub.get_dir(), "checkpoints", "relion_trace_amyloids")
        model_path = os.path.join(dest_dir, f"{name}.ckpt")
        model_path_gz = model_path + ".gz"
        completed_check_path = os.path.join(dest_dir, f"{name}_installed.txt")

        # Download file and install it if not already done
        if not os.path.isfile(completed_check_path):
            if verbose:
                print(f"Installing amyloid picker model ({name})...")
                os.makedirs(dest_dir, exist_ok=True)

            source_url = model_list[name][0]

            if verbose:
                print(f"Downloading model weights from:\n   {source_url}")

            import gzip, shutil
            torch.hub.download_url_to_file(source_url, model_path_gz, hash_prefix=model_list[name][1])
            with gzip.open(model_path_gz, 'rb') as f_in:
                with open(model_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                    os.remove(model_path_gz)

            with open(completed_check_path, "w") as f:
                f.write("Successfully downloaded model")

        if verbose:
            print(f"Amyloid picker model ({name}) successfully installed in {dest_dir}")

    else:
        model_path = name

    return model_path




###############################################################################
#                              MAIN
###############################################################################
def main():

    parser = argparse.ArgumentParser("Skeleton-based filament picking with instance merging, branch pruning, and ID labeling")
    parser.add_argument("-i", "--input", help="Input star file", default="None")
    parser.add_argument("-o", "--output", help="Output star file with coordinates", default="out")
    parser.add_argument("-r", "--radius", type=int, help="Minimum distance between ends of filaments (in A)", default=25)
    parser.add_argument("-l", "--minimum_length", type=int, default=250,
                        help="Remove final segments smaller than this size in Angstroms (default=250)")
    parser.add_argument("-p", "--psi_jump_threshold", default=45, type=float, help="Maximum difference in degrees of psi angles for consecutive points in skeletons (default=45)")
    parser.add_argument("-t", "--threshold", default=0.5, type=float, help="Foreground threshold (default=0.5)")
    parser.add_argument("-c", "--carbon", action='store_true', help="Ignore filaments on carbon")
    parser.add_argument("-ct", "--carbon_threshold", default=0.9, type=float, help="Carbon detection threshold (default=0.9)")
    parser.add_argument("-s", "--scale", type=float, default=1.0,
                        help="Multiply output coordinates by this factor (default=1.0)")
    parser.add_argument("-d", "--device", default="cuda:0", help="Which GPU device to use")
    parser.add_argument("-j", "--threads", type=int, default=1)
    parser.add_argument("-m", "--model_path", help="Path to PyTorch model checkpoint", default="amytracer-v2.0")
    parser.add_argument("-cm", "--carbon_model_path", help="Path to PyTorch model checkpoint for carbon detection", default="carbonpicker-v1.0")
    parser.add_argument("-a", "--abort", help="Abort if this file exists")
    parser.add_argument("-v", "--verb", type=int, help="Verbosity")
    parser.add_argument("--sample_step", type=int, default=100,
                        help="Sample a coordinate every 'sample_step' Angstroms along each path (default=50)")
    parser.add_argument("--plot", action='store_true',
                        help="If set, display/save a collage figure of the pipeline results.")
    args = parser.parse_args()

    # number of CPU threads allowed
    #torch.set_num_threads(args.threads)

    # ------------------- Load model ----------------------
    device = torch.device(args.device)

    # This is for when we call this program without any arguments upon installation of the relion package (in scripts/python_fetch_weights.in)
    if (args.input == "None"):
        my_model_path = install_model(args.model_path, verbose=True)
        my_carbon_model_path = install_model(args.carbon_model_path, verbose=True)
        sys.exit(0)
   
    # Load checkpoint file
    my_model_path = install_model(args.model_path)
    model = UNet.load_from_checkpoint(my_model_path, n_channels=3, n_classes=1)
    model.eval().to(device)
    if not args.carbon:
        if args.verb > 0:
            print("  + Not doing carbon detection")
        do_carbon = False
    else:
        do_carbon = True
        if args.verb > 0:
            print("  + Doing carbon detection")
        my_carbon_model_path = install_model(args.carbon_model_path)
        carbon_model = UNetCarbon.load_from_checkpoint(my_carbon_model_path, n_channels=1, n_classes=1)
        carbon_model.eval().to(device)

    if args.plot:
        print("  + Writing PNG file with intermediate results for every micrograph")
    input_star = starfile.read(args.input)

    @torch.compile
    def run_compiled_carbon_model(data):
        seg_logits = carbon_model(data)
        carbon_prob = torch.sigmoid(seg_logits.squeeze())
        return carbon_prob

    @torch.compile
    def run_compiled_model(data):
        seg_logits, embedding = model(data)
        seg_prob = torch.sigmoid(seg_logits.squeeze())
        return seg_prob

    # ------------------ Run Inference ----------------    
    iterator = tqdm(input_star.iterrows(), total=len(input_star), file=sys.stdout) if args.verb > 0 else input_star.iterrows()
    half_dtype = torch.bfloat16 if torch.cuda.is_bf16_supported() else torch.float16
    for i, row in iterator:

        if Path(args.abort).exists():
            print("Abort file exists. Aborting...")
            sys.exit()
            
	# ------------------- Load MRCs ----------------------
        fom_image, pixel_size, origin = load_mrc(row["rlnMicrographFomImage"])
        psi_image, _, _ = load_mrc(row["rlnMicrographPsiImage"])

        #--------------------- Carbon detection -----------------
        if do_carbon:
            #in_tensor = np.stack(mic_image, axis=0)
            #in_tensor = torch.from_numpy(in_tensor[None,None]).to(device=device, dtype=half_dtype)  # shape (1,1,H,W)
            mic_image, ori_pixel_size, _ = load_mrc(row["rlnMicrographName"])
            mic_image, resized_angpix = resize_and_norm(mic_image, 64, ori_pixel_size)
            in_tensor = torch.from_numpy(mic_image[None,None]).to(device=device, dtype=half_dtype)  # shape (1,1,H,W)
            with torch.no_grad():
                if device.type == "cuda":
                    with torch.autocast(device_type=device.type, dtype=half_dtype):
                        carbon_prob = run_compiled_carbon_model(in_tensor)
                else:
                    with torch.autocast(device_type="cpu", dtype=half_dtype):
                        carbon_prob = run_compiled_carbon_model(in_tensor)
                carbon_prob = carbon_prob.to(dtype=torch.float32, device="cpu").numpy()
                xsf, ysf= fom_image.shape
                if carbon_prob.shape != fom_image.shape:
                    carbon_prob  = cv2.resize(carbon_prob.astype(np.float32), (ysf, xsf), interpolation=cv2.INTER_AREA)
                carbon_mask = (carbon_prob>args.carbon_threshold)
                
        #--------------------- Filament segmentation -----------------
        outfile = row["rlnMicrographCoordinates"]
        directory = Path(outfile).parent
        # Check if the directory exists, and create it if not
        directory.mkdir(parents=True, exist_ok=True)
        
        H, W = fom_image.shape

        rad_psi = np.deg2rad(psi_image * 2)
        in_tensor = np.stack((fom_image, np.sin(rad_psi), np.cos(rad_psi)), axis=0)  # shape (3,H,W)
        in_tensor = torch.from_numpy(in_tensor[None]).to(device=device, dtype=half_dtype)  # shape (1,3,H,W)

        with torch.no_grad():
            if device.type == "cuda":
                with torch.autocast(device_type=device.type, dtype=half_dtype):
                    seg_prob = run_compiled_model(in_tensor)
            else:
                with torch.autocast(device_type="cpu", dtype=half_dtype):
                    seg_prob = run_compiled_model(in_tensor)
            foreground = seg_prob.to(dtype=torch.float32, device="cpu").numpy()

        # 2) Clean + threshold => remove small objects
        bin_foreground = (foreground > args.threshold)
        if do_carbon:
            bin_foreground = bin_foreground & (~carbon_mask)

        min_size = 100
        clean_foreground = remove_small_objects(bin_foreground, min_size=min_size)

        # 3) Crop edges if radius is given
        if args.radius is not None:
            d = round(args.radius / pixel_size)
            clean_foreground[:d, :] = False
            clean_foreground[-d:, :] = False
            clean_foreground[:, :d] = False
            clean_foreground[:, -d:] = False

        # 4) Skeletonize
        skeleton_img = skeletonize(clean_foreground)
        
        # if there is only a single non-zero pixel, Skeleton will give an error. Let's check for at least 10 non-zero pixels
        if np.count_nonzero(skeleton_img) > 10:

	    # 5) Prune skeleton based on Skan's branches: splits branches with jumps in psi-angles, remove any branch with length below args.final_min_size, and enforce minimal distance between end points
            pruned_skeleton, sampled_points = prune_skeleton(skeleton_img,
                                                             psi_image,
                                                             min_length=args.minimum_length / pixel_size,
                                                             min_distance=2 * args.radius / pixel_size,
                                                             step=int(args.sample_step / pixel_size),
                                                             psi_threshold=args.psi_jump_threshold)

            # 6) Sample coordinates every 'sample_step' from pruned skeleton, **including filament ID**.
            #sampled_points = sample_coordinates_with_ids(pruned_skeleton, step=int(args.sample_step / pixel_size))
            # sampled_points shape: (N, 3) => [x, y, filament_id]
            if sampled_points.size > 0:
        
                # 7) Scale, round, write star
                # scale only the x,y columns
                sampled_points[:, 0:2] = np.round(sampled_points[:, 0:2] * args.scale).astype(int)

                df = pd.DataFrame({
                    'rlnCoordinateX': sampled_points[:, 0],
                    'rlnCoordinateY': sampled_points[:, 1],
                    'rlnParticleSelectionType':  sampled_points[:, 2]
                })
                starfile.write(df, outfile, overwrite=True)
            else:
                Path(outfile).touch()
        else:
            Path(outfile).touch()
	    
        # ---------------------------------------------------
        # If plotting is requested, create a collage figure
        # ---------------------------------------------------
        if args.plot:

            fig, ax = plt.subplots(2, 4, figsize=(20,15))
            ax = ax.ravel()

            ax[0].imshow(fom_image, cmap='gray')
            ax[0].set_title("FOM Image")
            ax[0].axis('off')

            ax[1].imshow(psi_image, cmap='gray')
            ax[1].set_title("PSI Image")
            ax[1].axis('off')

            ax[2].imshow(color_foreground(foreground, args.threshold))
            ax[2].set_title("Foreground")
            ax[2].axis('off')

            ax[3].imshow(skeleton_img)
            ax[3].set_title("Raw Skeleton")
            ax[3].axis('off')

            if do_carbon:
                ax[4].imshow(mic_image, cmap='gray')
                ax[4].set_title("micrograph")
                ax[4].axis('off')

                ax[5].imshow(color_foreground(carbon_prob, args.carbon_threshold))
                ax[5].set_title("carbon detection")
                ax[5].axis('off')

            if skeleton_img.any() and np.count_nonzero(skeleton_img) > 10:
                ax[6].imshow(color_nonrandom(pruned_skeleton))
                ax[6].set_title("Pruned Skeleton")
                ax[6].axis('off')

                if sampled_points.size > 0:
                    ax[7].imshow(np.zeros_like(foreground), cmap='gray')  # blank background
                    ax[7].scatter(sampled_points[:,0]/args.scale, sampled_points[:,1]/args.scale, s=2, c='red')
                    ax[7].set_xlim([0, foreground.shape[1]])
                    ax[7].set_ylim([foreground.shape[0], 0])
                    ax[7].set_title("Final Sampled Picks (ID shown in STAR)")
                    ax[7].axis('off')

            plt.tight_layout()
            out_png = os.path.splitext(outfile)[0] + ".png"
            plt.savefig(out_png, dpi=200, bbox_inches="tight")
            plt.close()


if __name__ == "__main__":
    main()
