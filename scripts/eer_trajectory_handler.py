#!/bin/env python3

import argparse
from collections import OrderedDict
from math import floor
import numpy as np
import os
import sys

def load_star(filename):
    datasets = OrderedDict()
    current_data = None
    current_colnames = None
    
    in_loop = 0 # 0: outside 1: reading colnames 2: reading data

    for line in open(filename):
        line = line.strip()
        
        # remove comments
        comment_pos = line.find('#')
        if comment_pos > 0:
            line = line[:comment_pos]

        if line == "":
            if in_loop == 2:
                in_loop = 0
            continue

        if line.startswith("data_"):
            in_loop = 0

            data_name = line[5:]
            current_data = OrderedDict()
            datasets[data_name] = current_data

        elif line.startswith("loop_"):
            current_colnames = []
            in_loop = 1

        elif line.startswith("_"):
            if in_loop == 2:
                in_loop = 0

            elems = line[1:].split()
            if in_loop == 1:
                current_colnames.append(elems[0])
                current_data[elems[0]] = []
            else:
                current_data[elems[0]] = elems[1]

        elif in_loop > 0:
            in_loop = 2
            elems = line.split()

            assert len(elems) == len(current_colnames)
            for idx, e in enumerate(elems):
                current_data[current_colnames[idx]].append(e)        
        
    return datasets

def write_star(filename, datasets):
    f = open(filename, "w")

    for data_name, data in datasets.items():
        f.write( "\ndata_" + data_name + "\n\n")
        
        col_names = list(data.keys())
        need_loop = isinstance(data[col_names[0]], list)
        if need_loop:
            f.write("loop_\n")
            for idx, col_name in enumerate(col_names):
                f.write("_%s #%d\n" % (col_name, idx + 1))
            
            nrow = len(data[col_names[0]])
            for row in range(nrow):
                f.write("\t".join([data[x][row] for x in col_names]))
                f.write("\n")
        else:
            for col_name, value in data.items():
                f.write("_%s\t%s\n" % (col_name, value))
        
        f.write("\n")
    f.close()
    
def interpolate_trajectory(traj_star, eer_grouping, old_grouping):
    nz = int(traj_star['general']['rlnImageSizeZ'])
    if (old_grouping <= 0):
        if 'rlnEERGrouping' not in traj_star['general']:
            sys.stderr.write("ERROR: The trajectory STAR file does not contain rlnEERGrouping. You have to specify the old grouping as --old_group.\n")
            sys.exit(-1)
        old_grouping = float(traj_star['general']['rlnEERGrouping'])
    new_nz = int(floor(nz * old_grouping / eer_grouping))
    scale = eer_grouping / old_grouping

    traj_star['general']['rlnImageSizeZ'] = str(new_nz)
    traj_star['general']['rlnMicrographDoseRate'] = str(float(traj_star['general']['rlnMicrographDoseRate']) * scale)
    traj_star['general']['rlnEERGrouping'] = eer_grouping
 
    xs = np.array(traj_star['global_shift']['rlnMicrographShiftX'], dtype=np.float)
    ys = np.array(traj_star['global_shift']['rlnMicrographShiftY'], dtype=np.float)

    new_xs = np.zeros(new_nz)
    new_ys = np.zeros(new_nz)

    # This interpolation is not very accurate. We should take
    # the MIDDLE, not the start of a range, as an observation point.
    # However, such small error should be corrected in Polish anyway.
    for i in range(new_nz):
        src = i * scale

        src1 = int(floor(src))
        src2 = src1 + 1

        frac = src - src1
        #print(i, src, src1, src2)
        if src2 >= nz: # be lazy; don't extrapolate
            new_xs[i] = xs[nz - 1]
            new_ys[i] = ys[nz - 1]
        else:
            new_xs[i] = xs[src1] * (1 - frac) + xs[src2] * frac
            new_ys[i] = ys[src1] * (1 - frac) + ys[src2] * frac

    traj_star['global_shift']['rlnMicrographFrameNumber'] = list(np.linspace(1, new_nz, num=new_nz).astype(np.int).astype(np.str0))
    traj_star['global_shift']['rlnMicrographShiftX'] = list(new_xs.astype(np.str0))
    traj_star['global_shift']['rlnMicrographShiftY'] = list(new_ys.astype(np.str0))
 
    # z is not normalized, so have to be patched.
    if "local_motion_model" in traj_star:
        coeffs = np.array(traj_star['local_motion_model']['rlnMotionModelCoeff'], dtype=np.float)
        coeffs *= scale       # 1st-order in time(z)
        coeffs[1::3] *= scale # 2nd-order
        coeffs[2::3] *= scale # 3rd-order
        traj_star['local_motion_model']['rlnMotionModelCoeff'] = list(coeffs.astype(np.str0))

    return traj_star

def resample_image(traj_star, eer_upsampling):
    orig_size = int(traj_star['general']['rlnImageSizeX'])
    assert orig_size == int(traj_star['general']['rlnImageSizeY'])

    if (orig_size == 4096 and eer_upsampling == 2):
        scale = 2.0
    elif (orig_size == 8192 and eer_upsampling == 1):
        scale = 0.5
    else:
        raise "Illegal eer_upsampling"

    traj_star['general']['rlnImageSizeX'] = str(int(orig_size * scale))
    traj_star['general']['rlnImageSizeY'] = str(int(orig_size * scale))
    traj_star['general']['rlnMicrographBinning'] = str(eer_upsampling)
    traj_star['general']['rlnEERUpsampling'] = str(eer_upsampling)

    traj_star['general']['rlnMicrographOriginalPixelSize'] = str(float(traj_star['general']['rlnMicrographOriginalPixelSize']) / scale)

    xs = np.array(traj_star['global_shift']['rlnMicrographShiftX'], dtype=np.float) * scale
    ys = np.array(traj_star['global_shift']['rlnMicrographShiftY'], dtype=np.float) * scale
    traj_star['global_shift']['rlnMicrographShiftX'] = list(xs.astype(np.str0))
    traj_star['global_shift']['rlnMicrographShiftY'] = list(ys.astype(np.str0))

    # Hot pixels
    if 'hot_pixels' in traj_star:
        hot_xs = np.array(traj_star['hot_pixels']['rlnCoordinateX'], dtype=np.float)
        hot_ys = np.array(traj_star['hot_pixels']['rlnCoordinateY'], dtype=np.float)
        if scale == 2:
            hot_xs = np.hstack([2 * hot_xs, 2 * hot_xs, 2 * hot_xs + 1, 2 * hot_xs + 1])
            hot_ys = np.hstack([2 * hot_ys, 2 * hot_ys + 1, 2 * hot_ys, 2 * hot_ys + 1])
        elif scale == 0.5:
            tmp = np.floor(np.vstack([hot_xs, hot_ys]) / 2.0).astype(np.int)
            tmp = np.unique(tmp, axis = 1)
            hot_xs = tmp[0, :]
            hot_ys = tmp[1, :]

        traj_star['hot_pixels']['rlnCoordinateX'] = list(hot_xs.astype(np.str0))
        traj_star['hot_pixels']['rlnCoordinateY'] = list(hot_ys.astype(np.str0))

    return traj_star

def add_suffix(filename, suffix):
    tmp = os.path.splitext(filename)
    return "%s_%s%s" % (tmp[0], suffix, tmp[1])

parser = argparse.ArgumentParser(description='Tweak motion trajectory STAR files for EER movies')
parser.add_argument('--i', type=str, nargs='?', metavar='corrected_micrographs.star', required=True,
                    help='Motion correction STAR file')
parser.add_argument('--o', type=str, nargs='?', metavar='suffix', required=True,
                    help='Suffix for output files')
parser.add_argument('--old_group', type=int, nargs='?', metavar='group', default=0,
                    help='Old EER grouping (must be specified when not recorded in the STAR file)')
parser.add_argument('--regroup', type=int, nargs='?', metavar='group', default=0,
                    help='Regroup to this number of physical frames / fraction')
parser.add_argument('--resample', type=int, nargs='?', metavar='sampling', default=0,
                    help='Resample to this level. 1=4K, 2=8K (super-res)')

args = parser.parse_args()
#print(args)
fn_motioncorr_star = args.i
suffix = args.o

if (args.resample == 0 and args.regroup == 0):
    sys.stderr.write("Error: Nothing to do. Please specify --resample and/or --regroup.\n")
    sys.exit(-1)

motioncorr_star = load_star(fn_motioncorr_star)
print("Read %s" % fn_motioncorr_star)
print("Found %d movies" % len(motioncorr_star['micrographs']['rlnMicrographMetadata']))

for idx, fn_traj in enumerate(motioncorr_star['micrographs']['rlnMicrographMetadata']):
    fn_out = add_suffix(fn_traj, suffix)
    motioncorr_star['micrographs']['rlnMicrographMetadata'][idx] = fn_out
    print("Processing %s => %s" % (fn_traj, fn_out))

    traj_star = load_star(fn_traj)

    if (args.regroup > 0):
        interpolate_trajectory(traj_star, args.regroup, args.old_group)
    if (args.resample > 0):
        resample_image(traj_star, args.resample)

    # local_shift table is not updated, because it is not used by Polish.
    # To avoid confusion, delete it.
    if 'local_shift' in traj_star:
        del traj_star['local_shift']

    write_star(fn_out, traj_star)
    #break

fn_out = add_suffix(fn_motioncorr_star, suffix)
write_star(fn_out, motioncorr_star)
print("Written %s" % fn_out)
