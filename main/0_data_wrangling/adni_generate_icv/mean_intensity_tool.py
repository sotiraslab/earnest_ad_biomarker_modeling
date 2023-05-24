#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 20:42:02 2022

@author: earnestt1234
"""

import argparse
import os
import re

import nibabel as nib
import numpy as np
import pandas as pd

def find(directory, pattern):

    results = []
    for root, _, files in os.walk(directory):
        for file in files:
            if re.search(pattern, file):
                results.append(os.path.join(root, file))

    return results

def get_intesity_stats(nii, mask=None, binarize_mask=False, nonzero=True):

    """Returns mean intensity, number voxels, and voxel volume.  These calculations
    can be done in an ROI (`mask`) and with or without omitting zeros (`nonzero`)."""

    # read data
    data = nii.get_fdata().copy()

    if mask is not None:
        maskdata = mask.get_fdata().copy()
        if binarize_mask:
            maskdata = np.where(maskdata > 0, 1, 0)
        data = data * maskdata

    if nonzero:
        data[data <= 0] = np.nan

    # return 1: mean intensity
    mean_intensity = np.nanmean(data)

    # return 2: number voxels
    voxel_number = np.sum(~np.isnan(data))

    # return 3: voxel dimensions
    voxel_dims = (nii.header["pixdim"])[1:4]
    voxel_volume = np.prod(voxel_dims)

    return mean_intensity, voxel_number, voxel_volume

def parse():

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', help='Input directory to search', required=True)
    parser.add_argument('-p', '--pattern', help="Regex pattern to hook images", required=True)
    parser.add_argument('-o', '--output', help="Output CSV path, with extension", required=True)
    parser.add_argument('-m', '--mask', help="Calculate the volume within one or more ROIs", default=None, nargs='+')
    parser.add_argument('-b', '--binarizemask', help='Binarize any masks being passed (threshold at 0)', default=False, action='store_true')

    return parser.parse_args()

def main(infolder, pattern, output, nonzero=True, masks=None, binarize_masks=False):

    images= find(infolder, pattern)
    N = len(images)

    print()
    print(f"Found {N} images matching pattern '{pattern}'.")
    print(f"Using masks: {masks is not None}")

    if masks is not None:
        print("Masks:")
        for m in masks:
            print(f"    - {m}")

    rows = []

    print()
    for i, img in enumerate(images):

        # progress bar
        l = 15
        prop = (i+1) / N
        prog = int(round(prop * l, 0))
        s = '[' + '-' * prog + '*' + (l-prog) * ' ' + ']'

        name = os.path.basename(img)
        print(f"{s} ({i+1}/{N}) Image: {name}")

        nii = nib.load(img)

        if masks is None:
            mean_intensity, voxel_number, voxel_volume = get_intesity_stats(nii, nonzero=nonzero)
            row = {'Path':img,
                   'VoxelVolume': voxel_volume,
                   'MeanIntensity': mean_intensity,
                   'NumberVoxels': voxel_number}
        else:
            row = {"Path": img}
            for m in masks:

                mname = os.path.basename(m)
                mask = nib.load(m)
                mean_intensity, voxel_number, voxel_volume = get_intesity_stats(nii,
                                                                                mask=mask,
                                                                                binarize_mask=binarize_masks,
                                                                                nonzero=nonzero)

                if not row.get("VoxelVolume"):
                    row['VoxelVolume'] = voxel_volume

                row[f'MeanIntensity_{mname}'] = mean_intensity
                row[f'NumberVoxels_{mname}'] = voxel_number

        rows.append(row)

    df = pd.DataFrame(rows)

    if not df.empty:
        print()
        print(f"Saving output to {output}...")
        df.to_csv(output, index=False)
        print("Finished!")

if __name__ == '__main__':
    args = parse()
    main(infolder=args.input, pattern=args.pattern, output=args.output,
         masks=args.mask, binarize_masks=args.binarizemask)