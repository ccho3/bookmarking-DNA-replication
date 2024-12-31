#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script quantifies the maximal intensity of EGFP-Rpb3 (C2) using nuclear mask (C3) generated from mCherry-PCNA (C1).

Usage:
    python analysis6.py --input <input_directory> --output <output_directory>

"""

from analysis_functions import get_dirs, read_grouped_images, max_per_image
import pandas as pd
import os

# Parse arguments to get input and output directories
input_dir, output_dir = get_dirs()

# Set group labels
groups = ['water', 'XL413']

# Read images
imgs = read_grouped_images(input_dir, groups)

all_results=[]

for group, img_list in zip(groups, imgs):
    for replicate, img in enumerate(img_list, start=1):
        for frame in range(img.shape[0]):
            img_frame = img[frame,:,:]
            C1_max, C2_max = max_per_image(img_frame, probe_ch=(0, 1), mask_ch=2)
            result1 = (group, replicate, frame, 'Ch1', C1_max)
            all_results.append(result1)
            result2 = (group, replicate, frame, 'Ch2', C2_max)
            all_results.append(result2)
        print(f'Finished analyzing {group} replicate {replicate}')
     
all_results = pd.DataFrame(all_results, columns=['Group', "Replicate", "Frame", 'Channel',
                                        'Intensity_Max'])

results_file = os.path.join(output_dir, "PCNA_Rpb3_max_Cdc7i.csv")
all_results.to_csv(results_file, index=False)