#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script reads processed dataset 1 (movies of S-phase entry into NC11-14) from input directory, measures the mean and variance of 
intensities for mKate2-Brd4 (C1) and Cdc7-EGFP (C2) in the masked area (C3) across all the frames in each movie, and then
saves the result csv file in the output directory.

Usage:
    python analysis1.py --input <input_directory> --output <output_directory>
"""

from analysis_functions import get_dirs, read_grouped_images, mean_var_per_image
import pandas as pd
import os

# Parse arguments to get input and output directories
input_dir, output_dir = get_dirs()

# Set group labels
groups = ['M10', 'M11', 'M12', 'M13']

# Read images
imgs = read_grouped_images(input_dir, groups)

# Analyze data
all_results=[]

for group, img_list in zip(groups, imgs):
    for replicate, img in enumerate(img_list, start=1):
        for frame in range(7):
            img_frame = img[frame,:,:,:]
            C1_mean, C1_var, C2_mean, C2_var = mean_var_per_image(img_frame, probe_ch=(0, 1), mask_ch=2)
            result1 = (group, replicate, frame, 'Ch1', C1_mean, C1_var)
            all_results.append(result1)
            result2 = (group, replicate, frame, 'Ch2', C2_mean, C2_var)
            all_results.append(result2)
        print(f'Finished analyzing {group} replicate {replicate}')
        

# Export results
all_results = pd.DataFrame(all_results, columns=['Group', 'Replicate', "Frame", 'Channel',
                                        'Intensity_Mean', 'Intensity_Variance'])


results_file = os.path.join(output_dir, "Brd4_Cdc7_NC11-14.csv")
all_results.to_csv(results_file, index=False)