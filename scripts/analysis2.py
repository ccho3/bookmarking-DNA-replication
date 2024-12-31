#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script reads dataset2 (snapshots at S phase 0 sec from embryos with CBP RNAi or JQ1 injection) and measures the mean and 
variance of intensities for Brd4 (C1) and Cdc7 (C2) in each object (nucleus) in the masked area (C3) in each 2D image.

Usage:
    python analysis2.py --input <input_directory> --output <output_directory>

"""

from analysis_functions import get_dirs, read_grouped_images, mean_var_per_nucleus
import pandas as pd
import os

# Parse arguments to get input and output directories
input_dir, output_dir = get_dirs()

# Set group labels
groups = ['DMSO', 'JQ1', 'w', 'nej']

# Read images
imgs = read_grouped_images(input_dir, groups)

# Analyze data
all_results=[]

for group, img_list in zip(groups, imgs):
    for replicate, img in enumerate(img_list, start=1):
        results = mean_var_per_nucleus(img, probe_ch=(0,1), mask_ch=2)

        results_df = pd.DataFrame(results, columns=['Brd4_mean', 'Brd4_variance',
                                         'Cdc7_mean', 'Cdc7_variance'])
        results_df['Group'] = group
        results_df['Replicate'] = replicate

        all_results.append(results_df)
        print(f'Finished analyzing {group} replicate {replicate}')
        
all_results_df = pd.concat(all_results, ignore_index=True)


# Export results
results_file = os.path.join(output_dir, "Brd4_Cdc7_treatments.csv")
all_results_df.to_csv(results_file, index=False)