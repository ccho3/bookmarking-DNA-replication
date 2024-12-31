#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This scripts reads dataset3 (snapshots from embryos with control and Cdc7i injection) and measures the maximal intensity of Brd4 (C1) or Cdc7 (C2)
in each object (nucleus) in the masked area (C3) in each 2D image.

Usage:
    python analysis3.py --input <input_directory> --output <output_directory>

"""

from analysis_functions import get_dirs, read_grouped_images, max_per_nucleus
import pandas as pd
import os

# Parse arguments to get input and output directories
input_dir, output_dir = get_dirs()

# Set group labels
groups = ['water_0sec', 'XL413_0sec', 'water_80sec', 'XL413_80sec']

# Get group information
Treatments = [group.split('_')[0] for group in groups]
Times = [group.split('_')[1] for group in groups]

imgs = read_grouped_images(input_dir, groups)

# Analyze data
all_results=[]

for group, img_list, treatment, time in zip(groups, imgs, Treatments, Times):
    for replicate, img in enumerate(img_list, start=1):
        results = max_per_nucleus(img, probe_ch=(0,1), mask_ch=2)

        results_df = pd.DataFrame(results, columns=['Brd4_max', 'Cdc7_max'])
        results_df['Treatment'] = treatment
        results_df['Time'] = time
        results_df['Replicate'] = replicate

        all_results.append(results_df)
        print(f'Finished analyzing {group} replicate {replicate}')
        
all_results_df = pd.concat(all_results, ignore_index=True)

# Export results
results_file = os.path.join(output_dir, "Brd4_Cdc7_max_Cdc7i.csv")
all_results_df.to_csv(results_file, index=False)
