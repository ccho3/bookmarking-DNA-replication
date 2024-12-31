#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script performs segmentation of nuclei and MCP foci, tracks nuclei across frames, 
and quantifies MCP states in multi-channel microscopy movies.

The analysis involves the following steps:
1. Initialization (`--run initialize`): 
   - Reads TIFF files from the specified input directory.
   - Segments nuclei using the EGFP-Rpb3 channel (channel 2).
   - Tracks nuclei across frames.
   - Detects MCP foci using the MCP-mCherry channel (channel 1).
   - Visualizes and allows for manual editing of the segmentation and tracking results using napari.
   - Press 'u' in the napari viewer to saves intermediate data for further processing.
   - Repeat the above steps until all files have been processed.

2. Finalization (`--run finish`): 
   - Reads the intermediate data.
   - Assigns MCP states (on/off) to each tracked nucleus based on the edited MCP segmentation.
   - Summarizes the MCP state transitions and cumulative MCP status across frames.
   - Outputs the results to a CSV file.

Usage:
    python analysis4.py --input <input_directory> --output <output_directory> --run <initialize|finish>

    - `<input_directory>`: Path to the directory containing the input TIFF movies.
    - `<output_directory>`: Path to the directory where the results will be saved.
    - `<initialize|finish>`: Specify whether to initialize or finish the analysis.

Example:
    python analysis4.py --input "../data/dataset4/" --output "../results/" --run initialize
    python analysis4.py --input "../data/dataset4/" --output "../results/" --run finish


Dependencies:
- Python 3.12.4
- napari 0.5.1
- numpy 1.26.4 
- pandas 2.2.2 
- scipy 1.14.0
- scikit-image 0.24.0
- trackpy 0.6.4
"""

# Import packages
import os, glob, argparse, napari
import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage import io
from skimage.filters import gaussian, threshold_otsu
from skimage.measure import label, regionprops_table
from skimage.feature import peak_local_max, blob_log
from skimage.segmentation import watershed, clear_border
from skimage.morphology import binary_dilation, remove_small_objects
import trackpy as tp
from analysis_functions import read_single_image

# Define functions
def process_nuclei(img):
    """
    Generates nuclear masks from input movie using the EGFP-Rpb3 channel (channel 2).

    Parameters:
        img (numpy array): A 4D numpy array representing the multi-channel movie with dimensions (frames, channels, height, width). 
            Channel 1 is assumed to be MCP-mCherry, and channel 2 is EGFP-Rpb3.

    Returns:
        nuclear_masks: A 3D numpy array containing binary nuclear masks with dimensions (frames, height, width), where nuclei are marked as True (1) and the 
            background as False (0).
    """

    print("Getting nuclear masks")
    Rpb3 = img[:,1,:,:]
    nuclear_masks = np.zeros_like(Rpb3)

    for i in range(Rpb3.shape[0]):
        img = Rpb3[i, :, :]
        # Reduce bright HLBs
        img[img > 800] = 800
        img_blur = gaussian(img, sigma=4)
        # Otsu thresholding and dilation
        nuc_threshold = threshold_otsu(img_blur)
        img_mask = img_blur > nuc_threshold
        img_mask = binary_dilation(img_mask, footprint=np.ones((6, 6)))
        
        # Perform Watershed to separate nuclei
        distance = ndi.distance_transform_edt(img_mask)
        coords = peak_local_max(distance, min_distance=10)
        coord_mask = np.zeros(distance.shape, dtype=bool)
        coord_mask[tuple(coords.T)] = True
        markers = label(coord_mask)
        nuc_mask = watershed(-distance, markers, mask=img_mask, watershed_line=True)

        # Additional filtering
        nuc_mask = clear_border(nuc_mask)
        nuc_mask = remove_small_objects(nuc_mask, min_size=500)

        # Save the mask to the stack
        nuclear_masks[i,:,:] = nuc_mask>0

    return nuclear_masks


def track_nuclei(nuclear_masks):
    """
    Generates nuclear labels with unique IDs and tracks nuclei across frames.

    Parameters:
        nuclear_masks: A 3D numpy array containing binary nuclear masks.

    Returns:
        nuclear_labels: A 3D numpy array containing labeled nuclei. Each nucleus is assigned a unique ID.
       
        nuclear_table (pandas.DataFrame): A DataFrame containing the tracking results.
    """
    
    print("Labeling and tracking nuclei")
    
    frames = nuclear_masks.shape[0]
    
    nuclear_table = pd.DataFrame()
    nuclear_labels = np.zeros_like(nuclear_masks, dtype=int)
    for i in range(frames):
        nuc_mask = nuclear_masks[i,:,:]
        nuc_label = label(nuc_mask, connectivity=1)
        nuclear_labels[i,:,:] = nuc_label
        
        n = regionprops_table(nuc_label, properties=('centroid',))
        n = pd.DataFrame(n)
        n['frame'] = i
        n = n.rename(columns={"centroid-1": "x", "centroid-0": "y"})
        n = n[['frame', 'x', 'y']]
        nuclear_table = pd.concat([nuclear_table, n], ignore_index=True)
        
    # Tracking and filtering
    nuclear_table = tp.link_df(nuclear_table, search_range=20, memory=1)
    nuclear_table = tp.filter_stubs(nuclear_table, threshold=frames)
    nuclear_table = nuclear_table.reset_index(drop=True)

    return nuclear_labels,nuclear_table


def process_MCP(img, min_sigma=0.6, max_sigma=4, threshold=0.001):
    """
    Generates MCP masks using the movie's channel 1 using the Laplacian of Gaussian (LoG) method.

    Parameters:
        img : A 4D numpy array with dimensions (frames, channels, height, width) whose channel 1 corresponds to MCP-mCherry.
        min_sigma (float): The minimum standard deviation for Gaussian kernel used in the LoG blob detection.
        max_sigma (float): The maximum standard deviation for Gaussian kernel used in the LoG blob detection.
        threshold (float): The absolute lower bound for scale space maxima.

    Returns:
        MCP_masks: A 3D array with dimensions (frames, height, width) containing binary MCP masks. Each detected MCP 
            focus is marked as True (1).
    """
    
    print("Getting MCP masks")
    MCP = img[:,0,:,:]    
    MCP_masks = np.zeros_like(MCP, dtype=int)

    for i in range(MCP.shape[0]):
    #for i in range(8): # For testing
        img = MCP[i, :, :]
        img_blur = gaussian(img, sigma=2)
        blobs = blob_log(img_blur, min_sigma=min_sigma, max_sigma=max_sigma, threshold=threshold)
        print(f"Detected {blobs.shape[0]} blobs in frame {i}")
        # Create blob label
        blob_label = np.zeros_like(img, dtype=int)
        for j, blob in enumerate(blobs, start=1):
            y, x, r = blob
            # Create a circular mask for the blob
            rr, cc = np.ogrid[:img.shape[0], :img.shape[1]]
            mask = (rr - y)**2 + (cc - x)**2 <= r**2
            blob_label[mask] = j
        # Save the blob image to the stack
        MCP_masks[i,:,:] = blob_label>0

    return MCP_masks


def view_results(img, nuclear_masks, MCP_masks, nuclear_table):
    """Displays the movie, segmentation masks, and tracking data in a napari viewer."""
    
    Rpb3 = img[:,1,:,:]
    MCP = img[:,0,:,:]
    tracks_data = nuclear_table[['particle', 'frame', 'y', 'x']].values
    
    viewer = napari.Viewer()
    viewer.add_image(Rpb3, name='Rpb3', blending='additive', 
                     colormap='green', contrast_limits=[0,1500])
    viewer.add_image(nuclear_masks, name='nuclei', blending='additive',
                     opacity=0.50, colormap='blue')
    viewer.add_image(MCP, name='MCP', blending='additive', 
                     colormap='red', contrast_limits=[0,1500])
    viewer.add_labels(MCP_masks*9, name='MCP masks', blending='additive',
                     opacity=1.00)
    viewer.add_tracks(tracks_data, name='Tracks', tail_length=6, 
                      head_length=6, tail_width=4, colormap='hsv')
    viewer.layers['MCP masks'].brush_size=40
    
    return viewer


def update_MCP_masks(viewer):
    """
    Extracts the MCP masks that have been manually edited in the napari viewer.

    Parameters:
        viewer: The napari viewer instance where MCP masks have been edited.

    Returns:
        MCP_masks: The updated MCP masks.
    """

    MCP_masks = viewer.layers['MCP masks'].data
    print("MCP masks updated")

    return MCP_masks


def read_temp_data(file_path):
    """
    Reads temporary data files during "--run initialization" step

    Parameters:
        file_path (str): Path to the original movie file.

    Returns:
        
        nuclear_labels: A 3D numpy array representing labeled nuclei.

        nuclear_table: A DataFrame containing the tracking results.

        MCP_masks: The MCP masks as a 3D numpy array representing detected MCP foci.
    """
    temp_dir = os.path.join(os.path.dirname(file_path), "temp")
    filename = os.path.basename(file_path).split('.')[0]
    MCP_masks = io.imread(os.path.join(temp_dir, f'{filename}_MCP_masks.tif'))
    nuclear_labels = io.imread(os.path.join(temp_dir, f'{filename}_nuclear_labels.tif'))
    nuclear_table = pd.read_csv(os.path.join(temp_dir,f'{filename}_table.csv'))

    return nuclear_labels, nuclear_table, MCP_masks


def assign_MCP_state(nuclear_labels, MCP_masks, nuclear_table):
    """
    Uses the nuclear labels and MCP masks to determine the MCP state for each tracked nucleus. 
    It also calculates the transitions between MCP states (ON and OFF) for each nucleus over time.

    Parameters:
        nuclear_labels: A 3D numpy array containing labeled nuclei for each frame.
        
        MCP_masks: A 3D numpy array containing MCP masks.

        nuclear_table: A DataFrame containing tracking results.

    Returns:
        nuclear_table: The updated tracking DataFrame with additional columns:
                          - 'MCP_state': Binary indicator of MCP presence (1 if MCP is present, 0 otherwise).
                          - 'MCP_ON': Indicator of MCP state transition from OFF to ON compared to previous frame
                          - 'MCP_OFF': Indicator of MCP state transition from ON to OFF compared to previous frame
    """
    
    print("Assigning MCP state to tracked nuclei")
    # Convert data frame to numpy array
    t1 = np.asarray(nuclear_table)
    
    # Initialize lists to store the MCP data
    MCP_data = []

    # Iterate over each row in the tracking data array
    for row in t1:
        frame = int(row[0])
        x = int(row[1])
        y = int(row[2])
        nuclei = nuclear_labels[frame,:,:]
        MCP = MCP_masks[frame,:,:]
        nuc_id = nuclei[y,x]
        MCP_sum = MCP[nuclei==nuc_id].sum()
        MCP_data.append(MCP_sum)

    MCP_data = np.asarray(MCP_data)>0    
    nuclear_table['MCP_state'] = MCP_data.astype(int)

    # Add columns 'MCP_ON' and 'MCP_OFF' indicating the history of MCP state transitions.
    nuclear_table['MCP_ON'] = 0
    nuclear_table['MCP_OFF'] = 0
    
    # Group by particle to analyze each nucleus individually
    for particle, group in nuclear_table.groupby('particle'):
        group = group.sort_values(by='frame')
    
        # Find the first "on" frame (where MCP_state changes from 0 to 1)
        on_frames = group[group['MCP_state'].diff() == 1]['frame']
        if not on_frames.empty:
            on_frame = on_frames.iloc[0]
            # Mark all frames from the first "on" frame onward as MCP_ON = 1
            nuclear_table.loc[
                (nuclear_table['particle'] == particle) & 
                (nuclear_table['frame'] >= on_frame), 
                'MCP_ON'] = 1
    
        # Find the first "off" frame (where MCP_state changes from 0 to 1)
        off_frames = group[group['MCP_state'].diff() == -1]['frame']
        if not off_frames.empty:
            off_frame = off_frames.iloc[0]
            # Mark all frames from the first "off" frame onward as MCP_OFF = 1
            nuclear_table.loc[
                (nuclear_table['particle'] == particle) & 
                (nuclear_table['frame'] >= off_frame), 
                'MCP_OFF'] = 1

    return nuclear_table


def summarize(nuclear_table, file_path):
    """
    Calculates the number of active nuclei across frames and the cumulative counts of nuclei undergoing MCP state transitions
    up to the current frame. 

    Parameters:
        nuclear_table: A DataFrame containing tracking results.

        file_path (str): Path to the file used to generate the summary. The filename is used to extract
                         experimental metadata such as treatment, cycle, and replicate.

    Returns:
        summary: A DataFrame summarizing MCP presence and transitions.
    """
    print("Summarizing results")
    
    filename = os.path.basename(file_path).split('.')[0]

    summary = nuclear_table.groupby('frame').agg(
                        instantaneous_MCP_ON=('MCP_state', 'sum'),
                        cumulative_MCP_ON=('MCP_ON', 'sum'),
                        cumulative_MCP_off=('MCP_OFF', 'sum')).reset_index()
    
    summary['total_nuclei'] = nuclear_table.groupby('frame').size().values
    summary['treatment'] = filename.split('_')[0]
    summary['cycle'] = filename.split('_')[1]
    summary['replicate'] = filename.split('_')[2]
    
    return summary


def initialize_analysis(file_path):
    img, filename = read_single_image(file_path)
    input_dir = os.path.dirname(file_path)
    nuclear_masks = process_nuclei(img)
    nuclear_labels, nuclear_table = track_nuclei(nuclear_masks)
    MCP_masks = process_MCP(img, min_sigma=0.6, max_sigma=4, threshold=0.001)
    viewer = view_results(img, nuclear_masks, MCP_masks, nuclear_table)
    nuclear_labels = nuclear_labels.astype(np.uint16)


    def save_temp_data(viewer):
        temp_dir = os.path.join(input_dir, "temp")
        MCP_masks = update_MCP_masks(viewer)
        MCP_masks = MCP_masks.astype(np.uint16)

        io.imsave(os.path.join(temp_dir, f'{filename}_nuclear_labels.tif'), nuclear_labels, check_contrast=False)
        io.imsave(os.path.join(temp_dir, f'{filename}_MCP_masks.tif'), MCP_masks, check_contrast=False)
        nuclear_table.to_csv(os.path.join(temp_dir,f'{filename}_table.csv'), index=False)
        
        print("Temporary data saved and viewer closed.")
        viewer.close()

    viewer.bind_key('u', save_temp_data)

    return viewer


def finish_analysis(file_path):
    nuclear_labels, nuclear_table, MCP_masks = read_temp_data(file_path)
    nuclear_table = assign_MCP_state(nuclear_labels, MCP_masks, nuclear_table)
    summary = summarize(nuclear_table, file_path)
    return summary



# Analyze data

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Analyze Rpb3 MCP movies")
    parser.add_argument("--input", type=str, required=True, help="path to the data")
    parser.add_argument("--output", type=str, required=True, help="path to the output directory")
    parser.add_argument("--run", type=str, required=True, help="initialize or finish")

    args = parser.parse_args()

    input_dir = args.input
    output_dir = args.output

    
    file_path_list = glob.glob(os.path.join(input_dir, '*.tif'))


    if args.run == "initialize":

        temp_dir = os.path.join(input_dir, "temp")
        os.makedirs(temp_dir, exist_ok=True)

        for file_path in file_path_list:
            viewer = initialize_analysis(file_path)
            print('Press u after editing MCP masks')
            napari.run()

    elif args.run == "finish":
        summaries = []
        for file_path in file_path_list:
            summary = finish_analysis(file_path)
            summaries.append(summary)
        
        combined_summary = pd.concat(summaries, ignore_index=True)
        combined_summary.to_csv(os.path.join(output_dir, "Rpb3_MCP_summary.csv"), index=False)
        print("Analysis finished and summary saved.")

