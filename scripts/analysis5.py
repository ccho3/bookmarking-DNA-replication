#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script finds movies in input directory, performs the detection of MCP/PCP foci, track them across frames,
and then generates montages displaying cropped foci images for all tracks or the first 50 tracks in each movie.

Usage:
    python analysis5.py --input <input_directory> --output <output_directory>

    - `<input_directory>`: Path to the directory containing the input TIFF movies.
    - `<output_directory>`: Path to the directory where the montages will be saved.

Dependencies:
- Python 3.12.4
- napari 0.5.1
- numpy 1.26.4 
- pandas 2.2.2 
- scipy 1.14.0
- scikit-image 0.24.0
- trackpy 0.6.4
"""

# Import modules
import os, glob, napari, argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from skimage import io
from skimage.filters import gaussian
from skimage.feature import blob_dog
import trackpy as tp
from analysis_functions import read_single_image

# Define functions

def detect_blobs(img, min_sigma=0.1, max_sigma=8, threshold=0.0003):
    """
    Detects and labels transcriptional spots in a two-channel movie using the Difference of Gaussians (DoG) method.
    
    Parameters:
        img: A 4D numpy array with dimensions (frames, channels, height, width) whose channel 1 corresponds to PCP-mCherry and channel 2 
            corresponds to MCP-GFP.
        min_sigma (float): The minimum standard deviation for the Gaussian kernel.
        max_sigma (float): The maximum standard deviation for the Gaussian kernel.
        threshold (float): The threshold value for detecting blobs.

    Returns: 
        merge_blur: A numpy array for the merged and blurred images used to detect blobs.
        blob_labels: The labeled blobs for each frame.
        blobs_df: A pd DataFrame containing the detected blobs' information.
    """
    
    # Extract PCP and MCP movies, assuming that they are in the first and second channels
    PCP = img[:, 0, :, :]
    MCP = img[:, 1, :, :]
    
    # Merge the channels for tracking transcriptional loci
    merge = PCP + MCP
    
    # Initialize empty arrays and list to store results
    merge_blur = np.zeros_like(merge, dtype=np.float64)
    blob_labels = np.zeros_like(merge, dtype=int)
    blobs_all_frames = []
    
    # Detect blobs in each frame 
    for i in range(merge.shape[0]):
        # Apply Gaussian filter
        img_i = gaussian(merge[i, :, :], sigma=3)
        # Save the merged/blurred image
        merge_blur[i, :, :] = img_i
        
        # Detect blobs
        blobs = blob_dog(img_i, min_sigma=min_sigma, max_sigma=max_sigma, 
                         threshold=threshold)
        print("Detected", blobs.shape[0], "blobs in frame", i)
        
        # Create blob label image
        blob_label_i = np.zeros_like(img_i, dtype=int)
        for j, blob in enumerate(blobs, start=1):
            y, x, r = blob
            # Create a circular mask for the blob
            rr, cc = np.ogrid[:img_i.shape[0], :img_i.shape[1]]
            mask = (rr - y)**2 + (cc - x)**2 <= r**2
            blob_label_i[mask] = j
        # Save the blob image to the stack
        blob_labels[i,:,:] = blob_label_i
        
        # Save blob coordinates and frame information for tracking
        for blob in blobs:
            y, x, r = blob
            blobs_all_frames.append([i, x, y, r])
        
    # Convert blob data to DataFrame for tracking
    blobs_df = pd.DataFrame(blobs_all_frames, columns=['frame', 'x', 'y', 'radius'])
    
    return merge_blur, blob_labels, blobs_df


def measure_intensities(img, blobs_df, radius=8, bg_radius=15):
    """
    Measures and adds the integrated intensities of PCP and MCP channels to the blobs_df DataFrame.
    
    Parameters:
        img: A 4D numpy array with dimensions (frames, channels, height, width) whose channel 1 corresponds to PCP-mCherry and channel 2 
            corresponds to MCP-GFP.
        blobs_df (pandas.DataFrame): The DataFrame containing detected blobs from `detect_blobs`.
        radius (int): The radius of the region around each blob for measuring intensity.
        bg_radius (int): The radius of the region for measuring local background intensity.

    Returns:
        blobs_df: The input DataFrame with additional columns for PCP and MCP integrated intensities.
    """

    # Extract channels
    PCP = img[:, 0, :, :]
    MCP = img[:, 1, :, :]
    # Convert data frame to numpy array
    b1 = np.asarray(blobs_df)
    
    # Initialize lists to store the intensities
    PCP_data = []
    MCP_data = []

    # Iterate over each row in the tracking data array
    for row in b1:
        frame = int(row[0])
        x = int(row[1])
        y = int(row[2])
        
        # Define the region of interest (ROI) around the (x, y) coordinate
        y_min = max(y - radius, 0)
        y_max = min(y + radius + 1, PCP.shape[1])
        x_min = max(x - radius, 0)
        x_max = min(x + radius + 1, PCP.shape[2])
        
        # Calculate the integrated intensity for PCP and MCP channels
        PCP_intensity = np.mean(PCP[frame-1, y_min:y_max, x_min:x_max])
        MCP_intensity = np.mean(MCP[frame-1, y_min:y_max, x_min:x_max])
        
        
        # Calculate the local background intensity
        # Define the region of interest (ROI) around the (x, y) coordinate
        y1_min = max(y - 15, 0)
        y1_max = min(y + 15 + 1, PCP.shape[1])
        x1_min = max(x - 15, 0)
        x1_max = min(x + 15 + 1, PCP.shape[2])
        
        
        # Subtract the local background
        PCP_bg = np.median(PCP[frame-1, y1_min:y1_max, x1_min:x1_max])
        MCP_bg = np.median(MCP[frame-1, y1_min:y1_max, x1_min:x1_max])
        PCP_intensity = PCP_intensity - PCP_bg
        MCP_intensity = MCP_intensity - MCP_bg
        
        # Append the intensities to the lists
        PCP_data.append(PCP_intensity)
        MCP_data.append(MCP_intensity)
    
    # Add the intensities to the tracking DataFrame
    blobs_df['PCP_intensity'] = PCP_data
    blobs_df['MCP_intensity'] = MCP_data
    
    return blobs_df
    

def track_filter(img, blobs_df, search_range=20, memory=5, threshold=2, border=20):
    """
    Tracks particles across frames, filters out particles that touch image borders, and reorders tracks based on the first frame and the sum of PCP intensity.
    
    Parameters: 
        img: A 4D numpy array with dimensions (frames, channels, height, width) whose channel 1 corresponds to PCP-mCherry and channel 2 
            corresponds to MCP-GFP.
        blobs_df: The pd DataFrame containing detected blobs with measured intensities.
        search_range (int): The maximum distance particles can move between frames.
        memory (int): The maximum number of frames a particle can be lost before being removed.
        threshold (int): Minimum track length in frames to retain a particle.
        border (int): The distance from the image border within which particles are discarded.
    
    Returns: 
        pandas.DataFrame: The filtered and reordered tracking DataFrame.
    """

    # Perform tracking using trackpy 
    print("Performing tracking")
    t = tp.link_df(blobs_df, search_range=search_range, memory=memory)
    
    # Filter out short-lived tracks
    t = tp.filter_stubs(t, threshold=threshold)
    
    # Calculate first frame, outmost positions, and sum of PCP intensity for each tracked particle
    t_grouped = t.groupby('particle').agg(
        first_frame=('frame', 'min'),
        sum_PCP_intensity=('PCP_intensity', 'sum'),
        xmin=('x', 'min'),
        xmax=('x', 'max'),
        ymin=('y', 'min'),
        ymax=('y', 'max')
        ).reset_index()

    # Merge the calculated values back with the original DataFrame
    t_merge = t.merge(t_grouped, on='particle')

    # Keep particles that have not entered image borders
    t_merge = t_merge[(t_merge['xmin'] >= border) & 
                        (t_merge['xmax'] <= img.shape[3] - border -1) & 
                        (t_merge['ymin'] >= border) & 
                        (t_merge['ymax'] <= img.shape[2] - border -1)]

    # Sort by first frame and then by mean PCP intensity
    t_sorted = t_merge.sort_values(by=['first_frame', 'sum_PCP_intensity'], ascending=[True, False])

    # Reassign particle number so that they are pasted to the right position in montage
    unique_particles = list(dict.fromkeys(t_sorted['particle']))
    mapping = {value: i for i, value in enumerate(unique_particles)}
    t_sorted['particle'] = t_sorted['particle'].map(mapping)

    # Drop the temporary columns
    t_sorted = t_sorted.drop(columns=['first_frame', 'sum_PCP_intensity', 'xmin', 'xmax', 'ymin', 'ymax'])

    return t_sorted



def view_blobs_tracks(img, merge_blur, blob_labels, t):
    """
    Views the blob detection and tracking/filtering results in napari
    """
    # Prepare data for layers
    PCP = img[:, 0, :, :]
    MCP = img[:, 1, :, :]
    tracks_data = t[['particle', 'frame', 'y', 'x']].values
    # Create napari viewer
    viewer = napari.Viewer()
    viewer.add_image(PCP, name='PCP', contrast_limits=[0,250], 
                     colormap="magenta", blending="additive")
    viewer.add_image(MCP, name='MCP', contrast_limits=[0,1250], 
                     colormap="green", blending="additive")
    viewer.add_image(merge_blur, name='merge_blur', contrast_limits=[0,0.01], 
                     colormap="gray", blending="additive")
    viewer.add_image(blob_labels, name='blobs', contrast_limits=[0,1], 
                     colormap="yellow", blending="additive")
    viewer.add_tracks(tracks_data, name='Tracks', tail_length=100, 
                      head_length=100, tail_width=2, colormap='hsv')

    return viewer



def make_montage(img, t, radius=8, spacing=5, limit=999):
    """
    Generates a montage of cropped particle images for each track over time.

    Parameters:
        img: A 4D numpy array with dimensions (frames, channels, height, width) whose channel 1 corresponds to PCP-mCherry and channel 2 
            corresponds to MCP-GFP.        
        t: The filtered and reordered tracking DataFrame.
        radius (int): The radius around each particle to include in the montage.
        spacing (int): The distance between cropped tracks in the montage.
        limit (int): The maximum number of tracks to include in the montage.

    Returns:
        montage: The montage image as a 3D numpy array
    """

    r = radius
    d = 2*r + 1

    # Convert t to array and get numbers of frames and particles
    t_converted = np.asarray(t[['frame', 'x', 'y', 'particle']])
    n_frame = int(np.max(t_converted[:,0])) +1 # adding 1 to include frame 0
    n_particle = int(np.max(t_converted[:,3])) +1 # adding 1 to include particle 0

    # Check if particle number exceeds the specified limit
    n_particle = min(limit, n_particle)    
    t_converted = t_converted[t_converted[:, 3] < n_particle]
    

    # Initialize an empty array to store cropped images
    montage_data = np.zeros((2, n_frame*d, n_particle*(d+spacing)), dtype=np.uint16)

    # Loop through each row in t_converted, crop image, and then store in montage_data
    for row in t_converted:
        frame, x, y, particle = row
        frame = int(frame)
        x = int(x)
        y = int(y)
        particle = int(particle)
        
        # Calculate the bounding box for the spot image
        xmin = x - r
        xmax = x + r + 1
        ymin = y - r
        ymax = y + r + 1
        
        # Skip if the spot is at the boundary (i.e., if the bounding box exceeds image dimensions)
        if xmin < 0 or xmax > img.shape[3] or ymin < 0 or ymax > img.shape[2]:
            continue
        
        # Crop images
        PCP_spot_img = img[frame, 0, ymin:ymax, xmin:xmax]
        MCP_spot_img = img[frame, 1, ymin:ymax, xmin:xmax]
        
        # Calculate the bounding box in montage_data array
        row_start = frame * d
        row_end = row_start + d
        col_start = particle * (d + spacing)
        col_end = col_start + d
        
        # Paste cropped images to montage_data
        montage_data[0, row_start:row_end, col_start:col_end] = PCP_spot_img
        montage_data[1, row_start:row_end, col_start:col_end] = MCP_spot_img
    
    return montage_data


def main(file_path, output_dir, 
         min_sigma=0.1, max_sigma=8, threshold=0.00035, 
         radius=8, bg_radius=15, 
         search_range=20, memory=5, track_threshold=2, border=20,
         crop_radius=5, spacing=3, view=False):
    """
    Analyzes a single movie by detecting, tracking, and generating montages of transcriptional foci.
    
    Parameters:
        file_path (str): The path to the TIFF file to be analyzed.
        output_dir (str): The directory where the montages will be saved.
        min_sigma (float): Minimum sigma for blob detection.
        max_sigma (float): Maximum sigma for blob detection.
        threshold (float): Threshold for blob detection.
        radius (int): Radius for measuring foci intensity.
        bg_radius (int): Radius for measuring background intensity.
        search_range (int): Maximum distance for particle linking across frames.
        memory (int): Maximum frames a particle can disappear before being removed.
        track_threshold (int): Minimum length of a track to retain.
        border (int): Distance from image border to filter particles.
        crop_radius (int): Radius for cropping foci in the montage.
        spacing (int): Distance between cropped foci in the montage.
        view (bool): If True, displays the results using napari.

    Returns:
        None
    """

    img, filename = read_single_image(file_path)
    print(f'Start analyzing {filename}')
    merge_blur, blob_labels, blobs_df = detect_blobs(img, min_sigma=min_sigma, max_sigma=max_sigma, threshold=threshold)
    blobs_df = measure_intensities(img, blobs_df, radius=radius, bg_radius=bg_radius)
    t = track_filter(img, blobs_df, search_range=search_range, memory=memory, threshold=track_threshold, border=border)
    montage = make_montage(img, t, radius=crop_radius, spacing=spacing)
    io.imsave(os.path.join(output_dir, f'{filename}_all.tif'), montage, check_contrast=False)
    montage2 = make_montage(img, t, radius=crop_radius, spacing=spacing, limit=50)
    io.imsave(os.path.join(output_dir, f'{filename}_top50.tif'), montage2, check_contrast=False)
    print("Saving montages...")
    
    if view:
        view_blobs_tracks(img, merge_blur, blob_labels, t)
        napari.run()
        print('Close the napari viewer to continue analysis')


# Analyze data

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Analyze MCP PCPmovies")
    parser.add_argument("--input", type=str, required=True, help="path to the data")
    parser.add_argument("--output", type=str, required=True, help="path to the output directory")

    args = parser.parse_args()

    input_dir = args.input
    output_dir = args.output
    file_path_list = sorted(glob.glob(os.path.join(input_dir, '*.tif')))

    for file_path in file_path_list:
        main(file_path, output_dir, view=False)


#######################################################################################################

# Below is the ImageJ/FIJI macro to convert the output

# run("Stack to Hyperstack...", "order=xyczt(default) channels=2 slices=1 frames=1 display=Composite");
# Stack.setChannel(1);
# run("Magenta");
# setMinAndMax(5, 300);
# Stack.setChannel(2);
# setMinAndMax(66, 1500);
# run("Rotate 90 Degrees Left");
# run("Flip Vertically");
