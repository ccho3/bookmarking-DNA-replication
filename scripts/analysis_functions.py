"""Functions for reading images and simple quantifications"""

import os, glob, argparse
import numpy as np
from skimage import io, measure
from skimage.segmentation import clear_border
from skimage.morphology import remove_small_objects

def get_dirs():
    """Get input and output directories from command line inputs"""
    parser = argparse.ArgumentParser(description="Analyze Rpb3 MCP movies")
    parser.add_argument("--input", type=str, required=True, help="path to the data")
    parser.add_argument("--output", type=str, required=True, help="path to the output directory")
    args = parser.parse_args()
    input_dir = args.input
    output_dir = args.output
    return input_dir, output_dir


def read_ungrouped_images(file_path):
    """
    Read TIFF files from a directory without group labels.
    
    Parameters:
        file_path (str): The file path for TIFF files.

    Returns:
        img_list (list): A list of images stored as numpy arrays.
    """

    img_list=[]
    file_path_list = glob.glob(os.path.join(file_path, '*.tif'))

    for file in sorted(file_path_list):
        a = io.imread(file)
        img_list.append(a)
        filename = os.path.basename(file).split('.')[0]
        print(f'Found movie {filename}')

    return img_list


def read_grouped_images(file_path, groups):
    """
    Reads TIF files from the specified file path and groups them according to the provided labels.
    
    Parameters:
        file_path (str): The base file path  for TIF files, e.g., "path/to/data/", which contain TIF files with group headers.
        groups (list of str): A list of group labels to categorize the images.

    Returns:
        imgs (list of lists): A list where each element is a list of images belonging to a particular group.

    Example:
        imgs = read_grouped_images("data/", ["M10", "M11"])
        This will look for files like "data/M10*.tif" and "data/M11*.tif".    
    """

    imgs=[]

    for group in groups:
        group_path = os.path.join(file_path, f'{group}*.tif')
        group_files = sorted(glob.glob(group_path))

        group_imgs = []
        for file in group_files:
            filename = os.path.basename(file).split('.')[0]
            a = io.imread(file)
            group_imgs.append(a)
            print(f'Found movie {filename}')
        imgs.append(group_imgs)

    return imgs

def read_single_image(file_path):
    """
    Reads a TIFF file into a numpy array.

    Parameters:
        file_path (str): The path to the TIFF file.

    Returns:
        img (numpy array): The movie data as a numpy array.
        
        filename (str): The base name of the file (without extension).
    """

    img = io.imread(file_path)
    filename = os.path.basename(file_path).split('.')[0]
    print(f'Found movie {filename}')
    
    return img, filename


def mean_var_per_image(img, probe_ch, mask_ch):
    """
    Calculate the mean intensity and variance for two channels of an image within a masked region.

    Parameters:
        img: The image data as a 3D numpy array with the shape of (width, height, channels)
        probe_ch (list or tuple): A list or tuple containing two indices for the probe channels.
        mask_ch (int): The index of the channel to be used as the mask.

    Returns:
        C1_mean: Mean intensity of the first probe channel within the masked region, background subtracted.
        C1_var: Variance of the first probe channel within the masked region.
        C2_mean: Mean intensity of the second probe channel within the masked region, background subtracted.
        C2_var: Variance of the second probe channel within the masked region.
    """

    C1 = img[:,:,probe_ch[0]]
    C2 = img[:,:,probe_ch[1]]
    mask = img[:,:,mask_ch] >0

    C1_mask = np.mean(C1[mask])
    C1_bg = np.mean(C1[~mask])
    C1_mean = C1_mask - C1_bg
    C1_var = np.var(C1[mask])
    
    C2_mask = np.mean(C2[mask])   
    C2_bg = np.mean(C2[~mask])
    C2_mean = C2_mask - C2_bg
    C2_var = np.var(C2[mask])

    return C1_mean, C1_var, C2_mean, C2_var


def mean_var_per_nucleus(img, probe_ch, mask_ch):
    """
    Calculate the mean intensity and variance for two channels of individual nuclei in an image within a masked region.

    Parameters:
        img: The image data as a 3D numpy array with the shape of (width, height, channels)
        probe_ch (list or tuple): A list or tuple containing two indices for the probe channels.
        mask_ch (int): The index of the channel to be used as the mask.

    Returns:
        results (list): a list with the following measurements for all nuclei
            C1_mean: Mean intensity of the first probe channel within the masked region, background subtracted.
            C1_var: Variance of the first probe channel within the masked region.
            C2_mean: Mean intensity of the second probe channel within the masked region, background subtracted.
            C2_var: Variance of the second probe channel within the masked region.
    """
    C1 = img[:,:,probe_ch[0]]
    C2 = img[:,:,probe_ch[1]]    
    mask = img[:,:,mask_ch] > 0
    
    C1_bg = np.mean(C1[~mask])
    C2_bg = np.mean(C2[~mask])
    
    mask = clear_border(mask)
    mask = remove_small_objects(mask, min_size=300)

    mask_labels = measure.label(mask)
    label_id = np.unique(mask_labels)[1:]

    results = []    
    for i in label_id:
        C1_mean = C1[mask_labels == i].mean() - C1_bg
        C2_mean = C2[mask_labels == i].mean() - C2_bg
        C1_var = C1[mask_labels == i].var()
        C2_var = C2[mask_labels == i].var()
        result = (C1_mean, C1_var, C2_mean, C2_var)
        results.append(result)

    return results


def max_per_nucleus(img, probe_ch, mask_ch):
    """
    Calculate the maximal intensity for two channels of individual nuclei in an image within a masked region.

    Parameters:
        img: The image data as a 3D numpy array with the shape of (width, height, channels)
        probe_ch (list or tuple): A list or tuple containing two indices for the probe channels.
        mask_ch (int): The index of the channel to be used as the mask.

    Returns:
        results (list): a list with the measurement of maximal intensities in channels 1 and 2 for all nuclei
    """
    C1 = img[:,:,probe_ch[0]]
    C2 = img[:,:,probe_ch[1]]
    mask = img[:,:,mask_ch] >0

    mask = clear_border(mask)
    mask = remove_small_objects(mask, min_size=300)

    mask_labels = measure.label(mask)
    label_id = np.unique(mask_labels)[1:]

    results = []
    for i in label_id:
        C1_max = C1[mask_labels == i].max()
        C2_max = C2[mask_labels == i].max()
        result = (C1_max, C2_max)
        results.append(result)
        
    return results

def max_per_image(img, probe_ch, mask_ch):
    """
    Calculate the average of the maximum intensity in nuclei for two channels of an image.

    Parameters:
        img: The image data as a 3D numpy array with the shape of (width, height, channels)
        probe_ch (list or tuple): A list or tuple containing two indices for the probe channels.
        mask_ch (int): The index of the channel to be used as the mask.

    Returns:
        C1_max_average: The average of maximum intensity of the first probe channel across nuclei.
        C2_max_average: The average of maximum intensity of the second probe channel across nuclei
    """
    C1 = img[:,:,probe_ch[0]]
    C2 = img[:,:,probe_ch[1]]
    mask = img[:,:,mask_ch] >0

    mask = clear_border(mask)
    mask = remove_small_objects(mask, min_size=300)

    mask_labels = measure.label(mask)
    label_id = np.unique(mask_labels)[1:]
    
    C1_max_results = []
    C2_max_results = []
    for i in label_id:
        C1_max = C1[mask_labels == i].max()
        C2_max = C2[mask_labels == i].max()

        C1_max_results.append(C1_max)
        C2_max_results.append(C2_max)

        C1_max_average= np.asarray(C1_max_results).mean()
        C2_max_average= np.asarray(C2_max_results).mean()
    
    return C1_max_average, C2_max_average