# bookmarking-DNA-replication

This repository contains scripts and quantification results related to the manuscript **"A mitotic bookmark promotes early replication."** The scripts provide custom tools for analyzing fluorescent protein dynamics in early Drosophila embryos.

## Usage

The datasets generated for the manuscript and used for the following analyses are available at [Zenodo](https://zenodo.org/records/13761202).

To analyze the data, use the provided scripts by opening a terminal and running the corresponding commands. Replace `input_directory` and `output_directory` with the paths to the directories containing  TIFF files and the location where output should be saved, respectively. For more information, refer to the docstrings in each script.

### Analysis 1
This script analyzes `Dataset1` and quantifies the mean and variance of fluorescent intensities for mKate2-Brd4 and Cdc7-EGFP in movies of S-phase entry in nuclear cycles 11-14.

```bash
python analysis1.py --input <input_directory> --output <output_directory>
```

### Analysis 2

This script analyzes `Dataset2` and quantifies the mean and variance of fluorescent intensities for mKate2-Brd4 and Cdc7-EGFP in snapshots of embryos subjected to various CBP RNAi or JQ1 injection.

```bash
python analysis2.py --input <input_directory> --output <output_directory>
```

### Analysis 3

This script analyzes `Dataset3` and quantifies the maximum intensity of mKate2-Brd4 and Cdc7-EGFP in snapshots of Cdc7i-injected embryos.

```bash
python analysis1.py --input <input_directory> --output <output_directory>
```

### Analysis 4

This script analyzes `Dataset4` (movies of EGFP-Rpb3 and MCP-mCherry). It performs segmentation of nuclei and MCP foci, tracks nuclei across frames, and quantifies MCP state transitions.

The analysis involves two steps. 

1. Run the initilization command:

```bash
python analysis4.py --input <input_directory> --output <output_directory> --run initialize
```

2. Run the finishing command:
```bash
python analysis4.py --input <input_directory> --output <output_directory> --run finish
```

### Analysis 5

This script analyzes `Dataset5` (movies of MCP-GFP and mCherry-PCP in embryos carrying the hbP2-MS2-lacZ-PP7 reporter). It detects MCP/PCP foci, tracks them across frames, and then generates montages displaying cropped foci images for all tracks or the first 50 tracks in each movie.

```bash
python analysis5.py --input <input_directory> --output <output_directory>
```

### Analysis 6

This script analyzes `Dataset6` and quantifies the maximum intensity of EGFP-Rpb3 in control or Cdc7i-injected NC13 embryos.

```bash
python analysis6.py --input <input_directory> --output <output_directory>
```

### Data Visualization

The output from each analysis was used to create plots using the provided R scripts.

## Dependencies

- Python 3.12.4
- scipy 1.14.0
- numpy 1.26.4 
- pandas 2.2.2 
- scikit-image 0.24.0
- trackpy 0.6.4
- napari 0.5.1
- R 4.3.2
- tidyverse 2.0.0
