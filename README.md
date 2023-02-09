# Software package and datasets for "Towards quantitative super-resolution microscopy: Molecular maps with statistical guarantees"

## System requirements and code license

This software package contains examplary data, pre-computed MISCAT tests and
analysis scripts in Matlab (https://www.mathworks.com/products/matlab.html).
It has been successfully tested with Matlab R2020b.

The source code of the software is licensed under the MIT license (see file LICENSE).

The data is provided for academic and visualization purposes only.

## Setup

No specific installation is necessary. Execute "initialize.m" at least
once before running other scripts to include all paths.

## Examples

Make sure to run initialize.m at least once. 

Call examples/example_xxx.m to see the software estimating a hybrid segmentation,
  and confidence intervalls on the number of markers in each segment for experimental
  or simulated photon antibunching microscopy data.

## Usage

Make sure to run initialize.m at least once. 

The main routine is source/estimate_molecular_map.m which by itself
- performs a MISCAT segmentation test (source/miscat/*.m)
- estimates the background (quite simple)
- preprocesses the data (source/inference/preprocessing/preprocess_data.m)
- merges MISCAT results with another segmentation method
  (source/inference/boxes_to_segments/boxes_to_segments.m)
- estimate n,p for every obtained segment including confidence intervals
  (source\inference\estimate_np\estimate_np.m)

