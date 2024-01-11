# Elementary GIA Model with Spatially-Varying ELRA

## Overview
This script implements an Elementary Glacial Isostatic Adjustment (GIA) model with spatially-varying Effective Lithospheric Rigidity Anomalies (ELRA) as implemented in Kori-ULB based on Coulon et al. (2021). 
The script is written in MATLAB.

## Usage
1. **Input Data**: The script requires input data, including topography (H), bathymetry (B), and other parameters. Ensure you have the necessary input data or provide them as required. The given script is provided with the 'Antarctica25km.mat' input file as an example.

2. **Experiment Set-Up**: Modify the `forcing` variable to specify the type of experiment you want to simulate (uniform mass loss, uniform mass gain, or WAIS collapse). Adjust as needed. The WAIS collapse experiment is to be used with the 'Antarctica25km' input file.

3. **Constants and Parameters**: You can modify various constants and parameters, such as gravitational acceleration, densities, Poisson ratio, Earth's radius, and more, according to your need or application.

4. **Run the Script**: Execute the script in MATLAB or another compatible environment. It will iterate through time steps, calculate the geoid changes (if enabled), perform isostatic adjustments, and compute sea-level changes.

5. **Output**: The script generates several output variables, including geoid changes, volume above flotation, potential ocean volume, density corrections, and sea-level changes. These can be used for further analysis or visualization.

## Dependencies
This script relies on MATLAB and its built-in functions for numerical computations. Ensure you have MATLAB installed and set up properly.

## Reference
Coulon, V., Bulthuis, K., Whitehouse, P. L., Sun, S., Haubner, K., Zipf, L., & Pattyn, F. (2021). Contrasting response of West and East Antarctic ice sheets to glacial isostatic adjustment. Journal of Geophysical Research: Earth Surface, 126, e2020JF006003. https://doi.org/10.1029/2020JF006003

## Author
Violaine COULON
violaine.coulon@ulb.be

Feel free to contact the author with any questions or issues related to the script.
