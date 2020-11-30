# Use of neural networks for stable, accurate and physically consistent parameterization of subgrid atmospheric processes with good performance at reduced precision

Here we have the code and processed data from simulations and neural network parameterizations used in the manuscript ``Use of neural networks for stable, accurate and physically consistent parameterization of subgrid atmospheric processes with good performance at reduced precision''.

### code
The code is divided to three main directories:
1. sam_code_NN: Fortran code with all changes done to SAM used in the simulations (except for the reduced precision simulations).
   - The subdirectory sam_cases contains the namelist (prm), the main.f90 file and neural network subroutines (nn_convection_flux.f90, nn_diffusion.f90) for the runs used in the manuscript. The subdirectories in this directory have the naming convections:
    - run_files_x8_N_layers for data from the x8-NN simulation (96km grid spacing with NN parameterization) with neural-network parameterization (N layers) simulation
    - run_files_N_missing_bits_out_in_only for data from the x8-NN simulation (96km grid spacing with NN parameterization) with neural-network parameterization (5 layers) simulation with 23-N bits in the mantissa
2. NN_training: python code used for creating all Neural Networks used in the manuscript. There are two directories here.
   - run_training:  examples of the input files used to train the NNs (first the files starting with 'build' were run to create the train and test data sets, and later the files starting with 'run’ were run where we trained the neural networks).
   - src: python code to process the coarse-grained high-resolution data and to train the neural networks.
3. high_res_processing_code: matlab code used to calculate the coarse-grained and resolved tendencies, fluxes, diffusivity and input variables. This code uses the high-resolution data to calculate these quantities. 
The high-resolution simulation output and a readme.txt file describing the high-resolution data is found at [this google drive](https://drive.google.com/drive/folders/1TRPDL6JkcLjgTHJL9Ib_Z4XuPyvNVIyY).

### trained neural networks

All the neural networks used in this study are saved in the NNs directory.
The number of layers in each of the neural network is indicated in the file name (“...NN_layers*...”, where the * indicate the number of layers). 
 
### processed data from simulations
The processed data from simulations with neural-network parameterization is found at the data library - data_x8_x16_NN_log.
In this library there are different libraries for different simulations described in the manuscript:
- data_x8_5_layers - data from the x8-NN simulation (96km grid spacing) with neural-network parameterization (5 layers) simulation
- data_x8_4_layers - data from the x8-NN simulation (96km grid spacing) with neural-network parameterization (4 layers) simulation
- data_x8_3_layers - data from the x8-NN simulation (96km grid spacing) with neural-network parameterization (3 layers) simulation
- data_x8_2_layers - data from the x8-NN simulation (96km grid spacing) with neural-network parameterization (2 layers) simulation
- data_16_missing_bits_out_in_only - data from the x8-NN simulation (96km grid spacing) with neural-network parameterization (5 layers) simulation with reduced precision (7 bit in the mantissa)
- data_18_missing_bits_out_in_only - data from the x8-NN simulation (96km grid spacing) with neural-network parameterization (5 layers) simulation with reduced precision (5 bit in the mantissa)
- data_20_missing_bits_out_in_only - data from the x8-NN simulation (96km grid spacing) with neural-network parameterization (5 layers) simulation with reduced precision (3 bit in the mantissa)
- data_22_missing_bits_out_in_only - data from the x8-NN simulation (96km grid spacing) with neural-network parameterization (5 layers) simulation with reduced precision (1 bit in the mantissa)
- data_rf_x8_no_qp - data from the x8(96km grid spacing) with random forest parameterization
- data_x8_5_layers_solin - data from the x8-NN-solin simulation (96km grid spacing) with neural-network parameterization (5 layers) simulation with solar insolation as feature instead of distance from the equator
- data_x16_5_layers - data from the x16-NN (192km grid spacing) with neural-network parameterization (5 layers) simulation

Each folder contains a netcdf file with the following data:
- x - longitudinal coordinate (units:meters)
- y - latitudinal coordinate (units:meters)
- z(z) - vertical coordinate (units:meters)

The time and zonal average (taken from 3-hourly snapshot over 500 days):
- U(y, z) - zonal wind
- V(y, z) - meridional wind
- W(y, z) - vertical wind
- U2(y, z) - zonal wind squared
- V2(y, z) - meridional wind squared
- W2(y, z) - vertical wind squared
- T(y, z) - temperature
- QV(y, z) - water vapor
- QC(y, z) - cloud water
- QI(y, z) - cloud ice
- QP(y, z) - precipitable water
- RH(y, z) - PLEASE IGNORE
- P(z) - PLEASE IGNORE
- y_coarse - PLEASE IGNORE
- y_sim - PLEASE IGNORE

Precipitation averaged over 500 days:
- precip_avg_coarse - mean precipitation from the coarse-grained high-resolution data
- precip_avg_low - - mean precipitation from the low-resolution data (x8 for x8 simulations, x16 for x16 simulaitons)
- precip_avg_rf - mean precipitation from the low-resolution data with Neural network (OR random forest) parameterization (x8 for x8 simulations, x16 for x16 simulaitons)
- precip_xtreme_coarse -PLEASE IGNORE
- precip_xtreme_low - PLEASE IGNORE
- precip_xtreme_rf - PLEASE IGNORE

Precipitation frequency:

- bin_dim_full(bin_dim_full) - the bins used for calculating precipitation rate frequency (not used in the manuscript - bin every 1 mm/day)
- precip_dist_full(bin_dim_full) - the number of precipitation rate events in each bin for the neural network (random forest) simulation (not used in the manuscript- bin every 1 mm/day)
- precip_dist_trop_full(bin_dim_full) - the number of precipitation rate events (only in the tropics) in each bin for the neural network (random forest) simulation (not used in the manuscript- bin every 1 mm/day)
- precip_dist_ref(bin_dim_full) - the number of precipitation rate events in each bin for the low-resolution simulation (same resolution as the simulation with neural-network parameterization; not used in the manuscript- bin every 1 mm/day)
- precip_dist_trop_ref(bin_dim_full) - the number of precipitation rate events (only in the tropics) in each bin for the low-resolution simulation (same resolution as the simulation with neural-network parameterization; not used in the manuscript- bin every 1 mm/day)
- precip_dist_ref_high(bin_dim_full) - the number of precipitation rate events in each bin for the coarse-grained high-resolution simulation (not used in the manuscript- bin every 1 mm/day). The high resolution data was first coarse-grained to the same resolution that the simulations with neural-network parameterization were ran at. 
- precip_dist_trop_ref_high(bin_dim_full) - the number of precipitation rate events (only in the tropics) in each bin for the coarse-grained high-resolution simulation (not used in the manuscript- bin every 1 mm/day)
- bin_dim_full_log(bin_dim_full_log) - the bins used for calculating precipitation distribution (when precipitation>1mm/day we use bins equally spaced in the logarithm of precipitation rate).
- precip_dist_full_log(bin_dim_full_log) - the number of precipitation rate events in each bin for the neural network (random forest) simulation
- precip_dist_trop_full_log(bin_dim_full_log) - the number of precipitation rate events (only in the tropics) in each bin for the neural network (random fores) simulation
- precip_dist_ref_log(bin_dim_full_log) - the number of precipitation rate events in each bin for the low-resolution simulation (same resolution as the neural network simulation)
- precip_dist_trop_ref_log(bin_dim_full_log) - the number of prcipitation rate events (only in the tropics) in each bin for the low-resolution simulation (same resolution as the simulation with neural-network parameterization)
- precip_dist_ref_high_log(bin_dim_full_log) - the number of precipitation rate events in each bin for the coarse-grained high-resolution simulation
- precip_dist_trop_ref_high_log(bin_dim_full_log) - the number of precipitation rate events (only in the tropics) in each bin for the coarse-grained high-resolution simulation





