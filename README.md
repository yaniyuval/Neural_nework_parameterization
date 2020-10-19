# Use of neural networks for stable, accurate and physically consistent parameterization of subgrid atmospheric processes with good performance at reduced precision

The codes are divided to four main directories:
- sam_code_NN: Fortran code with all changes done to SAM used in the simulations (except for the reduced predcision simulations). 
- sam_code_NN_reduced_precision: Fortran code with all changes done to SAM when used reduced precision NNs. 
- NN_training: python code used for creating all Neural Networks used in the manuscript. There are two directories here. In the run_training direcory there are examples of the inputs used to train the NNs (first the files starting with 'build' were run, and later the files starting with 'run'). in src directory there are the python code. 
- high_res_processing_code: matlab code used to calculate the coarse-grained and resolved tendencies 
# Neural_nework_parameterization
