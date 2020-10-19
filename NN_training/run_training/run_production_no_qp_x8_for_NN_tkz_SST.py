import src.ml_train_nn as ml_train_nn
import src.ml_io_no_qp_flux_nn as ml_io_no_qp_flux_nn

import src.ml_load as ml_load
import numpy as np
import pdb
# run a RF

is_cheyenne = True
# Define preprocessor
f_ppi = {'name': 'StandardScaler'} # At the moment treated only the possibility to put here NoScalar (for RF it doesn't matter)
# f_ppi = {'name': 'StandardScaler'}

# o_ppi = {'name': 'SimpleO'}
o_ppi = {'name': 'StandardScaler'}

do_train = True
#max_z=20000
max_z=np.Inf
min_samples_leaf = 20
n_trees = 10#10
max_depth = 27

n_trn_exs = 16000000 #900000000
#n_trn_exs = 2000
use_rh = False
no_cos = True




flag_dict = dict()
only_plot = False # If I only put to plot....
scale_level = True #If true, each column is scaled seperately.
read_from_whole_data = 1 #In case we want to read the data from a pkl file that has all possible inputs

flag_dict['do_dqp'] = False
flag_dict['ver_adv_correct'] = False
flag_dict['do_hor_wind_input'] = True
flag_dict['do_ver_wind_input'] = False
flag_dict['do_z_diffusion'] = True #True
flag_dict['do_z_diffusion_correction'] = False #True

flag_dict['do_q_T_surf_fluxes'] = False #True  #Only relevant if vertical diffusion is included.
flag_dict['do_surf_wind']= True #True
flag_dict['do_q_surf_fluxes_out']= True #True

flag_dict['do_sedimentation'] = False # if I want to include the sedimentation tendencies (from the cloud scheme- Need to calculate it in matlab)
flag_dict['do_fall_tend'] = False
flag_dict['do_qp_as_var'] = False # If I want to run the simulation with qp as a prognostic parameter.

#Later I can consider doing the fluxes seperately from the radiation and all the micro tendencies and try running a NN with such
#conserving output.
flag_dict['do_radiation_output'] = False # If want to predict the radiation seperately.
flag_dict['rad_level'] = 0 #Should be 0 if no radiation should be used in hte RF

flag_dict['do_flux'] = False
flag_dict['do_hor_advection'] = False
flag_dict['do_hor_diffusion'] = False

# Input and output
flag_dict['Tin_feature'] = True
flag_dict['qin_feature'] = True
flag_dict['input_upper_lev'] = 15
flag_dict['Tin_z_grad_feature'] = False
flag_dict['qin_z_grad_feature'] = False
flag_dict['predict_tendencies'] = False #Usually True
flag_dict['do_qp_diff_corr_to_T']=False #Usually True
flag_dict['do_q_T_surf_fluxes_correction'] = True
flag_dict['do_t_strat_correction'] = False   #Usually True - I wanted not to correct nothing - include convection in stratosphere

flag_dict['output_precip'] = False
flag_dict['do_radiation_in_Tz'] = False #Usually true
flag_dict['calc_tkz_z'] = True
flag_dict['calc_tkz_z_correction'] = False

# flag_dict['calc_tkz_xy'] = False #This is a dummy at the moment

flag_dict['resolution'] = 8 #Usually true
flag_dict['tkz_data'] = True #change to true!! Check if I get the error
# flag_dict['do_dataframe'] = False

flag_dict['tkz_levels'] = 15

flag_dict['Tin_s_diff_feature'] = False
flag_dict['qin_s_diff_feature'] = False
flag_dict['dist_From_eq_in'] = False
flag_dict['T_instead_of_Tabs'] = False


flag_dict['tabs_resolved_init'] = True #Should be true - We only know the resolved tabs as far as I understand
flag_dict['qn_coarse_init'] = True #Should be true - want to have qn coarse from beginning of time step - to calculate qt
flag_dict['qn_resolved_as_var'] = False #Should be true - want to have qn coarse from beginning of time step - to calculate qt
flag_dict['do_zadv_sed_output'] = False
flag_dict['sed_level'] = 0
flag_dict['strat_corr_level'] = 0


flag_dict['do_surf_flux_hemispheric_symmetric_correction'] = True
flag_dict['Flux_for_TFULL'] = False

flag_dict['do_uv_surf_fluxes_out'] = False
flag_dict['do_momentum_output'] = False

flag_dict['solin_in'] = False
flag_dict['sst_in'] = True

flag_dict['exclusion_flag'] = False
flag_dict['ind1_exc'] = 0
flag_dict['ind2_exc'] = 0

flag_dict['batch_norm'] = False

dx = 12000*flag_dict['resolution']
dy = 12000*flag_dict['resolution']



data_specific_description = ml_io_no_qp_flux_nn.create_specific_data_string_desc(flag_dict)

training_expt1 = 'qobs'+data_specific_description
#training_expt2 = 'qobs4K'
do_wind_input = False #Yani added
do_diffution=False


input_vert_vars = ['Tin','qin','uin','vin_minusSH','usurf','SST'] #Dependent on the scenario we chose to model. This should give some flexibility to our coding.
output_vert_vars = ['tsurf','qsurf','tkz']
dim1 = 15
dim2 = dim1-1
input_vert_dim = [dim1,dim1,dim1,dim1,1,1] #Dependent on the scenario we chose to model. This should give some flexibility to our coding.
output_vert_dim = [1,1,dim1]


rewight_outputs = False #If I want to give more wight to certain features.
weight_list = [1,1]

ml_train_nn.train_wrapper(f_ppi, o_ppi, training_expt1, input_vert_dim, output_vert_dim,
                       input_vert_vars, output_vert_vars,flag_dict,
                       max_z=max_z, 
                       do_nn = True,
                       n_trn_exs=n_trn_exs, 
                       plot_training_results=False, 
                       do_train=do_train, 
                       use_rh=use_rh,
                       no_cos=no_cos,
                       do_wind_input=do_wind_input,
                       do_diffusion=do_diffution,
                       scale_level = scale_level,
                       rewight_outputs = rewight_outputs,
                       weight_list = weight_list,
                       is_cheyenne=is_cheyenne,
                       only_plot = only_plot,
                       n_in = sum(input_vert_dim), n_out = sum(output_vert_dim),n_layers=5, output_extreme_flag = True,
                       batch_norm = flag_dict['batch_norm'])


