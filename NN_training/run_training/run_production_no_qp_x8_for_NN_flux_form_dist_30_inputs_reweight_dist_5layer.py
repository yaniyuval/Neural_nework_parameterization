import src.ml_train_nn as ml_train_nn
import src.ml_io_no_qp_flux_nn as ml_io_no_qp_flux_nn

import src.ml_load as ml_load
import numpy as np
import pdb

is_cheyenne = True
f_ppi = {'name': 'StandardScaler'} 
o_ppi = {'name': 'StandardScaler'}

do_train = True
max_z=np.Inf
min_samples_leaf = 20
n_trees = 10
max_depth = 27

n_trn_exs = 16000000 
use_rh = False
no_cos = True




flag_dict = dict()
only_plot = False 
scale_level = True 
read_from_whole_data = 1 

flag_dict['do_dqp'] = True
flag_dict['ver_adv_correct'] = True
flag_dict['do_hor_wind_input'] = False
flag_dict['do_ver_wind_input'] = False
flag_dict['do_z_diffusion'] = False 
flag_dict['do_z_diffusion_correction'] = False 

flag_dict['do_q_T_surf_fluxes'] = False 
flag_dict['do_surf_wind']=False 
flag_dict['do_q_surf_fluxes_out']=False 

flag_dict['do_sedimentation'] = True 
flag_dict['do_fall_tend'] = False
flag_dict['do_qp_as_var'] = False 

flag_dict['do_radiation_output'] = False 
flag_dict['rad_level'] = 30 

flag_dict['do_flux'] = False
flag_dict['do_hor_advection'] = False
flag_dict['do_hor_diffusion'] = False

# Input and output
flag_dict['Tin_feature'] = True
flag_dict['qin_feature'] = True
flag_dict['input_upper_lev'] = 30
flag_dict['Tin_z_grad_feature'] = False
flag_dict['qin_z_grad_feature'] = False
flag_dict['predict_tendencies'] = True 
flag_dict['do_qp_diff_corr_to_T']=True 
flag_dict['do_q_T_surf_fluxes_correction'] = False
flag_dict['do_t_strat_correction'] = True   
flag_dict['output_precip'] = False
flag_dict['do_radiation_in_Tz'] = True 
flag_dict['calc_tkz_z'] = False
flag_dict['calc_tkz_z_correction'] = False


flag_dict['resolution'] = 8 
flag_dict['tkz_data'] = True 

flag_dict['tkz_levels'] = 0


flag_dict['Tin_s_diff_feature'] = False
flag_dict['qin_s_diff_feature'] = False
flag_dict['dist_From_eq_in'] = True
flag_dict['T_instead_of_Tabs'] = False


flag_dict['tabs_resolved_init'] = True 
flag_dict['qn_coarse_init'] = True 
flag_dict['qn_resolved_as_var'] = False 
flag_dict['do_zadv_sed_output'] = False
flag_dict['sed_level'] = 26
flag_dict['strat_corr_level'] = 99
flag_dict['qp_coarse_init'] = True

flag_dict['do_surf_flux_hemispheric_symmetric_correction'] = False
flag_dict['Flux_for_TFULL'] = False

flag_dict['do_uv_surf_fluxes_out'] = False
flag_dict['do_momentum_output'] = False
flag_dict['solin_in'] = False
flag_dict['sst_in'] = False


flag_dict['exclusion_flag'] = False
flag_dict['ind1_exc'] = 0
flag_dict['ind2_exc'] = 0

flag_dict['batch_norm'] = False

dx = 12000*flag_dict['resolution']
dy = 12000*flag_dict['resolution']



data_specific_description = ml_io_no_qp_flux_nn.create_specific_data_string_desc(flag_dict)

training_expt1 = 'qobs'+data_specific_description
do_wind_input = False 
do_diffution=False


input_vert_vars = ['Tin','qin','disteq'] 
output_vert_vars = ['Trad_rest','Tadv','qadv','qout','qsed_RESCALED_7epochs_no_drop']
dim1 = 30
dim2 = dim1-1
input_vert_dim = [dim1,dim1,1] 
output_vert_dim = [dim1,dim2,dim2,dim1,dim1]

rewight_outputs = True 
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
                       n_in = sum(input_vert_dim), n_out = sum(output_vert_dim),n_layers=5, output_extreme_flag = False,
                       batch_norm = flag_dict['batch_norm'])

