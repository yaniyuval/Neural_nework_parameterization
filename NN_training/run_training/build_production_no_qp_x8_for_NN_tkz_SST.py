import src.ml_io_no_qp_flux_nn as ml_io_no_qp_flux_nn
import numpy as np
import src.ml_io_separate_output as sep_out


training_expt1 = 'qobs'
training_expt2 = 'qobs4K'

start_time = 360000  + 1000*450 
end_time = start_time +2700 *450
interval = 450
train_size = 0.95
n_x_samp = 40
is_cheyenne =True

flag_dict = dict()


flag_dict['do_dqp'] = False
flag_dict['ver_adv_correct'] = False
flag_dict['do_hor_wind_input'] = True
flag_dict['do_ver_wind_input'] = False
flag_dict['do_z_diffusion'] = True 
flag_dict['do_z_diffusion_correction'] = False 

flag_dict['do_q_T_surf_fluxes'] = False 
flag_dict['do_surf_wind']= True 
flag_dict['do_q_surf_fluxes_out']= True 

flag_dict['do_sedimentation'] = False 
flag_dict['do_fall_tend'] = False
flag_dict['do_qp_as_var'] = False 

flag_dict['do_radiation_output'] = False 
flag_dict['rad_level'] = 0 

flag_dict['do_flux'] = False
flag_dict['do_hor_advection'] = False
flag_dict['do_hor_diffusion'] = False

flag_dict['Tin_feature'] = True
flag_dict['qin_feature'] = True
flag_dict['input_upper_lev'] = 15
flag_dict['Tin_z_grad_feature'] = False
flag_dict['qin_z_grad_feature'] = False
flag_dict['predict_tendencies'] = False 
flag_dict['do_qp_diff_corr_to_T']=False 
flag_dict['do_q_T_surf_fluxes_correction'] = True
flag_dict['do_t_strat_correction'] = False   
flag_dict['output_precip'] = False
flag_dict['do_radiation_in_Tz'] = False 
flag_dict['calc_tkz_z'] = True
flag_dict['calc_tkz_z_correction'] = False

flag_dict['resolution'] = 8 
flag_dict['tkz_data'] = True 

flag_dict['tkz_levels'] = 15

flag_dict['Tin_s_diff_feature'] = False
flag_dict['qin_s_diff_feature'] = False
flag_dict['dist_From_eq_in'] = False
flag_dict['T_instead_of_Tabs'] = False


flag_dict['tabs_resolved_init'] = True 
flag_dict['qn_coarse_init'] = True 
flag_dict['qn_resolved_as_var'] = False 
flag_dict['do_zadv_sed_output'] = False
flag_dict['sed_level'] = 0
flag_dict['strat_corr_level'] = 0 


flag_dict['do_surf_flux_hemispheric_symmetric_correction'] = True
flag_dict['Flux_for_TFULL'] = False

flag_dict['do_uv_surf_fluxes_out'] = False
flag_dict['do_momentum_output'] = False

flag_dict['solin_in'] = False
flag_dict['sst_in'] = True

rewight_outputs = False

dx = 12000*flag_dict['resolution']
dy = 12000*flag_dict['resolution']

do_shuffle = False

ml_io_no_qp_flux_nn.build_training_dataset(training_expt1, start_time, end_time, interval, n_x_samp=n_x_samp, train_size=train_size,
                             do_shuffle=do_shuffle, flag_dict =flag_dict,is_cheyenne=is_cheyenne  )

