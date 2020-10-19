import src.ml_io_no_qp_flux_nn as ml_io_no_qp_flux_nn
import numpy as np
import src.ml_io_separate_output as sep_out


training_expt1 = 'qobs'
training_expt2 = 'qobs4K'


start_time =  1732500 #669150
end_time = 2520000# 2262600

start_time = 360000  + 1000*450 #9000
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

rewight_outputs = False

dx = 12000*flag_dict['resolution']
dy = 12000*flag_dict['resolution']




do_shuffle = False # Yani added so I will not suffle between test and train - this simulates more the real situation when we run the model

ml_io_no_qp_flux_nn.build_training_dataset(training_expt1, start_time, end_time, interval, n_x_samp=n_x_samp, train_size=train_size,
                             do_shuffle=do_shuffle, flag_dict =flag_dict,is_cheyenne=is_cheyenne  )


# sep_out.build_training_dataset_separate_out(training_expt1, start_time, end_time, interval, n_x_samp=n_x_samp, train_size=train_size,
#                              do_shuffle=do_shuffle, flag_dict =flag_dict,is_cheyenne=is_cheyenne  )
#
# ll = ['Tin', 'qin', 'T_grad', 'q_grad','uin','vin','win','usurf','Tout','qout','Tout_adv','qout_adv','Tout_adv_x','qout_adv_x','Tout_adv_y','qout_adv_y','Tout_mic','qout_mic','Tout_dif','qout_dif','Tout_dif_x','qout_dif_x','Tout_dif_y','qout_dif_y','Tout_dif_surf_cor','qout_dif_surf_cor','Tout_strat_cor','qout_strat_cor','Tout_rad','Tout_cond','qSurfout','qSurfout_correction','precip_out','y', 'z', 'p', 'rho']


# ml_io.build_training_dataset(training_expt1, start_time, end_time, interval, n_x_samp=n_x_samp, train_size=train_size, do_wind_input=do_wind_input,
#                              do_shuffle=do_shuffle, do_z_diffusion=do_z_diffusion, do_q_T_surf_fluxes=do_q_T_surf_fluxes, do_surf_wind=do_surf_wind,
#                              do_q_surf_fluxes_out=do_q_surf_fluxes_out,
#                              do_sedimentation = do_sedimentation,
#                              do_radiation_output = do_radiation_output,
#                              do_fall_tend = do_fall_tend,
#                              do_qp_as_var = do_qp_as_var,
#                              Tin_feature=Tin_feature,
#                              qin_feature = qin_feature,
#                              Tin_z_grad_feature = Tin_z_grad_feature,
#                              qin_z_grad_feature = qin_z_grad_feature,
#                              predict_tendencies = predict_tendencies,
#                              do_flux = do_flux,
#                              do_hor_advection = do_hor_advection,
#                              do_hor_diffusion=do_hor_diff,
#                              do_qp_diff_corr_to_T=do_qp_diff_corr_to_T,
#                              do_q_T_surf_fluxes_correction=do_q_T_surf_fluxes_correction)
#ml_io.build_training_dataset(training_expt2, start_time, end_time, interval, n_x_samp=n_x_samp)
