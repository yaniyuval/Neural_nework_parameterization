import numpy as np
import matplotlib
import matplotlib.gridspec as gridspec

matplotlib.use('Agg')  # so figs just print to file. Needs to come before pyplot
import matplotlib.pyplot as plt
import scipy.stats
from sklearn import metrics
import pickle
import os
import sys
import src.ml_load as ml_load
import src.atmos_physics as atmos_physics
import warnings
from netCDF4 import Dataset

warnings.filterwarnings("ignore", category=DeprecationWarning)

# unpack routines
unpack_f = ml_load.unpack_f
unpack_f_extended = ml_load.unpack_f_extended

unpack_o = ml_load.unpack_o
pack_o = ml_load.pack_o

# plotting defaults
matplotlib.rcParams['agg.path.chunksize'] = 10000
matplotlib.rcParams.update({'font.size': 8})

# convert tendencies from K/s, kg/kg/s and kg/m^2/s to K/day, g/kg/day and mm/day for plots with unit labels
# note for precipitation kg/m^2 goes to mm
per_day = 3600 * 24
kg_to_g = 1000

# number of points for precipitation scatter plots
num_scatter = 10000


# ---   META PLOTTING SCRIPTS  --- #


def PlotAllFigs_nn(path_nn_stage2, est_str, datafile, do_nn, figpath, input_vert_vars, output_vert_vars, input_vert_dim, output_vert_dim,
                n_trn_exs=None,
                rain_only=False, no_cos=False, use_rh=False, wind_input=False, scale_per_column=False,
                rewight_outputs=False, weight_list=[1, 1], is_cheyenne=False, exclusion_flag=False, ind1_exc=0,
                ind2_exc=0):
    # Open the estimator and the preprocessing scheme
    base_dir = '/glade/scratch/janniy/'


    # Check how to upload my network...
    est_eval, _, errors, f_ppi, o_ppi, f_pp, o_pp, y, z, p, rho = \
        pickle.load(open(base_dir + 'mldata_tmp/regressors/' + est_str + '.pkl', 'rb'))


    # Load the data from the training/testing file
    f_scl, opred_scl, otrue_scl, f, opred, otrue = \
        ml_load.get_f_o_pred_true(est_str, datafile, max_z=max(z), input_vert_vars=input_vert_vars,
                                  output_vert_vars=output_vert_vars,
                                  input_vert_dim=input_vert_dim, output_vert_dim=output_vert_dim, n_trn_exs=n_trn_exs,
                                  rain_only=rain_only,
                                  no_cos=no_cos, use_rh=use_rh, wind_input=wind_input,
                                  scale_per_column=scale_per_column,
                                  rewight_outputs=rewight_outputs, weight_list=weight_list, is_cheyenne=is_cheyenne,do_nn = do_nn)

    print('Size of global dataset for plots after subsampling: ', f_scl.shape[0])

    # Do plotting
    print('Beginning to make plots and writing log files')

    # Create directory if it does not exist
    if not os.path.exists(figpath):
        os.makedirs(figpath)

    # text_file.close()

    ncfile = Dataset(figpath + 'data_test.nc', "w", format="NETCDF3_CLASSIC")
    ncfile.createDimension("tot_dim_column", sum(output_vert_dim))
    ncfile.createDimension("test_samples", None)
    nc_o_true = ncfile.createVariable("o_true", np.dtype('float32').char, ("test_samples", "tot_dim_column"))
    nc_o_pred = ncfile.createVariable("o_pred", np.dtype('float32').char, ("test_samples", "tot_dim_column"))
    nc_o_true[:] = otrue[:]
    nc_o_pred[:] = opred[:]
    ncfile.description = 'blabla-janniy'
    ncfile.close

    # Plot means and standard deviations
    plot_means_stds_generalized(otrue, opred, p, figpath, output_vert_vars, output_vert_dim)

    print('plotted means_stds')

    # Plot correlation coefficient, explained variance, and rmse
    plot_error_stats(otrue, opred, p, figpath, output_vert_vars, output_vert_dim)
    # Plot a scatter plot of true vs predicted precip
    print('plotted error_stats')
    if 'qout' in output_vert_vars:
        plot_scatter(f[:num_scatter], otrue[:num_scatter], opred[:num_scatter], z, p, rho, figpath, output_vert_vars,
                     output_vert_dim)
    # Plot residual from energy conservation
    print('plotted scatter precip')

    if ('Tin' in input_vert_vars and 'Tout' in output_vert_vars and 'qout' in output_vert_vars):
        plot_energy_conservation(f, otrue, opred, z, p, rho, figpath, input_vert_vars, input_vert_dim, output_vert_vars,
                                 output_vert_dim)

    print('Beginning to make y-z and y plots...')
    make_yz_plots(figpath, f_ppi, o_ppi, f_pp, o_pp, est_eval, y, z, p, rho, input_vert_vars, output_vert_vars,
                  input_vert_dim, output_vert_dim,
                  datafile, n_trn_exs, rain_only, no_cos, use_rh, wind_input=wind_input,
                  scale_per_column=scale_per_column,
                  rewight_outputs=rewight_outputs, weight_list=weight_list,do_nn=do_nn)
    print('Done with y-z and y plots.')
    print('Not making at the moment the importance plots...')



def plot_importance(figpath, est_eval, z, p, rho, use_rh, wind_input):
    # weighted averaged of decrease in variance
    importance = est_eval.feature_importances_

    # simple estimate of standard error of the mean of feature importance across trees
    estimators = est_eval.estimators_
    stderr = np.std([tree.feature_importances_ for tree in estimators], axis=0) / np.sqrt(len(estimators))

    # account for variations in interlevel spacing and rho
    rho_dz = atmos_physics.vertical_diff(rho, z)
    rescale = 1.0 / rho_dz

    fig = plt.figure(figsize=(3.0, 2.25))

    plt.plot(p * 0, p, 'k:', lw=0.5)

    plt.errorbar(rescale * unpack_f_extended(importance, 'T', axis=0, wind_input=wind_input), p,
                 xerr=rescale * unpack_f_extended(stderr, 'T', axis=0, wind_input=wind_input),
                 label='T')

    if use_rh:
        plt.errorbar(rescale * unpack_f_extended(importance, 'q', axis=0, wind_input=wind_input), p,
                     xerr=rescale * unpack_f_extended(stderr, 'q', axis=0, wind_input=wind_input),
                     label='rh')
    else:
        plt.errorbar(rescale * unpack_f_extended(importance, 'q', axis=0, wind_input=wind_input), p,
                     xerr=rescale * unpack_f_extended(stderr, 'q', axis=0, wind_input=wind_input),
                     label='q')

    plt.ylim([np.amax(p), np.amin(p)])
    plt.ylabel('Pressure (hPa)')
    plt.xlabel('Feature importance')
    plt.legend(loc="upper right")
    plt.legend(frameon=False)
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    plt.tight_layout()  # avoid overlap
    fig.savefig(figpath + 'feature_importance.eps', bbox_inches='tight')
    plt.close()

    text_file = open(figpath + 'feature_importance.txt', 'w')
    text_file.write(
        'sum of importance for q: %f \n' % np.sum(unpack_f_extended(importance, 'q', axis=0, wind_input=wind_input)))
    text_file.write(
        'sum of importance for T: %f \n' % np.sum(unpack_f_extended(importance, 'T', axis=0, wind_input=wind_input)))
    text_file.close()


def plot_importance_permute(figpath, f_scl, otrue_scl, est_eval, z, p, rho, use_rh, wind_input):
    acc = metrics.r2_score(otrue_scl, est_eval.predict(f_scl))
    importance = np.zeros(f_scl.shape[1])
    for i in range(f_scl.shape[1]):
        f_scl_shuff = f_scl.copy()
        np.random.shuffle(f_scl_shuff[:, i])
        shuff_acc = metrics.r2_score(otrue_scl, est_eval.predict(f_scl_shuff))
        importance[i] = (acc - shuff_acc) / acc

    # account for variations in interlevel spacing
    rescale = 1.0 / (rho * np.gradient(z))

    fig = plt.figure(figsize=(3.0, 2.25))

    plt.plot(p * 0, p, 'k:', lw=0.5)

    plt.plot(rescale * unpack_f_extended(importance, 'T', axis=0, wind_input=wind_input), p, label='T')

    if use_rh:
        plt.plot(rescale * unpack_f_extended(importance, 'q', axis=0, wind_input=wind_input), p, label='rh')
    else:
        plt.plot(rescale * unpack_f_extended(importance, 'q', axis=0, wind_input=wind_input), p, label='q')

    plt.ylim([np.amax(p), np.amin(p)])
    plt.ylabel('Pressure (hPa)')
    plt.xlabel('Feature importance')
    plt.legend(loc="upper right")
    plt.legend(frameon=False)
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    plt.tight_layout()  # avoid overlap
    fig.savefig(figpath + 'feature_importance_permute.eps', bbox_inches='tight')
    plt.close()


def plot_importance_precip_permute(figpath, f_scl, otrue_scl, est_eval, z, p, rho, f_ppi, f_pp, o_ppi, o_pp, use_rh,
                                   wind_input):
    otrue = ml_load.inverse_transform_data(o_ppi, o_pp, otrue_scl, z)
    precip_true = atmos_physics.calc_precip(unpack_o(otrue, 'q'),
                                            rho, z)

    opred_scl = est_eval.predict(f_scl)
    opred = ml_load.inverse_transform_data(o_ppi, o_pp, opred_scl, z)
    precip_pred = atmos_physics.calc_precip(unpack_o(opred, 'q'), rho, z)

    acc = metrics.r2_score(precip_true, precip_pred)

    importance = np.zeros(f_scl.shape[1])
    for i in range(f_scl.shape[1]):
        f_scl_shuff = f_scl.copy()
        np.random.shuffle(f_scl_shuff[:, i])
        o_scl_shuff = est_eval.predict(f_scl_shuff)
        o_shuff = ml_load.inverse_transform_data(o_ppi, o_pp, o_scl_shuff, z)
        precip_shuff = atmos_physics.calc_precip(unpack_o(o_shuff, 'q'),
                                                 rho, z)
        shuff_acc = metrics.r2_score(precip_true, precip_shuff)
        importance[i] = (acc - shuff_acc) / acc

    # account for variations in interlevel spacing
    rescale = 1.0 / (rho * np.gradient(z))

    # use sqrt to compare with linear_response

    fig = plt.figure(figsize=(3.0, 2.25))

    plt.plot(p * 0, p, 'k:', lw=0.5)

    plt.plot(rescale * unpack_f_extended(importance, 'T', axis=0, wind_input=wind_input), p, label='T')

    if use_rh:
        plt.plot(rescale * unpack_f_extended(importance, 'q', axis=0, wind_input=wind_input), p, label='rh')
    else:
        plt.plot(rescale * unpack_f_extended(importance, 'q', axis=0, wind_input=wind_input), p, label='q')

    plt.ylim([np.amax(p), np.amin(p)])
    plt.ylabel('Pressure (hPa)')
    plt.xlabel('Feature importance')
    plt.legend(loc="upper right")
    plt.legend(frameon=False)
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    plt.tight_layout()  # avoid overlap
    fig.savefig(figpath + 'feature_importance_precip_permute.eps', bbox_inches='tight')
    plt.close()


def plot_linear_response(figpath, f, est_eval, z, p, rho, f_ppi, f_pp, o_ppi, o_pp, use_rh=False, wind_input=False):
    # uses central difference of total width 1K, 1 g/kg
    if wind_input:
        q_t_size = int(48 * 2)
    else:
        q_t_size = f.shape[1]

    sensitivity = np.zeros(q_t_size)
    sensitivity_stderr = np.zeros(q_t_size)
    sensitivity_check_linear = np.zeros(q_t_size)
    sensitivity_stderr_check_linear = np.zeros(q_t_size)
    perturb = np.zeros(q_t_size)

    N = int(q_t_size / 2)
    perturb[0:N] = 1  # 1 K
    if use_rh:
        perturb[N:2 * N] = 0.05  # 5% RH
    else:
        perturb[N:2 * N] = 0.05  # 5% of qsat

    p_Pa = p * 100  # p is in hPa
    qsat = atmos_physics.sam_qsat(unpack_f_extended(f, 'T', wind_input=wind_input), p_Pa[None, :])

    for i in range(q_t_size):

        f_perturb_pos = f.copy()
        f_perturb_neg = f.copy()

        if i < N:
            f_perturb_pos[:, i] = f_perturb_pos[:, i] + 0.5 * perturb[i]
            f_perturb_neg[:, i] = f_perturb_neg[:, i] - 0.5 * perturb[i]
        else:
            f_perturb_pos[:, i] = f_perturb_pos[:, i] + 0.5 * perturb[i] * qsat[:, i - N]
            f_perturb_neg[:, i] = f_perturb_neg[:, i] - 0.5 * perturb[i] * qsat[:, i - N]

        f_perturb_pos_scl = ml_load.transform_data(f_ppi, f_pp, f_perturb_pos, z)
        f_perturb_neg_scl = ml_load.transform_data(f_ppi, f_pp, f_perturb_neg, z)
        o_perturb_pos_scl = est_eval.predict(f_perturb_pos_scl)
        o_perturb_neg_scl = est_eval.predict(f_perturb_neg_scl)
        o_perturb_pos = ml_load.inverse_transform_data(o_ppi, o_pp, o_perturb_pos_scl, z)
        o_perturb_neg = ml_load.inverse_transform_data(o_ppi, o_pp, o_perturb_neg_scl, z)
        precip_perturb_pos = atmos_physics.calc_precip(unpack_o(o_perturb_pos, 'q'), rho, z)
        precip_perturb_neg = atmos_physics.calc_precip(unpack_o(o_perturb_neg, 'q'), rho, z)
        sensitivity[i] = np.mean(precip_perturb_pos - precip_perturb_neg)
        sensitivity_stderr[i] = np.std(precip_perturb_pos - precip_perturb_neg) / np.sqrt(f.shape[0])
        # check linearity by using double the perturbation size and halving the resulting sensitivity
        f_perturb_pos = f.copy()
        f_perturb_neg = f.copy()
        if i < N:
            f_perturb_pos[:, i] = f_perturb_pos[:, i] + 1.0 * perturb[i]
            f_perturb_neg[:, i] = f_perturb_neg[:, i] - 1.0 * perturb[i]
        else:
            f_perturb_pos[:, i] = f_perturb_pos[:, i] + 1.0 * perturb[i] * qsat[:, i - N]
            f_perturb_neg[:, i] = f_perturb_neg[:, i] - 1.0 * perturb[i] * qsat[:, i - N]

        f_perturb_pos_scl = ml_load.transform_data(f_ppi, f_pp, f_perturb_pos, z)
        f_perturb_neg_scl = ml_load.transform_data(f_ppi, f_pp, f_perturb_neg, z)
        o_perturb_pos_scl = est_eval.predict(f_perturb_pos_scl)
        o_perturb_neg_scl = est_eval.predict(f_perturb_neg_scl)
        o_perturb_pos = ml_load.inverse_transform_data(o_ppi, o_pp, o_perturb_pos_scl, z)
        o_perturb_neg = ml_load.inverse_transform_data(o_ppi, o_pp, o_perturb_neg_scl, z)
        precip_perturb_pos = atmos_physics.calc_precip(unpack_o(o_perturb_pos, 'q'), rho, z)
        precip_perturb_neg = atmos_physics.calc_precip(unpack_o(o_perturb_neg, 'q'), rho, z)
        sensitivity_check_linear[i] = 0.5 * np.mean(precip_perturb_pos - precip_perturb_neg)
        sensitivity_stderr_check_linear[i] = 0.5 * np.std(precip_perturb_pos - precip_perturb_neg) / np.sqrt(f.shape[0])

        # account for variations in model level spacing
    dp = 50 * 100  # 50hPa
    rho_dz = atmos_physics.vertical_diff(rho, z)
    rescale = dp / (rho_dz * atmos_physics.g)

    # plot linear response on its own
    fig = plt.figure(figsize=(3.0, 2.25))

    plt.plot(p * 0, p, 'k:', lw=0.5)

    plt.plot(rescale * unpack_f_extended(sensitivity, 'T', axis=0, wind_input=wind_input) * per_day, p, 'C0',
             label='T:1K over dp=50hPa')
    plt.plot(rescale * unpack_f_extended(sensitivity_check_linear, 'T', axis=0, wind_input=wind_input) * per_day, p,
             'C0:')

    if use_rh:
        plt.plot(rescale * unpack_f_extended(sensitivity, 'q', axis=0, wind_input=wind_input) * per_day, p, 'C1',
                 label='rh:5% over dp=50hPa')
        plt.plot(rescale * unpack_f_extended(sensitivity_check_linear, 'q', axis=0, wind_input=wind_input) * per_day, p,
                 'C1:')
    else:
        plt.plot(rescale * unpack_f_extended(sensitivity, 'q', axis=0, wind_input=wind_input) * per_day, p, 'C1',
                 label='q:1g/kg over dp=50hPa')
        plt.plot(rescale * unpack_f_extended(sensitivity_check_linear, 'q', axis=0, wind_input=wind_input) * per_day, p,
                 'C1:')

    plt.ylim([np.amax(p), np.amin(p)])
    plt.ylabel('Pressure (hPa)')
    plt.xlabel('Sensitivity in mm $\mathregular{day^{-1}}$ per perturbation')
    plt.legend(loc="upper right")
    plt.legend(frameon=False)
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    plt.tight_layout()  # avoid overlap
    fig.savefig(figpath + 'linear_response.eps', bbof_inches='tight')
    plt.close()

    # plot linear response and feature importance in one graph with error bars for the paper
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6, 2.25))

    ax1.plot(p * 0, p, 'k:', lw=0.5)

    # temperature
    x_plot = rescale * unpack_f_extended(sensitivity, 'T', axis=0, wind_input=wind_input) * per_day
    x_error = rescale * unpack_f_extended(sensitivity_stderr, 'T', axis=0, wind_input=wind_input) * per_day
    ax1.plot(x_plot, p, label='Temperature')

    # water vapor
    x_plot = rescale * unpack_f_extended(sensitivity, 'q', axis=0, wind_input=wind_input) * per_day
    x_error = rescale * unpack_f_extended(sensitivity_stderr, 'q', axis=0, wind_input=wind_input) * per_day
    ax1.plot(x_plot, p, label='Humidity')

    ax1.set_ylim([np.amax(p), np.amin(p)])
    ax1.set_ylabel('Pressure (hPa)')
    ax1.set_xlabel('Sensitivity (mm $\mathregular{day^{-1}}$ per perturbation)')
    ax1.legend(loc="upper right")
    ax1.legend(frameon=False)
    ax1.spines['right'].set_color('none')
    ax1.spines['top'].set_color('none')

    # weighted averaged of decrease in variance
    importance = est_eval.feature_importances_

    # simple estimate of standard error of the mean of feature importance across trees (assumes trees are independent which is not correct)
    estimators = est_eval.estimators_
    importance_stderr = np.std([tree.feature_importances_ for tree in estimators], axis=0) / np.sqrt(len(estimators))

    # account for variations in interlevel spacing
    rescale_importance = 1.0 / rho_dz

    ax2.plot(p * 0, p, 'k:', lw=0.5)

    # temperature
    x_plot = rescale_importance * unpack_f_extended(importance, 'T', axis=0, wind_input=wind_input)
    x_error = rescale_importance * unpack_f_extended(importance_stderr, 'T', axis=0, wind_input=wind_input)
    ax2.plot(x_plot, p)

    # humidity
    x_plot = rescale_importance * unpack_f_extended(importance, 'q', axis=0, wind_input=wind_input)
    x_error = rescale_importance * unpack_f_extended(importance_stderr, 'q', axis=0, wind_input=wind_input)
    ax2.plot(x_plot, p)

    ax2.set_ylim([np.amax(p), np.amin(p)])
    ax2.set_xlabel('Importance (non-dimensional)')
    ax2.spines['right'].set_color('none')
    ax2.spines['top'].set_color('none')

    ax1.set_title('(a) Linear response', fontsize=8, loc="left")
    ax2.set_title('(b) Feature importance', fontsize=8, loc="left")

    plt.tight_layout()  # avoid overlap
    fig.savefig(figpath + 'linear_response_feature_importance.eps', bbox_inches='tight')
    plt.close()


def make_yz_plots(figpath, f_ppi, o_ppi, f_pp, o_pp, est_eval, y, z, p, rho, input_vert_vars, output_vert_vars,
                  input_vert_dim, output_vert_dim,
                  datafile, n_trn_exs, rain_only, no_cos, use_rh, wind_input=False, scale_per_column=False,
                  rewight_outputs=False, weight_list=[1, 1], do_nn=False):

    output_stat_dict = \
        ml_load.stats_by_yz(f_ppi, o_ppi, f_pp, o_pp, est_eval, y,
                            z, rho, datafile, n_trn_exs, input_vert_vars, output_vert_vars, input_vert_dim,
                            output_vert_dim, rain_only,
                            no_cos, use_rh, wind_input=wind_input, scale_per_column=scale_per_column,
                            rewight_outputs=rewight_outputs, weight_list=weight_list,do_nn=do_nn)

    # Make figs
    y_plot = (y - np.mean(y)) / 1e6  # in 1000 km
    # True means

    feature_list = ['_mean', '_bias', '_var', '_rmse', '_r', '_Rsq']

    ncfile = Dataset(figpath + 'data_test.nc', "a")
    ncfile.createDimension("lat_dim", output_stat_dict[next(iter(output_stat_dict))].shape[0])
    latitudes = ncfile.createVariable("lat", np.dtype('float32').char, ("lat_dim"))
    latitudes[:] = y.data[:]
    for key, value in output_stat_dict.items():
        if len(value.shape) > 1:
            temp_var = ncfile.createVariable(key, "f4", ("lat_dim", "press_dim"))
            temp_var[0:value.shape[0], 0:value.shape[1]] = value[:, :]
        else:
            temp_var = ncfile.createVariable(key, "f4", ("lat_dim"))
            temp_var[0:value.shape[0]] = value[:]
    ncfile.description = 'blabla-janniy'
    ncfile.close

    for out_feature, out_dim in zip(output_vert_vars, output_vert_dim):
        if (out_dim > 1 and out_feature != 'Tout' and out_feature != 'qout'):
            fig, ax1, ax2 = plot_contour(output_stat_dict[out_feature + feature_list[0]].T * per_day,
                                         output_stat_dict[out_feature + feature_list[1]].T * per_day, y_plot, z,
                                         p[0:out_dim],
                                         avg_hem=False)
            ax1.set_title(r'$\Delta$ ' + out_feature + ' True Mean')
            ax2.set_title(r'$\Delta$ ' + out_feature + ' bias')
            plt.tight_layout()  # avoid overlap
            fig.savefig(figpath + 'yz_' + out_feature + '_truemean_bias.eps', bbox_inches='tight')
            plt.close()

            fig, ax1, ax2 = plot_contour(output_stat_dict[out_feature + feature_list[2]].T * per_day,
                                         output_stat_dict[out_feature + feature_list[3]].T * per_day, y_plot, z,
                                         p[0:out_dim],
                                         avg_hem=False)
            ax1.set_title(r'$\Delta$ ' + out_feature + 'variance')
            ax2.set_title(r'$\Delta$ ' + out_feature + ' rmse')
            plt.tight_layout()  # avoid overlap
            fig.savefig(figpath + 'yz_' + out_feature + '_var_rmse.eps', bbox_inches='tight')
            plt.close()

            # mask where variance is very low
            mask_too_low = 0.01
            too_low = np.nonzero(output_stat_dict[out_feature + feature_list[2]] < mask_too_low * np.mean(
                output_stat_dict[out_feature + feature_list[2]]))
            Rsq_feature = output_stat_dict[out_feature + feature_list[5]]
            Rsq_feature[too_low] = np.nan

            fig, ax1, ax2 = plot_contour(output_stat_dict[out_feature + feature_list[4]].T,
                                         Rsq_feature.T, y_plot, z, p[0:out_dim],
                                         avg_hem=False)
            ax1.set_title(r'$\Delta$ ' + out_feature + 'r')
            ax2.set_title(r'$\Delta$ ' + out_feature + ' Rsq')
            plt.tight_layout()  # avoid overlap
            fig.savefig(figpath + 'yz_' + out_feature + '_r_Rsq.eps', bbox_inches='tight')
            plt.close()



        elif (
                out_dim == 1 and out_feature != 'Pmean_true' and out_feature != 'Pmean_pred' and out_feature != 'Pextreme_true' and out_feature != 'Pextreme_pred'):
            plt.figure(figsize=(3.0, 2))
            plt.plot(y_plot, output_stat_dict[out_feature + feature_list[0]] * per_day, label='True')
            plt.plot(y_plot, output_stat_dict[out_feature + feature_list[1]] * per_day, label='bias')
            plt.legend(loc="upper right")
            plt.legend(frameon=False)
            plt.xlabel('y (1000km)')
            plt.ylabel('mm $\mathregular{day^{-1}}$')
            ax = plt.gca()
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
            plt.tight_layout()  # avoid overlap
            plt.savefig(figpath + out_feature + 'Truemean_bias_y.eps', bbox_inches='tight')
            plt.close()

            plt.figure(figsize=(3.0, 2))
            plt.plot(y_plot, output_stat_dict[out_feature + feature_list[2]] * per_day, label='var')
            plt.plot(y_plot, output_stat_dict[out_feature + feature_list[3]] * per_day, label='rmse')
            plt.legend(loc="upper right")
            plt.legend(frameon=False)
            plt.xlabel('y (1000km)')
            plt.ylabel('mm $\mathregular{day^{-1}}$')
            ax = plt.gca()
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
            plt.tight_layout()  # avoid overlap
            plt.savefig(figpath + out_feature + 'var_rmse_y.eps', bbox_inches='tight')
            plt.close()

        plt.figure(figsize=(3.0, 2))
        plt.plot(y_plot, output_stat_dict['Pmean_true'] * per_day, label='True')
        plt.plot(y_plot, output_stat_dict['Pmean_pred'] * per_day, label='Pred')
        plt.legend(loc="upper right")
        plt.legend(frameon=False)
        plt.xlabel('y (1000km)')
        plt.ylabel('mm $\mathregular{day^{-1}}$')
        ax = plt.gca()
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        plt.tight_layout()  # avoid overlap
        plt.savefig(figpath + 'precip_mean.eps', bbox_inches='tight')
        plt.close()
        # Precipitation extremes
        plt.figure(figsize=(3.0, 2))
        plt.plot(y_plot, output_stat_dict['Pextreme_true'] * per_day, label='True')
        plt.plot(y_plot, output_stat_dict['Pextreme_pred'] * per_day, label='Pred')
        plt.legend(loc="upper right")
        plt.legend(frameon=False)
        plt.xlabel('y (1000km)')
        plt.ylabel('mm $\mathregular{day^{-1}}$')
        ax = plt.gca()
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        plt.tight_layout()  # avoid overlap
        plt.savefig(figpath + 'precip_extremes.eps', bbox_inches='tight')
        plt.close()

    if 'Tout_mean' in output_stat_dict.keys():
        fig, ax1, ax2 = plot_contour(output_stat_dict['Tout_mean'].T * per_day,
                                     output_stat_dict['qout_mean'].T * per_day * kg_to_g, y_plot, z, p, avg_hem=False)
        ax1.set_title(r'$\Delta$ Temp True Mean [K/day]')
        ax2.set_title(r'$\Delta$ Humid True Mean [g/kg/day]')
        plt.tight_layout()  # avoid overlap
        fig.savefig(figpath + 'yz_truemean.eps', bbox_inches='tight')
        plt.close()
        # Bias from true mean
        fig, ax1, ax2 = plot_contour(output_stat_dict['Tout_bias'].T * per_day,
                                     output_stat_dict['qout_bias'].T * per_day * kg_to_g, y_plot, z, p, avg_hem=False)
        ax1.set_title(r'$\Delta$ Temp Mean Bias [K/day]')
        ax2.set_title(r'$\Delta$ Humid Mean Bias [g/kg/day]')
        plt.tight_layout()  # avoid overlap
        fig.savefig(figpath + 'yz_bias.eps', bbox_inches='tight')
        plt.close()
        # Root mean squared error
        fig, ax1, ax2 = plot_contour(output_stat_dict['Tout_rmse'].T * per_day,
                                     output_stat_dict['qout_rmse'].T * per_day * kg_to_g, y_plot, z, p, avg_hem=False)
        ax1.set_title(r'$\Delta$ Temp RMSE [K/day]')
        ax2.set_title(r'$\Delta$ Humid RMSE [g/kg/day]')
        plt.tight_layout()  # avoid overlap
        fig.savefig(figpath + 'yz_rmse.eps', bbox_inches='tight')
        plt.close()
        # Pearson r Correlation Coefficient
        fig, ax1, ax2 = plot_contour_cc(output_stat_dict['Tout_r'].T, output_stat_dict['qout_r'].T, y_plot, z, p,
                                        avg_hem=False)
        ax1.set_title('(a) Temperature tendency', fontsize=8)
        ax2.set_title('(b) Humidity tendency', fontsize=8)
        #   plt.tight_layout() # avoid overlap
        fig.savefig(figpath + 'yz_corrcoeff.eps', bbox_inches='tight')
        plt.close()
        # Coefficient of determination
        # mask where variance is very low
        Rsq_T = output_stat_dict['Tout_Rsq']
        Rsq_q = output_stat_dict['qout_Rsq']

        mask_too_low = 0.01
        too_low = np.nonzero(output_stat_dict['Tout_var'] < mask_too_low * np.mean(output_stat_dict['Tout_var']))
        Rsq_T[too_low] = np.nan
        too_low = np.nonzero(output_stat_dict['qout_var'] < mask_too_low * np.mean(output_stat_dict['qout_var']))
        Rsq_q[too_low] = np.nan
        fig, ax1, ax2 = plot_contour_Rsq(Rsq_T.T, Rsq_q.T, y_plot, z, p, avg_hem=False)
        ax1.set_title('(a) Temperature tendency', fontsize=8)
        ax2.set_title('(b) Humidity tendency', fontsize=8)
        fig.savefig(figpath + 'yz_Rsq.eps', bbox_inches='tight')
        plt.close()


def plot_contour(T, q, y, z, p, avg_hem=False):
    if avg_hem:
        T, _ = ml_load.avg_hem(T, y, 1)
        q, y = ml_load.avg_hem(q, y, 1)
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    cax1 = ax1.contourf(y, p, T)
    ax1.set_ylim(np.amax(p), np.amin(p))
    ax1.set_ylabel('Pressure (hPa)')
    fig.colorbar(cax1, ax=ax1)
    cax2 = ax2.contourf(y, p, q)
    ax2.set_ylim(np.amax(p), np.amin(p))
    ax2.set_ylabel('Pressure (hPa)')
    fig.colorbar(cax2, ax=ax2)
    ax2.set_xlabel('y (10^3 km)')
    return fig, ax1, ax2


def plot_contour_cc(T, q, y, z, p, avg_hem=False):
    if avg_hem:
        T, _ = ml_load.avg_hem(T, y, 1)
        q, y = ml_load.avg_hem(q, y, 1)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6, 2.25))
    levels = np.linspace(0.2, 1, 9)

    cax1 = ax1.contourf(y, p, T, levels, extend='min')
    ax1.set_ylim(np.amax(p), np.amin(p))
    ax1.set_ylabel('Pressure (hPa)')
    ax1.set_xlabel('y (10^3 km)')

    cax2 = ax2.contourf(y, p, q, levels, extend='min')
    ax2.set_ylim(np.amax(p), np.amin(p))
    ax2.set_xlabel('y (10^3 km)')

    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.93, 0.15, 0.02, 0.7])
    fig.colorbar(cax2, cax=cbar_ax, ticks=np.arange(0.2, 1.1, 0.1))

    return fig, ax1, ax2



def plot_contour_Rsq(T, q, y, z, p, avg_hem=False):
    if avg_hem:
        T, _ = ml_load.avg_hem(T, y, 1)
        q, y = ml_load.avg_hem(q, y, 1)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6, 2))
    levels = np.linspace(0.2, 1, 9)

    cax1 = ax1.contourf(y, p, T, levels, extend='min')
    ax1.set_ylim(np.amax(p), np.amin(p))
    ax1.set_ylabel('Pressure (hPa)')
    ax1.set_xlabel('y (10^3 km)')

    cax2 = ax2.contourf(y, p, q, levels, extend='min')
    ax2.set_ylim(np.amax(p), np.amin(p))
    ax2.set_xlabel('y (10^3 km)')

    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.93, 0.15, 0.02, 0.7])
    fig.colorbar(cax2, cax=cbar_ax, ticks=np.arange(0.2, 1.1, 0.1))

    return fig, ax1, ax2


# Plot means and standard deviations
def plot_means_stds_generalized(y3_true, y3_pred, p, figpath, output_vert_vars, output_vert_dim):
    fig = plt.figure(figsize=(6.0, 4.5))
    dum_i = 1
    true_out_dict = ml_load.unpack_list(y3_true, output_vert_vars, output_vert_dim)
    pred_out_dict = ml_load.unpack_list(y3_pred, output_vert_vars, output_vert_dim)
    ncfile = Dataset(figpath + 'data_test.nc', "a")
    ncfile.createDimension("press_dim", p.shape[0])
    levels = ncfile.createVariable("press_levels", np.dtype('float32').char, ("press_dim"))
    levels[:] = p.data[:]
    ncfile.description = 'blabla-janniy'
    ncfile.close

    for key, dim in zip(output_vert_vars, output_vert_dim):
        # ncfile = Dataset(figpath + 'data_test.nc', "w", format="NETCDF3_CLASSIC")
        ncfile = Dataset(figpath + 'data_test.nc', "a")
        ncfile.createDimension(key + "_dim", dim)
        levels = ncfile.createVariable(key + "_level", np.dtype('float32').char, (key + "_dim"))
        levels[:] = p.data[0:dim]
        ncfile.description = 'blabla-janniy'
        ncfile.close

        if dim > 1:
            do_mean_or_std_y_generalized('mean', key, true_out_dict[key] * per_day, pred_out_dict[key] * per_day,
                                         p[0:dim], dum_i, len(output_vert_vars), figpath, dim)
            do_mean_or_std_y_generalized('std', key, true_out_dict[key] * per_day, pred_out_dict[key] * per_day,
                                         p[0:dim], dum_i + 1, len(output_vert_vars), figpath, dim)
            dum_i = dum_i + 2
        plt.tight_layout()  # avoid overlap
    fig.savefig(figpath + 'means_stds.eps', bbox_inches='tight')
    plt.close()




# Plot correlation coefficient, explained variance, and rmse
def plot_error_stats(y3_true, y3_pred, p, figpath, output_vert_vars, output_vert_dim):
    fig = plt.figure(figsize=(6.0, 4.5))
    true_out_dict = ml_load.unpack_list(y3_true, output_vert_vars, output_vert_dim)
    pred_out_dict = ml_load.unpack_list(y3_pred, output_vert_vars, output_vert_dim)
    plt.subplot(len(output_vert_vars), 2, 1)
    for key, dim in zip(output_vert_vars, output_vert_dim):
        if dim > 1:
            plot_pearsonr_generalized(true_out_dict[key], pred_out_dict[key], key, p[0:dim], label=key)
    plt.legend(loc="lower left")
    plt.legend(frameon=False)
    plt.subplot(len(output_vert_vars), 2, 2)
    for key, dim in zip(output_vert_vars, output_vert_dim):
        if dim > 1:
            plot_expl_var_generalized(true_out_dict[key], pred_out_dict[key], key, p[0:dim])
    fig_ind = 2
    for key, dim in zip(output_vert_vars, output_vert_dim):
        fig_ind = fig_ind + 1
        if len(output_vert_vars) == 1:
            plt.subplot(2, 2, fig_ind)
        else:
            plt.subplot(len(output_vert_vars), 2, fig_ind)
        if key == 'qout':
            factor1 = kg_to_g
        else:
            factor1 = per_day
        if dim > 1:
            plot_rmse_generalized(true_out_dict[key] * factor1, pred_out_dict[key] * factor1, key, p[0:dim])
    plt.tight_layout()  # avoid overlap
    fig.savefig(figpath + 'error_stats.eps', bbox_inches='tight')
    plt.close()


# Plot a scatter plot of true vs predicted
def plot_scatter(f, otrue, opred, z, p, rho, figpath, output_vert_vars, output_vert_dim):
    true_out_dict = ml_load.unpack_list(otrue, output_vert_vars, output_vert_dim)
    pred_out_dict = ml_load.unpack_list(opred, output_vert_vars, output_vert_dim)
    # Plot scatter of precipitation
    P_true = atmos_physics.calc_precip(true_out_dict['qout'], rho, z, output_vert_vars, true_out_dict)
    P_pred = atmos_physics.calc_precip(pred_out_dict['qout'], rho, z, output_vert_vars, pred_out_dict)

    ncfile = Dataset(figpath + 'data_test.nc', "a")
    P_true_var = ncfile.createVariable("P_true_samples", "f4", ('test_samples'))
    P_pred_var = ncfile.createVariable("P_pred_samples", "f4", ('test_samples'))
    P_true_var[:] = P_true[:]
    P_pred_var[:] = P_pred[:]
    ncfile.description = 'blabla-janniy'
    ncfile.close

    f = plt.figure(figsize=(3.0, 2))
    ax = plt.gca()
    _plot_scatter(ax, P_true * per_day, P_pred * per_day)
    ax.set_xlabel('SAM (mm $\mathregular{day^{-1}}$)')
    ax.set_ylabel('Random forest (mm $\mathregular{day^{-1}}$)')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    plt.tight_layout()  # avoid overlap
    f.savefig(figpath + 'precip_scatter.pdf', bbox_inches='tight')
    plt.close()

    Pneg = sum(P_pred < 0.0)
    Pnegpct = 100. * Pneg / len(P_pred)

    Pnegtrue = sum(P_true < 0.0)
    Pnegtruepct = 100. * Pnegtrue / len(P_true)

    Pzero = sum(P_pred == 0)
    Pzeropct = 100. * Pzero / len(P_pred)

    Pzerotrue = sum(P_true == 0)
    Pzerotruepct = 100. * Pzerotrue / len(P_pred)

    Psmall = sum(P_pred * per_day < 2.5)  # less than 2.5 mm/day
    Psmallpct = 100. * Psmall / len(P_pred)

    Psmalltrue = sum(P_true * per_day < 2.5)  #
    Psmalltruepct = 100. * Psmalltrue / len(P_true)

    Pcc, _ = scipy.stats.pearsonr(P_true, P_pred)
    PRsq = metrics.r2_score(P_true, P_pred)
    Pbias = np.mean(P_pred - P_true) * per_day

    text_file = open(figpath + 'precip.txt', 'w')
    text_file.write('Pred. P<0 %f percent of time \n' % Pnegpct)
    text_file.write('True. P<0 %f percent of time \n' % Pnegtruepct)
    text_file.write('Pred. P zero %f percent of time \n' % Pzeropct)
    text_file.write('True P zero %f percent of time \n' % Pzerotruepct)
    text_file.write('Pred. P<2.5 mm/day %f percent of time \n' % Psmallpct)
    text_file.write('True P<2.5 mm/day %f percent of time \n' % Psmalltruepct)
    text_file.write('Correlation coefficient for P: %f \n' % Pcc)
    text_file.write('Coefficient of determination for P: %f \n' % PRsq)
    text_file.write('Mean bias for P: %e \n' % Pbias)
    text_file.close()



def _plot_scatter(ax, true, pred, titstr=None):
    ax.scatter(true, pred, s=5, c='b', alpha=0.2)
    # Calculate mins and maxs and set axis bounds appropriately
    xmin = np.min(true)
    xmax = np.max(true)
    ymin = np.min(pred)
    ymax = np.max(pred)
    xymin = np.min([xmin, ymin])
    xymax = np.max([xmax, ymax])
    # Plot 1-1 line
    ax.plot([xymin, xymax], [xymin, xymax], color='k', ls='--')
    ax.set_xlabel('True')
    ax.set_ylabel('Predicted')
    if titstr is not None:
        ax.set_title(titstr)


# Plot energy conservation
def plot_energy_conservation(f, o_true, o_pred, z, p, rho, figpath, input_vert_vars, input_vert_dim, output_vert_vars,
                             output_vert_dim):
    true_out_dict = ml_load.unpack_list(o_true, output_vert_vars, output_vert_dim)
    pred_out_dict = ml_load.unpack_list(o_pred, output_vert_vars, output_vert_dim)
    feature_dict = ml_load.unpack_list(f, input_vert_vars, input_vert_dim)

    fig = plt.figure(figsize=(3.0, 4.5))
    plt.subplot(2, 1, 1)
    _plot_energy_conservation(feature_dict, true_out_dict, z, p, rho, output_vert_vars, figpath, label='true')
    plt.legend(loc="upper left")
    plt.legend(frameon=False)
    plt.subplot(2, 1, 2)
    _plot_energy_conservation(feature_dict, pred_out_dict, z, p, rho, output_vert_vars, figpath, label='predict')
    plt.legend(loc="upper left")
    plt.legend(frameon=False)
    plt.tight_layout()  # avoid overlap
    fig.savefig(figpath + 'energy_conservation.eps', bbox_inches='tight')
    plt.close()


# ----  PLOTTING SCRIPTS  ---- #
def do_mean_or_std_y_generalized(method, vari, true, pred, p, ind, plot_col, figpath, dim):
    methods = {'mean': np.mean, 'std': np.std}
    methods_ti = {'mean': 'Mean', 'std': 'Standard Deviation'}
    plt.subplot(plot_col, 2, ind)
    plt.plot(methods[method](true, axis=0).T, p, label='true')
    plt.plot(methods[method](pred, axis=0).T, p, label='pred')

    ncfile = Dataset(figpath + 'data_test.nc', "a")
    temp_true = ncfile.createVariable(vari + '_' + method + '_true', "f4", (vari + "_dim",))
    temp_pred = ncfile.createVariable(vari + '_' + method + '_pred', "f4", (vari + "_dim",))
    temp_true[:] = methods[method](true, axis=0).T[:]
    temp_pred[:] = methods[method](pred, axis=0).T[:]
    ncfile.description = 'blabla-janniy'
    ncfile.close

    plt.plot(0 * p, p, linewidth=0.3)
    plt.ylim(np.amax(p), np.amin(p))
    plt.ylabel('Pressure (hPa)')
    out_str_dict = {'T': 'K/day', 'q': 'g/kg/day'}
    # if ind > 2:
    #     plt.xlabel(out_str_dict[vari])
    plt.title(r'$\Delta$ ' + vari + " " + methods_ti[method])
    plt.legend()
    plt.legend(frameon=False)



def plot_pearsonr_generalized(o_true, o_pred, vari, p, label=None):
    r = np.empty(o_true.shape[1])
    prob = np.empty(o_true.shape[1])
    for i in range(o_true.shape[1]):
        r[i], prob[i] = scipy.stats.pearsonr(o_true[:, i], o_pred[:, i])
    plt.plot(r, p, label=label)
    plt.ylim([np.amax(p), np.amin(p)])
    # plt.xlim(0, 1)
    plt.ylabel('Pressure (hPa)')
    plt.title('Correlation Coefficient')




def plot_rmse_generalized(o_true, o_pred, key, p, label=None):
    rmse = np.sqrt(metrics.mean_squared_error(o_true, o_pred,
                                              multioutput='raw_values'))
    plt.plot(rmse, p, label=label)
    plt.ylim([np.amax(p), np.amin(p)])
    plt.xlim(xmin=0)
    plt.ylabel('Pressure (hPa)')
    plt.xlabel(key)
    plt.title('Root Mean Squared Error' + key)




def plot_expl_var_generalized(o_true, o_pred, vari, p, label=None):
    expl_var = metrics.explained_variance_score(o_true, o_pred,
                                                multioutput='raw_values')
    plt.plot(expl_var, p, label=label)
    plt.ylim([np.amax(p), np.amin(p)])
    plt.xlim(0, 1)
    plt.ylabel('Pressure (hPa)')
    plt.title('Explained Variance Regression Score')



def _plot_energy_conservation(f_dict, o_dict, z, p, rho, output_vert_vars, figpath, label=None):
    tend_residual = atmos_physics.energy_tendency_residual(f_dict['Tin'], o_dict['Tout'], o_dict['qout'], rho, z,
                                                           output_vert_vars, o_dict)
    ncfile = Dataset(figpath + 'data_test.nc', "a")
    energy_var = ncfile.createVariable("energy_tend_residual" + label, "f4", ('test_samples'))
    energy_var[:] = tend_residual[:]
    ncfile.description = 'blabla-janniy'
    ncfile.close
    n, bins, patches = plt.hist(tend_residual * per_day, 50, alpha=0.5, label=label)
    plt.title('RMS heating rate: {:1.2e}'.format(np.sqrt(np.mean(tend_residual ** 2)) * per_day))
    plt.xlabel('K/day over column')


def check_scaling_distribution(f, f_scl, o, o_scl, y, p,
                               figpath, use_rh=False, wind_input=False):
    # For input variables
    fig, ax = plt.subplots(2, 2)
    _plot_distribution(unpack_f_extended(f, 'T', wind_input=wind_input), y, p, fig, ax[0, 0], './figs/',
                       'T (unscaled) [K]', '')

    _plot_distribution(unpack_f_extended(f_scl, 'T', wind_input=wind_input), y, p, fig, ax[0, 1], './figs/',
                       'T (scaled) []', '')

    if use_rh:
        _plot_distribution(unpack_f_extended(f, 'q', wind_input=wind_input), y, p, fig, ax[1, 0], './figs/',
                           'rh (unscaled)', '')
    else:
        _plot_distribution(unpack_f_extended(f, 'q', wind_input=wind_input) * kg_to_g, y, p, fig, ax[1, 0], './figs/',
                           'q (unscaled) [g/kg]', '')

    _plot_distribution(unpack_f_extended(f_scl, 'q', wind_input=wind_input), y, p, fig, ax[1, 1], './figs/',
                       'q (scaled) []', '')

    plt.tight_layout()  # avoid overlap
    fig.savefig(figpath + 'input_scaling_check.png', bbox_inches='tight',
                dpi=450)
    plt.close()
    # For output variables
    fig, ax = plt.subplots(2, 2)
    _plot_distribution(unpack_o(o, 'T') * per_day, y, p, fig, ax[0, 0], './figs/', 'T tend (unscaled) [K/day]', '')

    _plot_distribution(unpack_o(o_scl, 'T'), y, p, fig, ax[0, 1], './figs/', 'T tend (scaled) []', '')

    _plot_distribution(unpack_o(o, 'q') * per_day * kg_to_g, y, p, fig, ax[1, 0], './figs/',
                       'q tend (unscaled) [g/kg/day]', '')

    _plot_distribution(unpack_o(o_scl, 'q'), y, p, fig, ax[1, 1], './figs/', 'q tend(scaled) []', '')

    plt.tight_layout()  # avoid overlap
    fig.savefig(figpath + 'output_scaling_check.png', bbox_inches='tight',
                dpi=450)
    plt.close()


def check_output_distribution(ot, ot_scl, op, op_scl, y, p,
                              figpath):
    # For unscaled variables
    fig, ax = plt.subplots(2, 2)
    x1, x2, bins = _plot_distribution(unpack_o(ot, 'T') * per_day, y, p, fig,
                                      ax[0, 0], './figs/',
                                      r'$\Delta$T true [K/day]', '')
    _plot_distribution(unpack_o(op, 'T') * per_day, y, p, fig,
                       ax[0, 1], './figs/', r'$\Delta$T pred [K/day]', '', x1,
                       x2, bins)
    x1, x2, bins = _plot_distribution(unpack_o(ot, 'q') * per_day * kg_to_g, y, p, fig,
                                      ax[1, 0], './figs/',
                                      r'$\Delta$q true [g/kg/day]', '')
    _plot_distribution(unpack_o(op, 'q') * per_day * kg_to_g, y, p, fig, ax[1, 1],
                       './figs/', r'$\Delta$q pred [g/kg/day]', '', x1, x2,
                       bins)
    plt.tight_layout()  # avoid overlap
    fig.savefig(figpath + 'output_compare_true_pred_unscaled.png',
                bbox_inches='tight', dpi=450)
    plt.close()
    # For scaled variables
    fig, ax = plt.subplots(2, 2)
    x1, x2, bins = _plot_distribution(unpack_o(ot_scl, 'T'), y, p, fig,
                                      ax[0, 0], './figs/',
                                      r'$\Delta$T true (scld) []', '')
    _plot_distribution(unpack_o(op_scl, 'T'), y, p, fig, ax[0, 1], './figs/',
                       r'$\Delta$T pred (scld) []', '', x1, x2, bins)
    x1, x2, bins = _plot_distribution(unpack_o(ot_scl, 'q'), y, p, fig,
                                      ax[1, 0], './figs/',
                                      r'$\Delta$q true (scld) []', '')
    _plot_distribution(unpack_o(op_scl, 'q'), y, p, fig, ax[1, 1], './figs/',
                       r'$\Delta$q pred (scld) []', '', x1, x2, bins)
    plt.tight_layout()  # avoid overlap
    fig.savefig(figpath + 'output_compare_true_pred_scaled.png',
                bbox_inches='tight', dpi=450)
    plt.close()


def _plot_distribution(z, y, p, fig, ax, figpath, titlestr, xstr, xl=None,
                       xu=None, bins=None):
    """Plots a stack of histograms of log10(data) at all levels"""
    # Initialize the bins and the frequency
    num_bins = 100
    if bins is None:
        bins = np.linspace(np.percentile(z, .02), np.percentile(z, 99.98),
                           num_bins + 1)
    n = np.zeros((num_bins, p.size))
    # Calculate distribution at each level
    for i in range(p.size):
        n[:, i], _ = np.histogram(z[:, i], bins=bins)
    # Take a logarithm and deal with case where we take log of 0
    n = np.log10(n)
    n_small = np.amin(n[np.isfinite(n)])
    n[np.isinf(n)] = n_small
    # Plot histogram
    ca = ax.contourf(bins[:-1], p, n.T)
    ax.set_ylim(np.amax(p), np.amin(p))
    if xl is not None:
        ax.set_xlim(xl, xu)
    plt.colorbar(ca, ax=ax)
    ax.set_xlabel(xstr)
    ax.set_ylabel('Pressure (hPa)')
    ax.set_title(titlestr)
    xl, xr = ax.set_xlim()
    return xl, xr, bins


def plot_sample_profiles(num_prof, f, otrue, opred, p, figpath, samp=None, wind_input=False):
    # Make directory if one does not exist
    if not os.path.exists(figpath + '/samples/'):
        os.makedirs(figpath + '/samples/')
    # Plot some number of sample profiles
    for i in range(num_prof):
        if samp is None:
            samp = np.random.randint(0, f.shape[0])
        sample_filename = figpath + '/samples/' + str(samp) + '.eps'
        plot_sample_profile(f[samp, :], otrue[samp, :], opred[samp, :],
                            p, filename=sample_filename, wind_input=wind_input)
        # Need to reset sample number for next iteration
        samp = None


def plot_sample_profile(f, o_true, o_pred, p, filename=None, pflag=False, wind_input=False):
    """Plots the vertical profiles of input T & q and predicted and true
    output tendencies"""
    fig, (ax1, ax3) = plt.subplots(1, 2, figsize=(7.5, 5))
    T = unpack_f_extended(f, 'T', axis=0, wind_input=wind_input)
    q = unpack_f_extended(f, 'q', axis=0, wind_input=wind_input)
    # Plot input temperature profile
    ax1.plot(T, p, label=r'$T$')
    ax3.set_ylim(np.amax(p), np.amin(p))
    ax1.set_title('Input Profiles')
    ax1.grid(True)
    ax1.legend(loc='upper left')
    ax1.legend(frameon=False)
    cp = atmos_physics.cp
    L = atmos_physics.L
    kJ_scale = 0.001
    k_per_day = 3600 * 24
    ax3.plot(cp * ml_load.unpack_o(o_true, 'T', axis=0) * per_day * kJ_scale, p, color='red',
             ls='-', label=r'$\Delta$T true')
    ax3.plot(cp * ml_load.unpack_o(o_pred, 'T', axis=0) * per_day * kJ_scale, p, color='red',
             ls='--', label=r'$\Delta$T pred')
    ax3.plot(L * ml_load.unpack_o(o_true, 'q', axis=0) * per_day * kJ_scale, p, color='blue',
             ls='-', label=r'$\Delta$q true')
    ax3.plot(L * ml_load.unpack_o(o_pred, 'q', axis=0) * per_day * kJ_scale, p, color='blue',
             ls='--', label=r'$\Delta$q pred')
    ax3.set_ylim(np.amax(p), np.amin(p))
    ax3.set_xlabel('Cp*T or L*q [kJ/day/kg]')
    ax1.set_ylabel('Pressure [hPa]')
    ax3.set_title('Output Tendencies')
    ax3.legend(loc="upper left")
    ax3.legend(frameon=False)
    ax3.grid(True)
    fig.tight_layout()
    # Save file if requested
    if filename is not None:
        fig.savefig(filename, bbox_inches='tight')
        plt.close()


def plot_model_error_over_time(errors, est_str, fig_dir):
    x = np.arange(errors.shape[0])
    ytix = [.5e-3, 1e-3, 2e-3, 5e-3, 10e-3, 20e-3, 50e-3, 100e-3]
    # Plot error rate vs. iteration number
    fig = plt.figure(figsize=(3.0, 2.25))
    # Plot training errors from cost function
    plt.semilogy(x, np.squeeze(errors[:, 0]), alpha=0.5, color='blue',
                 label='Training (cost function)')
    plt.yticks(ytix, ytix)
    plt.ylim((np.nanmin(errors), np.nanmax(errors)))
    # Plot training errors that are not associated with cost function
    plt.semilogy(x, np.squeeze(errors[:, 4]), alpha=0.5, color='red',
                 label='Training')
    # Plot cross-validation errors
    plt.semilogy(x, np.squeeze(errors[:, 2]), alpha=0.5, color='green',
                 label='Cross-Val')
    plt.legend()
    plt.legend(frameon=False)
    plt.title('Error for ' + est_str)
    plt.xlabel('Iteration Number')
    fig.savefig(fig_dir + 'error_history.png', bbox_inches='tight', dpi=450)
    plt.close()


def test_estimator(est_eval, f_ppi, o_ppi, f_pp, o_pp, z):
    # test f90 version works for an example case
    # rf:
    f = np.array([[271.2355, 270.9278, 269.5019, 269.3549, 267.3458, 264.9370, 262.6389, 259.4129, 255.3006, 251.5926,
                   246.6970, 243.2816, 239.1782, 234.0406, 229.5529, 224.2106, 218.8746, 214.7380, 208.3731, 204.8103,
                   200.6235, 197.5632, 196.3092, 195.5485, 193.2881, 189.4108, 185.2819, 186.3155, 182.0662, 182.9386,
                   179.2453, 180.2408, 183.9875, 183.4066, 179.1729, 179.6529, 182.1947, 180.6949, 180.8699, 186.4952,
                   190.7054, 2.1817051E-03, 2.1621622E-03, 1.4375404E-03, 1.3706952E-03, 1.3715313E-03, 1.3734254E-03,
                   6.2456832E-04, 6.5306388E-04, 4.8721873E-04, 4.4362518E-04, 3.7724295E-04, 3.6783234E-04,
                   1.5334306E-04, 1.0800199E-04, 7.0741400E-05, 6.1005565E-05, 4.3920740E-05, 1.4823842E-05,
                   1.0503937E-05, 4.5622778E-06, 2.9868793E-06, 1.9148797E-06, 1.8365218E-06, 1.6058746E-06,
                   1.5466301E-06, 1.6219693E-06, 1.7039781E-06, 1.8224653E-06, 1.8860034E-06, 2.2180418E-06,
                   2.2034710E-06, 2.5527524E-06, 3.0400804E-06, 3.1567427E-06, 3.1452885E-06, 3.2680146E-06,
                   3.3647266E-06, 3.3535289E-06, 3.3496365E-06, 3.1498450E-06, 2.9508053E-06]])
    f_scl = ml_load.transform_data(f_ppi, f_pp, f, z)
    o_scl = est_eval.predict(f_scl)
    o = ml_load.inverse_transform_data(o_ppi, o_pp, o_scl, z)
    print(o)
    sys.exit()

