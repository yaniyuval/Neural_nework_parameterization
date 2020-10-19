import numpy as np
# import sknn_jgd.mlp
import time
from sklearn.ensemble import RandomForestRegressor
# from src.ml_io as  create_specific_data_string_desc
from src.ml_io import write_netcdf_nn
import src.ml_io as ml_io
import src.ml_load as ml_load
import src.ml_train as ml_train


import pickle
# import src.ml_plot as ml_plot
import os
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import cross_val_score
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.datasets import load_digits
from sklearn.model_selection import learning_curve
from sklearn.model_selection import ShuffleSplit


# from xgboost.sklearn import XGBRegressor

def evaluate(model, test_features, test_labels):
    predictions = model.predict(test_features)
    errors = abs(predictions - test_labels)
    mape = 100 * np.mean(errors / abs(test_labels / 2 + predictions / 2))
    accuracy = 100 - mape
    print('Model Performance')
    print('Average Error: {:0.4f} degrees.'.format(np.mean(errors)))
    print('Accuracy = {:0.2f}%.'.format(accuracy))
    return accuracy


def PreprocessData(f_ppi, f, o_ppi, o, pp_str, n_trn_exs, z):
    """Transform data according to input preprocessor requirements and make
    make preprocessor string for saving"""
    f_pp = ml_load.init_pp(f_ppi, f)
    f = ml_load.transform_data(f_ppi, f_pp, f, z)
    o_pp = ml_load.init_pp(o_ppi, o)
    o = ml_load.transform_data(o_ppi, o_pp, o, z)
    # Make preprocessor string for saving
    pp_str = pp_str + 'F-' + f_ppi['name'] + '_'
    pp_str = pp_str + 'O-' + o_ppi['name'] + '_'
    # Add number of training examples to string
    pp_str = pp_str + 'Ntrnex' + str(n_trn_exs) + '_'
    return f_pp, f, o_pp, o, pp_str


def plot_learning_curve(estimator, title, X, y, ylim=None, cv=None, n_jobs=None, train_sizes=np.linspace(.1, 1.0, 5)):
    """
    Generate a simple plot of the test and training learning curve.

    Parameters
    ----------
    estimator : object type that implements the "fit" and "predict" methods
        An object of that type which is cloned for each validation.

    title : string
        Title for the chart.

    X : array-like, shape (n_samples, n_features)
        Training vector, where n_samples is the number of samples and
        n_features is the number of features.

    y : array-like, shape (n_samples) or (n_samples, n_features), optional
        Target relative to X for classification or regression;
        None for unsupervised learning.

    ylim : tuple, shape (ymin, ymax), optional
        Defines minimum and maximum yvalues plotted.

    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:
          - None, to use the default 3-fold cross-validation,
          - integer, to specify the number of folds.
          - :term:`CV splitter`,
          - An iterable yielding (train, test) splits as arrays of indices.

        For integer/None inputs, if ``y`` is binary or multiclass,
        :class:`StratifiedKFold` used. If the estimator is not a classifier
        or if ``y`` is neither binary nor multiclass, :class:`KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validators that can be used here.

    n_jobs : int or None, optional (default=None)
        Number of jobs to run in parallel.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    train_sizes : array-like, shape (n_ticks,), dtype float or int
        Relative or absolute numbers of training examples that will be used to
        generate the learning curve. If the dtype is float, it is regarded as a
        fraction of the maximum size of the training set (that is determined
        by the selected validation method), i.e. it has to be within (0, 1].
        Otherwise it is interpreted as absolute sizes of the training sets.
        Note that for classification the number of samples usually have to
        be big enough to contain at least one sample from each class.
        (default: np.linspace(0.1, 1.0, 5))
    """
    plt.figure()
    plt.title(title)
    if ylim is not None:
        plt.ylim(*ylim)
    plt.xlabel("Training examples")
    plt.ylabel("Score")
    train_sizes, train_scores, test_scores = learning_curve(estimator, X, y, cv=cv, n_jobs=n_jobs, train_sizes=train_sizes)
    train_scores_mean = np.mean(train_scores, axis=1)
    train_scores_std = np.std(train_scores, axis=1)
    test_scores_mean = np.mean(test_scores, axis=1)
    test_scores_std = np.std(test_scores, axis=1)
    plt.grid()
    plt.fill_between(train_sizes, train_scores_mean - train_scores_std, train_scores_mean + train_scores_std, alpha=0.1,color="r")
    plt.fill_between(train_sizes, test_scores_mean - test_scores_std,test_scores_mean + test_scores_std, alpha=0.1, color="g")
    plt.plot(train_sizes, train_scores_mean, 'o-', color="r",label="Training score")
    plt.plot(train_sizes, test_scores_mean, 'o-', color="g",label="Cross-validation score")
    plt.legend(loc="best")
    return plt


do_train = True
#max_z=20000
max_z=np.Inf
min_samples_leaf = 3
n_trees = 10#10
max_depth = 27

n_trn_exs = 100000000#990000000
#n_trn_exs = 2000
use_rh = False
no_cos = True


flag_dict = dict()
only_plot = False # If I only put to plot....
scale_level = False #If true, each column is scaled seperately.
read_from_whole_data = 1 #In case we want to read the data from a pkl file that has all possible inputs

flag_dict['do_dqp'] = False
flag_dict['ver_adv_correct'] = False
flag_dict['do_hor_wind_input'] = True
flag_dict['do_ver_wind_input'] = False
flag_dict['do_z_diffusion'] = True
flag_dict['do_z_diffusion_correction'] = False  # True

flag_dict['do_q_T_surf_fluxes'] = False  # True  #Only relevant if vertical diffusion is included.
flag_dict['do_surf_wind'] = True  # True
flag_dict['do_q_surf_fluxes_out'] = True

flag_dict[
    'do_sedimentation'] = False  # if I want to include the sedimentation tendencies (from the cloud scheme- Need to calculate it in matlab)
flag_dict['do_fall_tend'] = False
flag_dict['do_qp_as_var'] = False  # If I want to run the simulation with qp as a prognostic parameter.

# Later I can consider doing the fluxes seperately from the radiation and all the micro tendencies and try running a NN with such
# conserving output.
flag_dict['do_radiation_output'] = False  # If want to predict the radiation seperately.
flag_dict['rad_level'] = 0  # Should be 0 if no radiation should be used in hte RF

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
flag_dict['dist_From_eq_in'] = True
flag_dict['T_instead_of_Tabs'] = False


flag_dict['tabs_resolved_init'] = True #Should be true - We only know the resolved tabs as far as I understand
flag_dict['qn_coarse_init'] = True #Should be true - want to have qn coarse from beginning of time step - to calculate qt
flag_dict['qn_resolved_as_var'] = False #Should be true - want to have qn coarse from beginning of time step - to calculate qt
flag_dict['do_zadv_sed_output'] = False

flag_dict['sed_level'] = 0
flag_dict['strat_corr_level'] = 0


dx = 12000*flag_dict['resolution']
dy = 12000*flag_dict['resolution']



data_specific_description = ml_io.create_specific_data_string_desc(flag_dict)

training_expt1 = 'qobs'+data_specific_description
#training_expt2 = 'qobs4K'
do_wind_input = False #Yani added
do_diffution=False

input_vert_vars = ['Tin','qin','uin','vin_minusSH','usurf','latin'] #Dependent on the scenario we chose to model. This should give some flexibility to our coding.
output_vert_vars = ['tsurf','qsurf','tkz']
input_vert_dim = [15,15,15,15,1,1] #Dependent on the scenario we chose to model. This should give some flexibility to our coding.
output_vert_dim = [1,1,15]


is_cheyenne = True
# Define preprocessor
f_ppi = {'name': 'NoScaler'} # At the moment treated only the possibility to put here NoScalar (for RF it doesn't matter)
# f_ppi = {'name': 'StandardScaler'}

# o_ppi = {'name': 'SimpleO'}
o_ppi = {'name': 'StandardScaler'}
rewight_outputs = False #If I want to give more wight to certain features.
weight_list = [1,1]
rain_only = False


datadir, trainfile, testfile, pp_str = ml_load.GetDataPath(training_expt1,is_cheyenne = is_cheyenne)

f, o, y, z, rho, p = ml_load.LoadData(trainfile, max_z, input_vert_vars=input_vert_vars,
                                      output_vert_vars=output_vert_vars,
                                      rain_only=rain_only, n_trn_exs=n_trn_exs, no_cos=no_cos, use_rh=use_rh,
                                      wind_input=do_wind_input)
print('read train data')

# load test data
tf, to, ty, tz, trho, tp = ml_load.LoadData(testfile, max_z, input_vert_vars=input_vert_vars,
                                            output_vert_vars=output_vert_vars,
                                            rain_only=rain_only, n_trn_exs=n_trn_exs, no_cos=no_cos, use_rh=use_rh)
print('read test data')

# Scale data (both train and test)
f_pp, f_scl, tf_scl, o_pp, o_scl, to_scl, pp_str = ml_train.PreprocessData_tr_ts(f_ppi, f, tf, o_ppi, o, to, pp_str,
                                                                        n_trn_exs, z, input_vert_dim, input_vert_vars,
                                                                        output_vert_dim, output_vert_vars, scale_level,
                                                                        rewight_outputs=rewight_outputs,
                                                                        weight_list=weight_list)  # Yani TO DO!!!



# X_train = np.concatenate((f_scl, tf_scl))
# Y_train = np.concatenate((o_scl, to_scl))
# per_test = 0.03
# per_train = 1 - per_test
# max_ind = np.int(np.ceil(X_train.shape[0] * per_train))
# f_scl = X_train[0:max_ind, :]
# tf_scl = X_train[max_ind:, :]
# o_scl = Y_train[0:max_ind, :]
# to_scl = Y_train[max_ind:, :]



n_estimators = 10
max_features = 1.0/3.0  # ['sqrt']
max_depth = 27
# Minimum number of samples required to split a node
# min_samples_split = [10]
# Minimum number of samples required at each leaf node
min_samples_leaf = 20
# Method of selecting samples for training each tree
bootstrap = True
# Create the random grid
# num_train = 1000#f_scl.shape[0]
# print("the number of training samples are",num_train)
# f_scl = f_scl[0:num_train]
# o_scl = o_scl[0:num_train]

# for i in n_estimators:
#     print(i)
#     n_estimators_tmp = i
#     for j in max_depth:
#         max_depth_tmp = j
#         for k in min_samples_leaf:
#             min_samples_leaf_tmp = k
#             # print("n estimators: ", n_estimators_tmp, "max_depth", max_depth_tmp, "min_samples_leaf", min_samples_leaf_tmp )
#             start = time.time()
#             rf = RandomForestRegressor(max_depth=max_depth_tmp, n_estimators=n_estimators_tmp, min_samples_leaf=min_samples_leaf, n_jobs=11)
            # rf.fit(f_scl, o_scl)
            # scores = cross_val_score(estimator=rf,
            #                 X=f_scl,
            #                y=o_scl,
            #               cv=3,   ## This is the split size of the train data (should be larger for smaller data - to reduce bias)
            #              n_jobs=1)
        # print('CV accuracy scores: %s' % scores)
        # print('CV accuracy: %.3f +/- %.3f' % (np.mean(scores), np.std(scores)))
    # for item in scores:
    #     file.write("scores: %s\n" %item)
    # file.write('CV accuracy: %.3f +/- %.3f\n' % (np.mean(scores), np.std(scores)))

    # train_score = rf.score(f_scl, o_scl)
    # test_score = rf.score(tf_scl, to_scl)
    # print("train score:", train_score, "test score: ", test_score)
    # file.write('train_score : %.3f\n' % (train_score))
    # file.write('test_score : %.3f\n' % (test_score))

    # end = time.time()
    # print("time to run the model ({:.1f} seconds)".format( end-start))
    # file.write('time of run : %.3f\n' % (end- start))
do_nn = False
rf = RandomForestRegressor(max_depth=max_depth, n_estimators=n_estimators, min_samples_leaf=min_samples_leaf, n_jobs=11)
title = "Learning Curve RF"
cv = ShuffleSplit(n_splits=1, test_size=0.15, random_state=0)
# plot_learning_curve(rf, title, f_scl, o_scl, ylim=(0.0, 1.01), cv=cv, n_jobs=4)
ylim = (0.0, 1.01)
fig = plt.figure()
plt.title(title)
plt.ylim(ylim)
plt.xlabel("Training examples")
plt.ylabel("Score")
# data_part = [0.001,0.003, 0.01,0.03,0.1,0.3,0.5713,0.99]
data_part = [5000, 10000, 30000, 100000, 300000, 1000000, 5000000,f_scl.shape[0]]
number_of_data_points = len(data_part)

train_score = np.zeros([number_of_data_points])
test_score = np.zeros([number_of_data_points])
train_sizes= np.zeros([number_of_data_points])
Output_ind1 = 2
Output_ind2 = 17


for batch_i in range(number_of_data_points):
    # ind_data = int(np.ceil(data_part[batch_i] *f_scl.shape[0]))
    ind_data = int(np.ceil(data_part[batch_i]))
    # print(ind_data)
    train_sizes[batch_i] = ind_data
    print('number of points is', ind_data)
    # ind_test= min(ind_data*5,tf_scl.shape[0])
    ind_test = np.int(tf_scl.shape[0] / 2)  # take only CV
    est, est_errors, train_score[batch_i], test_score[batch_i] = ml_train.train_est(rf, 'rf_str', f_scl[0:ind_data, Output_ind1:Output_ind2], o_scl[0:ind_data, Output_ind1:Output_ind2], tf_scl[0:ind_test, Output_ind1:Output_ind2], to_scl[0:ind_test,  Output_ind1:Output_ind2], do_nn)



figpath = '/glade/scratch/janniy/offline_figs/qobsTTFFFFFTT26TTTFTF48TFFFFFTFFFFF80FFFFTTF2626_F-NoSc_O-Stan_Ntr5000000_Nte972360_F_Tin_qin_qpin_O_Tout_qout_qpout_qrad_RF_NTr10_MinS3max_d27_maxzinf_nocos_te70_tr79/'
figpath = '/glade/scratch/janniy/offline_figs/'
pickle.dump([train_sizes,train_score, test_score], open(figpath + 'Learning_curve_RFdiff' + str(Output_ind1) + '_' + str(Output_ind2) + '.pkl', 'wb'))
# a,b,c = pickle.load(open(figpath +'Learning_curve'+ '.pkl', 'rb'))
