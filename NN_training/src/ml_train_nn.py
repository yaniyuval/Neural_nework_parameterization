import numpy as np
# import sknn_jgd.mlp #uncomment if I want to run NN
import time
from sklearn.ensemble import RandomForestRegressor
from src.ml_io import write_netcdf_rf
from src.ml_io import write_netcdf_nn
import src.ml_load as ml_load
import pickle
import src.ml_plot_nn as ml_plot_nn
import os
import math

# For NN pytorch training:
import torch
from torch.autograd import Variable
import torch.nn.functional as F
import torch.utils.data as Data
import torchvision
from torch import nn, optim
#from time import time
from netCDF4 import Dataset
import netCDF4
from sklearn.metrics import r2_score


# ---  build random forest or neural net  ---
def train_wrapper(f_ppi, o_ppi, training_expt, input_vert_dim, output_vert_dim,
                  input_vert_vars, output_vert_vars, flag_dict,
                  do_nn=False, n_iter=None, do_train=True, 
                  no_cos=True, use_rh=False,
                  max_z=40000.0,  
                  rain_only=False, 
                  n_trn_exs=None,  
                  plot_training_results=False,
                  n_layers=2, batch_norm = True, do_wind_input = False, do_diffusion=True,
                  scale_level=False,rewight_outputs = False, weight_list = [1,1],is_cheyenne=False,only_plot = False,
                  n_in = 91, n_out = 210, dropoff= 0.02, output_extreme_flag = False):
    """Loads training data and trains and stores estimator

    Args:
        f_ppi (dict): The type of preprocessing to do to the features (inputs)
        o_ppi (dict): The type of preprocessing to do to the targets (outputs)
        n_iter (int): Number of iterations
        n_stable (int): Number of iterations after stability reached
        max_z (float): Don't train on data above this level
        weight_decay (float): Regularization strength. 0 is no regularization
        rain_only (bool): Only train on precipitating examples
        n_trn_exs (int): Number of training examples to learn on
        do_nn (bool): Use an ANN instead of a random forest 
        no_cos (bool): If true, don't weight by cosine(latitude)
        min_samples_leaf (int): minimum samples per leaf
        plot_training_results (bool): Whether to also plot the model on training data
        use_rh (bool): use generalized relative humidity instead of total non-precip water as feature
        do_train (bool): whether to train (just plot the results if false)
    Returns:
        str: String id of trained NN
    """
    # Load data (note LoadData seeds the random number generator)

    if not only_plot:
        datadir, trainfile, testfile, pp_str = ml_load.GetDataPath_nn(training_expt, wind_input=do_wind_input,
                                                                   is_cheyenne=is_cheyenne)

        f, o, y, z, rho, p, weight_list = ml_load.LoadData(trainfile, max_z, input_vert_vars=input_vert_vars, output_vert_vars=output_vert_vars,
                                              rain_only=rain_only, n_trn_exs=n_trn_exs, no_cos=no_cos, use_rh=use_rh,wind_input = do_wind_input,
                                              exclusion_flag=flag_dict['exclusion_flag'],ind1_exc=flag_dict['ind1_exc'],ind2_exc=flag_dict['ind2_exc'],rewight_outputs=rewight_outputs)

        #load test data
        tf, to, ty, tz, trho, tp, tweight_list = ml_load.LoadData(testfile, max_z, input_vert_vars=input_vert_vars, output_vert_vars=output_vert_vars,
                                                    rain_only=rain_only, n_trn_exs=n_trn_exs, no_cos=no_cos, use_rh=use_rh,rewight_outputs=rewight_outputs)

        # Scale data (both train and test)
        f_pp, f_scl, tf_scl, o_pp, o_scl, to_scl, pp_str = PreprocessData_tr_ts(f_ppi, f, tf, o_ppi, o, to, pp_str,
                                                                                n_trn_exs, z, input_vert_dim, input_vert_vars,
                                                                                output_vert_dim, output_vert_vars,scale_level,
                                                                                rewight_outputs=rewight_outputs,weight_list=weight_list) #Yani TO DO!!!




        # Either build a random forest or build a neural netowrk
        if do_nn:
            est, est_str = BuildNN(pp_str, n_in, n_out, n_layers,dropoff, batch_norm) # Should include in

        if output_extreme_flag:
            print('Removing extremes - the plots are not done on the correct data set I think - to take care if choosing to remove extremes... ')
            f_scl, tf_scl, o_scl, to_scl, indices_rm_tr, indices_rm_test = remove_extremes(f_scl, tf_scl, o_scl, to_scl)
            est_str = est_str + '_rm_ex'

        # Print details about the ML algorithm we are using
        print(est_str + ' Using ' + str(f.shape[0]) + ' training examples with ' +
              str(f.shape[1]) + ' input features and ' + str(o.shape[1]) +
              ' output targets')


        # Train the estimator
        if do_nn:
            start1 = time.time()
            est, path_nn_stage2, est_str = train_nn(est, est_str, f_scl, o_scl, tf_scl, to_scl, output_vert_dim, output_vert_vars)
            end1 = time.time()
            print("The training time in seconds is:")
            print(end1 - start1)
            #exit()
            save_nn(est, est_str, n_layers, o_pp, f_pp, f_ppi, o_ppi, y, z, p, rho, batch_norm=batch_norm)
        # else:
        #     if do_train:
        #         est, est_errors, train_score, test_score = train_est(est, est_str, f_scl, o_scl, tf_scl, to_scl, do_nn)
        #
        #         est_str = est_str + 'te' + str(int(str(test_score)[2:4])) + '_tr' + str(int(str(train_score)[2:4]))
        #     if flag_dict['exclusion_flag']:
        #         est_str = est_str + str(flag_dict['ind1_exc']) + str(flag_dict['ind2_exc'])
        #   # Save the estimator to access it later
        #     save_est(est, est_str, est_errors, f_ppi, o_ppi, f_pp, o_pp, y, z, p, rho, train_score, test_score,is_cheyenne,exclusion_flag=flag_dict['exclusion_flag'],ind1_exc=flag_dict['ind1_exc'],ind2_exc=flag_dict['ind2_exc'])
        # # Write a netcdf file for the gcm
        # if do_nn:
        #     write_netcdf_nn(est_str, trainfile, rain_only, no_cos, use_rh,is_cheyenne)
        # else:
        #     write_netcdf_rf(est_str, trainfile, output_vert_vars, output_vert_dim, rain_only, no_cos, use_rh,scale_level,
        #                      rewight_outputs=rewight_outputs,weight_list=weight_list,is_cheyenne=is_cheyenne)#,exclusion_flag=flag_dict['exclusion_flag'],ind1_exc=flag_dict['ind1_exc'],ind2_exc=flag_dict['ind2_exc'])



        # Plot figures with testing data using all of it
    if only_plot:
        trainfile = '/glade/scratch/janniy/mldata_tmp/training_data/qobsTTFFFFFTF26TTTFTF48TFFFFFTFTFFF40FFTFTTF4848_training_x_no_subsampling.pkl'
        testfile = '/glade/scratch/janniy/mldata_tmp/training_data/qobsTTFFFFFTF26TTTFTF48TFFFFFTFTFFF40FFTFTTF4848_testing_x_no_subsampling.pkl'
        est_str = 'qobsTTFFFFFTF26TTTFTF48TFFTFTFFF40FFTFTTF4848FF_X01_F-NoSc_O-Stan_Ntr5000000_Nte972360_F_Tin_qin_qpin_latin_O_Tout_qout_qpout_RF_NTr10_MinS20max_d27_maxzinf_nocos_te50_tr54'
        #Momentum
        trainfile = '/glade/scratch/janniy/mldata_tmp/training_data/qobsFFTTFFTFF0TFTFTF48FFFFFFFFFFFF1648FFFFTTF00_training_x_no_subsampling.pkl'
        testfile = '/glade/scratch/janniy/mldata_tmp/training_data/qobsFFTTFFTFF0TFTFTF48FFFFFFFFFFFF1648FFFFTTF00_testing_x_no_subsampling.pkl'
        est_str = 'qobsFFTTFFTFF0FFTFTF48FFFFFFFFF1648FFFFTTF00FF_X01_F-NoSc_O-Stan_Ntr5000000_Nte607770_F_Tin_qin_qpin_uin_vin_minusSH_win_surf_wind_O_u_corr_RF_NTr10_MinS20max_d27_maxzinf_nocos_te13_tr19'

    if only_plot:
        figpath = '/glade/scratch/janniy/figs_tmp_xy' + est_str + '/'
    else:
        figpath = './figs/' + est_str + '/'

    figpath = '/glade/scratch/janniy/NN_offline_figs/' + est_str + '/'



    if do_nn:
        print('To DO later')
        ml_plot_nn.PlotAllFigs_nn(path_nn_stage2, est_str, testfile, do_nn, figpath, input_vert_vars, output_vert_vars,input_vert_dim,output_vert_dim, rain_only=rain_only,
                           n_trn_exs=n_trn_exs, no_cos=no_cos, use_rh=use_rh, wind_input = do_wind_input,scale_per_column=scale_level,
                            rewight_outputs=rewight_outputs,weight_list=weight_list,is_cheyenne=is_cheyenne,exclusion_flag=flag_dict['exclusion_flag'],ind1_exc=flag_dict['ind1_exc'],ind2_exc=flag_dict['ind2_exc'])

    if plot_training_results: # note use n_trn_exs here as training data
        if do_nn:
            print('To DO later')
        else:
            figpath = figpath + 'training_data/'
            ml_plot.PlotAllFigs(est_str, trainfile, do_nn, figpath, input_vert_vars, output_vert_vars,input_vert_dim,output_vert_dim,
                                rain_only=rain_only, n_trn_exs=n_trn_exs,
                                no_cos=no_cos, use_rh=use_rh, wind_input = do_wind_input,
                                rewight_outputs=rewight_outputs,weight_list=weight_list,is_cheyenne=is_cheyenne,exclusion_flag=flag_dict['exclusion_flag'],ind1_exc=flag_dict['ind1_exc'],ind2_exc=flag_dict['ind2_exc'])
    return est_str



def PreprocessData_tr_ts(f_ppi, f, test_f, o_ppi, o, test_o, pp_str, n_trn_exs, z, input_vert_dim,
                         input_vert_vars, output_vert_dim, output_vert_vars,scale_per_column,
                         rewight_outputs=False,weight_list=[1,1]):
    """Transform data according to input preprocessor requirements and make
    make preprocessor string for saving"""

    #Preprocessing train features
    f_dict = ml_load.unpack_list(f, input_vert_vars, input_vert_dim)
    f_pp_dict = ml_load.init_pp_generalized(f_ppi,f_dict,input_vert_vars,scale_per_column)
    f_dict = ml_load.transform_data_generalized(f_ppi, f_pp_dict, f_dict, input_vert_vars, z,scale_per_column,rewight_outputs=False) #For random forest this is not necessary
    f = ml_load.pack_list(f_dict,input_vert_vars)

    #Preprocessing test features
    t_f_dict = ml_load.unpack_list(test_f, input_vert_vars, input_vert_dim,scale_per_column)
    t_f_dict = ml_load.transform_data_generalized(f_ppi, f_pp_dict, t_f_dict, input_vert_vars, z,scale_per_column,rewight_outputs=False) #For random forest this is not necessary
    t_f = ml_load.pack_list(t_f_dict, input_vert_vars)

    #Preprocessing train output
    o_dict = ml_load.unpack_list(o, output_vert_vars, output_vert_dim)
    o_pp_dict = ml_load.init_pp_generalized(o_ppi,o_dict,output_vert_vars,scale_per_column=False)
    if rewight_outputs:
        for ind, name in enumerate(output_vert_vars,start=0):
        #     list_res = [1.0,2.25,4.0,4.0,1.0]
            o_pp_dict[name].var_ = o_pp_dict[name].var_/(weight_list[ind]**2)
            print('rescaling output!!!')

    print('Note - for outputs - no scale per column')
    # rewight_outputs =True
    # print('RESCALING OUT111!')
    o_dict = ml_load.transform_data_generalized(o_ppi, o_pp_dict, o_dict, output_vert_vars, z,scale_per_column=False,rewight_outputs=rewight_outputs,weight_list=weight_list) #For random forest this is not necessary
    o = ml_load.pack_list(o_dict, output_vert_vars)

    #Preprocessing test output
    t_o_dict = ml_load.unpack_list(test_o, output_vert_vars, output_vert_dim,scale_per_column)
    t_o_dict = ml_load.transform_data_generalized(o_ppi, o_pp_dict, t_o_dict, output_vert_vars, z,scale_per_column,
                                                  rewight_outputs=rewight_outputs,weight_list=weight_list) #For random forest this is not necessary
    t_o = ml_load.pack_list(t_o_dict, output_vert_vars)

    #output string
    pp_str =pp_str + 'F-' + f_ppi['name'][0:4] + '_'
    pp_str = pp_str + 'O-' + o_ppi['name'][0:4] + '_'
    # Add number of training examples to string
    pp_str = pp_str + 'Ntr' + str(f.shape[0]) + '_'
    pp_str = pp_str + 'Nte' + str(t_f.shape[0]) + '_'
    pp_str = pp_str + 'F_'
    for i in range(len(input_vert_dim)):
        pp_str = pp_str + input_vert_vars[i] + '_'
    pp_str = pp_str + 'O_'
    for i in range(len(output_vert_dim)):
        pp_str = pp_str  + output_vert_vars[i] + '_'


    return f_pp_dict,f, t_f, o_pp_dict,o, t_o, pp_str


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

# def CatchRegularization(weight_decay):
#     """scikit-neuralnetwork seems to have a bug if regularization is set to zero"""
#     if weight_decay > 0.0:
#         regularize = 'L2'
#     else:
#         regularize = None
#     return regularize

def UpdateName(no_cos, use_rh, rain_only, est_str):
    if no_cos:
        est_str = est_str + '_nocos_'
    if use_rh:
        est_str = est_str + '_rh'
    if rain_only:
        est_str = est_str + '_rain'

    return est_str

def save_est(est, est_str, est_errors, f_ppi, o_ppi, f_pp, o_pp, y, z, p, rho, train_score, test_score,is_cheyenne=False,exclusion_flag=False,ind1_exc=0,ind2_exc=0):
    """Save estimator"""
    if is_cheyenne == False:  # On aimsir/esker
        base_dir = '/net/aimsir/archive1/janniy/'
    else:
        base_dir = '/glade/scratch/janniy/'

    if not os.path.exists(base_dir + 'mldata_tmp/regressors/'):
        os.makedirs(base_dir + 'mldata_tmp/regressors/')
    # if exclusion_flag:
    #     pickle.dump([est, est_str, est_errors, f_ppi, o_ppi, f_pp, o_pp, y, z, p, rho], open(base_dir + 'mldata_tmp/regressors/' + est_str+ str(ind1_exc)+ str(ind2_exc) + '.pkl', 'wb'))
    # else:
    pickle.dump([est, est_str, est_errors, f_ppi, o_ppi, f_pp, o_pp, y, z, p, rho],
                open(base_dir + 'mldata_tmp/regressors/' + est_str + '.pkl', 'wb'))

def store_stats(i, avg_train_error, best_train_error, avg_valid_error,
                best_valid_error, avg_train_obj_error, best_train_obj_error,
                **_):
    if i == 1:
        global errors_stored
        errors_stored = []
    errors_stored.append((avg_train_error, best_train_error,
                          avg_valid_error, best_valid_error,
                          avg_train_obj_error, best_train_obj_error))


def BuildNN(pp_str, n_in, n_out, n_layers =2, dropoff=0.0, batch_norm = True):

    """Builds an NN using pytorch
    Currently only two options - 2 layers or 5 layers
    """
    if batch_norm == True:
        print('Using batch normalization NNs')
        if n_layers == 2:
            est =  Net_ANN(n_in, n_out)
        elif n_layers == 5:
            est = Net_ANN_5(n_in, n_out)
        else:
            raise Exception('Can only deal with DNN with 2 or 5 layers (with BN)')
    else:
        print('No batch normalization')
        if n_layers == 2:
            est =  Net_ANN_no_BN(n_in, n_out, neurons=4096)
        elif n_layers == 3:
            est = Net_ANN_3_no_BN(n_in, n_out, neurons=128)
        elif n_layers == 4:
            est = Net_ANN_4_no_BN(n_in, n_out, neurons=128)
        elif n_layers == 5:
            est = Net_ANN_5_no_BN(n_in, n_out, neurons=128)
        elif n_layers == 6:
            est = Net_ANN_6_no_BN(n_in, n_out, neurons=128)
        else:
            raise Exception('Can only deal with DNN with 2 or 5 layers (without BN)')
    # Construct name
    est_str = pp_str
    # Add the number of iterations too
    est_str = est_str + 'NN_layers' + str(n_layers) + 'in'+str(n_in) + 'out' + str(n_out)+ '_BN_' + str(batch_norm)[0]
    return est, est_str


def Load_NN(pp_str, n_in, n_out, n_layers =2, dropoff=0.0, batch_norm = True):
    BuildNN(pp_str, n_in, n_out, n_layers, dropoff, batch_norm)


def BuildRandomForest(max_z, n_trees, min_samples_leaf, pp_str, max_depth, do_diffusion):

    est = RandomForestRegressor(n_estimators=n_trees, min_samples_leaf = min_samples_leaf, max_features = 1.0/3.0, max_depth=max_depth, n_jobs=10, random_state = 123)
    est_str = pp_str + 'RF_NTr' + str(n_trees) + '_MinS' + str(min_samples_leaf) + 'max_d' + str(max_depth)  + '_maxz' + str(max_z)

    return est, est_str


def train_est(est, est_str, f_scl, o_scl, tf_scl, to_scl, do_nn):
    """Train estimator"""

    # Initialize
    start = time.time()

    # Train the model using training data
    est.fit(f_scl, o_scl)
    train_int_score = max(int(np.ceil(f_scl.shape[0] / 10)), 10000)
    train_score = est.score(f_scl[0:train_int_score, :], o_scl[0:train_int_score, :])
    test_score = est.score(tf_scl, to_scl)

    end = time.time()
    print("Training Score: {:.4f} for Model {:s} ({:.1f} seconds)".format(
                                              train_score, est_str, end-start))
    print("Test  Score: {:.4f} for Model {:s} ({:.1f} seconds)".format(
                                              test_score, est_str, end - start))
    if do_nn:
     # This is an n_iter x 4 array...see score_stats
     errors = np.asarray(errors_stored)
    else:
     errors = np.empty(0)

    # Return the fitted models and the scores
    return est, errors, train_score, test_score

# Yani adding for local predictions
def train_est_local(est, est_str, f_scl, o_scl, tf_scl, to_scl, input_vert_vars, output_vert_vars, train_lev_num = 5):
    """Train estimator locally - using only limited number of levels"""
    f_scl_local, o_scl_local = convert_local(f_scl, o_scl,train_lev_num, input_vert_vars, output_vert_vars) # Converting to local data (each sample has the number of levels and a single output.
    tf_scl_local, to_scl_local = convert_local(tf_scl, to_scl, train_lev_num, input_vert_vars, output_vert_vars)


#Yani adding for nn training
def train_nn(net, est_str, f_scl, o_scl, tf_scl,to_scl, output_vert_dim, output_vert_vars,epochs =7, min_lr = 2e-4, max_lr = 2e-3, step_size_up=4000,batch_size = 1024):



    y_train_small_py = torch.from_numpy(o_scl.reshape(-1, o_scl.shape[1])).float()
    # y_train_val_py = torch.from_numpy(y_train_val.reshape(-1, y_train_small.shape[1])).float()

    X_norm_py = torch.from_numpy(f_scl.reshape(-1, f_scl.shape[1])).float()
    X_train_val_norm_py = torch.from_numpy(tf_scl.reshape(-1, tf_scl.shape[1])).float()


    torch_dataset = Data.TensorDataset(X_norm_py, y_train_small_py)
    loader = Data.DataLoader(
        dataset=torch_dataset,
        batch_size=batch_size,
        shuffle=True, num_workers=4,)



    optimizer = optim.Adam(net.parameters(), lr=1e-7)
    loss_func = torch.nn.MSELoss()
    scheduler = optim.lr_scheduler.CyclicLR(optimizer, base_lr=min_lr, max_lr= max_lr, step_size_up=step_size_up, cycle_momentum=False)
    torch.set_num_threads(10)  # This is related to multiple processors but not what I really need.
    for epoch in range(epochs):  # Can use only 4 maybe
        train_model_cyclic(net, loss_func, loader, optimizer, scheduler)
        test_model(net, X_train_val_norm_py, to_scl, output_vert_dim, output_vert_vars)

    PATH = '/glade/u/home/janniy/convection_parametrization/paul_codes/ML-convection_sam_flex_io/NN_saved/' + est_str +  'stage0.pth'
    torch.save(net.state_dict(), PATH)

    #Run a few epochs with lower learning rate.
    scheduler = optim.lr_scheduler.CyclicLR(optimizer, base_lr=min_lr/10, max_lr= max_lr/10, step_size_up=step_size_up, cycle_momentum=False)
    for epoch in range(epochs-2):
        train_model_cyclic(net, loss_func, loader, optimizer, scheduler)
        test_score= test_model(net, X_train_val_norm_py, to_scl, output_vert_dim, output_vert_vars)
    test_score = test_model(net, X_train_val_norm_py, to_scl, output_vert_dim, output_vert_vars)
    train_score = test_model(net, X_norm_py[0:500000,:], o_scl[0:500000,:], output_vert_dim, output_vert_vars)

    # train_score = test_model(net, X_train_val_norm_py, y_train_val)
    est_str = est_str +  '_te' + str(int(str(test_score)[2:4])) + '_tr' + str(int(str(train_score)[2:4]))
    PATH = '/glade/u/home/janniy/convection_parametrization/paul_codes/ML-convection_sam_flex_io/NN_saved/' + est_str +  'stage2.pth'
    torch.save(net.state_dict(), PATH)

    return net, PATH, est_str

def rmse(x, y): return math.sqrt(((x - y) ** 2).mean())

def train_model(net,criterion,trainloader,optimizer,batchsize):
    net.train()
    test_loss = 0
    for step, (batch_x, batch_y) in enumerate(trainloader):  # for each training step
        b_x = Variable(batch_x)
        b_y = Variable(batch_y)
        prediction = net(b_x)  # input x and predict based on x
        loss = criterion(prediction, b_y)  # must be (1. nn output, 2. target)
        optimizer.zero_grad()  # clear gradients for next train
        loss.backward()  # backpropagation, compute gradients
        optimizer.step()  # apply gradients
        test_loss = test_loss + loss.data.numpy()
    print('the loss in this Epoch',test_loss)

def test_model(net,X_train_val_norm_py,y_train_val, output_vert_dim, output_vert_vars):
    net.eval()
    pred_val = net(Variable(X_train_val_norm_py))
    print('RMSE: ',rmse(pred_val.data.numpy(),y_train_val), ' R2:' ,r2_score(y_train_val[:,:], pred_val.data.numpy()[:,:],multioutput='variance_weighted'))
    idim_now = 0
    for dim1,name in zip(output_vert_dim,output_vert_vars):
        print(name + 'R2:',r2_score(y_train_val[:, idim_now:idim_now+dim1], pred_val.data.numpy()[:, idim_now:idim_now+dim1], multioutput='variance_weighted'))
        idim_now = idim_now  + dim1
    return r2_score(y_train_val[:,:], pred_val.data.numpy()[:,:],multioutput='variance_weighted')

def test_model_per_var(net,X_train_val_norm_py,y_train_val, output_vert_dim):
    net.eval()
    pred_val = net(Variable(X_train_val_norm_py))
    j = 0
    for i in output_vert_dim:
        print('R2:' ,r2_score(y_train_val[:,j:j+i], pred_val.data.numpy()[:,j:j+i],multioutput='variance_weighted'))
        j = j + i
    return r2_score(y_train_val[:,:], pred_val.data.numpy()[:,:],multioutput='variance_weighted')


def train_model_cyclic(net,criterion,trainloader,optimizer,scheduler):
#     amp.initialize(model, optimizer, opt_level)
    net.train()
    test_loss = 0
    for step, (batch_x, batch_y) in enumerate(trainloader):  # for each training step
        b_x = Variable(batch_x)
        b_y = Variable(batch_y)
        prediction = net(b_x)  # input x and predict based on x
        loss = criterion(prediction, b_y)  # must be (1. nn output, 2. target)
        optimizer.zero_grad()  # clear gradients for next train
        loss.backward()  # backpropagation, compute gradients
        optimizer.step()  # apply gradients
        optimizer.zero_grad()
        test_loss = test_loss + loss.data.numpy()
        scheduler.step() #defined I think on the epoches...
    print('the loss in this Epoch',test_loss)





def save_nn(net, est_str, n_layers, o_pp, f_pp, f_ppi, o_ppi, y, z, p, rho, batch_norm = True):
    '''Write nn as nc file - only equipped to deal with 2 and 5 layers at the moment'''
    # Set output filename
    base_dir = '/glade/scratch/janniy/'

    output_filename = base_dir + 'mldata_tmp/gcm_regressors/'+est_str+'.nc'
    net2 = net
    net2.eval()

    in_dim = net2.linear1.weight.T.shape[0]
    X_mean = np.zeros(in_dim)
    X_std = np.zeros(in_dim)
    ind_now = 0
    print('check that I am iterating correctly over variables')
    for key, value in f_pp.items():  # Iterate over the different features mean and std. Note that I need to verify that the order of items is what I think it is....
        ind_tmp = ind_now + value.mean_.shape[0]
        X_mean[ind_now:ind_tmp] = value.mean_
        X_std[ind_now:ind_tmp] = np.sqrt(value.var_)  # To get the std I am taking the sqrt of variance
        ind_now = ind_tmp

    # if n_layers == 2:
    #     out_dim = net2.linear2.weight.T.shape[1]
    # elif  n_layers == 5:
    #     out_dim = net2.linear5.weight.T.shape[1]

    Y_mean = np.zeros(len(o_pp))
    Y_std = np.zeros(len(o_pp))
    ind_now = 0
    print('check that I am iterating correctly over outputs')
    for key, value in o_pp.items():  # Iterate over the different features mean and std.
        # ind_tmp = ind_now + value.mean_.shape[0]
        Y_mean[ind_now] = value.mean_
        Y_std[ind_now] = np.sqrt(value.var_)  # To get the std I am taking the sqrt of variance
        ind_now = ind_now +1

    if n_layers == 2:
        ncfile = Dataset(output_filename, 'w', format="NETCDF3_CLASSIC")
        ncfile.createDimension('single', 1)
        ncfile.createDimension('N_in', net2.linear1.weight.T.shape[0])
        ncfile.createDimension('N_h1', net2.linear1.weight.T.shape[1])
        ncfile.createDimension('N_out', net2.linear2.weight.T.shape[1])
        ncfile.createDimension('N_out_dim', len(o_pp))
        # Create variable entries in the file
        nc_w1 = ncfile.createVariable('w1', np.dtype('float32').char,
                                      ('N_h1', 'N_in'))  # Reverse dims
        nc_w2 = ncfile.createVariable('w2', np.dtype('float32').char,
                                      ('N_out', 'N_h1'))
        nc_b1 = ncfile.createVariable('b1', np.dtype('float32').char,
                                      ('N_h1'))
        nc_b2 = ncfile.createVariable('b2', np.dtype('float32').char,
                                      ('N_out'))

        if batch_norm:
            nc_batch_mean = ncfile.createVariable('batch_mean',
                                                  np.dtype('float32').char, ('N_h1'))
            nc_batch_stnd = ncfile.createVariable('batch_stnd',
                                                  np.dtype('float32').char, ('N_h1'))

            nc_batch_weight = ncfile.createVariable('batch_weight',
                                                    np.dtype('float32').char, ('N_h1'))
            nc_batch_bias = ncfile.createVariable('batch_bias',
                                              np.dtype('float32').char, ('N_h1'))

        nc_oscale_mean = ncfile.createVariable('oscale_mean',
                                               np.dtype('float32').char, ('N_out_dim'))
        nc_oscale_stnd = ncfile.createVariable('oscale_stnd',
                                               np.dtype('float32').char, ('N_out_dim'))

        nc_fscale_mean = ncfile.createVariable('fscale_mean',
                                               np.dtype('float32').char, ('N_in'))
        nc_fscale_stnd = ncfile.createVariable('fscale_stnd',
                                               np.dtype('float32').char, ('N_in'))

        nc_w1[:] = net2.linear1.weight.data.numpy()
        nc_w2[:] = net2.linear2.weight.data.numpy()
        nc_b1[:] = net2.linear1.bias.data.numpy()
        nc_b2[:] = net2.linear2.bias.data.numpy()

        if batch_norm:
            print('saving NC NN with BN')
            nc_batch_mean[:] = net2.dense1_bn.running_mean.data.numpy()
            nc_batch_stnd[:] = torch.sqrt(net2.dense1_bn.running_var).data.numpy()
            nc_batch_weight[:] = net2.dense1_bn.weight.data.numpy()
            nc_batch_bias[:] = net2.dense1_bn.bias.data.numpy()

        nc_oscale_mean[:] = Y_mean
        nc_oscale_stnd[:] = Y_std

        nc_fscale_mean[:] = X_mean
        nc_fscale_stnd[:] = X_std

        ncfile.description = 'NN flux Created with ml_train_nn'
        ncfile.close()

    elif n_layers == 3:
        ncfile = Dataset(output_filename, 'w', format="NETCDF3_CLASSIC")
        ncfile.createDimension('single', 1)
        ncfile.createDimension('N_in', net2.linear1.weight.T.shape[0])
        ncfile.createDimension('N_h1', net2.linear1.weight.T.shape[1])
        ncfile.createDimension('N_h2', net2.linear2.weight.T.shape[1])
        ncfile.createDimension('N_out', net2.linear3.weight.T.shape[1])
        ncfile.createDimension('N_out_dim', len(o_pp))
        # Create variable entries in the file
        nc_w1 = ncfile.createVariable('w1', np.dtype('float32').char,
                                      ('N_h1', 'N_in'))  # Reverse dims

        nc_w2 = ncfile.createVariable('w2', np.dtype('float32').char,
                                      ('N_h2', 'N_h1'))  # Reverse dims
        nc_w3 = ncfile.createVariable('w3', np.dtype('float32').char,
                                      ('N_out', 'N_h2'))  # Reverse dims

        nc_b1 = ncfile.createVariable('b1', np.dtype('float32').char,
                                      ('N_h1'))
        nc_b2 = ncfile.createVariable('b2', np.dtype('float32').char,
                                      ('N_h2'))
        nc_b3 = ncfile.createVariable('b3', np.dtype('float32').char,
                                      ('N_out'))

        if batch_norm:
            nc_batch_mean = ncfile.createVariable('batch_mean',
                                                  np.dtype('float32').char, ('N_h1'))
            nc_batch_stnd = ncfile.createVariable('batch_stnd',
                                                  np.dtype('float32').char, ('N_h1'))

            nc_batch_weight = ncfile.createVariable('batch_weight',
                                                    np.dtype('float32').char, ('N_h1'))
            nc_batch_bias = ncfile.createVariable('batch_bias',
                                                  np.dtype('float32').char, ('N_h1'))

            nc_batch_mean2 = ncfile.createVariable('batch_mean2',
                                                   np.dtype('float32').char, ('N_h2'))
            nc_batch_stnd2 = ncfile.createVariable('batch_stnd2',
                                                   np.dtype('float32').char, ('N_h2'))
            nc_batch_weight2 = ncfile.createVariable('batch_weight2',
                                                     np.dtype('float32').char, ('N_h2'))
            nc_batch_bias2 = ncfile.createVariable('batch_bias2',
                                                   np.dtype('float32').char, ('N_h2'))

        nc_oscale_mean = ncfile.createVariable('oscale_mean',
                                               np.dtype('float32').char, ('N_out_dim'))
        nc_oscale_stnd = ncfile.createVariable('oscale_stnd',
                                               np.dtype('float32').char, ('N_out_dim'))

        nc_fscale_mean = ncfile.createVariable('fscale_mean',
                                               np.dtype('float32').char, ('N_in'))
        nc_fscale_stnd = ncfile.createVariable('fscale_stnd',
                                               np.dtype('float32').char, ('N_in'))

        #     nc_w1[:] = net2.linear1.weight.data.numpy().T
        #     nc_w2[:] = net2.linear2.weight.data.numpy().T
        nc_w1[:] = net2.linear1.weight.data.numpy()
        nc_w2[:] = net2.linear2.weight.data.numpy()
        nc_w3[:] = net2.linear3.weight.data.numpy()

        nc_b1[:] = net2.linear1.bias.data.numpy()
        nc_b2[:] = net2.linear2.bias.data.numpy()
        nc_b3[:] = net2.linear3.bias.data.numpy()

        if batch_norm:
            print('saving NC NN with BN')
            nc_batch_mean[:] = net2.dense1_bn.running_mean.data.numpy()
            nc_batch_stnd[:] = torch.sqrt(net2.dense1_bn.running_var).data.numpy()
            nc_batch_weight[:] = net2.dense1_bn.weight.data.numpy()
            nc_batch_bias[:] = net2.dense1_bn.bias.data.numpy()

            nc_batch_mean2[:] = net2.dense2_bn.running_mean.data.numpy()
            nc_batch_stnd2[:] = torch.sqrt(net2.dense2_bn.running_var).data.numpy()
            nc_batch_weight2[:] = net2.dense2_bn.weight.data.numpy()
            nc_batch_bias2[:] = net2.dense2_bn.bias.data.numpy()

        nc_oscale_mean[:] = Y_mean
        nc_oscale_stnd[:] = Y_std

        nc_fscale_mean[:] = X_mean
        nc_fscale_stnd[:] = X_std

        ncfile.description = 'NN flux Created with ml_train_nn'
        ncfile.close()


    #######
    elif n_layers == 4:
        ncfile = Dataset(output_filename, 'w', format="NETCDF3_CLASSIC")
        ncfile.createDimension('single', 1)
        ncfile.createDimension('N_in', net2.linear1.weight.T.shape[0])
        ncfile.createDimension('N_h1', net2.linear1.weight.T.shape[1])
        ncfile.createDimension('N_h2', net2.linear2.weight.T.shape[1])
        ncfile.createDimension('N_h3', net2.linear3.weight.T.shape[1])
        ncfile.createDimension('N_out', net2.linear4.weight.T.shape[1])
        ncfile.createDimension('N_out_dim', len(o_pp))
        # Create variable entries in the file
        nc_w1 = ncfile.createVariable('w1', np.dtype('float32').char,
                                      ('N_h1', 'N_in'))  # Reverse dims

        nc_w2 = ncfile.createVariable('w2', np.dtype('float32').char,
                                      ('N_h2', 'N_h1'))  # Reverse dims
        nc_w3 = ncfile.createVariable('w3', np.dtype('float32').char,
                                      ('N_h3', 'N_h2'))  # Reverse dims

        nc_w4 = ncfile.createVariable('w4', np.dtype('float32').char,
                                      ('N_out', 'N_h3'))  # Reverse dims

        nc_b1 = ncfile.createVariable('b1', np.dtype('float32').char,
                                      ('N_h1'))
        nc_b2 = ncfile.createVariable('b2', np.dtype('float32').char,
                                      ('N_h2'))
        nc_b3 = ncfile.createVariable('b3', np.dtype('float32').char,
                                      ('N_h3'))
        nc_b4 = ncfile.createVariable('b4', np.dtype('float32').char,
                                      ('N_out'))
        if batch_norm:
            nc_batch_mean = ncfile.createVariable('batch_mean',
                                                  np.dtype('float32').char, ('N_h1'))
            nc_batch_stnd = ncfile.createVariable('batch_stnd',
                                                  np.dtype('float32').char, ('N_h1'))

            nc_batch_weight = ncfile.createVariable('batch_weight',
                                                    np.dtype('float32').char, ('N_h1'))
            nc_batch_bias = ncfile.createVariable('batch_bias',
                                                  np.dtype('float32').char, ('N_h1'))

            nc_batch_mean2 = ncfile.createVariable('batch_mean2',
                                                   np.dtype('float32').char, ('N_h2'))
            nc_batch_stnd2 = ncfile.createVariable('batch_stnd2',
                                                   np.dtype('float32').char, ('N_h2'))
            nc_batch_weight2 = ncfile.createVariable('batch_weight2',
                                                     np.dtype('float32').char, ('N_h2'))
            nc_batch_bias2 = ncfile.createVariable('batch_bias2',
                                                   np.dtype('float32').char, ('N_h2'))

            nc_batch_mean3 = ncfile.createVariable('batch_mean3',
                                                   np.dtype('float32').char, ('N_h3'))
            nc_batch_stnd3 = ncfile.createVariable('batch_stnd3',
                                                   np.dtype('float32').char, ('N_h3'))
            nc_batch_weight3 = ncfile.createVariable('batch_weight3',
                                                     np.dtype('float32').char, ('N_h3'))
            nc_batch_bias3 = ncfile.createVariable('batch_bias3',
                                                   np.dtype('float32').char, ('N_h3'))

        nc_oscale_mean = ncfile.createVariable('oscale_mean',
                                               np.dtype('float32').char, ('N_out_dim'))
        nc_oscale_stnd = ncfile.createVariable('oscale_stnd',
                                               np.dtype('float32').char, ('N_out_dim'))


        nc_fscale_mean = ncfile.createVariable('fscale_mean',
                                               np.dtype('float32').char, ('N_in'))
        nc_fscale_stnd = ncfile.createVariable('fscale_stnd',
                                               np.dtype('float32').char, ('N_in'))

        #     nc_w1[:] = net2.linear1.weight.data.numpy().T
        #     nc_w2[:] = net2.linear2.weight.data.numpy().T
        nc_w1[:] = net2.linear1.weight.data.numpy()
        nc_w2[:] = net2.linear2.weight.data.numpy()
        nc_w3[:] = net2.linear3.weight.data.numpy()
        nc_w4[:] = net2.linear4.weight.data.numpy()

        nc_b1[:] = net2.linear1.bias.data.numpy()
        nc_b2[:] = net2.linear2.bias.data.numpy()
        nc_b3[:] = net2.linear3.bias.data.numpy()
        nc_b4[:] = net2.linear4.bias.data.numpy()

        if batch_norm:
            print('saving NC NN with BN')
            nc_batch_mean[:] = net2.dense1_bn.running_mean.data.numpy()
            nc_batch_stnd[:] = torch.sqrt(net2.dense1_bn.running_var).data.numpy()
            nc_batch_weight[:] = net2.dense1_bn.weight.data.numpy()
            nc_batch_bias[:] = net2.dense1_bn.bias.data.numpy()

            nc_batch_mean2[:] = net2.dense2_bn.running_mean.data.numpy()
            nc_batch_stnd2[:] = torch.sqrt(net2.dense2_bn.running_var).data.numpy()
            nc_batch_weight2[:] = net2.dense2_bn.weight.data.numpy()
            nc_batch_bias2[:] = net2.dense2_bn.bias.data.numpy()

            nc_batch_mean3[:] = net2.dense3_bn.running_mean.data.numpy()
            nc_batch_stnd3[:] = torch.sqrt(net2.dense3_bn.running_var).data.numpy()
            nc_batch_weight3[:] = net2.dense3_bn.weight.data.numpy()
            nc_batch_bias3[:] = net2.dense3_bn.bias.data.numpy()

        nc_oscale_mean[:] = Y_mean
        nc_oscale_stnd[:] = Y_std

        nc_fscale_mean[:] = X_mean
        nc_fscale_stnd[:] = X_std

        ncfile.description = 'NN flux Created with ml_train_nn'
        ncfile.close()

#########
    elif n_layers == 5:
        ncfile = Dataset(output_filename, 'w', format="NETCDF3_CLASSIC")
        ncfile.createDimension('single', 1)
        ncfile.createDimension('N_in', net2.linear1.weight.T.shape[0])
        ncfile.createDimension('N_h1', net2.linear1.weight.T.shape[1])
        ncfile.createDimension('N_h2', net2.linear2.weight.T.shape[1])
        ncfile.createDimension('N_h3', net2.linear3.weight.T.shape[1])
        ncfile.createDimension('N_h4', net2.linear4.weight.T.shape[1])
        ncfile.createDimension('N_out', net2.linear5.weight.T.shape[1])
        ncfile.createDimension('N_out_dim', len(o_pp))
        # Create variable entries in the file
        nc_w1 = ncfile.createVariable('w1', np.dtype('float32').char,
                                      ('N_h1', 'N_in'))  # Reverse dims

        nc_w2 = ncfile.createVariable('w2', np.dtype('float32').char,
                                      ('N_h2', 'N_h1'))  # Reverse dims
        nc_w3 = ncfile.createVariable('w3', np.dtype('float32').char,
                                      ('N_h3', 'N_h2'))  # Reverse dims
        nc_w4 = ncfile.createVariable('w4', np.dtype('float32').char,
                                      ('N_h4', 'N_h3'))  # Reverse dims
        nc_w5 = ncfile.createVariable('w5', np.dtype('float32').char,
                                      ('N_out', 'N_h4'))  # Reverse dims

        nc_b1 = ncfile.createVariable('b1', np.dtype('float32').char,
                                      ('N_h1'))
        nc_b2 = ncfile.createVariable('b2', np.dtype('float32').char,
                                      ('N_h2'))
        nc_b3 = ncfile.createVariable('b3', np.dtype('float32').char,
                                      ('N_h3'))
        nc_b4 = ncfile.createVariable('b4', np.dtype('float32').char,
                                      ('N_h4'))
        nc_b5 = ncfile.createVariable('b5', np.dtype('float32').char,
                                      ('N_out'))
        if batch_norm:
            nc_batch_mean = ncfile.createVariable('batch_mean',
                                                  np.dtype('float32').char, ('N_h1'))
            nc_batch_stnd = ncfile.createVariable('batch_stnd',
                                                  np.dtype('float32').char, ('N_h1'))

            nc_batch_weight = ncfile.createVariable('batch_weight',
                                                    np.dtype('float32').char, ('N_h1'))
            nc_batch_bias = ncfile.createVariable('batch_bias',
                                                  np.dtype('float32').char, ('N_h1'))

            nc_batch_mean2 = ncfile.createVariable('batch_mean2',
                                                   np.dtype('float32').char, ('N_h2'))
            nc_batch_stnd2 = ncfile.createVariable('batch_stnd2',
                                                   np.dtype('float32').char, ('N_h2'))
            nc_batch_weight2 = ncfile.createVariable('batch_weight2',
                                                     np.dtype('float32').char, ('N_h2'))
            nc_batch_bias2 = ncfile.createVariable('batch_bias2',
                                                   np.dtype('float32').char, ('N_h2'))

            nc_batch_mean3 = ncfile.createVariable('batch_mean3',
                                                   np.dtype('float32').char, ('N_h3'))
            nc_batch_stnd3 = ncfile.createVariable('batch_stnd3',
                                                   np.dtype('float32').char, ('N_h3'))
            nc_batch_weight3 = ncfile.createVariable('batch_weight3',
                                                     np.dtype('float32').char, ('N_h3'))
            nc_batch_bias3 = ncfile.createVariable('batch_bias3',
                                                   np.dtype('float32').char, ('N_h3'))

            nc_batch_mean4 = ncfile.createVariable('batch_mean4',
                                                   np.dtype('float32').char, ('N_h4'))
            nc_batch_stnd4 = ncfile.createVariable('batch_stnd4',
                                                   np.dtype('float32').char, ('N_h4'))
            nc_batch_weight4 = ncfile.createVariable('batch_weight4',
                                                     np.dtype('float32').char, ('N_h4'))
            nc_batch_bias4 = ncfile.createVariable('batch_bias4',
                                                   np.dtype('float32').char, ('N_h4'))

        nc_oscale_mean = ncfile.createVariable('oscale_mean',
                                               np.dtype('float32').char, ('N_out_dim'))
        nc_oscale_stnd = ncfile.createVariable('oscale_stnd',
                                               np.dtype('float32').char, ('N_out_dim'))


        nc_fscale_mean = ncfile.createVariable('fscale_mean',
                                               np.dtype('float32').char, ('N_in'))
        nc_fscale_stnd = ncfile.createVariable('fscale_stnd',
                                               np.dtype('float32').char, ('N_in'))

        #     nc_w1[:] = net2.linear1.weight.data.numpy().T
        #     nc_w2[:] = net2.linear2.weight.data.numpy().T
        nc_w1[:] = net2.linear1.weight.data.numpy()
        nc_w2[:] = net2.linear2.weight.data.numpy()
        nc_w3[:] = net2.linear3.weight.data.numpy()
        nc_w4[:] = net2.linear4.weight.data.numpy()
        nc_w5[:] = net2.linear5.weight.data.numpy()

        nc_b1[:] = net2.linear1.bias.data.numpy()
        nc_b2[:] = net2.linear2.bias.data.numpy()
        nc_b3[:] = net2.linear3.bias.data.numpy()
        nc_b4[:] = net2.linear4.bias.data.numpy()
        nc_b5[:] = net2.linear5.bias.data.numpy()

        if batch_norm:
            print('saving NC NN with BN')
            nc_batch_mean[:] = net2.dense1_bn.running_mean.data.numpy()
            nc_batch_stnd[:] = torch.sqrt(net2.dense1_bn.running_var).data.numpy()
            nc_batch_weight[:] = net2.dense1_bn.weight.data.numpy()
            nc_batch_bias[:] = net2.dense1_bn.bias.data.numpy()

            nc_batch_mean2[:] = net2.dense2_bn.running_mean.data.numpy()
            nc_batch_stnd2[:] = torch.sqrt(net2.dense2_bn.running_var).data.numpy()
            nc_batch_weight2[:] = net2.dense2_bn.weight.data.numpy()
            nc_batch_bias2[:] = net2.dense2_bn.bias.data.numpy()

            nc_batch_mean3[:] = net2.dense3_bn.running_mean.data.numpy()
            nc_batch_stnd3[:] = torch.sqrt(net2.dense3_bn.running_var).data.numpy()
            nc_batch_weight3[:] = net2.dense3_bn.weight.data.numpy()
            nc_batch_bias3[:] = net2.dense3_bn.bias.data.numpy()

            nc_batch_mean4[:] = net2.dense4_bn.running_mean.data.numpy()
            nc_batch_stnd4[:] = torch.sqrt(net2.dense4_bn.running_var).data.numpy()
            nc_batch_weight4[:] = net2.dense4_bn.weight.data.numpy()
            nc_batch_bias4[:] = net2.dense4_bn.bias.data.numpy()

        nc_oscale_mean[:] = Y_mean
        nc_oscale_stnd[:] = Y_std

        nc_fscale_mean[:] = X_mean
        nc_fscale_stnd[:] = X_std

        ncfile.description = 'NN flux Created with ml_train_nn'
        ncfile.close()

    elif n_layers == 6:
        ncfile = Dataset(output_filename, 'w', format="NETCDF3_CLASSIC")
        ncfile.createDimension('single', 1)
        ncfile.createDimension('N_in', net2.linear1.weight.T.shape[0])
        ncfile.createDimension('N_h1', net2.linear1.weight.T.shape[1])
        ncfile.createDimension('N_h2', net2.linear2.weight.T.shape[1])
        ncfile.createDimension('N_h3', net2.linear3.weight.T.shape[1])
        ncfile.createDimension('N_h4', net2.linear4.weight.T.shape[1])
        ncfile.createDimension('N_h5', net2.linear5.weight.T.shape[1])
        ncfile.createDimension('N_out', net2.linear6.weight.T.shape[1])
        ncfile.createDimension('N_out_dim', len(o_pp))
        # Create variable entries in the file
        nc_w1 = ncfile.createVariable('w1', np.dtype('float32').char,
                                      ('N_h1', 'N_in'))  # Reverse dims

        nc_w2 = ncfile.createVariable('w2', np.dtype('float32').char,
                                      ('N_h2', 'N_h1'))  # Reverse dims
        nc_w3 = ncfile.createVariable('w3', np.dtype('float32').char,
                                      ('N_h3', 'N_h2'))  # Reverse dims
        nc_w4 = ncfile.createVariable('w4', np.dtype('float32').char,
                                      ('N_h4', 'N_h3'))  # Reverse dims
        nc_w5 = ncfile.createVariable('w5', np.dtype('float32').char,
                                      ('N_h5', 'N_h4'))  # Reverse dims
        nc_w6 = ncfile.createVariable('w6', np.dtype('float32').char,
                                      ('N_out', 'N_h5'))  # Reverse dims

        nc_b1 = ncfile.createVariable('b1', np.dtype('float32').char,
                                      ('N_h1'))
        nc_b2 = ncfile.createVariable('b2', np.dtype('float32').char,
                                      ('N_h2'))
        nc_b3 = ncfile.createVariable('b3', np.dtype('float32').char,
                                      ('N_h3'))
        nc_b4 = ncfile.createVariable('b4', np.dtype('float32').char,
                                      ('N_h4'))
        nc_b5 = ncfile.createVariable('b5', np.dtype('float32').char,
                                      ('N_h5'))
        nc_b6 = ncfile.createVariable('b6', np.dtype('float32').char,
                                      ('N_out'))
        if batch_norm:
            raise Exception('No BN with 6 layers!')

        nc_oscale_mean = ncfile.createVariable('oscale_mean',
                                               np.dtype('float32').char, ('N_out_dim'))
        nc_oscale_stnd = ncfile.createVariable('oscale_stnd',
                                               np.dtype('float32').char, ('N_out_dim'))


        nc_fscale_mean = ncfile.createVariable('fscale_mean',
                                               np.dtype('float32').char, ('N_in'))
        nc_fscale_stnd = ncfile.createVariable('fscale_stnd',
                                               np.dtype('float32').char, ('N_in'))

        #     nc_w1[:] = net2.linear1.weight.data.numpy().T
        #     nc_w2[:] = net2.linear2.weight.data.numpy().T
        nc_w1[:] = net2.linear1.weight.data.numpy()
        nc_w2[:] = net2.linear2.weight.data.numpy()
        nc_w3[:] = net2.linear3.weight.data.numpy()
        nc_w4[:] = net2.linear4.weight.data.numpy()
        nc_w5[:] = net2.linear5.weight.data.numpy()
        nc_w6[:] = net2.linear6.weight.data.numpy()

        nc_b1[:] = net2.linear1.bias.data.numpy()
        nc_b2[:] = net2.linear2.bias.data.numpy()
        nc_b3[:] = net2.linear3.bias.data.numpy()
        nc_b4[:] = net2.linear4.bias.data.numpy()
        nc_b5[:] = net2.linear5.bias.data.numpy()
        nc_b6[:] = net2.linear6.bias.data.numpy()


        nc_oscale_mean[:] = Y_mean
        nc_oscale_stnd[:] = Y_std

        nc_fscale_mean[:] = X_mean
        nc_fscale_stnd[:] = X_std

        ncfile.description = 'NN flux Created with ml_train_nn'
        ncfile.close()

    else:
        raise Exception('Can only save DNN with 2 or 5 layers')


    """Save dictionary for nn rescaling and other properies"""

    if not os.path.exists(base_dir + 'mldata_tmp/regressors/'):
        os.makedirs(base_dir + 'mldata_tmp/regressors/')
    est_errors = 0 # This is dummy not to change how we dump.
    pickle.dump([net, est_str, est_errors, f_ppi, o_ppi, f_pp, o_pp, y, z, p, rho],
                open(base_dir + 'mldata_tmp/regressors/' + est_str + '.pkl', 'wb'))

def remove_extremes(f_scl, tf_scl, o_scl,to_scl, drop_criterion_std = 60):
    '''Remove samples with output which is unrealistically large
    I am not sure yet if this is a good idea, but it is possible that some samples are somehow corrupted'''

    o_scl_mean = np.mean(o_scl, axis=0, dtype=np.float64)
    o_scl_std = np.std(o_scl, axis=0, dtype=np.float64)

    o_scl_temp = (o_scl - o_scl_mean) / (o_scl_std + 0.000000000001)
    to_scl_temp = (to_scl - o_scl_mean) / (o_scl_std + 0.000000000001)

    indices_to_drop_y = np.unique(np.nonzero(o_scl_temp[:, :] > drop_criterion_std)[0])
    indices_to_drop_y_val = np.unique(np.nonzero(to_scl_temp[:, :] > drop_criterion_std)[0])

    tot_train_ind_drop = np.unique(indices_to_drop_y)
    tot_val_ind_drop = np.unique(indices_to_drop_y_val)
    print('size of indices Y train', indices_to_drop_y.shape)
    print('size of indices Y val', indices_to_drop_y_val.shape)

    # indices_to_drop_X = np.unique(np.nonzero(f_scl[:, :] > drop_criterion_std)[0])
    # print('size of indices X train', indices_to_drop_X.shape)
    # indices_to_drop_X_val = np.unique(np.nonzero(tf_scl[:, :] > drop_criterion_std)[0])
    # print('size of indices X val', indices_to_drop_X_val.shape)
    #
    # tot_train_ind_drop = np.concatenate((indices_to_drop_X, indices_to_drop_y), axis=0)
    # tot_train_ind_drop = np.unique(tot_train_ind_drop)
    #
    # tot_val_ind_drop = np.concatenate((indices_to_drop_X_val, indices_to_drop_y_val), axis=0)
    # tot_val_ind_drop = np.unique(tot_val_ind_drop)


    # For train:
    f_scl = np.delete(f_scl, tot_train_ind_drop, axis=0)
    o_scl = np.delete(o_scl, tot_train_ind_drop, axis=0)

    # For validation:
    tf_scl = np.delete(tf_scl, tot_val_ind_drop, axis=0)
    to_scl = np.delete(to_scl, tot_val_ind_drop, axis=0)

    return f_scl, tf_scl, o_scl, to_scl, tot_train_ind_drop, tot_val_ind_drop


class Net_ANN(nn.Module):
    def __init__(self,n_in, n_out,dropoff=0.0):
        super(Net_ANN, self).__init__()
        self.linear1 = nn.Linear(n_in, 256)
        self.linear2 = nn.Linear(256, n_out)
        self.dense1_bn = nn.BatchNorm1d(256)
        self.lin_drop = nn.Dropout(dropoff)  # regularization method to prevent overfitting.

    def forward(self, x):
        x = self.dense1_bn(F.relu(self.linear1(x)))
        x = self.lin_drop(x)
        x = self.linear2(x)
        return x

class Net_ANN_no_BN(nn.Module):
    def __init__(self,n_in, n_out, neurons = 128, dropoff=0.0):
        super(Net_ANN_no_BN, self).__init__()
        self.linear1 = nn.Linear(n_in, neurons)
        self.linear2 = nn.Linear(neurons, n_out)
        self.lin_drop = nn.Dropout(dropoff)  # regularization method to prevent overfitting.

    def forward(self, x):
        x = F.relu(self.linear1(x))
        x = self.lin_drop(x)
        x = self.linear2(x)
        return x

class Net_ANN_5(nn.Module):
    def __init__(self,n_in, n_out,dropoff=0.0):
        super(Net_ANN_5, self).__init__()
        self.linear1 = nn.Linear(n_in, 256)
        self.linear2 = nn.Linear(256, 256)
        self.linear3 = nn.Linear(256, 256)
        self.linear4 = nn.Linear(256, 256)
        self.linear5 = nn.Linear(256, n_out)

        self.dense1_bn = nn.BatchNorm1d(256)
        self.dense2_bn = nn.BatchNorm1d(256)
        self.dense3_bn = nn.BatchNorm1d(256)
        self.dense4_bn = nn.BatchNorm1d(256)

        self.lin_drop = nn.Dropout(dropoff)  # regularization method to prevent overfitting.
        self.lin_drop2 = nn.Dropout(dropoff)
        self.lin_drop3 = nn.Dropout(dropoff)
        self.lin_drop4 = nn.Dropout(dropoff)

    def forward(self, x):
        x = self.dense1_bn(F.relu(self.linear1(x)))
        x = self.lin_drop(x)
        x = self.dense2_bn(F.relu(self.linear2(x)))
        x = self.lin_drop2(x)
        x = self.dense3_bn(F.relu(self.linear3(x)))
        x = self.lin_drop3(x)
        x = self.dense4_bn(F.relu(self.linear4(x)))
        x = self.lin_drop4(x)
        x = self.linear5(x)
        return x


class Net_ANN_3_no_BN(nn.Module):
    def __init__(self,n_in, n_out, neurons = 128, dropoff=0.0):
        super(Net_ANN_3_no_BN, self).__init__()
        self.linear1 = nn.Linear(n_in, neurons)
        self.linear2 = nn.Linear(neurons, neurons)
        self.linear3 = nn.Linear(neurons, n_out)

        self.lin_drop = nn.Dropout(dropoff)  # regularization method to prevent overfitting.
        self.lin_drop2 = nn.Dropout(dropoff)

    def forward(self, x):
        x = F.relu(self.linear1(x))
        x = self.lin_drop(x)
        x = F.relu(self.linear2(x))
        x = self.lin_drop2(x)
        x = self.linear3(x)
        return x


class Net_ANN_4_no_BN(nn.Module):
    def __init__(self,n_in, n_out, neurons = 128, dropoff=0.0):
        super(Net_ANN_4_no_BN, self).__init__()
        self.linear1 = nn.Linear(n_in, neurons)
        self.linear2 = nn.Linear(neurons, neurons)
        self.linear3 = nn.Linear(neurons, neurons)
        self.linear4 = nn.Linear(neurons, n_out)

        self.lin_drop = nn.Dropout(dropoff)  # regularization method to prevent overfitting.
        self.lin_drop2 = nn.Dropout(dropoff)
        self.lin_drop3 = nn.Dropout(dropoff)

    def forward(self, x):
        x = F.relu(self.linear1(x))
        x = self.lin_drop(x)
        x = F.relu(self.linear2(x))
        x = self.lin_drop2(x)
        x = F.relu(self.linear3(x))
        x = self.lin_drop3(x)
        x = self.linear4(x)
        return x

class Net_ANN_5_no_BN(nn.Module):
    def __init__(self,n_in, n_out, neurons = 128, dropoff=0.0):
        super(Net_ANN_5_no_BN, self).__init__()
        self.linear1 = nn.Linear(n_in, neurons)
        self.linear2 = nn.Linear(neurons, neurons)
        self.linear3 = nn.Linear(neurons, neurons)
        self.linear4 = nn.Linear(neurons, neurons)
        self.linear5 = nn.Linear(neurons, n_out)

        self.lin_drop = nn.Dropout(dropoff)  # regularization method to prevent overfitting.
        self.lin_drop2 = nn.Dropout(dropoff)
        self.lin_drop3 = nn.Dropout(dropoff)
        self.lin_drop4 = nn.Dropout(dropoff)

    def forward(self, x):
        x = F.relu(self.linear1(x))
        x = self.lin_drop(x)
        x = F.relu(self.linear2(x))
        x = self.lin_drop2(x)
        x = F.relu(self.linear3(x))
        x = self.lin_drop3(x)
        x = F.relu(self.linear4(x))
        x = self.lin_drop4(x)
        x = self.linear5(x)
        return x


class Net_ANN_6_no_BN(nn.Module):
    def __init__(self,n_in, n_out, neurons = 128, dropoff=0.0):
        super(Net_ANN_6_no_BN, self).__init__()
        self.linear1 = nn.Linear(n_in, neurons)
        self.linear2 = nn.Linear(neurons, neurons)
        self.linear3 = nn.Linear(neurons, neurons)
        self.linear4 = nn.Linear(neurons, neurons)
        self.linear5 = nn.Linear(neurons, neurons)
        self.linear6 = nn.Linear(neurons, n_out)

        self.lin_drop = nn.Dropout(dropoff)  # regularization method to prevent overfitting.
        self.lin_drop2 = nn.Dropout(dropoff)
        self.lin_drop3 = nn.Dropout(dropoff)
        self.lin_drop4 = nn.Dropout(dropoff)
        self.lin_drop5 = nn.Dropout(dropoff)

    def forward(self, x):
        x = F.relu(self.linear1(x))
        x = self.lin_drop(x)
        x = F.relu(self.linear2(x))
        x = self.lin_drop2(x)
        x = F.relu(self.linear3(x))
        x = self.lin_drop3(x)
        x = F.relu(self.linear4(x))
        x = self.lin_drop4(x)
        x = F.relu(self.linear5(x))
        x = self.lin_drop5(x)
        x = self.linear6(x)
        return x

class Net_CNN_3_channels(nn.Module):
    def __init__(self,n_in, n_out):
        super(Net_CNN_3_channels, self).__init__()
        self.conv1 = nn.Sequential(
            nn.Conv1d(in_channels=3, out_channels=10, kernel_size=5, padding=2),
            nn.BatchNorm1d(10),
            nn.ReLU(),
            nn.MaxPool1d(2)
        )
        self.conv2 = nn.Sequential(
            nn.Conv1d(in_channels=10, out_channels=20, kernel_size=5, padding=2),
            nn.BatchNorm1d(20),
            nn.ReLU(),
            #             nn.MaxPool1d(2)
        )

        self.linear1 = nn.Linear(20 * (n_in-1)/3/2, 256)
        self.linear2 = nn.Linear(256, n_out)
        self.dense0_bn = nn.BatchNorm1d(20 * (n_in-1)/3/2)
        self.dense1_bn = nn.BatchNorm1d(256)
        self.lin_drop = nn.Dropout(0.02)  # regularization method to prevent overfitting.
        if (n_in-1)/3/2% 1 != 0:
            raise Exception('The number of inputs after convolutional layers is not an integer and not difined well {}'.format((n_in-1)/3/2))

    def forward(self, x):
        x1 = torch.reshape(x[:, :], [-1, 3, (n_in-1)/3])
        x1 = self.conv1(x1)
        x1 = self.conv2(x1)
        x4 = x[:,-1].unsqueeze(1)
        x1 = x1.view(x1.size(0), -1)
        x = torch.cat((x1, x4), dim=1)
                # print('I think I should have bn here as well')
        x = self.dense0_bn(x)
        x = self.dense1_bn(F.relu(self.linear1(x)))
        x = self.lin_drop(x)
        x = self.linear2(x)
        return x

