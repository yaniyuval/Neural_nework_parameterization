import numpy as np
from sklearn import preprocessing, metrics
import scipy.stats
import pickle
import warnings
import src.atmos_physics as atmos_physics
import pandas as pd
from netCDF4 import Dataset
import torch
from torch.autograd import Variable
import torch.nn.functional as F
import torch.utils.data as Data
import torchvision
from torch import nn, optim

def LoadData(filename, max_z, input_vert_vars, output_vert_vars, all_ys=True, ind_y=None, n_trn_exs=None,
             rain_only=False, no_cos=True, verbose=False, use_rh=False, wind_input = False, exclusion_flag=False,
             ind1_exc = 0, ind2_exc=0, rewight_outputs = False):
    """v2 of the script to load data. See prep_convection_output.py for how
       the input filename is generated.

    Args:
      filename:  The file to be loaded.
      max_z:    The topmost model level for which to load data. 
      all_ys:  Logical value for whether to load data from all y's
      ind_y:    If all_ys is false, give the index value for the
                 y at which to load data.
      n_trn_exs: Number of training examples to load. If set to None, or
                 if requested number exceeds max available will load all.
      rain_only:  If true, only return training examples of when it is raining
      no_cos:   If true, don't use cos(lat) weighting for loading training examples
      verbose:   If true, prints some basic stats about training set

    Returns:
      f       : 2-d numpy array of input features (m_training examples x
                n_input features). 
      o       : 2-d numpy array of output targets (m_traning examples x
                n_output targets).
    """
    # Data to read in is n_lev n_y (SH & NH) n_samples
    # Samples are quasi indpendent with only a few (e.g. 10) from each y
    data_l = pickle.load(open(filename, 'rb'))

    varis = input_vert_vars + output_vert_vars
    v = dict()

    for ind,var in enumerate(varis,start=0):
        v[var] = np.float32(data_l[ind]) # I put it as float 32 since some of the python functions that I use have a problem with float16.
        if ind==len(varis)-1:
            y = data_l[ind+1]
            z = data_l[ind + 2]
            p = data_l[ind + 3]
            rho = data_l[ind + 4]
            if rewight_outputs:
                print('weights to rescale outputs.')
                weight_list = data_l[ind + 5]
            else:
                weight_list = [1,1] # Not sure this will be used...

    # Added this to
    if exclusion_flag:
        exclution_lat_list_1 = list(range(90 - ind2_exc, 90 - ind1_exc))  # list(range(0, 90))
        exclution_lat_list_2 = list(range(90 + ind1_exc, 90 + ind2_exc))
        exclution_lat_list = exclution_lat_list_1 + exclution_lat_list_2
        if len(exclution_lat_list)>0:
            y_tot_len = y.shape[0]
            # y2 = np.delete(y, exclution_lat_list, axis=0)
            for ind,var in enumerate(varis,start=0):
                dim_of_lat = v[var].shape.index(y_tot_len)
                v[var] = np.delete(v[var], exclution_lat_list, axis=dim_of_lat)
            print('I chose a subset of y indices for generalization tests.  ')




    # Limit levels to those specified
    ind_z = np.less_equal(z, max_z)
    z = z[ind_z]
    p = p[ind_z]
    rho = rho[ind_z]


    rank_of_vars = len(v[var].shape)
    # Reshape the arrays
    for var in varis:
        # Change shape of data to be n_samp n_z
        if (v[var].shape[0] > 1 and len(v[var].shape) == 3 and v[var].shape[0]!=360 and v[var].shape[0]!=45 and v[var].shape[0]!=90): 
            if all_ys:
                if no_cos:
                    v[var] = reshape_all_ys(v[var], ind_z[0:v[var].shape[0]])
                else:
                    v[var] = reshape_cos_ys(v[var], ind_z[0:v[var].shape[0]], y)
            else:
                if ind_y is not None:
                    v[var] = reshape_one_y(v[var], ind_z[0:v[var].shape[0]], ind_y)
                else:
                    raise TypeError('Need to set an index value for ind_y')
        elif len(v[var].shape) == 2: 
            if all_ys:
                v[var] = v[var].swapaxes(0, 1)
                v[var] = v[var].reshape(-1,1) 
            else:
                if ind_y is not None:
                    v[var] = reshape_one_y(v[var], ind_z,ind_y)
                else:
                    raise TypeError('Need to set an index value for ind_y')
        elif (v[var].shape[0] > 1 and len(v[var].shape) == 4): 
            if all_ys:
                if no_cos:
                    v[var] = reshape_all_ys_4d(v[var], ind_z[0:v[var].shape[0]])
                else:
                    TypeError('Need to set an index value for ind_y 4 for with cosine option')
            else:
                if ind_y is not None:
                    v[var] = reshape_one_y(v[var], ind_z[0:v[var].shape[0]], ind_y)
                else:
                    raise TypeError('Need to set an index value for ind_y 4')

        elif ((v[var].shape[0] == 360 or v[var].shape[0] == 45 or v[var].shape[0] == 90) and len(v[var].shape) == 3):  
            if all_ys:
                if no_cos:
                    v[var] = reshape_all_ys_4d(v[var], ind_z[0:1])
                else:
                    TypeError('Need to set an index value for ind_y 4 for with cosine option')
            else:
                if ind_y is not None:
                    v[var] = reshape_one_y(v[var], ind_z[0:1], ind_y)
                else:
                    raise TypeError('Need to set an index value for ind_y 4')


        else:
            raise TypeError('There is a variable that has larger dimentions than 2 but it is not treated properly') 

    # Use relative humidity as a feature
    if use_rh: 
      # define generalized rh as qt/qsat based on coarse-grained qt and T
      p_Pa = p*100 # p is in hPa
      v['qin'] = v['qin']/atmos_physics.sam_qsat(v['Tin'],p_Pa[None,:])

     
    # Randomize the order of these events to ensure different y's are mixed up
    np.random.seed(123) # seed random number generator so always get same set of data when call this function (e.g., if call again to make plots)
    m = v[input_vert_vars[0]].shape[0]

    if rank_of_vars == 4: 
        print('We have 4D variables - xy structure')
    else:
        randind = np.random.permutation(m)
        for var in varis:
            v[var] = v[var][randind, :]

    # Concatenate feature and output variables together

    f = pack_list(v, input_vert_vars)
    o = pack_list(v, output_vert_vars)

    if rain_only:
       raise ValueError('rain_only not implemented')

    # Limit to only certain events if requested
    if n_trn_exs is not None:
        if n_trn_exs > o.shape[0]:
            warnings.warn('Requested more samples than available. Using the ' +
                          'maximum number available')
            n_trn_exs = o.shape[0]
        ind = np.arange(n_trn_exs)
        f = f[ind, :]
        o = o[ind, :]


    return f, o, y, z, rho, p , weight_list


def reshape_cos_ys(z, ind_z, y, is_sfc=False):
    if is_sfc:
        z = z.swapaxes(0, 1)
        z2 = np.empty((0))
    else:
        z = z[ind_z, :, :]
        z = z.swapaxes(0, 2)
        z2 = np.empty((0, sum(ind_z)))
    n_ex = z.shape[0]
    for i, yval in enumerate(y):
        # cosine of pseudo latitude
        Ninds = int(n_ex * np.cos((yval-np.mean(y))/6.37122e6))
        if is_sfc:
            z2 = np.concatenate((z2, z[0: Ninds, i]), axis=0)
        else:
            z2 = np.concatenate((z2, z[0:Ninds, i, :]), axis=0)
    return z2


def     reshape_all_ys(z, ind_z):
    # Expects data to be n_z n_y n_samples and returns
    # (n_y*n_samp n_z)
    z = z[ind_z, :, :]
    z = z.swapaxes(0, 2)
    return np.reshape(z, (-1, sum(ind_z)))

def reshape_all_ys_4d(z, ind_z):
    # Expects data to be n_z n_y n_samples and returns
    # (n_y*n_samp n_z)
    if len(z.shape) == 4:
        z = z[ind_z, :, :, :]
        z = z.swapaxes(0, 3)
    elif len(z.shape) == 3:
        z = z[:, :, :]
        z = np.transpose(z, axes=(2, 0, 1))
    else:
        TypeError('Cannot reshape this becuse dealing only with 3 and 4D arrays')

    return np.reshape(z, (-1, sum(ind_z)))


def reshape_one_y(z, ind_z, ind_y):
    # Expects data to be (n_z n_y n_samples) and returns (n_samp n_z)
    if len(z.shape) == 3 and ind_z.shape[0] > 1:
        z = z[ind_z, ind_y, :]
        z = z.swapaxes(0, 1)
    elif len(z.shape) == 3 and ind_z.shape[0] == 1:
        z = z[ind_y, :, :]
        # z = z.swapaxes(0, 1)
        z = np.reshape(z,(-1,1))
    elif len(z.shape) == 2:
        z = z[ind_y, :]
        z = np.reshape(z,(z.shape[0],1))
        # z = z.swapaxes(0, 1)
    elif len(z.shape) == 4:
        z = z[ind_z, ind_y,:,  :]
        z = np.reshape(z, (z.shape[0], -1))
        z = z.swapaxes(0, 1)
    else:
        raise TypeError('number of dimensions is unexpected')
    return z

def pack_f(T, q, axis=1):
    """Combines input profiles"""
    return np.concatenate((T, q), axis=axis)

def pack_list(v, vars_list ,axis=1):
    """gets a dictionary and makes it a large array"""
    inp_array = v[vars_list[0]] #initialize the array
    for var in vars_list[1:]:
        inp_array = np.concatenate((inp_array, v[var]), axis)
    return inp_array

def unpack_list(l_array, vars_list, vars_z_size ,axis=1):
    """Takes a large array, and give back a dictionary with the relevant fields"""
    v = dict()
    curr_dim = 0
    if sum(vars_z_size) >1:
        for name, dim in zip(vars_list, vars_z_size):
            v[name] = l_array[:, curr_dim:dim + curr_dim]
            curr_dim = curr_dim + dim
    else: #The case I only have one dimention....
        v[vars_list[0]] = l_array[:,None]
    return v


def pack_f_extended(T, q, u, v, w, axis=1):
    """Combines input profiles"""
    return np.concatenate((T, q, u, v, w), axis=axis)

def unpack_f(data, vari, axis=1):
    """Reverse pack operation"""
    N = int(data.shape[axis]/2)
    varipos = {'T': np.arange(N), 'q': np.arange(N,2*N)}
    out = np.take(data, varipos[vari], axis=axis)
    return out

def unpack_f_extended(data, vari, axis=1, wind_input=False):
    """Reverse pack operation"""
    if wind_input:
        Num_vars = int(data.shape[axis]/48)
        N = int(data.shape[axis]/Num_vars)
    else:
        N = int(data.shape[axis] / 2)

    varipos = {'T': np.arange(N), 'q': np.arange(N,2*N)}
    out = np.take(data, varipos[vari], axis=axis)
    return out

def pack_o(d1, d2, axis=1):
    """Combines T & q profiles"""
    return np.concatenate((d1, d2), axis=axis)


def choose_output_from_dic():
    """Gets an output from dictionary of outputs"""

def unpack_o(data, vari, axis=1):
    """Reverse pack operation"""
    N = int(data.shape[axis]/2)
    varipos = {'T': np.arange(N), 'q': np.arange(N, 2*N)}
    out = np.take(data, varipos[vari], axis=axis)
    return out

# Initialize & fit scaler Modified by Yani to fit for each generalized feature together
def init_pp_generalized(ppi, dict_data, input_vert_vars,scale_per_column):
    # Initialize list of scaler objects
    pp_dict = dict()
    for name in input_vert_vars:
        if ppi['name'] == 'MinMax':
            pp_dict[name] = preprocessing.MinMaxScaler(feature_range=(-1.0, 1.0))
            pp_dict[name].fit(np.reshape(dict_data[name],(-1,1)))
        elif ppi['name'] == 'MaxAbs':
            pp_dict[name] = preprocessing.MaxAbsScaler()
            pp_dict[name].fit(np.reshape(dict_data[name],(-1,1)))
        elif ppi['name'] == 'StandardScaler':
            pp_dict[name] = preprocessing.StandardScaler()
            if scale_per_column: #If yes it should scale every feature differently!
                pp_dict[name].fit(dict_data[name])
            else:
                pp_dict[name].fit(np.reshape(dict_data[name],(-1,1)))
        elif  ppi['name'] == 'F_stscl_add':
            pp_dict[name] = preprocessing.StandardScaler()
            if scale_per_column:  # Should scle each column seperately - to verify!
                pp_dict[name].fit(dict_data[name])
                std_add = 0.00001
                X_std = np.std(dict_data[name], axis=0, dtype=np.float64) + std_add
                pp_dict[name].mean_ = np.mean(dict_data[name], axis=0, dtype=np.float64)
                pp_dict[name].var_ = X_std*X_std
            else:
                raise TypeError('Choosing F_stscl_add was coded to assume we scale features for each column1')

        elif ppi['name'] == 'RobustScaler':
            pp_dict[name] = preprocessing.RobustScaler()
            pp_dict[name].fit(np.reshape(dict_data[name],(-1,1)))
        elif ppi['name'] == 'SimpleO':
            if len(input_vert_vars) !=2:
                print('Note that all variables but the first two are not normalized with 1!')
                # raise ValueError('Incorrect scaler name')
            pp_dict[name] = [atmos_physics.cp, atmos_physics.L]
            for i in range(len(input_vert_vars) - 2):
                pp_dict[name].append(1)
        elif ppi['name'] == 'SimpleO_expz':
            if len(input_vert_vars) !=2:
                # raise ValueError('Incorrect scaler name')
                print('Note that all variables but the first two are not normalized with 1!')
            pp_dict[name] = [atmos_physics.cp, atmos_physics.L]
            for i in range(len(input_vert_vars) - 2):
                pp_dict[name].append(1)
            else:
                pp_dict[name] = [atmos_physics.cp, atmos_physics.L]
        elif ppi['name'] == 'NoScaler':
            pp_dict[name] = []
        else:
            raise ValueError('Incorrect scaler name')

    return pp_dict



# Initialize & fit scaler
def init_pp(ppi, raw_data):
    # Initialize list of scaler objects
    if ppi['name'] == 'MinMax':
        pp = preprocessing.MinMaxScaler(feature_range=(-1.0, 1.0))
        pp.fit(raw_data)
    elif ppi['name'] == 'MaxAbs':
        pp = preprocessing.MaxAbsScaler() 
        pp.fit(raw_data)
    elif ppi['name'] == 'StandardScaler':
        pp = preprocessing.StandardScaler() 
        pp.fit(raw_data)
    elif ppi['name'] == 'RobustScaler':
        pp = preprocessing.RobustScaler()
        pp.fit(raw_data)
    elif ppi['name'] == 'SimpleO':
        pp = [atmos_physics.cp, atmos_physics.L]  
    elif ppi['name'] == 'SimpleO_expz':
        pp = [atmos_physics.cp, atmos_physics.L]  
    elif ppi['name'] == 'NoScaler':
        pp = []
    else:
        raise ValueError('Incorrect scaler name')

    return pp


# Transform data using initialized scaler
def transform_data_generalized(ppi, f_pp_dict, f_dict, input_vert_vars, z,scale_per_column=False,rewight_outputs=False,weight_list=[1,1]):
    if ppi['name'] == 'SimpleO':
        trans_data_dic = dict()
        for (index, name) in enumerate(input_vert_vars):
            trans_data_dic[name]= f_dict[name]*f_pp_dict[name][index]
    elif ppi['name'] == 'SimpleO_expz':
        trans_data_dic = dict()
        for (index, name) in enumerate(input_vert_vars):
            trans_data_dic[name]= f_dict[name]*f_pp_dict[name][index]*np.exp(-z/7000.0)
    elif ppi['name'] == 'NoScaler':
        trans_data_dic = f_dict
    elif ppi['name'] == 'F_stscl_add': 
        trans_data_dic = dict()
        for name in input_vert_vars:
            if scale_per_column:  
                trans_data_dic[name] = (f_dict[name] - f_pp_dict[name].mean_)/np.sqrt(f_pp_dict[name].var_)
            else:
                raise TypeError('Choosing F_stscl_add was coded to assume we scale features for each column')

    else: #Using standard scalar to renormalize
        trans_data_dic = dict()
        for name in input_vert_vars:
            if scale_per_column: 
                trans_data_dic[name] = f_pp_dict[name].transform(f_dict[name])
            else: 
                trans_data_dic[name] = np.reshape(f_pp_dict[name].transform(np.reshape(f_dict[name],(-1,1))),(f_dict[name].shape[0],f_dict[name].shape[1]))

    if rewight_outputs: 
        print('rescaling outputs')
        print('length of the weight list is:', len(weight_list))
        for ind, name in enumerate(input_vert_vars,start=0):
            trans_data_dic[name] = trans_data_dic[name]*weight_list[ind]

    return trans_data_dic


def inverse_transform_data_generalized(ppi, f_pp_dict, f_dict, input_vert_vars,
                                       z,scale_per_column=False,rewight_outputs=False,weight_list=[1,1]):

    if rewight_outputs: 
        for ind, name in enumerate(input_vert_vars,start=0):
            f_dict[name] = f_dict[name]/weight_list[ind]

    if ppi['name'] == 'SimpleO':
        trans_data_dic = dict()
        for (index, name) in enumerate(input_vert_vars):
            trans_data_dic[name]= f_dict[name]/f_pp_dict[name][index]
    elif ppi['name'] == 'SimpleO_expz':
        trans_data_dic = dict()
        for (index, name) in enumerate(input_vert_vars):
            trans_data_dic[name]= f_dict[name]/f_pp_dict[name][index]/np.exp(-z/7000.0)
    elif ppi['name'] == 'NoScaler':
        trans_data_dic = f_dict
    else:
        trans_data_dic = dict()
        for name in input_vert_vars:
            if scale_per_column: 
                trans_data_dic[name] = f_pp_dict[name].inverse_transform(f_dict[name])
            else: 
                trans_data_dic[name] = np.reshape(f_pp_dict[name].inverse_transform(np.reshape(f_dict[name],(-1,1))),(f_dict[name].shape[0],f_dict[name].shape[1]))
    return_data = pack_list(trans_data_dic,input_vert_vars)
    # Return a numpy array of the transformed data output
    return return_data





# Transform data using initialized scaler
def transform_data(ppi, pp, raw_data, z):
    if ppi['name'] == 'SimpleO':
        T_data = unpack_o(raw_data, 'T')*pp[0]
        q_data = unpack_o(raw_data, 'q')*pp[1]
        return_data = pack_o(T_data, q_data)
    elif ppi['name'] == 'SimpleO_expz':
        T_data = unpack_o(raw_data, 'T')*pp[0]*np.exp(-z/7000.0)
        q_data = unpack_o(raw_data, 'q')*pp[1]*np.exp(-z/7000.0)
        return_data = pack_o(T_data, q_data)
    elif ppi['name'] == 'NoScaler':
        return_data = raw_data
    else:
        return_data = pp.transform(raw_data)

    # Return single transformed array as output
    return return_data 


# Apply inverse transformation to unscale data
def inverse_transform_data(ppi, pp, trans_data, z):
    if ppi['name'] == 'SimpleO':
        T_data = unpack_o(trans_data, 'T')/pp[0]
        q_data = unpack_o(trans_data, 'q')/pp[1]
        return_data = pack_o(T_data, q_data)
    elif ppi['name'] == 'SimpleO_expz':
        T_data = unpack_o(trans_data, 'T')/pp[0]*np.exp(z/7000.0)
        q_data = unpack_o(trans_data, 'q')/pp[1]*np.exp(z/7000.0)
        return_data = pack_o(T_data, q_data)
    elif ppi['name'] == 'NoScaler':
        return_data = trans_data
    else:
        return_data = pp.inverse_transform(trans_data)
    return return_data



def load_one_y(f_ppi, o_ppi, f_pp, o_pp, est, ind_y, datafile, max_z, input_vert_vars, output_vert_vars, input_vert_dim, output_vert_dim,
                 n_trn_exs, rain_only, no_cos, use_rh, wind_input = False,scale_per_column=False,
                 rewight_outputs=False,weight_list=[1,1],do_nn=False):
    """Returns n_samples 2*n_z array of true and predicted values
       at a given y"""
    # Load data
    f, o, y, z, rho, p, weight_list = \
        LoadData(datafile, max_z, input_vert_vars, output_vert_vars, all_ys=False, ind_y=ind_y,
                 verbose=False, n_trn_exs=None, rain_only=rain_only, 
                 no_cos=no_cos, use_rh=use_rh, wind_input = wind_input, rewight_outputs =rewight_outputs )
    # Calculate predicted output

    f_dict = unpack_list(f, input_vert_vars,input_vert_dim)
    f_scl_dict = transform_data_generalized(f_ppi, f_pp, f_dict, input_vert_vars, z,scale_per_column, rewight_outputs=False)
    # f_scl = transform_data(f_ppi, f_pp, f, z)
    f_scl = pack_list(f_scl_dict, input_vert_vars)

    if do_nn:
        tmp_f_scl = torch.from_numpy(f_scl)
        est.eval()
        o_pred_scl = est(tmp_f_scl.float()) 
        o_pred_scl = o_pred_scl.detach().numpy()
    else:
        o_pred_scl = est.predict(f_scl)
    o_pred_scl_dict = unpack_list(o_pred_scl, output_vert_vars,output_vert_dim)
    o_pred = inverse_transform_data_generalized(o_ppi, o_pp, o_pred_scl_dict,output_vert_vars, z,scale_per_column
                                                ,rewight_outputs=rewight_outputs,weight_list=weight_list)
    o_pred_dict = unpack_list(o_pred, output_vert_vars,output_vert_dim)

    o_dict = unpack_list(o, output_vert_vars,output_vert_dim)


    return o_dict, o_pred_dict

def stats_by_yz(f_ppi, o_ppi, f_pp, o_pp, est, y, z, rho, datafile, n_trn_exs, input_vert_vars, output_vert_vars,
                input_vert_dim, output_vert_dim, rain_only, no_cos, use_rh, wind_input = False,scale_per_column=False,
                rewight_outputs=False,weight_list=[1,1],do_nn=False):
    # Initialize
    output_stat_dict = dict()
    feature_list = ['_mean','_var','_bias','_rmse','_r','_Rsq']
    for output_name,z_dim in zip(output_vert_vars,output_vert_dim):
        for feature in feature_list:
            output_stat_dict[output_name+feature] = np.zeros((len(y), z_dim))

    output_stat_dict['Pmean_true'] = np.zeros((len(y)))
    output_stat_dict['Pmean_pred']= np.zeros((len(y)))
    output_stat_dict['Pextreme_true']= np.zeros((len(y)))
    output_stat_dict['Pextreme_pred']= np.zeros((len(y)))
    #
    for i in range(len(y)):
        o_true_dict, o_pred_dict = \
            load_one_y(f_ppi, o_ppi, f_pp, o_pp, est, i, datafile,
                         np.max(z), input_vert_vars, output_vert_vars, input_vert_dim, output_vert_dim, n_trn_exs, rain_only,
                         no_cos, use_rh, wind_input = wind_input,scale_per_column = scale_per_column,
                         rewight_outputs=rewight_outputs,weight_list=weight_list,do_nn=do_nn)

        if i==0:
         print('size of test dataset for a given y and level', o_true_dict[output_vert_vars[0]].shape[0])

        for output_name,z_dim in zip(output_vert_vars,output_vert_dim):
            output_stat_dict[output_name+'_mean'][i,:] = np.mean(o_true_dict[output_name],axis=0)
            output_stat_dict[output_name+'_var'][i,:] = np.var(o_true_dict[output_name],axis=0)
            output_stat_dict[output_name+'_bias'][i,:] = np.mean(o_pred_dict[output_name],axis=0) - output_stat_dict[output_name+'_mean'][i,:]
            output_stat_dict[output_name+'_rmse'][i,:] = np.sqrt(
                metrics.mean_squared_error(
                o_true_dict[output_name], o_pred_dict[output_name],
                                           multioutput='raw_values'))
            for j in range(z_dim):
                if np.sum(o_true_dict[output_name][:, j]==0) >  o_true_dict[output_name][:, j].shape[0]*0.99 and output_name!='qpout': 
                    output_stat_dict[output_name + '_Rsq'][i, j] = np.nan
                    continue
                output_stat_dict[output_name +'_r'][i,j] = scipy.stats.pearsonr(
                    o_true_dict[output_name][:, j], o_pred_dict[output_name][:, j])[0]
                output_stat_dict[output_name + '_Rsq'][i,j] = metrics.r2_score(o_true_dict[output_name][:, j], o_pred_dict[output_name][:, j])
                if output_stat_dict[output_name + '_Rsq'][i, j] < -10:
                    output_stat_dict[output_name + '_Rsq'][i, j] = -10
            if output_name == 'qout':
                P_true = atmos_physics.calc_precip(o_true_dict['qout'], rho, z,output_vert_vars,o_true_dict)
                P_pred = atmos_physics.calc_precip(o_pred_dict['qout'], rho, z,output_vert_vars,o_pred_dict)
                output_stat_dict['Pmean_true'][i] = np.mean(P_true)
                output_stat_dict['Pmean_pred'][i] = np.mean(P_pred)
                output_stat_dict['Pextreme_true'][i] = np.percentile(P_true, 99.9)
                output_stat_dict['Pextreme_pred'][i] = np.percentile(P_pred, 99.9)




    return output_stat_dict

def GetDataPath(training_expt, wind_input = False,is_cheyenne=False,full_data_separate=False):

    if is_cheyenne == False:  # On aimsir/esker
        base_dir = '/net/aimsir/archive1/janniy/'
    else:
        base_dir = '/glade/scratch/janniy/'

    if wind_input:
        datadir = base_dir + 'mldata/training_data_tmp/'
    else:
        datadir = base_dir + 'mldata_tmp/training_data/'


    # practice_flag = False
    if full_data_separate:
        trainfile = datadir + training_expt + '_training_short.pkl'
        testfile = datadir + training_expt + '_testing_short.pkl'
    else:
        trainfile = datadir + training_expt + '_training.pkl'
        testfile = datadir + training_expt + '_testing.pkl'

    pp_str = training_expt + '_'

    print(trainfile)
    print(testfile)
    return datadir, trainfile, testfile, pp_str



def get_f_o_pred_true(est_str, training_file, max_z, input_vert_vars, output_vert_vars,input_vert_dim,output_vert_dim,
                      all_ys=True, ind_y=None, 
                      n_trn_exs=None, rain_only=False,  
                      no_cos=False, use_rh=False, wind_input = False, scale_per_column=False,
                      rewight_outputs=False,weight_list=[1,1],is_cheyenne=False, do_nn =False):
    # Load model and preprocessors

    if is_cheyenne == False:  # On aimsir/esker
        base_dir = '/net/aimsir/archive1/janniy/'
    else:
        base_dir = '/glade/scratch/janniy/'

    est, _, errors, f_ppi, o_ppi, f_pp, o_pp, y, z, _, _ = \
        pickle.load(open(base_dir + 'mldata_tmp/regressors/' + est_str + '.pkl', 'rb'))
    # Load raw data from file
    f, otrue, _, _, _, _, weight_list = \
        LoadData(training_file, max_z=max_z, input_vert_vars=input_vert_vars, output_vert_vars=output_vert_vars, all_ys=all_ys, ind_y=ind_y, n_trn_exs=n_trn_exs, rain_only=rain_only, no_cos=no_cos, use_rh=use_rh, wind_input = wind_input, rewight_outputs =rewight_outputs )
    print('JY - added weight list to plot - need to think if necessary')
    # Scale true values
    otrue_dict = unpack_list(otrue,output_vert_vars,output_vert_dim)
    otrue_scl_dict = transform_data_generalized(o_ppi, o_pp, otrue_dict,output_vert_vars, z,scale_per_column, rewight_outputs=rewight_outputs,weight_list=weight_list)
    otrue_scl= pack_list(otrue_scl_dict, output_vert_vars)


    # Apply f preprocessing to scale f-data and predict output
    #
    f_dict = unpack_list(f, input_vert_vars, input_vert_dim)
    f_scl_dict = transform_data_generalized(f_ppi, f_pp, f_dict, input_vert_vars, z,scale_per_column,rewight_outputs=False)
    f_scl = pack_list(f_scl_dict, input_vert_vars)

    if do_nn:
        tmp_f_scl = torch.from_numpy(f_scl)
        est.eval()
        opred_scl = est(tmp_f_scl.float()) 
        opred_scl=opred_scl.detach().numpy()
    else:
        opred_scl = est.predict(f_scl) 
    opred_scl_dict = unpack_list(opred_scl, output_vert_vars, output_vert_dim)
    opred = inverse_transform_data_generalized(o_ppi, o_pp, opred_scl_dict,output_vert_vars, z, scale_per_column,
                                               rewight_outputs=rewight_outputs,weight_list=weight_list)


    return f_scl, opred_scl, otrue_scl, f, opred, otrue


def load_error_history(est_str,is_cheyenne=False):
    if is_cheyenne == False:  # On aimsir/esker
        base_dir = '/net/aimsir/archive1/janniy/'
    else:
        base_dir = '/glade/scratch/janniy/'
    _, _, err, _, _, _, _, _, _, _ = pickle.load(open(base_dir + 'mldata_tmp/regressors/' +
                                                      est_str, + 'pkl', 'rb'))
    return err



def GetDataPath_nn(training_expt, wind_input = False,is_cheyenne=False,full_data_separate=False):
    if is_cheyenne == False:  # On aimsir/esker
        base_dir = '/net/aimsir/archive1/janniy/'
    else:
        base_dir = '/glade/scratch/janniy/'

    if wind_input:
        datadir = base_dir + 'mldata/training_data_tmp/'
    else:
        datadir = base_dir + 'mldata_tmp/training_data/'

    if full_data_separate:
        trainfile = datadir + training_expt + '_training_short.pkl'
        testfile = datadir + training_expt + '_testing_short.pkl'
    else:
        trainfile = datadir + training_expt + '_training.pkl'  
        testfile = datadir + training_expt + '_testing.pkl' 

    pp_str = training_expt + '_'

    print(trainfile)
    print(testfile)
    return datadir, trainfile, testfile, pp_str
