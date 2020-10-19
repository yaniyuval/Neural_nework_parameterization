import numpy as np
import scipy.stats
import matplotlib
matplotlib.use('Agg')  # so figs just print to file. 
import matplotlib.pyplot as plt
import src.ml_train as ml_train
import src.ml_load as ml_load
import pickle
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score


# n_cv-fold cross validation
n_cv = 10
n_jobs = 4

# note cv used elsewhere to mean 'convection active'

def plot_error_vs_n_trees(f_ppi, o_ppi, max_z=20000.0,
                     training_expt = 'abs1.0_norf_ras',
                     min_samples_leaf = 10,  
                     n_trn_exs=None, rain_only=False, 
                     no_cos = True,
                     use_rh=False,
                     load_results=True,
                     save_results=False):

    if load_results:
     print('loading results')
     n_trees, cv_error, cv_error_std =pickle.load(open('figs_errors/error_vs_n_trees.pkl', 'rb'))
    else:
     datadir, trainfile, _, _ = ml_load.GetDataPath(training_expt)

     f, o, y, z, rho, p  = ml_load.LoadData(trainfile, max_z, n_trn_exs=n_trn_exs, rain_only=rain_only, no_cos=no_cos, use_rh=use_rh)


     # scale data
     f_pp = ml_load.init_pp(f_ppi, f)
     o_pp = ml_load.init_pp(o_ppi, o)
     f_scl = ml_load.transform_data(f_ppi, f_pp, f, z)
     o_scl = ml_load.transform_data(o_ppi, o_pp, o, z)

     rf = RandomForestRegressor(min_samples_leaf = min_samples_leaf, max_features = 1.0/3.0, random_state = 123, warm_start = False)

     min_n_trees = 1
     max_n_trees = 21

     n_trees = range(min_n_trees, max_n_trees + 1,2)
     cv_error = np.zeros(len(n_trees))
     cv_error_std = np.zeros(len(n_trees)) # standard deviation across folds

     for i in range(len(n_trees)):
        rf.set_params(n_estimators=n_trees[i])
        scores = cross_val_score(rf, f_scl, o_scl, cv=n_cv, n_jobs=n_jobs)
        cv_error[i] = 1-scores.mean()
        cv_error_std[i] = scores.std()
        print(str(n_trees[i]) + ': ' + str(cv_error[i]))

     if save_results:
      print('saving results')
      pickle.dump([n_trees, cv_error, cv_error_std], open('figs_errors/error_vs_n_trees.pkl', 'wb'))

    print(list(n_trees))
    print(cv_error)
    # some confusion in literature as to whether should include 1/sqrt(n_cv) in standard error
    print(np.max(cv_error_std/np.sqrt(n_cv))) # max of standard error
    fig = plt.figure(figsize=(3.0,2.25))
    plt.plot(n_trees, cv_error, 'o-')
    plt.xlim(0, 22.8)
    plt.ylim(0,0.51) 
    plt.xlabel("Number of trees")
    plt.ylabel("Error")
    #plt.legend(loc="upper right")
    #plt.legend(frameon=False)
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    plt.tight_layout() # avoid overlap
    plt.show()
    fig.savefig('figs_errors/error_vs_ntrees.eps', bbox_inches='tight')
    plt.close()


##############

def plot_error_vs_min_samples_leaf(f_ppi, o_ppi, max_z=20000.0,
                                   training_expt = 'abs1.0_norf_spinup',
                                   n_trees=10, 
                                   n_trn_exs=None, rain_only=False, 
                                   no_cos = True, use_rh=False,
                                   load_results=True,
                                   save_results=False):

    if load_results:
     print('loading results')
     min_samples_leaf, cv_error, cv_error_std =pickle.load(open('figs_errors/error_vs_min_samples_leaf.pkl', 'rb'))
    else:

     datadir, trainfile, _, _ = ml_load.GetDataPath(training_expt)

     f, o, _, z, rho, _  = ml_load.LoadData(trainfile, max_z, n_trn_exs=n_trn_exs, rain_only=rain_only, no_cos=no_cos, use_rh=use_rh)

     # scale data
     f_pp = ml_load.init_pp(f_ppi, f)
     o_pp = ml_load.init_pp(o_ppi, o)
     f_scl = ml_load.transform_data(f_ppi, f_pp, f, z)
     o_scl = ml_load.transform_data(o_ppi, o_pp, o, z)

     rf = RandomForestRegressor(n_estimators = n_trees, random_state = 123, max_features = 1.0/3.0, warm_start = False)

     min_min_samples_leaf = 1
     max_min_samples_leaf = 16
     step_min_samples_leaf = 3

     min_samples_leaf = range(min_min_samples_leaf, max_min_samples_leaf + 1, step_min_samples_leaf)

     cv_error = np.zeros(len(min_samples_leaf))
     cv_error_std = np.zeros(len(min_samples_leaf))

     for i in range(len(min_samples_leaf)):
        print(min_samples_leaf[i])
        rf.set_params(min_samples_leaf=min_samples_leaf[i])
        scores = cross_val_score(rf, f_scl, o_scl, cv=n_cv, n_jobs=n_jobs)
        cv_error[i] = 1-scores.mean()
        cv_error_std[i] = scores.std()
        print(str(min_samples_leaf[i]) + ': ' + str(cv_error[i]))

     if save_results:
      print('saving results')
      pickle.dump([min_samples_leaf, cv_error, cv_error_std], open('figs_errors/error_vs_min_samples_leaf.pkl', 'wb'))

    print(list(min_samples_leaf))
    print(cv_error)
    print(np.max(cv_error_std/np.sqrt(n_cv)))
    fig = plt.figure(figsize=(3.0,2.25))
    plt.plot(min_samples_leaf, cv_error, '-o')
    plt.xlim(0, 16.5)
    plt.ylim(0,0.51) 
    plt.xlabel("Minimum sample size for a leaf")
    plt.ylabel("Error")
    #plt.legend(loc="upper right")
    #plt.legend(frameon=False)
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    plt.tight_layout() # avoid overlap
    plt.show()
    fig.savefig('figs_errors/error_vs_min_samples_leaf.eps', bbox_inches='tight')
    plt.close()

##############

def plot_error_vs_n_trn_exs(f_ppi, o_ppi, max_z=40000.0,
                            training_expt = 'abs1.0_norf_spinup',
                            n_trees=10, min_samples_leaf=10,
                            rain_only=False, 
                            no_cos = True, 
                            use_rh=False,
                            load_results=True,
                            save_results=False):
   

    if load_results:
     print('loading results')
     n_trn_exs, cv_error, cv_error_std =pickle.load(open('figs_errors/error_vs_n_trn_exs.pkl', 'rb'))
    else:
     datadir, trainfile, _, _ = ml_load.GetDataPath(training_expt)

     rf = RandomForestRegressor(n_estimators = n_trees, min_samples_leaf = min_samples_leaf, max_features = 1.0/3.0, random_state = 123, warm_start = False)

     n_trn_exs = np.array([1, 5, 10, 30, 50, 70, 80, 90])*10000

     cv_error = np.zeros(len(n_trn_exs)) 
     cv_error_std = np.zeros(len(n_trn_exs)) # standard deviation of error estimate across folds

     for i in range(len(n_trn_exs)):
        f, o, _, z, rho, _  = ml_load.LoadData(trainfile, max_z, n_trn_exs=n_trn_exs[i], rain_only=rain_only, no_cos=no_cos, use_rh=use_rh)
        f_pp = ml_load.init_pp(f_ppi, f)
        o_pp = ml_load.init_pp(o_ppi, o)
        f_scl = ml_load.transform_data(f_ppi, f_pp, f, z)
        o_scl = ml_load.transform_data(o_ppi, o_pp, o, z)
        scores = cross_val_score(rf, f_scl, o_scl, cv=n_cv, n_jobs=n_jobs)
        cv_error[i] = 1-scores.mean()
        cv_error_std[i] = scores.std()
        print(str(n_trn_exs[i]) + ': ' + str(cv_error[i]))

     if save_results:
      print('saving results')
      pickle.dump([n_trn_exs, cv_error, cv_error_std], open('figs_errors/error_vs_n_trn_exs.pkl', 'wb'))


    print(n_trn_exs)
    print(cv_error)
    print(np.max(cv_error_std/np.sqrt(n_cv)))

    fig = plt.figure(figsize=(3.0,2.25))
    fscale = 100000.0
    plt.plot(n_trn_exs/fscale, cv_error, '-o')
    plt.xlim(-0.15, 9.5)
    plt.ylim(0,0.51) 
    plt.xlabel("Number of training examples ($10^5$)")
    plt.ylabel("Error")
    #plt.legend(loc="upper right")
    #plt.legend(frameon=False)
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    plt.tight_layout() # avoid overlap
    plt.show()
    fig.savefig('figs_errors/error_vs_n_trn_exs.eps', bbox_inches='tight')
    plt.close()

def plot_train_test_error_nocv_vs_n_trn_exs(f_ppi, o_ppi, max_z=40000.0,
                            training_expt = 'abs1.0_norf_spinup',
                            n_trees=10, min_samples_leaf=10,
                            rain_only=False,
                            no_cos=True,
                            use_rh=False):

    datadir, trainfile, testfile, _ = ml_load.GetDataPath(training_expt)

    f_test, o_test, _, z, _, _ = ml_load.LoadData(testfile, max_z, n_trn_exs=None, rain_only=rain_only, no_cos=no_cos, use_rh=use_rh)

    rf = RandomForestRegressor(n_estimators = n_trees, min_samples_leaf = min_samples_leaf, random_state = 123, warm_start = False)

    n_trn_exs = np.array([1, 5, 10, 30, 50, 70, 80, 90])*10000

    test_error = np.zeros(len(n_trn_exs))
    train_error = np.zeros(len(n_trn_exs))

    for i in range(len(n_trn_exs)):
       f, o, _, z, rho, _  = ml_load.LoadData(trainfile, max_z, n_trn_exs=n_trn_exs[i], rain_only=rain_only, no_cos=no_cos, use_rh=use_rh)
       f_pp = ml_load.init_pp(f_ppi, f)
       o_pp = ml_load.init_pp(o_ppi, o)
       f_scl = ml_load.transform_data(f_ppi, f_pp, f, z)
       o_scl = ml_load.transform_data(o_ppi, o_pp, o, z)
       rf.fit(f_scl, o_scl)
       f_test_scl = ml_load.transform_data(f_ppi, f_pp, f_test, z)
       o_test_scl = ml_load.transform_data(o_ppi, o_pp, o_test, z)
       test_error[i] = 1.0-rf.score(f_test_scl,o_test_scl)
       train_error[i] = 1.0-rf.score(f_scl,o_scl)
       print(str(n_trn_exs[i]) + ': ' + str(test_error[i]))
       print(str(n_trn_exs[i]) + ': ' + str(train_error[i]))

    print(test_error)
    print(train_error)
    fig = plt.figure()
    fscale = 100000.0
    plt.plot(n_trn_exs/fscale, test_error, '-o', label='test')
    plt.plot(n_trn_exs/fscale, train_error, '-o', label='train')
    plt.xlim(-0.15, 9.5)
    plt.ylim(0,0.51) 
    plt.xlabel("n_trn_exs")
    plt.ylabel("Error")
    plt.legend(loc="upper right")
    plt.legend(frameon=False)
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    plt.tight_layout() # avoid overlap
    plt.show()
    fig.savefig('figs_errors/error_test_train_nocv_vs_n_trn_exs.eps', bbox_inches='tight')
    plt.close()

