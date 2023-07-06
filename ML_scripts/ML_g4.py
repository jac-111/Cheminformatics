# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 20:56:19 2022

This code is forms the bulk of the ML. 

Scikit-learn estimators are used to predict individual QM properties 
given in the csv file.

This code has an associated input script and paramter library and WILL NOT
function without them. The input file describes the estimators to use, hyperparameter
choice, data split etc and must be supplied at execution. the /lib/ directory contains
the hyperparameters for each esitmator in the grid or random regime.

Histogram sampling is employed to generate reperesentative sampling.

Output from the program are put into the /CrossValidation/ directory.

@author: jacob
"""

import sys, os
import pandas as pd
import numpy as np
import random
import csv

from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RandomizedSearchCV
from sklearn.neighbors import KNeighborsRegressor
from sklearn.cross_decomposition import PLSRegression
from sklearn.linear_model import SGDRegressor
from sklearn.kernel_ridge import KernelRidge
from sklearn.neural_network import MLPRegressor
from sklearn.svm import SVR


from scipy.stats import expon
from scipy.stats import uniform
from scipy.stats import loguniform

#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                   ~       SETUP DATA - std.inp       ~
#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if len(sys.argv) > 2:                                                                       # Check for too many sys.args
    print("\n**** ERROR ****\n\nToo many arguements passed to script.\n")
    print("Only pass the .txt file with the ML input keys and values\n")
    print("See /Examples directory for correct formatting.\n\n")
    quit()
elif len(sys.argv) == 1:                                                                    # Check for too few sys.args
    print("\n**** ERROR ****\n\nPass .txt with ML input keys and values in command line")
    print("e.g.\npython ML_g1.py ML_input.txt")
    print("See /Examples directory for correct formatting.\n\n")
    quit()
elif len(sys.argv) == 2:                                                                    # Good execution
    print("working with %s" % sys.argv[1])


inp_prams = {}                                                                              # Dictionary to store ML input parameters

read_ml_inp = open(sys.argv[1], 'r')                                                        # Read ML_input.txt and assign values
for i, line in enumerate(read_ml_inp):                              
    if line.split():
        if line.split()[0] != '#':            
            inp_prams[line.split()[0]] = float(line.split()[1])
            

default_keys = ['TMQM_portion', 'training_split', 'num_racs', 'num_data_bins',\
                'CV_folds', 'CV_grid', 'CV_random', 'CV_random_iter', 'KNeighboursRegression',\
                'PLSRegression','StochasticGradient', 'KRR', 'SVR', 'MLPRegression', 'CrossValidation',\
                ]                                                                           # Default 'keys' for the input file/dict

default_vals = [0.5, 0.2, 151, 200, 10, 5, 1, 0, 0, 0 ,0, 0, 0, 0, 0]                                      # Default 'values' for the input file/dict

estimators = ['KNeighboursRegression', 'PLSRegression', 'StochasticGradient', 'KRR', 'SVR',\
              'MLPRegression']                                                               # Names of estimators

estimator_types = {"KNeighboursRegression": KNeighborsRegressor(), "PLSRegression" : PLSRegression(), \
                   "StochasticGradient" : SGDRegressor(), "KRR" : KernelRidge(), "SVR" : SVR(), \
                    "MLPRegression" : MLPRegressor()}                                       # Dict of estimator funcs
    
bad_inp = 0
if len(inp_prams) != len(default_keys):                                 # Check number of keys given in file
    print("###\tERROR\t###\n\nMISSING OR TOO MANY INPUT VALUES:")
    bad_inp = 1
for key in inp_prams:                                                   # Check key names
    if key not in default_keys:
        print("UNKNOWN INPUT : %s" % key)
        bad_inp = 1
        
if bad_inp == 1:                                                        # Quit if bad input keys/values
    print("\n\n###\tERROR\t###\n\nPLEASE REASSESS INPUT FILE")
    print("SHOULD HAVE THESE VALUES:")
    for key in default_keys:
        print(key)
    quit()

for key in inp_prams:                                                   # Quick read of input params
    print(key, inp_prams[key])

    
num_bins = int(inp_prams['num_data_bins'])                              # num bins for training data histogram 

if not os.path.exists('./CrossValidation'):
    os.makedirs('./CrossValidation')

#path_to_data = '../Processed_data/tmqm_RACs_part10.csv'
#path_to_data = '../Processed_data/tmqm_RACs-151.csv'                        # sample path to data
path_to_data = '../Processed_data/tmqm_n_DB-Heterocyclic_RAC-151_EDis.csv'                    # path to data

tmqm_data = pd.read_csv(path_to_data)                                   # data as pandas DF

print("\nUsing %s for working data\n\nis this correct?" % path_to_data) # check .csv file
check = input('yes or no?\n')

if check not in ("yes", "y"):
    
    print("\nNOTE:\nthe default tmqm_RAC data is found in the Processed_data directory\n")
    print("please edit the path_to_data variable in this program to the .csv file of choice")
    quit()


total_num_data = (len(tmqm_data))                                       # check data type and number of columns

if inp_prams['num_racs']:
    num_RACs = int(inp_prams['num_racs'])
else:
    print("PLEASE CHECK ML INPUT FILE:\nit appears that no value for the number of racs is given.")
    quit()
    

num_columns = len(tmqm_data.columns)                                    # total number of columns inc qm

tmqm_titles = []                                                        # all column titles

for col in tmqm_data.columns:
    tmqm_titles.append(col)

QM_names = (tmqm_titles[2 : len(tmqm_titles) - (num_RACs)])             # QM column titles only
num_QM = len(QM_names)
print("These are the QM properties to be considered:")
print(QM_names)



if inp_prams['CV_grid'] == 1 and inp_prams['CV_random'] == 0:
    CV_type = './lib/grid/'
elif inp_prams['CV_grid'] == 0 and inp_prams['CV_random'] == 1:
    CV_type = './lib/random/'
else:
    print("\n\nERROR\n\Check CV Settings in ML_input file\n\n")
    quit()

num_training_vals = int(len(tmqm_data.iloc[:,2]) * inp_prams['training_split']) 


#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                   ~     MachineLearning function      ~
#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


def run_ml(estimator, QM_id, parameters):
    
    print(estimator, QM_id, parameters)
    

                                                                                                    # HISTOGRAM DATA TO GENERATE REPRESENTATIVE SAMPLE
    hist1 = pd.cut(np.array(tmqm_data.iloc[:,(QM_id + 2)]), bins=num_bins, labels=False)                  # makes list of histogram 'positions'
                                                                                                        # i.e. list of bin position for each element
    hist_labels = pd.cut(tmqm_data.iloc[:,(QM_id + 2)], bins=num_bins, labels=False, retbins=True)        # Bin intervals
    
    arr_binned_vals = [[] for x in range(num_bins)]                                             # list for list containing values in each bin
    arr_binned_portion = [[] for x in range(num_bins)]                                          # list for list containing values in each bin
    arr_binned_training = [[] for x in range(num_bins)]                                         # list for list containing values in each bin
                                                                                                # with training partition from random sampling 
    
    for x in range(num_bins):                                                                   # collect values for list of list
        for y in range(len(hist1)):                                                             # Each nested list is a bin. Fill it with index of element in df
            if hist1[y] == x:
                arr_binned_vals[x].append(y)                                            
    
    my_sum = 0      
    
    # for each bin, randomly sample from whole dataset at given percentage. = arr_binned_portion
    # ^ this is for smaller dataset analysis.
    # from these nested lists, randomly sample again to produce training set.
    for x in range(len(arr_binned_vals)):
        my_sum += int(len(arr_binned_vals[x]) * inp_prams['training_split'])                    # count number of values being used
        arr_binned_portion[x] = sorted(random.sample(arr_binned_vals[x], (int(len(arr_binned_vals[x]) * inp_prams['TMQM_portion'])))) # randomly select from binned data
        arr_binned_training[x] = sorted(random.sample(arr_binned_portion[x], (int(len(arr_binned_portion[x]) * inp_prams['training_split'])))) # randomly select from binned data
    
    
    
    portion_ids = [j for i in arr_binned_portion for j in i]            # get index of all utilsed values
    training_ids = [j for i in arr_binned_training for j in i]          # get training indices
    test_ids = [x for x in portion_ids if x not in training_ids]        # get test indices
    
    # ASSIGNMENT OF DATAFRAMES FOR ML
    
    train_vals = tmqm_data.iloc[:,QM_id + 2][training_ids]                                  # Training QM data - individual basis
    
    train_racs = tmqm_data.iloc[:, (num_columns-num_RACs):num_columns].iloc[training_ids]   # Training RACs
    
    test_vals = tmqm_data.iloc[:,QM_id + 2].iloc[test_ids]                                  # Test QM data
    
    test_racs = tmqm_data.iloc[:, (num_columns-num_RACs):num_columns].iloc[test_ids]        # Test RACs
    
    train_concat = pd.concat([train_vals,train_racs],axis=1)                                # join training data into one df
    test_concat = pd.concat([test_vals,test_racs], axis=1)                                  # join test data into one df
    
    train_center = pd.DataFrame.mean(train_concat)                                          # training mean
    train_scale =  pd.DataFrame.std(train_concat)                                           # training standard deviation
    
    
    
    train_scaled = (train_concat - train_center) / train_scale                              # normalisation
    test_scaled = (test_concat - train_center) / train_scale
    
    if inp_prams['CV_grid'] == 1 and inp_prams['CV_random'] == 0:                           # Check for random or grid search
        
        # GRID CROSSVALIDATION
        cv = GridSearchCV(
            estimator_types[estimator],
            parameters,
            cv = int(inp_prams['CV_folds']),
            verbose = 1)
        
    elif inp_prams['CV_grid'] == 0 and inp_prams['CV_random'] == 1:
        
        # RANDOM CROSSVALIDATION
        cv = RandomizedSearchCV(
            estimator_types[estimator],
            parameters,
            cv = int(inp_prams['CV_folds']),
            n_iter = int(inp_prams['CV_random_iter']),
            verbose = 1)

    
    cv.fit(train_scaled.iloc[:,1:], train_scaled.iloc[:,0])                 # Execute crossvalidation
    
    results_dict = cv.cv_results_                                           # All results from crossvalidation
  
    best_params = cv.best_params_                                           # list of best parameters
  
    train_score = cv.best_score_                                            # training score
  
    test_score = cv.score(test_scaled.iloc[:,1:], test_scaled.iloc[:,0])
    
    
    return(results_dict, best_params, train_score, test_score)
    

#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                   ~     Function to parse random
#                       and grid crossvalidation params ~
#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This is a utility function to get the hyper parameters selected and organised.
def read_params(ml_type):
    ml_path = CV_type + ml_type + "_data.txt"
    param_file = open(ml_path, 'r')

    param_dict = {}
    for x, line in enumerate(param_file):
        if line.split():
            if line.split()[0] == '#':
                pass
            else:
                
                line_list = line.split()
                if line_list[1] == "uniform" or line_list[1] == "expon" or line_list[1] == "log":
                    print(line_list[1])
                    param_type = line_list.pop(0)
                    if line_list[0] == "uniform":
                        param_dict[param_type] = uniform(float(line_list[1]), float(line_list[2]))
                    if line_list[0] == "expon":
                        param_dict[param_type] = expon(float(line_list[1]), float(line_list[2]))
                    if line_list[0] == "log":
                        param_dict[param_type] = loguniform(float(line_list[1]), float(line_list[2]))
                    
                else:
                        
                    param_type = line_list.pop(0)
                    new_list = []
                    for element in line_list:
                        try:
                            elem = int(element)
                        except ValueError:
                            elem = element
                        if type(elem) != int:
                            try:
                                elem = float(element)
                            except ValueError:
                                elem = element
                        if type(elem) != float and type(elem) != int:
                            try:
                                elem = tuple(map(int, element.split(',')))
                            except ValueError:
                                elem = element
                        
                        new_list.append(elem)
                        
                    param_dict[param_type] = new_list
    
    return param_dict



#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                   ~     Main parse declare vars      ~
#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


all_params_list = np.array([{} for y in range(len(QM_names))]).astype(object)       # all results from cv
best_params_dict = np.array([{} for y in range(len(QM_names))])                     # best parameters for each QM and estimator
train_score_list = np.array([0 for y in range(len(QM_names))]).astype(float)        # all training scores from CV for QM
test_score_list  = np.array([0 for y in range(len(QM_names))]).astype(float)        # all test scores from CV for QM

results_df_list = [[] for y in range(len(estimator_types))]                         

file_loc_results = './CrossValidation/Results.csv'
file_loc_param = './CrossValidation/bestParams.csv'

csv_macro_results = open(file_loc_results, 'w')                 # Initialise results .csv file
write_result = csv.writer(csv_macro_results)
csv_macro_results.write(',')
write_result.writerow(QM_names)

csv_macro_param = open(file_loc_param, 'w')                     # Initialise best parameter .csv file
write_param = csv.writer(csv_macro_param)
csv_macro_param.write(',')
write_param.writerow(QM_names)




engaged_estimators = []

for x in estimators:
    if x in inp_prams:
        
        if inp_prams[x] == 1:
            engaged_estimators.append(x)
        
print("\n\nEstimators to use:\n")
print(engaged_estimators)

#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                   ~     Main loop for ML runs      ~
#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for estim in (engaged_estimators):
    params = read_params(estim)
    
    csv_macro_results.write(estim + "\n")
    csv_macro_param.write(estim + "\n")
    
    print("esimator being used:")
    print(estim)
    for qm in range(len(QM_names)):
        
        
        all_params_list[qm], best_params_dict[qm], \
            train_score_list[qm], test_score_list[qm] = \
                run_ml(estim, qm, params)
         
    
    print("End of %s Estimator run" % estim)
    
    result_df = pd.DataFrame({'Best_Train': train_score_list, 'Best_Test':test_score_list}).T
    result_df.to_csv(csv_macro_results,header=False)
    param_df = pd.DataFrame({'Best_parms': best_params_dict}).T
    param_df.to_csv(csv_macro_param,header=False)
    
    print(result_df)
    
    

