# -*- coding: utf-8 -*-
"""
Created on Tue May  3 21:38:18 2022

@author: jacob
"""

import os
import pandas as pd

ligand_spe =                                                # the Final SPE for the ligand that has been sliced e.g. float -924.925739523055
ligand_name = ''                                       # name for the .csv file with E_dissociation data


tmqm_df = pd.read_csv('../Processed_data/tmqm_RACs-151.csv')       # location of all tmqm_data

dir_path = os.path.dirname(os.path.realpath(__file__))                  # check current working directory

split_dir = dir_path.split('/')                                         # isolate directory from name

if split_dir[-1] == 'orca_jobs':                                        # check working directory name
    print("working directory = " + split_dir[-1])                       # SIMPLE CHECK
    check = 'yes'
else:
    print("working directory = " + split_dir[-1])                       # if not in 'correct' directory
    print("\nWill now use all_final_spe.txt files within:\n")           # 
    print("og_orca_strucs\nsliced_strucs_orca\n")
    print("to generate new tmqm dataset with Energy of Dissociation values\n\n")
    check = input("Does this directory contain your orca job directories? (yes or no)\n")   # user input request for yes or no to check in subdirectories are present
    print(check)
    
if check != 'yes':
    quit()
    
og_orca_path = dir_path + "/og_strucs_orca/all_final_spe.txt"           # path to og_strucs_orca/all_final_spe.txt
sliced_orca_path = dir_path + "/sliced_strucs_orca/all_final_spe.txt"   # path to sliced_strucs_orca/all_final_spe.txt

og_spe_file = open(og_orca_path, 'r')                   # open final spe files
sliced_spe_file = open(sliced_orca_path, 'r')

og_spe_vals = {}
slice_spe_vals = {}

for x, line in enumerate(og_spe_file):
    og_spe_vals[int(line.split()[0].split('_')[0]) + 1] = float(line.split()[1])        # collect csd index and spe value - full struc

for x, line in enumerate(sliced_spe_file):
    slice_spe_vals[int(line.split()[0].split('_')[0]) +1] = float(line.split()[1])      # collect csd index and spe value - sliced struc


E_dis_dict = {}

for keys in og_spe_vals:
    E_dis_dict[keys] = (((slice_spe_vals[keys] + ligand_spe) - og_spe_vals[keys]) * 627.51) # calculate E_Dis for each structure


isolated_tmqm = (tmqm_df.loc[tmqm_df['runs'].isin(og_spe_vals)])                        # get df of tmqm and rac data for
                                                                                        # E_Dis molecules only
avail_Racs = isolated_tmqm.runs.to_list()
                                         
avail_edis = {}
for x in avail_Racs:
    avail_edis[x] = E_dis_dict[x]

e_dis_vals = []

for keys, values in avail_edis.items():
    e_dis_vals.append(values)

                                  
isolated_tmqm.insert(loc=9, column='E_dis', value = e_dis_vals)                              # Append E_Dis to df

isolated_tmqm.to_csv(ligand_name+'_E_dis.csv', index=False)                             # print new small df

print(isolated_tmqm['E_dis'])
    