#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 17:18:14 2022

This script grabs the index and CSD of each individual .xyz from the tmqm dataset. 
Following this, it will convert it to the molSimplify 3dmol format and generate 
the RACs for the molecule. The RACs will then be put into a .csv file along with
the respective QM data. - this dataset forms the foundation of the ML aspect of 
this project.

This code has been been adapted from Heather Kuliks group for the production of
RACs for the QM9 dataset:
http://hjkgrp.mit.edu/content/qm9-kernel-models-using-molsimplify-racs-and-r-part-1
I have made some modifications for my preferred formatting and for the tmqm dataset

@author: jacob
"""

import glob, os, sys
import random
import math
import time
from molSimplify.Informatics.autocorrelation import*
from molSimplify.Informatics.misc_descriptors import*
from molSimplify.Informatics.graph_analyze import*
import pandas as pd
import inspect


percentage = 20
decim = percentage/100

expected_num_descrip = 151

autocorrel_list = ["full_complex_autocorrelation", "all_ligand_autocorrelation",\
                           "all_ligand_delta_autocorrelation", "metal_centric_autocorrelation",\
                           "metal_delta"]
failed_struc_count = [0, 0, 0, 0, 0]
start_time = time.time()
print("--- %s seconds ---" % (time.time() - start_time))

class quickrun:                                     # A Class is a separate entity from the main code
    def __init__(self,name,mol,prop):               # Errors do not drop out in std.err
        self.name =  name
        self.mol = mol
        ## descriptors
        self.set_desc = False
        self.descriptors =  list()
        self.descriptor_names =  list()
        self.prop = prop
        
        
    def check_all_desc(self, loud=False):
        
        check_arr = list()
        try:
                
            check_arr.append(generate_full_complex_autocorrelations(self.mol,depth=3,loud=False))
            check_arr.append(generate_all_ligand_autocorrelations(self.mol,depth=3,loud=False))
            check_arr.append(generate_all_ligand_deltametrics(self.mol,depth=3,loud=False))
            check_arr.append(generate_metal_autocorrelations(self.mol,depth=3,loud=False))
            check_arr.append(generate_metal_deltametrics(self.mol,depth=3,loud=False))
            for x in range(len(check_arr)):
                for key, value in check_arr[x].items() :
                    for i in value:
                       # print((i))
                        if len(i) != 4:
                            print("\n\n FAILED DUE TO ERROR IN ")
                            print("\n" + autocorrel_list[x])
                            failed_struc_count[x] += 1
                            return 1
                        break
        except:
            return 1
            
                #print (key,value)
            #print(value)
            #print(len(value[0]))
            
        return 0
        
    def get_descriptor_vector(self,loud=False,name=True):          # Get the RAC descriptions
        #print("ping_start")
        results_dict_all = generate_full_complex_autocorrelations(self.mol,depth=3,loud=loud) # molsimplify function
        #for key, value in results_dict_all.items() :
        #    print (key,value)
        #print("ping_full_complex")
        
        results_dict_All_lig = generate_all_ligand_autocorrelations(self.mol,depth=3,loud=loud)
        results_dict_all_lig_delta = generate_all_ligand_deltametrics(self.mol,depth=3,loud=loud) # molsimplify function
        #print("results_dict_all_lig_delta")
        #print(results_dict_all_lig_delta)
        #for key, value in results_dict_All_lig.items() :
        #    print (key)
        #    print (value)
        #print("ping_all_lig")
        
        #results_dict_atmOnly = generate_atomonly_autocorrelations(self.mol,depth=4,loud=loud)
        #results_dict_coul = generate_full_complex_coulomb_autocorrelations(self.mol,depth=3,loud=loud)
        #for key, value in results_dict_coul.items() :
        #    print (key)
        #print("ping")
        #results_dict_mc_eq_ax = generate_mc_eq_ax_autocorrelation(self.mol,depth=3,loud=loud)
        #for key, value in results_dict_mc_eq_ax.items() :
        #    print (key)
        #print("ping")
        results_dict_metal = generate_metal_autocorrelations(self.mol,depth=3,loud=loud)
        results_dict_metal_delta = generate_metal_deltametrics(self.mol,depth=3,loud=loud)
        #print ("results_dict_metal_delta")
        #print(results_dict_metal_delta)
        #for key, value in results_dict_metal.items() :
        #    print (key,value)
        #print("ping_metal_only")
        
        def I_0_tear(I_0_tear_dict):
            x_track = 0
            for x in I_0_tear_dict:
                
                for y in range(len(I_0_tear_dict[x])):    
                    I_0_tear_dict[x][y] = list(I_0_tear_dict[x][y])
                
                del I_0_tear_dict[x][2][0] 
                
            #print(I_0_tear_dict)
            return (I_0_tear_dict)
        
        def I01_S0_tear(I01_S0_tear_dict):
            for x in I01_S0_tear_dict:
                for y in range(len(I01_S0_tear_dict[x])): 
                    #print(I01_S0_tear_dict[x][y])
                    I01_S0_tear_dict[x][y] = list(I01_S0_tear_dict[x][y])
                    
                    #print(I01_S0_tear_dict[x][y])
                #print("POG")
                #print(I01_S0_tear_dict[x])
                del I01_S0_tear_dict[x][2][0] 
                del I01_S0_tear_dict[x][2][0]
                del I01_S0_tear_dict[x][4][0]
                #print(I01_S0_tear_dict[x])
            return (I01_S0_tear_dict)
        
        
        
        def delta_tear(delta_dict):         # ~~~~~~~ NOTE ~~~~~~~~~~~
                                            # This function will remove empty RACs from the dictionary printed by molsimplify
                                            # This is purely for convenience for further processing - the data is empty always
                                            # and is never included in the Kulik datasets. Essentially a truncation.
                                            # Chi_0, Z_0, I_0, I_1, I_2, I_3, T_0, S_0 are torn.
            for x in delta_dict:
                for y in range(len(delta_dict[x])):
                    delta_dict[x][y] = list(delta_dict[x][y])
                    del delta_dict[x][y][0]
                del delta_dict[x][2]
            return delta_dict
        
        
        
        
        f_all_lig = {key:results_dict_All_lig[key].copy() for key in ['colnames', 'result_ax_full', 'result_eq_full']}
        lc_all_lig = {key:results_dict_All_lig[key].copy() for key in ['colnames', 'result_ax_con', 'result_eq_con']}
        
        trunc_lc_all_lig = I_0_tear(lc_all_lig)
        trunc_results_dict_metal = I01_S0_tear(results_dict_metal)
        
        trunc_results_dict_metal_delta = delta_tear(results_dict_metal_delta)
        trunc_results_dict_all_lig_delta = delta_tear(results_dict_all_lig_delta)
        
        
        self.append_descriptors(f_all_lig['colnames'],f_all_lig['result_ax_full'],'f','ax')
        self.append_descriptors(f_all_lig['colnames'],f_all_lig['result_eq_full'],'f','eq')
        self.append_descriptors(trunc_lc_all_lig['colnames'],trunc_lc_all_lig['result_ax_con'],'lc','ax')
        self.append_descriptors(trunc_lc_all_lig['colnames'],trunc_lc_all_lig['result_eq_con'],'lc','eq')
        self.append_descriptors(trunc_results_dict_metal['colnames'],trunc_results_dict_metal['results'],'f','mc')
        self.append_descriptors(results_dict_all['colnames'],results_dict_all['results'],'f','all')
        self.append_descriptors(trunc_results_dict_all_lig_delta['colnames'],trunc_results_dict_all_lig_delta['result_ax_con'], 'f', 'ax_delta')
        self.append_descriptors(trunc_results_dict_all_lig_delta['colnames'],trunc_results_dict_all_lig_delta['result_eq_con'], 'f', 'eq_delta')
        self.append_descriptors(trunc_results_dict_metal_delta['colnames'],trunc_results_dict_metal_delta['results'],'f','mc_delta')
        #self.append_descriptors(results_dict_coul['colnames'],results_dict_coul['results'],'f','all_coulomb')
        #self.append_descriptors(results_dict_atmOnly['colnames'],results_dict_atmOnly['results'],'f','atom_only')
        #self.append_descriptors(results_dict_mc_eq_ax['colnames'],results_dict_mc_eq_ax['result_mc_ax_ac'],'f','mc_ax_ac')
        #self.append_descriptors(results_dict_mc_eq_ax['colnames'],results_dict_mc_eq_ax['result_mc_eq_ac'],'f','mc_eq_ac')
        #self.append_descriptors(results_dict_metal['colnames'],results_dict_metal['results'],'f','metal')
        #self.append_descriptors(results_dict_metOx['colnames'],results_dict_metOx['results'],'f','metal_ox')
        #self.append_descriptors(results_dict_metOxEff['colnames'],results_dict_metOxEff['results'],'f','metal_ox_eff')
        #self.append_descriptors(results_dict_multiAtm['colnames'],results_dict_multiAtm['results'],'f','multi_atm')
        
        self.set_desc = True
        

    def append_descriptors(self,list_of_names,list_of_props,prefix,suffix): # Puts the descriptions in csv - obsolete from Kulik code
        for names in list_of_names:
            if hasattr(names, '__iter__'):
                names = ["-".join([prefix,str(i),suffix]) for i in names]
                self.descriptor_names += names
            else:
                names = "-".join([prefix,str(names),suffix])
                self.descriptor_names.append(names)
                
        for values in list_of_props:
            if hasattr(values, '__iter__'):
                self.descriptors.extend(values)
            else:
                self.descriptors.append(values)
        #print(self.descriptor_names)
        #print(self.descriptors)


def write_data_csv(list_of_runs, check):
    with open(racs_csv,'a') as f: 
        if check != 1:
            f.write('runs,'+",".join(prop_strings)+',')
            n_cols = len(list_of_runs[0].descriptor_names)
        
            for i,names in enumerate(list_of_runs[0].descriptor_names): # writing the RAC descriptor names
                if i<(n_cols-1):
                    f.write(names+',')
                else:
                    f.write(names+'\n')
        for runs in list_of_runs:
            try:
                f.write(runs.name)  # Writes the molecule index and CSD_code
                for properties in runs.prop:
                    f.write(','+str(properties))
                for descrip in runs.descriptors:
                    f.write(','+str(descrip))
                f.write('\n')
            except:
                pass
#def check_descrip_arr_len():
    
def write_failed_csv(job_path, failed_index):
    with open("../Proecessed_data/%s/failed_strucs_%s.csv" % (percentage, percentage), 'a') as f:       
        print()

print('Finding file names..')                                               # Searches the directory tmqm_indiv for the file paths
target_paths=sorted(glob.glob('../Processed_data/tmqm_indiv/*.xyz'))        # to the individual .xyzs for parsing
#print(target_paths)
idicies = random.sample(range(len(target_paths)), math.floor(len(target_paths) * decim))                          # sort them numerically into a list        
#print(idicies)
sampled_paths=[target_paths[i] for i in idicies]
print('found ' + str(len(target_paths)) + ' molecules to read')             # Describe number of mols to read
print('Sample: ' + target_paths[0])



Raw_QM = pd.read_csv('../Raw_data/tmQM_y.csv',sep=';')                      # collect tmqm data
 
prop_strings = ['CSD_code', 'Electronic_E', 'Dispersion_E', 'Dipole_M',\
                'Metal_q', 'HL_Gap',\
                'HOMO_Energy', 'LUMO_Energy', 'Polarizability']             # Column title decleration for tmqm qm

racs_csv = '../Processed_data/tmqm_RACs-151.csv'
failed_csv = '../Processed_data/tmqm_RACs-151_FAILED_STRUCS.csv'
    
myCsv = open(racs_csv, 'w')                        # open .csv for QM and rac printing
my_failed = open(failed_csv, 'w')         # store CSD of structures that failed to be read by molSimplify
with open(failed_csv, 'w') as f:
    f.write("Index,CSD\n")

myCsv.close()                                                               # close file (wipe .csv if already exists)
my_failed.close()
    
list_of_runs = list()
count = 0                                                                   # count number of geometries for memory dump
tick = 0                                                                    # track position of memory dump


for geoms in target_paths:                                                  # get each .xyz path
    count += 1                                                              # track index of .xyz
    name = geoms.split("_")[2]                                              # get CSD
    name1 = str(count)                                                      # .xyz index as string for 'quickrun'
    
    geomfile =  open(geoms, 'r')
    
    mymol = mol3D()                                                         # instance MolSimplify mol3d
    mymol.readfromxyz(geoms)                                                # load xyz data to instance
    #print(mymol.printxyz())
    #print("BOO")
    mymol.findMetal(transition_metals_only=False)[0]
    prop = (Raw_QM[prop_strings].iloc[count-1])                             # get QM data for given index
    this_run = quickrun(name1,mymol,prop)                                   # create quickrun class with xyz data
    fault_check = this_run.check_all_desc()
    
    #if len(this_run.descriptor_names) != expected_num_descrip:
    #    fault_check = 1
    
    
    if fault_check != 1:
        
        this_run.get_descriptor_vector()                                        # Generate RACs for this structure
        #print(len(this_run.descriptor_names))
        #print(fault_check)
        
        list_of_runs.append(this_run)                                           # append this_run.class to a list (LARGE)
        sys.stdout.write('\r number of molecules read = '+str(count) + "/"+str(len(target_paths))) # track mols read
        sys.stdout.flush()                                                      # keep updated progress track
        
    else:
        print("\nFailed Structure = %s\n" % geoms)
        with open(failed_csv,'a') as f: 
            
            f.write(geoms.split("/")[3].split("_")[0])
            f.write("," + geoms.split("/")[3].split("_")[1] + "\n")
            
            
            
    if  count % 50 == 0:                                                    # data dump
        write_data_csv(list_of_runs, tick)                                  # write qm and RAC data to csv
        del list_of_runs                                                    # list_of_runs must be cleared to stop memory crash
        list_of_runs = list()
        tick = 1
        print("\ndistribution of structure failures:\n")
        print(autocorrel_list)
        print(failed_struc_count)
        print("\n\n--- %s seconds ---" % (time.time() - start_time) + "\n")

print("\n\nMETA ANALYSIS")

write_data_csv(list_of_runs, tick)
myCsv.close()

all_desc_pd = pd.read_csv(racs_csv)

print("number of successfully read structures:\n")
print(all_desc_pd.shape[0])
all_failed_pd = pd.read_csv(failed_csv)
print("number of failed structures:\n")
print(all_failed_pd.shape[0])
print("\nnumber of descriptors used:")
print(all_desc_pd.iloc[:,10:].shape[1])




