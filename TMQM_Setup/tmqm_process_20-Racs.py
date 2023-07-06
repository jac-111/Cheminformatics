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

This script is solely for the prodcution of the 20 RACs describing the whole molecule
as described in the above tutorial. tmqm_process_151-Racs.py will generate the 151 RACS
described in the Kulik paper:
    "Resolving Transition Metal Chemical Space: Feature Selection for Machine Learning
    and Structure-Property Relationships"

@author: jacob
"""

#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                   ~       Declare Libraries          ~
#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import glob, os, sys
from molSimplify.Informatics.autocorrelation import*
from molSimplify.Informatics.misc_descriptors import*
from molSimplify.Informatics.graph_analyze import*
import pandas as pd


#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                   ~       molSimplify Class          ~
#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


class quickrun:                                     # A Class is a separate entity from the main code
    def __init__(self,name,mol,prop):               # Errors do not drop out in std.err
        self.name =  name
        self.mol = mol
        ## descriptors
        self.set_desc = False
        self.descriptors =  list()
        self.descriptor_names =  list()
        self.prop = prop

    def get_descriptor_vector(self,loud=False,name=False):                                      # Get the RAC descriptions
        
        results_dictionary = generate_full_complex_autocorrelations(self.mol,depth=4,loud=loud) # molsimplify function for 20 RACs
        
        
        self.append_descriptors(results_dictionary['colnames'],results_dictionary['results'],'f','all') # Append descriptors for
                                                                                                        # final .csv
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



def write_data_csv(list_of_runs, check):                                # Simple write to .csv function
    with open('../Processed_data/tmqm_RACs-20.csv','a') as f: 
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
        
        

print('Finding file names..')                                               # Searches the directory tmqm_indiv for the file paths
target_paths=sorted(glob.glob('../Processed_data/tmqm_indiv/*.xyz'))        # to the individual .xyzs for parsing
                                                                            # sort them numerically into a list        
print('found ' + str(len(target_paths)) + ' molecules to read')             # Describe number of mols to read
print('Sample: ' + target_paths[0])



Raw_QM = pd.read_csv('../Raw_data/tmQM_y.csv',sep=';')                      # collect tmqm data
 
prop_strings = ['CSD_code', 'Electronic_E', 'Dispersion_E', 'Dipole_M',\
                'Metal_q', 'HL_Gap',\
                'HOMO_Energy', 'LUMO_Energy', 'Polarizability']             # Column title decleration for tmqm qm

myCsv = open('../Processed_data/tmqm_RACs-20.csv', 'w')                        # open .csv for QM and rac printing
myCsv.close()                                                               # close file (wipe .csv if already exists)
    
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
    prop = (Raw_QM[prop_strings].iloc[count-1])                             # get QM data for given index
    this_run = quickrun(name1,mymol,prop)                                   # create quickrun class with xyz data
    this_run.get_descriptor_vector()                                        # Generate RACs for this structure
    list_of_runs.append(this_run)                                           # append this_run.class to a list (LARGE)
    sys.stdout.write('\r number of molecules read = '+str(count) + "/"+str(len(target_paths))) # track mols read
    sys.stdout.flush()                                                      # keep updated progress track
    if  count % 50 == 0:                                                    # data dump
        write_data_csv(list_of_runs, tick)                                  # write qm and RAC data to csv
        del list_of_runs                                                    # list_of_runs must be cleared to stop memory crash
        list_of_runs = list()
        tick = 1

write_data_csv(list_of_runs, tick)





