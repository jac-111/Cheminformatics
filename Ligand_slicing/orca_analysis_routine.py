# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 23:33:42 2022

@author: jacob
"""

import glob, os, sys

print('Finding file names..')                                   # Searches the directory tmqm_indiv for the file paths
if len(sorted(glob.glob('./jobs/*'))) == 0:
    print("\n**** ERROR ****\n\nNo files found to be read\n")
    print("Have you executed this code within the /orca_jobs/* directory?\n")
    print("This code looks for a directory named '/jobs'")
    
    
orca_dir_paths=sorted(glob.glob('./jobs/*'))                        # find all CSD directories
print('found ' + str(len(orca_dir_paths)) + ' orca directories')    # print located directories
print('Sample: ' + orca_dir_paths[0])                               # print first CSD for a check
final_spe = [None] * len(orca_dir_paths)                            # set up list to collect final SPE from output file

count = 0 

for dirs in orca_dir_paths:                                         # loop over collected CSD paths 
    file_name = dirs.split('/')[2] + '.out'                         # get CSD output file name
    out_path = dirs + '/' + file_name                               # get path to output file
    read_out = open(out_path, 'r')                                  # read the output file
    for i, line in enumerate(read_out):                             # loop over each line in output file
        
        if line.split():                                            # if line has text in it:
            
            if line.split()[0] == "FINAL":                          # if "FINAL" is the first word
                final_spe[count] = float(line.split()[4])           # Get 5th word/number ~ THE FINAL SPE
    count += 1                                                      # count for each CSD to put into list of SPEs

macro_results = open("all_final_spe.txt", 'w')                      # results file

for i in range(len(orca_dir_paths)):  
    macro_results.writelines("%s\t%s\n" % (orca_dir_paths[i].split('/')[2], final_spe[i]))  # print output file with all SPEs

