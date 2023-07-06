# -*- coding: utf-8 -*-
"""
Spyder Editor

Program designed to separate the TMQM xyz data into individual .xyz files
for further comprehension. The original TMQM xyz data contains all xyzs separated into
2 files. I have made this script to split these files into individual .xyz.

POSSIBLE EXTENSIONS:
    - make more generic in handling of input data e.g. pass $PATH at
      commandline execution
    - less tmqm specificity
    - these loops can be easily reduced into a single function by
      putting mol_data and mol_data2 into a list and looping over it
      ^ superficial improvement ~ better coding practice
"""

                
import os, sys                                                # Libraries for python to interface with commandline

if not os.path.isdir('../Processed_data/tmqm_indiv'):           # Make individual xyz directory if
        os.makedirs('../Processed_data/tmqm_indiv')             # one doesnt already exist

file_ext = '../Processed_data/tmqm_indiv/'                      # path to printing site of indiv .xyz file
mol_data = open('../Raw_data/tmQM_X1.xyz', 'r')                 # load path 1
mol_data2 = open('../Raw_data/tmQM_X2.xyz', 'r')                # load path 2



tick = 0                                                        # used to count the appearance of new xyz
X1 = str(42375)                                                 # num of xyzs in X1 - used to track progress
X2 = str(44290)                                                 # ^
track = 0                                                       # used to restart search for line_0 and line_1 of xyz
my_str=''                                                       # blank str appended with .xyz data

for i, line in enumerate(mol_data):                             # enumerate = parse through all lines

    if track == 1:                                              # if line 1 of xyz:
        CSD_temp = line.split()[2] + '_'                        # collect CSD
        
    if line != '\n':                                            # if not at end of xyz
        my_str += line                                          # append line to xyz str
        track += 1                                              # track .xyz num lines
        
    if line == '\n' or line == '':                              # if end of xyz
        num_str = str(tick)                                     # store index of .xyz
        temp_num = num_str.zfill(6) + '_'                       # apply number to 6 digits i.e. 1 --> 000001
        file_name = file_ext + temp_num + CSD_temp + ".xyz"     # path to write .xyz
        temp_file = open(file_name, 'w')                        # open .xyz 
        temp_file.write(my_str)                                 # write .xyz file
        temp_file.close()                                       # close .xyz
        my_str = ''                                             # reset str to store xyz data
        tick += 1                                               # count xyz index
        track = 0                                               # reset track to identify line 1 for csd code   
        
        sys.stdout.write('\r number of molecules read = '+str(tick) + "/" + X1)     # read out progress
        sys.stdout.flush()                                      # flush screen to update progress



for i, line1 in enumerate(mol_data2):
    if track == 1:
        CSD_temp = line1.split()[2] + '_'
    if line1 != '\n':
        my_str += line1
        track += 1
    if line1 == '\n':
        num_str = str(tick)
        temp_num = num_str.zfill(6) + '_'
        file_name = file_ext + temp_num + CSD_temp + ".xyz"
        #print(file_name)
        temp_file = open(file_name, 'w')
        temp_file.write(my_str)
        temp_file.close()
        my_str = ''
        tick += 1
        track = 0
        sys.stdout.write('\r number of molecules read = '+str(tick - 42375) + "/" + X2)
        sys.stdout.flush()
        

print("Final tick = ", tick)
