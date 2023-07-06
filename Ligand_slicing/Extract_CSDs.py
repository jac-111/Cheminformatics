# -*- coding: utf-8 -*-
"""
Spyder Editor

A simple script to extract the desired CSDs QM and RAC data from the whole TMQM
QM RAC .csv. 

POSSIBLE EXTENSIONS:
    - have the code read a keyword, given at execution, that relates the desired 
      structure to the .csv output ~ making the code more interchangeable for 
      different ligands.
"""

import pandas as pd
import os
import sys



# CHECK STANDARD INPUT
if len(sys.argv) > 2:
    print("\n**** ERROR ****\n\nToo many arguements passed to script.\n")
    print("Only pass the .csv file with the CSDs you want to isolate\n")
    print("See /Examples directory for correct formatting.")
    quit()
elif len(sys.argv) == 1:
    print("\n**** ERROR ****\n\nPass .csv with CSDs to isolate to this script in command line")
    print("e.g.\npython Extract_CSDs.py H2O.csv")
    print("See /Examples directory for correct formatting.")
    quit()
elif len(sys.argv) == 2:
    print("working with %s" % sys.argv[1])
    CSD_csv = sys.argv[1]
    
if not os.path.exists('./ligand_indiv'):                                        # make local directory to store
    os.makedirs('./ligand_indiv')                                               # desired CSD .xyz files

tmqm_df = pd.read_csv('../Processed_data/tmqm_RACs-20.csv',header=None,dtype=str)  # Read tmqm_RAC data as dataframe
ligand_CSDs = pd.read_csv(CSD_csv, header=None,dtype=str)                       # Read list of CSDs for slicing

tmqm_csd = tmqm_df.iloc[:,1]                                                    # get all TMQM CSDs in list
ligand_CSDs = ligand_CSDs.values.tolist()                                       # Get all desired CSDs in list
ligand_CSDs = [item for sublist in ligand_CSDs for item in sublist]             # Clean list
appearance = tmqm_csd.isin(ligand_CSDs).astype(int)                             # Find match in TMQM CSDs
locat = appearance.index[appearance == 1].tolist()                              # Extract match CSD


ligand_CSDs_all = tmqm_df.iloc[locat]                                           # Get matched CSDs QM and RAC data
print(ligand_CSDs_all.iloc[:,1])                                                # print matched CSDs

csv_name = sys.argv[1].split('/')[-1].split('.')[0]
tmqm_csv_name = csv_name + "_tmqm.csv"

ligand_CSDs_all.to_csv(tmqm_csv_name,header=None, index=False)                  # write matched QM and RAC

for csd in ligand_CSDs_all.iloc[:,1]:                                           # for all matched CSDs
    var = str(csd)
    var_str = "../Processed_data/tmqm_indiv/*" + var + "*"                      # copy the matched CSD xyz to
    os.system("cp %s ./ligand_indiv/" % var_str)                                # local directory

