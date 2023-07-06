SOME NOTES ON LIGAND SLICING:

To start, you will need a csv with your CCSD CSD codes. e.g.
ABEMAX
ACIVEQ
ACIVIU
AGOJUC
AGUTOO
AGUTUU
AGUXIM
AGUXOS
AGUXUY
...

Then execute:
python ./Extract_CSDs.py Your_CSV.csv

This will pull the .xyz files from the ../Processed directory/tmqm_indiv/ directory
if they are present. The associated QM data will also be pulled.

From here, you will then need to provide the .xyz for the ligand you would like to slice from 
the structures. I would reccomend taking the structure from one of the xyz files that have
been identified as containing the ligand.

Then execute:
python ./slice_ligand.py ligand_xyz.xyz

You will then have a directory of sliced structures and a directory of faulty structures.

This part must be done by hand.
 - go into the faulty structures and fix them. Thereasons for their failure should be printed
   and an explaination for this can be found in the original TMQM-ML_Doc.txt

Once you've fixed all of these structures and added them to the correctly sliced structure directory,
you can now generate the orca jobs for processing.
prior to execution, check your system requirements for a sequential execution of orca jobs. I can't
predict this part, so I have left the settings I used as default.
This shouldn't be too difficult, just edit the strings within the program to suit your needs/hardware.

Then execute:
python ./make_orca_inps.py

You will then be able to pass the orca_jobs directory to the HPC or machine you can execute orca on.
From there, execute the bash script within orca_jobs.

Once they've all run and you're happy with the results, you can put the orca_analysis_routine.py script
in the orca_jobs directory to pull the final energies of the whole and sliced structures.

You are now ready to compile the Energy of dissociation. 

First you will need to open and edit the make_EDis_tmqm_dataset.py script to have the spe for the single
ligand at line 11 and choose a name for the final csv file at line 12.

Now, put the make_EDis_tmqm_dataset.py into the orca_jobs directory and execute to calculate the Energy
of dissociation and put it into a format the ML_scripts can understand.

If all has gone to plan, you should have correctly formatted data for ML observation as covered in the ML_scripts
section.

Good-Luck.