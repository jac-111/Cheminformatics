     _____ ___  ___ _____ ___  ___        ___  ___ _     
    |_   _||  \/  ||  _  ||  \/  |        |  \/  || |    
      | |  | .  . || | | || .  . | ______ | .  . || |    
      | |  | |\/| || | | || |\/| ||______|| |\/| || |    
      | |  | |  | |\ \/' /| |  | |        | |  | || |____
      \_/  \_|  |_/ \_/\_\\_|  |_/        \_|  |_/\_____/
                                                         
                                                         
														 
TMQM_Analysis_doc written by JE 22/07/2022


I Finished working on this project with some ideas for the future, please see the final words of this document for a few starting ideas for picking up this project!

PROJECT OVERVIEW

This package has been primarily developed to handle and analyse the TMQM dataset through machine learning methods. these data can be found at:

"tmQM Dataset—Quantum Geometries and Properties of 86k Transition Metal Complexes", 
David Balcells and Bastian Bjerkem Skjelstad
Journal of Chemical Information and Modelling, 2020, 60, 12, 6135–6146

TMQM contains Quantum Data (QM) and xyz data relating to a selection 
of mono-nucleaic Transition metal (TM) complexes. Each structure has been taken from Cambridge Crystallography Data Center (CCDC) along with the CSD code. This package enables the xyz data to be used to generate attributes for the prediction of the QM data through machine learning algorithms. To do this, we utilise the Revised Auto-corellations (RACs) 
developed by the Kulik group to describe the molecular structure with a number of terms. These terms are then used to predict the QM data. We have principally used molSimplify for this:

### MOLSIMPLIFY README PAGE
https://molsimplify.readthedocs.io/en/latest/index.html
 
### KULIK TUTORIALS PAGE 
http://hjkgrp.mit.edu/tutorials/

This document provides an overview for each program.
These scripts are coded in python and depend on a number of libraries.
It is recommended that conda be installed on the local system, be that full
conda or miniconda. Conda can be used for the installation of the RDKit and
Molsimplify libraries, both of which are essential. For information regarding
the setup of the environment, please read the Prerequisite section of this doc.
I have written this section to be as generic as possible for a linux
environment. Note I have used ubuntu Debian type.

This code is focussed on processing data to interface with the python ML library scikit. However, I have also developed a set of codes to identify and remove specific ligands from the TM complexes. We have done this to see how far we push the application of RACs for QM prediction. These codes are able to generate orca jobs for the selected TM structures with and without the seleceted ligand. The orca jobs can be executed to calculate the energy of dissociation via a rudimentary Bash script. The energy of Dissociation is then employed as a predictable molecular feature.

Some of the code has been optimised from code that was developed in
Google Colab, where each terminal is instanced, meaning after each reconnection the native packages are reset. Therfore conda and the required libraries must be installed each time. This will become more important as the teaching aspect is implemented.

This package has been written in context of a local linux distribution. I planned to write a colab package of the same type although currently some of the necessary libraries (molSimplify) are not available in colab to date.


# PACKAGE BREAKDOWN


The contents of this package are split into 5 directories. Overall, there are
4 key areas different areas to handle for TMQM and 4 other directories.

TMQM HANDLING:

 - Raw_data/       - This contaisn the TMQM .xyz and QM .csv files only.

 - TMQM_Setup/ 	   - This contains the programs that parse and separate the TMQM data sets and store them in the /Processed_data/ directory in the necessary format.

 - Processed_data/ - Processed data is held here ~ Both individiual .xyz files and the .csv file containing the QM and RAC data.

 - Ligand_slicing/ - All functions required to handle the removal of ligands, their specific QMs and RACs and the setup of the Orca jobs required for the calculation of the energy of dissociation.

## OTHER ML PROGRAMS:

 - ML_scripts/     - This directory contains a number of different ML scripts.

 - Legacy/	   - Contains old files which may have some useful code parts within them. - These directories are prevailent within all directories.

 - Clustering/ 	   - Analysing the TMQM data through ML clustering methods. Including Kmean method.

 - Kulik_data/     - Some analysis and preprocessing routines that use the data provided by the Kulik Research group.

--------------------------------------------------------------------------------------


This document will walk through the use of this package IN ORDER OF NECESSARY EXECUTION, from system requirements through to the ML description.


It is important that some programs be executed before others. For example, for ligand separation, individual .xyzs must be generated from the /TMQM_Setup/ scripts in the processed_data directory before execution. This is also true for the ML scripts.

Each code should be executed within the directory heirarchy as given.
When certain scripts are executed, they will be looking for data in:

e.g. ../processed_data/ 	~ If this is not available, the codes will not work.

All nuance in execution and personalisation are covered here. 


	****				TLDR				****

	****	THIS DOCUMENT GOES IN ORDER OF NECESSARY EXECUTION	****

	****	MAINTAIN DIRECTORY POSITIONING FOR PROPER EXECUTION	****




# Local Prerequisites

I have installed a fresh version of Ubuntu 22.04 to install Conda and associated libraries.

CONDA
To install Conda, navigate to the conda page and download the bash installer of choice. then install via:
bash Miniconda3-py39_4.12.0-Linux-x86-64.sh

RDKIT
To install the peripheral libraries:
conda install -y cmake cairo pillow eigen pkg-config
conda install -y boost-cpp boost py-boost
conda install -y gxx_linux-64
conda install -c rdkit rdkit

OPENBABEL
conda install -c conda-forge openbabel

MOLSIMPLIFY
conda install -c hjkgroup molsimplify

For my fresh install, these commands bring in all libraries that are called for this package. I have also found that the whole package will work after these sections have been installed.


# Execution

All codes have been developed with python 3 in an ubuntu environment. The heirarchy of directories is intentional and all codes should be executed within their directory as such:

python ./tmqm_process.py

Of course, these codes can be modified to work on an individual basis, however to work within this package. Extension should be fairly straight forward on an ad hoc basis. I have annotated the codes and tried to keep them as short as possible.




# TMQM_Setup

Codes designed to handle the TMQM dataset and produce RACs


## /TMQM_Setup/tmqm_xyz_split.py

This program parses the TMQM xyz data to produce individual .xyz files required for further structure comprehension in the /Processed_data/ directory.

The default xyz files provided by the TMQM xyz data comprises of 2 xyz files containing all of the xyz structures of the molecules within the dataset. To discretise these to allow for individual xyz to be grabbed when needed, this script breaks the xyzs into single files. File names relate to their
position in the original master xyz file and the associated csd 'code'.

With some small modifications, a different dataset could be managed but pay attention to the breakers between xyz files. In the case of TMQM, this is denoted by a blank string ''. For context, openbabel uses '$$$$' as
an xyz breaker in .xyz files containing more than one structure.


## /TMQM_Setup/tmqm_process_20-Racs.py & /TMQM_Setup/tmqm_process_151-Racs.py

These codes generates the RAC descriptors for our tmqm .xyz files. To do so, MolSimplify must be installed into the environment. Refer to the prerequesites or the MolSimplify website for installation instructions.

These codes iterate through all of the .xyz files within the /Processed_data/tmqm_indiv/ directory and determines the RAC through the molsimplify.generate_*****_autocorrelations() method within the  
MolSimplify library. The program will then map these values to the corresponding QM values and collate them in a .csv file as /Processed_data/tmqm_RACs.csv. This file forms the basis for many ML cycles.

Both codes produce a .csv with the same formatting.

***		tmqm_process_20-Racs.py 		***
Will produce the 20 RACs discussed in the Kulik tutorial:
http://hjkgrp.mit.edu/content/qm9-kernel-models-using-molsimplify-racs-and-r-part-1
These are whole molecule descriptors.

***		tmqm_process_151-Racs.py 		***
This code is slightly more involded than the 20 RACs counter part. As TMQM contains crystal structures, for some structures, molSimplify fails to comprehend the .xyz correctly, so some error management has
been implemented. 

Also, the 151 RACs are cut down from the original 180 produced by molSimplify. These cut values are either 0 or constant and are removed in the Kulik paper:
"Resolving Transition Metal Chemical Space: Feature Selection for MAchine Learning and structure-Property Relationships" 


# ML_scripts

This directory contains a single .py script that handles the tmqm RACs data produced by the tmqm_process***.py. A number of options are available to this code including fraction of original dataset to use, training split, active esimators etc. These choices are expressed in the input file - ML_input.txt.

This input file must be provided at execution. The program will ignore all blank lines and those containging '#' at the first character. parameter order does not matter. unrecognised parameters with throw error.

EXECUTION
python ./ML_g4.py ML_input.txt

OVERVIEW
ML methods are implemented with the scikit-learn library. By default, the available esimators are:
KNeighbour regression
Kernel Ridge Regression
MLP Regression
PLS Regression
Stochastic Gradient
Support Vector Machine

This code will execute CrossValidation for each estimator and print the results into the /CrossValidation directory. 

INPUT FILE (ML_input.txt)
This file contains inforation relating to the parameters of the ML runs i.e. CV style, estimators, training split etc. This file must be provided at execution. Edit this file to meet your specifications. Here is the input file format with explainations of each term infront of # symbol
```
#	Here are the input parameters for the ML script to parse
#	White space does not matter

#	Splitting TMQM dataset

TMQM_portion	1				# FLOAT # amount of dataset to use. 1 = 100%
training_split	0.6				# FLOAT # training split. 1 = 100%
num_racs	151				# INT	# given number of RACs
num_data_bins	5				# INT	# number of bins for Histogram

#	Cross validation search settings

CV_folds	5				# INT	# number of crossvalidation folds
CV_grid		0				# INT	# 1 = USE GRID CV
CV_random	1				# INT	# 1 = USE RANDOM CV
CV_random_iter	100				# INT	# number of iterations for random CV

#	Which ML algorithms to activate
KNeighboursRegression	0			# INT	# 1 = INCLUDE ML METHOD
PLSRegression		0
StochasticGradient	0
KRR			0
SVR			0
MLPRegression		1

#
CrossValidation 	1
#
```


LIB DIRECTORY
Lib/Grid/:
This set of txt files describe the parameters for a grid search based crossvalidation.
Lib/random/:
These txt files describe parameters for a random CV search.
NOTE: These files require descriptions of distribution limits in the python.scipy.stats.expon format

FORMATTING
This can handle any number of RACs, assuming that the correct number is given in the ML_input.txt file. There can also be any number QM properties. 
HOWEVER, the first 2 columns will always be disregarded as they were originally used to describe the molecule index and CSD code.
Here is the original column setup:
ID, CSD, Electronic_E, Dispersion_E, Dipole_M, Metal_q, HL_Gap, HOMO_Energy, LUMO_Energy, Polarizability, RACS

OUTPUT
This script will produce a number of .csv files in the /CrossValidation/ directory.
The Results.csv file describes the coefficient of determination for each QM property for each estimator. In the same format, the bestParams.csv file gives the best parameters from the Crossvalidation step for the training score.



# Ligand_slicing

Codes for handling the identification, removal and processing of Ligands from the TM complexes. These codes also produce the orca files for the calculation of the energy of dissociation.
## PREWORD:
This code relies on the user chosing a ligand and finding the tm complexes that contain this ligand in the CCSD database prior to execution. The list of CSDs that contain these ligands must be stored in a .csv file as shown in /Ligand_slicing/n-heterocyclic_CSD.csv example. This will allow the
   program to identify which QM, RAC and xyz data to handle for the production of the orca jobs.  I have been able to do this easily through the CCSD program called Conquest. When all of the CSDs have been identified, the .csv produced by default by Conquest can be used in these codes. This set of codes could be extended to screen each tmqm structure. However, this would be time consuming and also runs the risk of missed structure comprehsion as we know the tmqm dataset is taken from the ccsd.

 ###	EXAMPLE DATA	
This section of the package has been written in context of n_heterocyclic Carbene ligand. We beleived it to be a good choice for the development of this program. I have kept the outputs from these slicings within the /Ligand_slicing/Examples/ directory for reference.
   

## /Ligand_slicing/Extract_CSDs.py

This is a relatively simple code that cherry picks the associated .xyzs and QM data for 
CSDs that appear in a .csv.

This code is used to collect the .xyz files that contain the desired ligand within the TMQM (QM-RAC) dataset. It will then compile and print the QM data into a separate .csv file (e.g. TMQM_n-heterocyclic.csv) for further processing. The code will also identify their .xyz file from the collection of individual TMQM .xyzs and copy them to a local directory (ligand_indiv). 

This file is looking at files in the Processed_data/ directory. 


THE SCRIPT MUST BE EXECUTED WITH THE .csv FILE CONTAINING THE CSDs TO ISOLATE AS FOLLOWS:

EXECUTION:
python Extract_CSDs.py CSDs.csv

See the /Example/ directory for an example use case.




## /Ligand_slicing/slice_ligand.py

This code reads the .xyz of the ligand to extract, identifies its presence in the .xyz structure and removes it. These 'sliced' structures are then collected in a separate directory.

This code requires the individual ligand .xyz file to be passed at execution:

EXECUTION:
python ./slice_ligand.py ligand.xyz


**** 	THIS CODE IS NOT BULLET PROOF 	****

I have tried to catch as many errors and faults as possible. This code relies
on the structure comprehension of openbabel and RDKit ~ which can both be problematic. This is not at the fault of the programs, these are very powerful and sophisticated systems. However, TM complexes are very diverse and difficult to manange so, inevitably; the simple rules for organic character 
these programs are built on tends to fail, producing a variety of errors.

Here is the flowchart for  this script in the context of a single structure. I have
outlined the errors I have encountered and their nature. Unfortunately, to the best of my ability, they cannot be avoided.


## STAGE 1 : HANDLE LIGAND STRUCTURE (easy)
 - OpenBabel converts the LIGAND .xyz file structure to .mol file
 - RDKit reads ligand .mol file into Mol in memory within code. Mol is an RDKit object.
 - RDKit converts ligand to Smarts 
 - RDKit converts the Smarts back to Mol
	~ This makes the structure more generic and permits for conformational difference in the TM complex structure.
	~ I have found this to be entirely necessary for the extraction of ligands. There may be issues with complex ligands, however, n_heterocyclic is relatively complex and has been handled very well, though this is not a guarentee. I have made the code print the Mol > Smarts > Mol structure for inspection. It is saved with the same file name + _simple.mol.

## STAGE 2 : HANDLE TM COMPLEX STRUCTURE (hard)
 - OpenBabel converts the individual TMQM .xyz file to .mol format and stores in /temp_folder/
 - This is done as RDKit DOES NOT READ .xyz file format.
 	~ Due to the nature of atom coordination in TM complexes, Nitrogen or other trivalent atoms may be perceived as having 4 bonds if coordinating to the metal. Openbabel will raise an error but write the .mol file regardless with 4 bonds where there 'shouldn't be' - this is problematic for RDKit so accomodations must be made (see next paragraph). The 4 bonds associated to the Nitrogen atom can seen in the .mol file. Errors may also arrise due to the Kekulization of aromatic groups however RDKit is able to overlook these issues. (sporadically).

### RDKit reads .mol into Mol in memory
 - RDKit may simply refuse to read the .mol if the valence is incorrect (see above paragraph). These misread structures are stored in /faulty_strucs/failed_rdkit_descrip. Misreads are typically due to the error mentioned above. From my own test, around 5% of the database are unreadable in this way, so ~ 4.3K structures total. New errors may arise. This is a known issue in cheminformatics. it is hard to make a fix all, especially for TM complexes.
 - these faulty structures must be manipulated by hand. I do not know of another way around this.

### RDKit identifies the presence of the ligand within the TM complex
 - If all has gone to plan, every structure that is being observed should have the ligand to be removed. The code will test for the presence of the ligand. iIf not detected, it will pass the structure to the directory /faulty_strucs/no_match_substruc/ for user consideration. NOTE, there is a certain likelyhood that this will proc even though the ligand is infact present. This is again due to the complex structural nature of TM complexes.
- These no_match_substruc/ conformations are infact very very likely to have the ligand, however, RDKit is failing to identify. The user will need to look at these on an individual basis. it is likely you will find the structure is infact there, but in a peculiar arrangement. These situations must be managed on an ad hoc basis.

### RDKit slices the ligand from the complex

- At this point, RDKit can identify and remove 1 or more ligands. if there are multiple present in the structure, RDKit will remove them all. To the best of my knowledge, there is no way to stop this,  although there might be. If this condition is met, the multiple ligand complex is pass to  /faulty_strucs/multiple_substruc_removed for user consideration.
 - The appearance of multiple ligands is very common and an unfortunate problem. RDKit cannot decide which of the multiple appearances to remove - this must be done by the user.

## STAGE 3 : STRUCTURE RECONSTITUTION (hard)
- RDKit adds Hydrogens to the structure
	~ The RDKit Mol format in memory holds the structure of the given complex with all Hydrogens stripped. This is a necessity of the program and cannot be sidestepped. Fortunately, they can be replaced. However, this process is done without direct knowledge of previous Hydrogen position and is based entirely on valence rules ~ it is predicting where the Hydrogens need to go. THIS IS PROBLEMATIC FOR TM COMPLEXES. For instance, if the TM complex has a Hydride ligand, it is generally never replaced properly. To catch this, I have the program count the number of atoms of the whole molecule and the sliced structure + ligand. If they match, they are assumed OK. However, if they do not, then the structure must be considered by hand and is passed to the directory /faulty_strucs/non_matching_hs. 
- Sliced structures are compiled into the directory /sliced_strucs_h as .xyzs.


## STAGE 4 : MANUAL SLICING ()

###		USER TRAULING 		
###		  READ THIS		

 - Use your prefered visualiser to remove the ligands from the structures in /faulty_strucs by hand
	~ This section can be quite tedious but cannot be avoided. For the reasons mentioned above, some systems cannot be comprehended correctly by these software.
 - Once you've corrected the struc, Introduce the manually handled structures into the /sliced_strucs_h directory
	~ This allows the follow up programs to work with all of the structures.

As mentioned there are a number of issues with this code. However, these problems are a product of the complex structure of TM complexes. RDKit and openbabel are built to be generic and struggle to perfectly charcterise complexes. I beleive I have accounted for all of the errors I have encountered
with this code, however, if applied to more diverse systems, new errors may arise that I have not accounted for. For this I reccomend refering to the openbabel and RDKit manuals
https://openbabel.org/wiki/Main_Page
https://www.rdkit.org/docs/index.html

They both have worked examples and are a very useful resource.

I should also note that openbabel is also able to do substructure searches as a standalone program. This may be superior to the RDKit method, however due to my familiarity with RDKit from earlier on in my time with this project, I used RDKit. - This could be a site for program extension.

This is an unfortunate and annoying problem, but I have not been able to devise a more sophisticated solution.


## /Ligand_slicing/make_orca_inps.py

This is a relatively short utility program that is designed to generate the orca input files and bash script for the sequential execution of orca jobs on the linux machines used by Alex Hamiltons research group. As there are likely to be many jobs, we are using the B97-3c level of theory by default.
This program will generate a directory containing the full structures (og_strucs_orca) and the sliced structures (sliced_strucs_orca) and a bash script within both for mass orca execution.

B97-3c is a relatively low level of theory for the optimisation of the structure. Each orca job is contained within its own directory. I found that producing a single orca input file with multiple inputs fell over as the structures were not concordent with one another. This may be a site for program optimisation.

The header for the orca input file is generic and located at line 77.

The string for the bash script to execute the orca jobs is held within this code at line 47 and must be altered if the machine is not formatted in the same way. This is relatively straight-forward. Simply edit the string containing the orca command to your own system. Or, use the bash template to generate
your own bash script.

Upon Execution of the Bash script, it should give you an idea of progress. This is arbitrary as the job time is dependent on the number and nature of the atoms in each job. 

###	LIKELY ERRORS	

At this point some orca jobs may fail. This is due to the imperfect nature of the ligand slicing process. Mostly, I have seen errors resulting from the implicit Hydrogen replacement. For my n-heterocyclic example, I have seen the Carbon atoms of the benzene ring may also be coordinating. This means the structure has been recognised, but the number of Hydrogens is inconsistent with the spin and charge. These jobs must be addressed by hand unfortunately. The Bash script should pick up on the jobs that fail
and place the path to the failed directory into a .txt file. At the end of the execution, it should put all faulty jobs into a separate folder /err_jobs/ where each structure can be inspected.


## /Ligand_slicing/orca_analysis_routine.py

Another relatively short utility code. This code is to be used after the orca jobs have been executed. 

THIS PROGRAM SHOULD BE EXECUTED WITHIN THE /orca_jobs/ DIRECTORY WHERE THE COMPLETED ORCA JOBS ARE

All of the final Single Point Energies are taken from the  respective orca output files and compiles them into a single .txt file with the CSD in tact. 

This program should be executed within the 'sliced' and 'og' directory where the /jobs/ directory exists that holds the individual directories for orca jobs.

This program only uses basic libraries so should be useable with a basic python 3 install. - i.e. can be executed on the orca machine.


### /Ligand_slicing/make_EDis_tmqm_dataset.py

Another utility program to compile the all_final_spe.txt data found in the sliced_struc_orca and og_strucs_orca directories into the Energy of Dissociation in the same format as the .csv files found in the Processed_data
home directory.

This program should be executed in the orca_jobs directory

## IMPORTANT NOTE

THIS PROGRAM DOES NOT KNOW THE ENERGY OF YOUR LIGAND! YOU MUST ADD THIS TO THE CODE!

-  At line line 11, add the final spe for the ligand which has been sliced.
- This job needs to be done separately.
- This would be a site of extension for this package.

You must also provide the path to the tmqm_racs .csv file at line 15. This will work for both 20 and 151 racs.

At line 12, you must also choose the name of your output .csv file.

The final .csv output file can be plugged straight into the ML scripts for prediction.




# Clustering

this is the least developed section within this package. Here, we were looking at methods to cluster the data using principle component analysis (PCA) and Kmean clustering. 


### CODE FLOWCHART
The code is currently setup to store the RACs for each molecule. Following this, the Transition metal associated with each molecule are parsed from each .xyz molecule. 2 dimensional PCA analysis is then carried out on a scaled version of the RAC data. This is a dimension compression technique. K means
clustering is an ML technique that is then employed to cluster these PCA values into a given number of groups.

The associated TM Elements are then mapped to the clusters and printed to a graph via matplotlib. This code has also been used to  look at clustering of d-block groups and rows.

### NOTES ON THE CODE:
The parsing of .xyz files to identify the TM elements is highly time consuming. As a work around, I have implemented a pickling  'save' system. On your first run of this program, edit line 34 and 35 to your preferred fraction of the TMQM dataset and edit line 35
'bypass = 1'
to
'bypass = 0'

This will disengage the pickle loop. When disengaged, the code will follow another loop that will write pickled lists of Transition metal occurance in each cluster. This effectively acts as a save feature, enabling different Kmeans and PCA parameters to be used on the same dataset without parsing for the TM element everytime. This will also produce a .csv with the RAC and QM data.

After your first run, edit line 35:
bypass = 0
to 
bypass = 1

this will engage the save feature.

There are a number of sample images generated with this program in /legacy_images/.




## Kulik_data

This directory contains the data from the Kulik group papers and some preprocessing and analysis scripts:
"tmQM Dataset—Quantum Geometries and Properties of 86k Transition Metal Complexes"
,
"predicting electronic structure properties of transition metal complexes with neural networks"
and the QM9 RAC tutorial.

I have mostly worked with the "predicting electronic structure..." structural and spin splitting data. From their work, the Kulik group are able to predict the Spin splitting Energy for the given structures to within ~ 1 Kcal/mol. I have been able to replicate this.

I have also studied the effect of adding some extra information to the RACs, including ligand denticity, oxidation and spin. My numbers suggest that I get the best response is produced by only including the RAC data. 

This section of code was effectively used to benchmark my understanding of ML methods and my implementation of Molsimplify etc - I have tried to replicate the data produced in "Resolving Transition Metal Chemical Space: Feature Selection for Machine Learning  and structure-Property Relationships" at section 4.1.

This paper provides the structures of the studied molecules, but not the RACs. 

The first step is calculating the RACs.

## /Kulik_data/Setup/paper1_geom_Rac_conv.py

PREPROCESSING!

This progam reads all of the structures provided by the supplementary info in:
"predicting electronic structure properties of transition metal complexes with neural networks"


The supplied geometries are optimised structures with different levels of exchanage correlation included during the dft optimisation at high or low spin state. Here is a breakdown of the structure and file name:

5_cr_2_3_ax_h2o_eq_h2o_00.xyz

Here, from left to right, the molecular structure is given in the file name:
5 = molecule ID
cr = Core Transition metal Element
2 = Oxidation State
3 = Spin
ax_h2o = Axial molecule is h2o
eq_h2o = Equitorial molecule is h2o
00 = fraction of exchange correlation included in dft geom opt.

CODE BREAKDOWN
This code is able to generate either the 20 full_molecule RACs only, or the 151-RACs discussed within the paper.
Each RAC is determined in the same way as previously covered in /TMQM_Setup/tmqm_process_151-Racs.py.

Switching between 20 and 151 RACs can be done at line 22:
rac_20 = 1 or 
rac_20 = 0 for 151 RACs

Most of these structures are able to be read by molSimplfy as the structures were generated by the molecule construction utility of molSimplify. However, I have observed that ~10 structures fail to be correctly read. This may be due to updates to the cheminformatics utility of molsimplify.
Although, I am not certain of this.

## /Kulik_data/Setup/Assign_property_1.py

This is a small utility program to assign the Splitting energy provided by the supplementary information to the associated RAC. Within the .csv printed by /Kulik_data/Setup/paper1_geom_Rac_conv.py, each structure has a high and low spin state. This code will firstly give each RAC the splitting energy within the dataframe. This means each splitting energy appears twice in the dataframe, once for low and once for high spin state. 

The second part of this code splits this dataframe into high and low spin structures. To do so, the script will parse through the dataframe and find the same molecule data in high and low spin state. This is done by looking for
the metal, oxidating state, HF_exchange, Axial ligand and equitorial ligand. These are identical in the high and low spin state. The lower spin state is then put into a 'low-spin' dataframe, and the high spin state in the 'high-spin' dataframe. However, as this parses through the dataframe, each molecule is passed twice, once for low and once for high spin. This means each dataframe has each molecule duplicated. To manage this, pandas has a utility to remove all duplicates .drop_duplicates().



## /Kulik_data/ML/paper1_analysis_dent.py

This is the ML script that is capable of replicating the results produced in "Resolving Transition Metal Chemical Space: Feature Selection for MAchine Learning and structure-Property Relationships" at section 4.1. By default, this code will also test for the best combination of parameters:
OXIDATION : SPIN : AXIAL DENTICITY : EQUITORIAL DENTICITY : HARTREE-FOCK EXCHANGE : RACs

There are 32 unique combinations. Kernel Ridge is also engaged for this code and can be implemented with random search or grid search crossvalidation for hyper parameter tuning.

see line 163

CODE BREAKDOWN

Setup Denticity dictionary:
by default, the data included in "predicting electronic structure properties of transition metal complexes with neural networks" contains axial and equitorial molecule type information. To replicate the data given in "Resolving Transition metal Chemical..." These must be converted into the denticity values of each molecule type. This section contains the conversion chart/dict for this.

Determine all ax and eq denticity:
This follows up on the previous section and determines the ax and eq denticity and handles the dataframes




## /Kulik_data/QM_9

This is a small demo where I have carried out the QM9 RAC tutorial given by the Kulik group at:
http://hjkgrp.mit.edu/tutorials/2018-02-20-qm9-kernel-models-using-molsimplify-racs-and-r-part-1/

this is small code that incorporates some conversion from R to python. See the full tutorial for a breakdown.


#Closing words



Here are a few ideas for extension and next steps.

The Energy of dissociation should be an intrinsic value for each structure, so it is quite surprising that we have not been able to predict it. To run some more test on improving these values:

As mentioned, the structures used for the Kulik data have been generated by Molsimplify, meaning they are highly regular and perfect. Dataset quality is a huge part of ML and if I was able to continue with this project, I would like to see if the the Energy of dissociation could be predicted from well defined regular TM structures generated by molsimplify. I would also include the test of including peripheral data alongside the RACs as covered in the /Kulik_data/ML/paper1_analysis_dent.py section.

For the Kulik paper, the number of structurally distinct molecules is ~150, so getting a similar number of molecules shouldnt be too difficult.

I think this would be a method to deterimine whether the energy of dissociation is within the scope of the RACs.

Following this, there are a number of ways that the codes and regime may be improved.

The values for the energy of dissociation that we have determined are within feasible limits, however the low level of theory could be improved upon. This is unlikely to improve things drastically.

Better methods for hyperparameter optimisation could be used. The Kulik group have a highly sophisticated set of methods to approach this topic.

Certain pieces of code could be improved up. They have been marked as inefficient in the code.



I hope you find this code useful for the future of this project. 
Best of luck.
