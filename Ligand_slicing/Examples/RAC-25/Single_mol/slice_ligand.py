#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 16:02:01 2022

@author: jacob
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import openbabel
import glob, os, sys
from rdkit.Chem import Draw

print('Finding file names..')                                       # Searches the directory ligand_indiv for the file paths
og_xyz_paths=sorted(glob.glob('ligand_indiv/*.xyz'))                # sort all indiv .xyz files into a list
print('found ' + str(len(og_xyz_paths)) + ' molecules to read')     # print num of files found
print('Sample: ' + og_xyz_paths[0])                                 # print first file path as example

if len(sys.argv) > 2:                                                                   # check commandline arguments are present
    print("\n**** ERROR ****\n\nToo many arguements passed to script.\n")               # ERROR print out
    print("Only pass the ligand .xyz file with the CSDs you want to slice\n")           # ERROR print out
    print("See /Examples directory for correct formatting.")                            # ERROR print out
    quit()                                                                              # Exit program
elif len(sys.argv) == 1:                                                                # Check CL
    print("\n**** ERROR ****\n\nPass ligand .xyz file to this script in command line")
    print("e.g.\npython Extract_CSDs.py ligand.xyz")
    print("See /Examples directory for correct formatting.")
    quit()
elif len(sys.argv) == 2:
    print("working with %s" % sys.argv[1])
    CSD_csv = sys.argv[1]

substruc_path = sys.argv[1]                                         # Set CL argument to variable

if not os.path.exists('./temp_folder'):                             # Setup folder for temporary files    
    os.makedirs('./temp_folder')

if not os.path.exists('./temp_folder/indiv_mols'):                  # Temp_folder individual .mol files
    os.makedirs('./temp_folder/indiv_mols')

if not os.path.exists('./faulty_strucs'):                           # Setup folder for faulty structures
    os.makedirs('./faulty_strucs')
  
if not os.path.exists('./faulty_strucs/failed_rdkit_descrip'):      # folder for structures that failed to be read by RDkit
    os.makedirs('./faulty_strucs/failed_rdkit_descrip')

if not os.path.exists('./faulty_strucs/multiple_substruc_removed'): # folder for structures with multiple ligands present
    os.makedirs('./faulty_strucs/multiple_substruc_removed')

if not os.path.exists('./faulty_strucs/non_matching_hs'):           # folder for structures where number of atoms do not match
    os.makedirs('./faulty_strucs/non_matching_hs')
  
if not os.path.exists('./faulty_strucs/no_match_substruc'):         # folder for structures where RDKit cannot identify ligand
    os.makedirs('./faulty_strucs/no_match_substruc')    

if not os.path.exists('./temp_folder/hold_folder'):
    os.makedirs('./temp_folder/hold_folder')
    
if not os.path.exists('./sliced_strucs_h'):                         # folder for sliced structures
    os.makedirs('./sliced_strucs_h')




with open(substruc_path, 'r') as f:                                 # read ligand substructure from .xyz
     substruc_num_atoms = int(f.readline())                         # get number of atoms in ligand to exract
f.close()

substruc_name = substruc_path.split('/')[-1].split(".")[0]          # get filename without suffix
substruc_mol = './' + substruc_name + '.mol'                        # make string for .mol path
substruc_simple = substruc_name + "_simple.mol"                     # make string for simplified substructure

obConversion = openbabel.OBConversion()                             # make openbabel conversion object
obConversion.SetInAndOutFormats("xyz", "mol")                       # make object convert xyz to mol
mol = openbabel.OBMol()                                             # make openbabel mol object

obConversion.ReadFile(mol, substruc_path)                           # read substructure .xyz file
obConversion.WriteFile(mol, substruc_mol)                           # write substructre as .mol file


substruc_mol1 = Chem.MolFromMolFile(substruc_mol)                   # RDKit read substructure from .mol
substruc_smarts = Chem.MolToSmarts(substruc_mol1)                   # RDKit convert Mol to Smarts
print(substruc_smarts)
substruc_conv = Chem.MolFromSmarts(substruc_smarts)                 # RDKit Mol from Smarts

Chem.MolToMolFile(substruc_conv, substruc_simple)                   # RDKit write .mol file of simplifed substructure


def make_mols(whole_xyz):                                           # CONVERT .xyz TO .mol
    
    mol_name = whole_xyz.split(".")[0].split("/")[1]                # get filename without suffix
    mol_path = "./temp_folder/indiv_mols/" + mol_name + ".mol"      # make .mol file path
    obConversion = openbabel.OBConversion()                         # make openbabel conversion object
    obConversion.SetInAndOutFormats("xyz", "mol")                   # conversion object xyz to mol
    mol = openbabel.OBMol()                                         # make openbabel mol object
    
    obConversion.ReadFile(mol, whole_xyz)                           # read .xyz as mol
    
    obConversion.WriteFile(mol, mol_path)                           # write .mol from .xyz
    


for struc in og_xyz_paths:                                          # make all indiv .xyz to .mol
    make_mols(struc)
    
    
print('Finding file names..')                                   # Searches the directory temp_folder/indiv_mols for the file paths
target_mols=sorted(glob.glob('temp_folder/indiv_mols/*.mol'))   # list all .mols in temp_folder/indiv_mols
print('found ' + str(len(target_mols)) + ' molecules to read')  # print num of .mols
print('Sample: ' + target_mols[0])                              # print first .mol as example                   

count = 0                                                       # index of each .mol
check_conv = 0                                                  # check all structures are attempted to be converted


failed_strucs = []                                              # instance list for failed misread structures
failed_strucs_h = []                                            # instance list for non-matching number of atoms
no_match = []                                                   # instance list for no ligand match structures

for struc in target_mols:                                       # For all paths in sorted list of mols
    mol = Chem.MolFromMolFile(struc)                            # make RDKit Mol object from indiv .mol
    print(mol)
    Draw.MolToFile(mol, 'mol_img.png', kekulize=False, wedgeBonds=True )
    mol_smart = Chem.MolToSmarts(mol)
    print(mol_smart)
    
    if not mol:                                                 # if making Mol fails:
        print(str(struc) + " has failed")                       # RDKit could not read this structure
        os.system("cp " + str(struc) + " ./faulty_strucs/failed_rdkit_descrip/")    # copy structure to faulty directory
        check_conv += 1                                         # +1 for check all structures are attempted
        failed_strucs.append(struc.strip('/')[1])               # add filename to list of failed
        
    else:
        if mol.HasSubstructMatch(substruc_conv):                # if structure has substructure ligand present
            
            sliced_struc = AllChem.DeleteSubstructs(mol, substruc_conv) # delete ligand from structure
            
            ''' THESE LINES PREPARE THE PATH EXTENSIONS FOR CORRECT ALLOCATION DEPENDING ON NATURE OF LIGAND SLICING'''
            
            sliced_h = './temp_folder/hold_folder/' + struc.split(".")[0].split("/")[2] + 'hold.mol'
            
            sliced_h_xyz = './temp_folder/hold_folder/' + struc.split(".")[0].split("/")[2] + 'h_slice.xyz'
            
            sliced_h_xyz_pass = './sliced_strucs_h/' + struc.split(".")[0].split("/")[2] + 'sliced_strucs_h.xyz'
            
            multi_substruc_remov = './faulty_strucs/multiple_substruc_removed/' \
                + struc.split(".")[0].split("/")[2] + '_multipleSubstruc.xyz'
            
            non_matching_hs = './faulty_strucs/non_matching_hs/' \
                + struc.split(".")[0].split("/")[2] + '_nonMatch.xyz'
            
            
            ''' WRITING AND RE-READING THE SLICED STRUCTURE '''
            ''' THIS CLEANS THE WRITING OF THE .mol FILE '''

            Chem.MolToMolFile(sliced_struc, sliced_h)       # Mol to .mol file
            mol1 = Chem.MolFromMolFile(sliced_h)            # Same .mol Re-converted to Mol
            Chem.MolToMolFile(mol1, sliced_h)               # Rewrting Mol to .mol    
            
            os.system("obabel %s -O %s -h" % (sliced_h, sliced_h_xyz)) # Use Openbabel to convert .mol to .xyz
            
            with open(og_xyz_paths[count]) as f:            # get full molecule number of atoms
                mol_num_atoms = int(f.readline())
            f.close()
            
            with open(sliced_h_xyz, 'r') as f:              # get sliced structure number of atoms
                sliced_num_atoms = int(f.readline())
            f.close()
            
            # Check number of atoms in ligand + sliced structure = whole molecule number of atoms
            multi_slice_check = ((mol_num_atoms - sliced_num_atoms) / substruc_num_atoms) 
            
            if multi_slice_check > 1:               # If this ^ check is False, check why.
               print("NUM ATOMS DO NOT MATCH")
               
               if multi_slice_check.is_integer():   # detect if multiple substructure ligands have been removed
                   print("MULTIPLE SUBSTRUC REMOVAL WARNING : " + struc.split(".")[0].split("/")[2])
                   os.system("obabel %s -O %s -h" % (sliced_h, multi_substruc_remov))
                   
                   os.system("cp " + target_mols[count] + ' ./faulty_strucs/multiple_substruc_removed/')
                   check_conv += 1
                   
               elif multi_slice_check.is_integer() == False:    # detect if Hydrogen is missing
                   print("SUSPECT NUM HYDROGENS IS INCORRECT : " + struc.split(".")[0].split("/")[2])
                   os.system("obabel %s -O %s -h" % (sliced_h, non_matching_hs))
                   
                   os.system("cp " + target_mols[count] + ' ./faulty_strucs/non_matching_hs/')
                   check_conv += 1
                   
            else:
                print("NUMBER OF ATOMS CONSERVED AND SINGLE SUBSTRUC EXTRACTED : " + struc.split(".")[0].split("/")[2])
                os.system("obabel %s -O %s -h" % (sliced_h, sliced_h_xyz_pass))
                check_conv += 1
                
        else:
            os.system("cp " + str(struc) + " ./faulty_strucs/no_match_substruc/")
            no_match.append(struc.strip('/')[1])
            check_conv += 1
    
    count +=1
    
print("\n\nNum files to convert = %s\nNum files converted = %s" % (str(len(target_mols)), str(check_conv)))


