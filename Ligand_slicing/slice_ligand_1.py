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


print('Finding file names..')                                   # Searches the directory tmqm_indiv for the file paths
og_xyz_paths=sorted(glob.glob('ligand_indiv/*.xyz'))
print('found ' + str(len(og_xyz_paths)) + ' molecules to read')
print('Sample: ' + og_xyz_paths[0])

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


if not os.path.exists('./temp_folder'):
    os.makedirs('./temp_folder')

if not os.path.exists('./temp_folder/indiv_mols'):
    os.makedirs('./temp_folder/indiv_mols')

if not os.path.exists('./faulty_strucs'):
    os.makedirs('./faulty_strucs')
  
if not os.path.exists('./faulty_strucs/failed_rdkit_descrip'):
    os.makedirs('./faulty_strucs/failed_rdkit_descrip')

if not os.path.exists('./faulty_strucs/multiple_substruc_removed'):
    os.makedirs('./faulty_strucs/multiple_substruc_removed')

if not os.path.exists('./faulty_strucs/non_matching_hs'):
    os.makedirs('./faulty_strucs/non_matching_hs')
  
if not os.path.exists('./faulty_strucs/no_match_substruc'):
    os.makedirs('./faulty_strucs/no_match_substruc')    

if not os.path.exists('./temp_folder/hold_folder'):
    os.makedirs('./temp_folder/hold_folder')
    
if not os.path.exists('./sliced_strucs_h'):
    os.makedirs('./sliced_strucs_h')




with open(substruc_path, 'r') as f:
     substruc_num_atoms = int(f.readline())
f.close()

substruc_name = substruc_path.split(".")[0]
substruc_mol = './' + substruc_name + '.mol'


obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats("xyz", "mol")
mol = openbabel.OBMol()

obConversion.ReadFile(mol, substruc_path) 
obConversion.WriteFile(mol, substruc_mol)


substruc_mol1 = Chem.MolFromMolFile(substruc_mol)
substruc_smarts = Chem.MolToSmarts(substruc_mol1)
substruc_conv = Chem.MolFromSmarts(substruc_smarts)


def make_mols(whole_xyz):
    
    mol_name = whole_xyz.split(".")[0].split("/")[1]
    mol_path = "./temp_folder/indiv_mols/" + mol_name + ".mol"
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("xyz", "mol")
    mol = openbabel.OBMol()
    
    obConversion.ReadFile(mol, whole_xyz) 
    
    
    obConversion.WriteFile(mol, mol_path)
    


for struc in og_xyz_paths:
    make_mols(struc)
    
    
print('Finding file names..')                                   # Searches the directory tmqm_indiv for the file paths
target_mols=sorted(glob.glob('temp_folder/indiv_mols/*.mol'))
print('found ' + str(len(target_mols)) + ' molecules to read')
print('Sample: ' + target_mols[0])

count = 0
check_conv = 0


failed_strucs = []
failed_strucs_h = []
no_match = []

for struc in target_mols:
    mol = Chem.MolFromMolFile(struc)
    
    if not mol:
        print(str(struc) + " has failed")
        os.system("cp " + str(struc) + " ./faulty_strucs/failed_rdkit_descrip/")
        check_conv += 1
        failed_strucs.append(struc.strip('/')[1])
        
    else:
        if mol.HasSubstructMatch(substruc_conv):
            
            sliced_struc = AllChem.DeleteSubstructs(mol, substruc_conv)
            print(struc)
            sliced_h = './temp_folder/hold_folder/' + struc.split(".")[0].split("/")[2] + 'hold.mol'
            
            sliced_h_xyz = './temp_folder/hold_folder/' + struc.split(".")[0].split("/")[2] + 'h_slice.xyz'
            
            sliced_h_xyz_pass = './sliced_strucs_h/' + struc.split(".")[0].split("/")[2] + 'sliced_strucs_h.xyz'
            
            
            multi_substruc_remov = './faulty_strucs/multiple_substruc_removed/' \
                + struc.split(".")[0].split("/")[2] + '_multipleSubstruc.xyz'
            
            non_matching_hs = './faulty_strucs/non_matching_hs/' \
                + struc.split(".")[0].split("/")[2] + '_nonMatch.xyz'
            
            
            Chem.MolToMolFile(sliced_struc, sliced_h)
            
            mol1 = Chem.MolFromMolFile(sliced_h)
            Chem.MolToMolFile(mol1, sliced_h)
            
            os.system("obabel %s -O %s -h" % (sliced_h, sliced_h_xyz))
            
            with open(og_xyz_paths[count]) as f:
                mol_num_atoms = int(f.readline())
            f.close()
            
            with open(sliced_h_xyz, 'r') as f:
                sliced_num_atoms = int(f.readline())
            f.close()
            
            multi_slice_check = ((mol_num_atoms - sliced_num_atoms) / substruc_num_atoms)
            
            if multi_slice_check > 1:
               print("NUM ATOMS DO NOT MATCH")
               
               if multi_slice_check.is_integer():
                   print("MULTIPLE SUBSTRUC REMOVAL WARNING : " + struc.split(".")[0].split("/")[2])
                   os.system("obabel %s -O %s -h" % (sliced_h, multi_substruc_remov))
                   
                   os.system("cp " + target_mols[count] + ' ./faulty_strucs/multiple_substruc_removed/')
                   check_conv += 1
               elif multi_slice_check.is_integer() == False:
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


