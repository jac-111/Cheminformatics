# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 20:04:00 2022

@author: jacob
"""


import glob, os, sys

if not os.path.exists('orca_jobs'):
        os.makedirs('orca_jobs')

def mk_dirs(paths, dir_):
    
    if not os.path.exists(dir_):
        os.makedirs(dir_)
    dir_2 = dir_ + "/jobs"
    if not os.path.exists(dir_2):
        os.makedirs(dir_2)
    for struc in paths:
        print(struc)
        job_name = struc.split(".")[0].split("/")[1]
        if not os.path.exists(dir_2 + "/" + job_name):
           os.makedirs(dir_2 + "/" + job_name)
        
        orca_inp = open("./" + dir_2 + "/" + job_name + "/" + job_name + ".inp", 'w')
        
        orca_inp.write(orca_header1)
        raw_struc = open(struc, 'r')
        raw_struc_str = ''
        
        for x, line in enumerate(raw_struc):
            if x == 1:
                char = line.split('|')[1].split('=')[1]
            if x > 1:
                raw_struc_str += line
        
        orca_header2 =  ('"%s"\n\n*xyz %s 1\n' % (job_name, char))
        
        orca_inp.write(orca_header2)
        orca_inp.write(raw_struc_str + "*")
        raw_struc.close()
        orca_inp.close()
    
    bash_script = open(dir_ + "/exe_jobs.sh", 'w')
    bash_script.write("#!/usr/bin/env bash\necho -n './err_jobs.txt'\nfor file in ./jobs/*\ndo\n\tvar1=$file/*\n\tarrIN=(${file//// })\
                      \n\tfile_name=${arrIN[2]}\
                      \n\tcd $file\
                      \n\tif /opt/orca4/orca $file_name.inp > $file_name.out; then\
                      \n\t\techo $file\
                      \n\t\techo success\
                      \n\telse\
                      \n\t\techo $file\
                      \n\t\techo failed\
                      \n\t\techo $file > '../../err_jobs.txt'\
                      \n\tfi\
                      \n\tcd ../../\
                      \ndone\
                      \n\nmkdir err_jobs\
                      \nwhile read line; do\
                      \n\techo '$line'\
                      \n\tcp -r '$line' ./err_jobs/\
                      \ndone <err_jobs.txt")
                          
                      
print('Finding file names..')                                   # Searches the directory tmqm_indiv for the file paths
sliced_xyz_paths=sorted(glob.glob('sliced_strucs_h/*.xyz'))
print('found ' + str(len(sliced_xyz_paths)) + ' xyz to read to orca inp')
print('Sample: ' + sliced_xyz_paths[0])

print('Finding file names..')                                   # Searches the directory tmqm_indiv for the file paths
og_xyz_paths=sorted(glob.glob('ligand_indiv/*.xyz'))
print('found ' + str(len(og_xyz_paths)) + ' molecules to read')
print('Sample: ' + og_xyz_paths[0])

orca_header1 = '%pal nprocs 20 end\n!B97-3c opt slowconv \n%base '

dir_name = "orca_jobs/sliced_strucs_orca"

mk_dirs(sliced_xyz_paths, dir_name)

dir_name = "orca_jobs/og_strucs_orca"

mk_dirs(og_xyz_paths, dir_name)


