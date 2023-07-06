# -*- coding: utf-8 -*-
"""
Created on Wed May 25 13:48:56 2022

This is a small utility program to split the processed tmqm data into smaller
datasets. I did this to work with smaller subset for ML work.


@author: jacob
"""


import pandas as pd


tmqm_df = pd.read_csv("../Processed_data/tmqm_RACs.csv")
print(tmqm_df)
'''
part_50 = tmqm_df.sample(frac = 0.5)
part_50.to_csv("../Processed_data/tmqm_RACs_part50.csv", index=False)
part_25 = tmqm_df.sample(frac = 0.25)
part_25.to_csv("../Processed_data/tmqm_RACs_part25.csv", index=False)
part_10 = tmqm_df.sample(frac = 0.1)
part_10.to_csv("../Processed_data/tmqm_RACs_part10.csv", index=False)
part_5 = tmqm_df.sample(frac = 0.05)
part_5.to_csv("../Processed_data/tmqm_RACs_part5.csv", index=False)
part_1 = tmqm_df.sample(frac = 0.01)
part_1.to_csv("../Processed_data/tmqm_RACs_part1.csv", index=False)
part_05 = tmqm_df.sample(frac = 0.005)
part_05.to_csv("../Processed_data/tmqm_RACs_part5.csv", index=False)
'''
part_01 = tmqm_df.sample(frac = 0.001)
part_01.to_csv("../Processed_data/tmqm_RACs_part01.csv", index=False)


    