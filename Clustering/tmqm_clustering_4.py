# -*- coding: utf-8 -*-
"""
Created on Wed May 18 15:57:23 2022

@author: jacob
"""


import os, shutil, glob
import time
import pandas as pd
import csv
import string
import pickle

from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split

import matplotlib.pyplot as plt
import numpy as np
from random import randint



#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                   ~       SETUP DATA                 ~
#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



fraction = 0.1
bypass = 0
num_clusters = 20                                           # Number of Kmeans Clusters

folder = './temp_folder'

if not os.path.isdir(folder):           # Make individual xyz directory if doesn't exist
        os.makedirs(folder) 

for filename in os.listdir(folder):     # clear folder of previous run files
    file_path = os.path.join(folder, filename)
    try:
        if os.path.isfile(file_path) or os.path.islink(file_path):
            os.unlink(file_path)
        elif os.path.isdir(file_path):
            shutil.rmtree(file_path)
    except Exception as e:
        print('Failed to delete %s. Reason: %s' % (file_path, e))

if bypass == 0:
    
    tmqm_df = pd.read_csv("../Processed_data/tmqm_RACs-20.csv")
    
    tmqm_df = tmqm_df.sample(frac=fraction)
    
    tmqm_df.to_csv("./tmqm_cluster_sample.csv", header=True, index=False)

elif bypass == 1:
    tmqm_df = pd.read_csv('./tmqm_cluster_sample.csv')
    

tmqm_racs = tmqm_df.iloc[:,10:]

scaler = StandardScaler()
scaler.fit(tmqm_racs)
tmqm_racs_scaled = pd.DataFrame(scaler.transform(tmqm_racs))

tmqm_index = tmqm_df.iloc[:,0]


#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                   ~       PCA Analysis of RACs        ~
#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pca = PCA(n_components=2)
pca.fit(tmqm_racs_scaled)

exp_vari = pca.explained_variance_ratio_
print("Retained variance: cumulativeSum()")
print(exp_vari.cumsum())
pca = pca.fit_transform(tmqm_racs_scaled)

#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                   ~       Kmeans Analysis           ~
#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



kmeans = KMeans(init="k-means++", n_clusters=num_clusters, n_init=4) 


tmqm_cluster_index = (kmeans.fit_predict(pca))              # list containing which cluster each index belongs to

cluster_lists = [[] for x in range(num_clusters)]               # list of lists of 'in program' index clustering
                                                                # i.e. each list is a cluster with 'in program' index
                                                                # of each TM in it. - NOT CSD INDEX.
                                                                
cluster_csd_index = [[] for x in range(num_clusters)]           # list of lists of csd index in clusters
                                                                # i.e each list is a cluster with CSD associated
                                                                # number.

for x in range(num_clusters):
    for y in range(len(tmqm_cluster_index)):
        if tmqm_cluster_index[y] == x:
            cluster_lists[x].append(y)                          # fill cluster list with in program index
            cluster_csd_index[x].append(tmqm_index[y] - 1)      # fill cluster list with CSD numbers
            
print("CLUSTER_LISTS")
print(cluster_lists)
print(cluster_csd_index)
print("Cluster centers")
print(kmeans.cluster_centers_)


#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                   ~ Def functions and periodic table  ~
#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

main_elms = ['B', 'C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl', 'Se', 'Br', 'I']
main_track = [[0] * len(main_elms) for x in range(num_clusters)]

tm_elms = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', \
           'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', \
           'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', \
           'Ac', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn']

tm_rows = [['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'],
         ['Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd'],
         ['La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg'],
         ['Ac', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn']]
    
tm_track = [[0] * len(tm_elms) for x in range(num_clusters)]

tm_store = [[] for x in range(len(tm_elms))]                    # list to store in program index
                                                                # associated with each element


cluster_sizes = [[] for x in range(len(cluster_lists))]         # list of lengths of clusters

for x in range(len(cluster_sizes)):                             
    cluster_sizes[x] = [len(cluster_lists[x])]


with open("cluster_occup.csv", "w", newline="") as f:           # clusters with indexes of molecules
    writer = csv.writer(f)
    writer.writerows(cluster_sizes)

with open("cluster_indexes.csv", "w", newline="") as f:         # length of lindexes (num mols)
    writer = csv.writer(f)
    writer.writerows(cluster_lists)

centroids = kmeans.cluster_centers_

with open("Kmean_Centroids.csv", "w", newline="") as f:         # position of centroids
    writer = csv.writer(f)
    writer.writerows(centroids)



#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                   ~ Main Parse loop for TM identification  ~
#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start_time = time.time()
print("--- %s seconds ---" % (time.time() - start_time))

if bypass == 0:   
    tick = 0
    
    for x in range(len(cluster_csd_index)):
        for y in range(len(cluster_csd_index[x])):
            tick += 1
            xyz_num = str(cluster_csd_index[x][y]).zfill(6)                         # make CSD number in print format
            
            xyz_path = glob.glob('../Processed_data/tmqm_indiv/%s*.xyz' % xyz_num)  # get path of xyz_num to xyz file
            
            mol_data = open(xyz_path[0], 'r')           # open .xyz file
            for z, line in enumerate(mol_data): 
                if z > 1:
                    elem = line.split()[0]              # Get Element string
                    
                    if elem in tm_elms:
                        
                        tm_posit = tm_elms.index(elem)  # Get index of TM Element in list of TM elements (tm_elms)
                        tm_track[x][tm_posit] += 1      # Track occurance of TM in each cluster
    
                        tm_store[tm_posit].append(cluster_lists[x][y])  # Assign list for each TM Element the in program index
                        
                        break
                    
                    
            mol_data.close()     
            if tick % 100 == 0 :                # Simple check section
                print("--- %s seconds ---" % (time.time() - start_time))
                print(str(tick), "/", str(len(tmqm_cluster_index)))
                
    
    with open("tm_track.pic", "wb") as fp:   #Pickling
        pickle.dump(tm_track, fp)
    
    with open("tm_store.pic", "wb") as fp:   #Pickling
        pickle.dump(tm_store, fp)


#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                   ~ Main Parse ByPass with Pickle  ~
#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

if bypass == 1:    
    with open("tm_track.pic", "rb") as fp:   # unPickling
        tm_track = pickle.load(fp)
    
    with open("tm_store.pic", "rb") as fp:   # unPickling
        tm_store = pickle.load(fp)



#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                   ~ MatPlotLib section for plotting voro  ~
#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


h = 0.02  # point in the mesh [x_min, x_max]x[y_min, y_max].

# Plot the decision boundary. For that, we will assign a color to each
x_min, x_max = pca[:, 0].min() - 1, pca[:, 0].max() + 1
y_min, y_max = pca[:, 1].min() - 1, pca[:, 1].max() + 1
xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))

# Obtain labels for each point in mesh. Use last trained model.
Z = kmeans.predict(np.c_[xx.ravel(), yy.ravel()])

# Put the result into a color plot
Z = Z.reshape(xx.shape)
plt.figure(1)
plt.clf()
plt.imshow(
    Z,
    interpolation="nearest",
    extent=(xx.min(), xx.max(), yy.min(), yy.max()),
    cmap=plt.cm.Paired,
    aspect="auto",
    origin="lower",
)


#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                   ~ Parsing Cluster lists to collect
#                       grouping data i.e. rows etc     ~
#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tm_dict = dict.fromkeys(tm_elms)            # Make dict with TM Elements for keys
tm_row_dict = {0: [],                       # Dictionary with 4 rows for rows in TM block
               1: [],
               2: [],
               3: []}

tm_row_cluster_index = {0: [],                       # Dictionary with 4 rows for rows in TM block
               1: [],
               2: [],
               3: []}

tm_groups_dict = {0: [],                    # Dict for groups
               1: [],
               2: [],
               3: [],
               4: [],
               5: [],
               6: [],
               7: [],
               8: [],
               9: [],
               10: []}
               


for x in tm_dict:           # make each Key have a list value
    tm_dict[x] = []

for x in range(len(tm_store)):

    for y in range(len(tm_store[x])):                                           # add PCA coordinates to each Element key in dictionary
        tm_dict[tm_elms[x]].append((pca[tm_store[x][y]].tolist()))
        
        for z in range(len(tm_rows)):                                           # Assign rows PCA coordinates of associated Elements
            if tm_elms[x] in tm_rows[z]:
                tm_row_dict[z].append((pca[tm_store[x][y]].tolist()))

        for v in range(len(tm_rows)):                                           # Assign groups PCA coordinates of associated Elements
            for b in range(len(tm_rows[0])):
                if tm_elms[x] == tm_rows[v][b]:
                    tm_groups_dict[b].append((pca[tm_store[x][y]].tolist()))


for x in range(len(tm_store)):
    for y in range(len(tm_store[x])):
        for z in range(len(cluster_lists)):
            if tm_store[x][y] in cluster_lists[z]:
                tm_row_cluster_index[np.floor(x/10)].append(z)                  # Assign [[row1], [row2] ...] to have cluster index for
                                                                                # all occupants in each row




#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                   ~ Converting dictionary of lists to
#                         dictionary of np.arrays  ~
#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for x in tm_dict:                       # Convert tm_dict lists into np.array
    tm_dict[x] = np.array(tm_dict[x])

tm_dict_occ_keys = []                   # list for parsing occupied keys:
                                        # required for graphing later
for x in tm_dict.keys():
    if len(tm_dict[x]) > 0:
        tm_dict_occ_keys.append(x)      # identify and collect occupied keys


for x in tm_row_dict:
    tm_row_dict[x] = np.array(tm_row_dict[x])
tm_row_dict_occ_keys = []
for x in tm_row_dict.keys():
    if len(tm_row_dict[x]) > 0:
        tm_row_dict_occ_keys.append(x)

for x in tm_groups_dict:
    tm_groups_dict[x] = np.array(tm_groups_dict[x])
tm_groups_dict_occ_keys = []
for x in tm_groups_dict.keys():
    if len(tm_groups_dict[x]) > 0:
        tm_groups_dict_occ_keys.append(x)


#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                   ~ Colour coding the PCA coordinates  ~
#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

colors = []
for i in range(len(tm_elms)):
    colors.append('#%06X' % randint(0, 0xFFFFFF))


for data_list in tm_row_dict.values():
    
    if len(data_list) > 0:
        x = data_list[:,0]
        y = data_list[:,1]
        plt.scatter(x,y,s=25)#,color=colors.pop())


plt.legend(tm_row_dict_occ_keys,bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, ncol=2)




#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                   ~ Adding the Centroid crosses to graph  ~
#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# coordinates for centroids of groupings
centroids = kmeans.cluster_centers_

tick2 = 0
for x, y in zip(centroids[:, 0], centroids[:, 1],):
    plt.text(x, y, tick2, color="black", fontsize=12)
    tick2 += 1

plt.title(
    "K-means clustering on the tmqm RAC dep 5 data (PCA-reduced data)\n"
    "Centroids are marked with white cross"
)
plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
plt.xticks(())
plt.yticks(())
plt.savefig('MY_FIGURE_VORO', bbox_inches='tight')


#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                   ~ Creating Stacked bar chart  ~
#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sum1 = 0
for x in range(len(tm_row_dict)):
    sum1 += np.sum(tm_row_dict[x])
    print(len(tm_row_dict[x]))
print(sum1)
fig, ax = plt.subplots()

labels = list(range(num_clusters))
row_cluster_freq = [[] for x in range(len(tm_rows))]

for x in range(len(tm_row_cluster_index)):
    for y in range(num_clusters):
        row_cluster_freq[x].append(tm_row_cluster_index[x].count(y))


sum2 = 0
for x in range(len(row_cluster_freq)):
    sum2 += np.sum(row_cluster_freq[x])

ax.bar(labels, row_cluster_freq[0], label="Row_0")
ax.bar(labels, row_cluster_freq[1], label="Row_1", bottom=row_cluster_freq[0])
bottom1 = np.array(row_cluster_freq[0]) + np.array(row_cluster_freq[1])
ax.bar(labels, row_cluster_freq[2], label="Row_2", bottom=bottom1)

ax.set_ylabel('frequency')
ax.set_title('frequency of each row appearance within each Kmeans cluster')
ax.legend()

plt.xticks(labels)

plt.savefig('MY_FIGURE_BAR', bbox_inches='tight')


#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                   ~ Classification ML Algo  ~
#                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cluster_row_lis = []
for x in range(len(tm_row_cluster_index)):
    for y in range(len(tm_row_cluster_index[x])):
        temp_clus_row = [tm_row_cluster_index[x][y], x]
        cluster_row_lis.append(temp_clus_row)

cluster_row_arr = np.array(cluster_row_lis)

X_train, X_test, y_train, y_test = train_test_split(cluster_row_arr[:,0].reshape(-1, 1), \
                                                    cluster_row_arr[:,1].reshape(-1, 1), random_state=0)

print(X_train, X_test, y_train, y_test)

sc = StandardScaler()
sc.fit(X_train)
X_train_std = sc.transform(X_train)
X_test_std  = sc.transform(X_test)

svm = SVC(kernel='rbf', random_state=0, gamma=0.7, C=1.0)
svm.fit(X_train_std, y_train)
print("test set accuracy: {:.2f}".format(svm.score(X_test_std, y_test)))


'''


