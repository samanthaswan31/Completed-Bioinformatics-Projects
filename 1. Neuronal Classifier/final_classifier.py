#!/usr/bin/env python
# coding: utf-8

# In[1]:


#import libraries
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA


# In[3]:


#initialize 'final_score_overall' to store the final scores for all subjects
final_score_overall = []

#loop through each subject (20 subjects in total, subj_6 data incomplete)
for subj_num in range(1, 21):
    #handle file path formatting for subjects with single, double-digit IDs
    if subj_num < 10:
        f = h5py.File('Data/GTH_s0' + str(subj_num) + '_decision_power_struct_nobs.mat', 'r')
    else:
        f = h5py.File('Data/GTH_s' + str(subj_num) + '_decision_power_struct_nobs.mat', 'r')
    
    #import the dataset, convert to array
    dataset = np.array(f["power_struct"]["highgamma"]["powspctrm"])
    print(f"Import Shape - {dataset.shape}")  #print the shape of the imported dataset

    #transpose the dataset to rearrange dimensions
    dataset = np.transpose(dataset, (2, 0, 1))
    print(f"Transposed Shape - {dataset.shape}")  #print the new shape after transposition

    #import the labels (indicating gamble v not gamble)
    labelsset = np.array(f["power_struct"]["beh"]["gambles"])
    f.close()  # Close the file after reading data

    #initialize 'final_score_subject' to store scores for the current subject
    final_score_subject = []

    #unpack the dimensions of the dataset for trials, time, and channels
    trial, time, channel = dataset.shape
    print(f"Trial = {trial}, Time = {time}, channel = {channel}")  #print dataset dimensions

    #loop through different numbers of principal components (1 to 10)
    for n_component in range(1, 11):  
        #apply PCA with the current number of components
        pca_dataset = PCA(n_components=n_component)
        dataset_mid = np.reshape(dataset, (time * trial, channel))  #flatten dataset for PCA
        pca_initial = pca_dataset.fit_transform(dataset_mid)  #perform PCA
        pca_initial = np.reshape(pca_initial, (n_component, time, trial))  # reshape output

        #initialize empty lists for separating 'safe' and 'not safe' trials
        notsafe_presorting = list()
        safe_presorting = list()

        #sort trials based on labels into 'safe' and 'not safe'
        for postsorting in range(0, pca_initial.shape[2]):
            if labelsset[0, postsorting] == 0:  #label 0 == 'safe'
                safe_presorting.append(pca_initial[:, :, postsorting])

        for postsorting in range(0, pca_initial.shape[2]):
            if labelsset[0, postsorting] == 1:  # label 1 == 'not safe'
                notsafe_presorting.append(pca_initial[:, :, postsorting])

        #compute mean across sorted 'safe' and 'not safe' trials
        mean_notsafe = np.mean(notsafe_presorting, axis=0)
        mean_safe = np.mean(safe_presorting, axis=0)

        #transpose PCA dataset for comparison operations
        pc_dataset = np.transpose(pca_initial, (2, 0, 1))

        #initialize 'trials_postcomp' to store trial classifications
        trials_postcomp = []

        #classify each trial by comparing distances from "safe" and "not safe" means
        for trials in pc_dataset:
            if np.mean(sum(abs(mean_safe - trials))) > np.mean(sum(abs(mean_notsafe - trials))):
                trials_postcomp.append(1)  # classify as "not safe"
            else: 
                trials_postcomp.append(0)  # classify as "safe"

        #convert trial classifications to array
        comparison = np.array(trials_postcomp)

        #compute accuracy by comparing classifications to original labels
        accuracy_score = (labelsset == comparison)
        final_score_subject.append(((np.count_nonzero(accuracy_score)) / (len(accuracy_score[0]))) * 100)

    #append subject-specific scores to the overall scores list
    final_score_overall.append(final_score_subject)

#plot the final scores for analysis
for subj_num in range(0, 1):  # results for the first subject
    y1 = final_score_overall[subj_num]  #scores for the subject
    x1 = np.arange(1, 11, 1)  # X-axis values (number of components)
    plt.plot(x1, y1)  #generate plot

#add plot title and axis labels for clarity, show plot
plt.title("Final Accuracy Score Per Subject")
plt.xlabel("Principal Components")
plt.ylabel("Final Accuracy Score")
plt.show()  

