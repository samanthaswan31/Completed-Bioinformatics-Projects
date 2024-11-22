#!/usr/bin/env python
# coding: utf-8

# **Developing Custom Functions for Mean and Covariance Matrix Calculation and Visualization Using the Iris Dataset**
# 
# Let X be an N × M data matrix with M variables (features) and N experiments (samples, measurements).
# 
# Write the following two functions:
# 
# • mean: The mean function should return the average values of the data matrix variables. The average is taken over different experiments (samples).
# 
# • cov: The cov method that calculates a covariance matrix C for X. It assumes that the different variables (features) of the data will be aligned along columns and that different experiments will be aligned along rows. So each row of X represents a different experiment, consisting of columns corresponding to the different features of a single experiment.

# In[6]:


import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sps
import pandas as pd
import seaborn as sns
import sklearn as skl
sns.set()


# In[8]:


#create function mean() 
def mean(input_matrix):

    #print purpose of the function
    print("This function returns the average values of the data matrix variables.")

    #check if the first element of input_matrix is a list, if not print statement, return input_matrix
    if not isinstance(input_matrix[0],list):
        print("This matrix is a single experiment")
        return input_matrix

    else:
        #format input_matrix so if there are missing measurements it can still calculate mean
        #find the longest row
        longest_row = max(len(row) for row in input_matrix)

        #append 'None' value so that all experiments have same number of measurements
        for row in input_matrix:
            while len(row) < longest_row:
                row.append(None)
        
        #format matrix for finding the mean for each list (feature/variable)
        number_features = len(input_matrix[0])
        calc_mean = [0 for _ in range(len(input_matrix[0]))]
        format_matrix = [[row[i] for row in input_matrix] for i in range(number_features)]

        #after formatting, now remove None for calculation
        filter_format = [[x for x in row if x is not None] for row in format_matrix]
        
        #iterate through each variable (row) in input_matrix
        for i, variable in enumerate(filter_format):
            #calculate the variable mean
            variable_mean = sum(variable) / len(variable)
            #add it to the formatted matrix using indexing for correct placement
            calc_mean[i] = variable_mean
            
        #return the calculated means
    return calc_mean


# In[10]:


mean()


# In[12]:


#create function cov to find the covariance matrix
def cov(input_matrix):
    def mean(matrix):

    #check if the first element of input_matrix is a list, if not print statement, return input_matrix
        if not isinstance(matrix[0],list):
            print("This matrix is a single experiment")
            return input_matrix

        else:
            #format input_matrix so if there are missing measurements it can still calculate mean
            #find the longest row
            longest_row = max(len(row) for row in matrix)

            #append 'None' value so that all experiments have same number of measurements
            for row in matrix:
                while len(row) < longest_row:
                    row.append(None)
        
            #format matrix for finding the mean for each list (feature/variable)
            number_features = len(matrix[0])
            calc_mean = [0 for _ in range(len(matrix[0]))]
            format_matrix = [[row[i] for row in matrix] for i in range(number_features)]

            #after formatting, now remove None for calculation
            filter_format = [[x for x in row if x is not None] for row in format_matrix]
        
            #iterate through each variable (row) in input_matrix
            for i, variable in enumerate(filter_format):
                #calculate the variable mean
                variable_mean = sum(variable) / len(variable)
                #add it to the formatted matrix using indexing for correct placement
                calc_mean[i] = variable_mean
            
            #return the calculated means
        return calc_mean
        
    #print purpose of the function
    print("This function returns the covariance matrix C for an inputted matrix.")
    
    #check if matrix is list, if not convert to list format
    if not isinstance(input_matrix,list):
        input_matrix = input_matrix.tolist()

    #standardize input_matrix list values so that all lengths are the same
    #find the longest row
    long_row = max(len(row) for row in input_matrix)

            #append 'None' value so that all experiments have same number of measurements
    for row in input_matrix:
        while len(row) < long_row:
            row.append(None)
    
    #create empty formatted matrix to store covariable values
    num_features = len(input_matrix[0])
    cov_matrix = [[0 for _ in range(num_features)] for _ in range(num_features)]

    #calculate N as the number of experiments (or rows)
    N = float(len(input_matrix)) #number of rows

    #create variable_means to store the means from each column
    feature_means = mean(input_matrix)

    #calculate covariance
    for i in range(len(input_matrix[0])): #for row/experiment
            for j in range(len(input_matrix[0])): #for row/experiment
                covariance_sum = 0 #rest covariance sum once done looping through row
                valid_entry = 0 #checks to make sure a valid covariance sum was computed
                for k in range(len(input_matrix)): #for column/variables
                    if input_matrix[k][i] is not None and input_matrix[k][j] is not None:
                    #want to handle unspecified number of features so input_matrix indexing compares column by column
                    #use feature_means calculated in means function, use i and j for each column
                        covariance_sum += (input_matrix[k][i] - feature_means[i]) * (input_matrix[k][j] - feature_means[j])
                        cov_matrix[i][j] = covariance_sum / len(input_matrix)
                    else:
                        cov_matrix[i][j] = None
                
    #return calculated cov_matrix  
    return cov_matrix


# In[14]:


cov()


# • Load the iris data from scikit-learn using load iris (first install sklearn (https://scikit-learn.org)
# if you haven’t yet, using pip).
# 
# • Without using pandas or any built-in functions, calculate the covariance matrix for the iris data using
# your cov function from problem 1.
# 
# • Visualize, show the heatmap of the iris data covariance matrix using imshow. Add a colorbar, and
# label the matrix axis ticks with the feature names from the data.
# 
# • Using pandas, create a DataFrame for the iris data. Print (show) it in the notebook.
# 
# • Calculate the covariance matrix using pandas (use ddof=0). Compare to the result of your function
# cov, are they the same?
# 
# • Show the heatmap of the covariance matrix using seaborn. Is it the same as the one you created with
# imshow? (You do not need to match the colormaps)

# In[16]:


#load iris data from scikit-learn
from sklearn.datasets import load_iris
iris_entire = load_iris()
#can subset so have just data without labels
iris_data = iris_entire['data']
print(iris_entire)


# In[18]:


#calculate covariance matrix for iris data using cov function, converts to list format within the cov function
cov_iris = cov(iris_data)


# In[20]:


#create a DataFrame for iris data, print
iris_df = pd.DataFrame(iris.data, columns=iris.feature_names)
print(iris_df)


# In[ ]:


#calculate covariance matrix with pandas, ddof = 0
iris_covmatrix = iris_df.cov(ddof=0)
print(iris_covmatrix)


# Comparing the results, they are the same. This is especially easy to see when I print my own result row by row and get an exact overlaid match numerically.

# In[1332]:


#my calculated C matrix
for row in cov_iris:
    print(row)

#pandas calculated C matrix
print(iris_covmatrix)


# In[1180]:


#show heatmap of iris data covariance matrix using imshow

#want cov_iris generated by my function cov
plt.figure(figsize=(8, 8))
plot = plt.imshow(cov_iris, cmap='YlOrBr')
#colorbar corolates to relative covariance between features
plt.colorbar(label='Relative Covariance Between Features')

#set up tick marks to correspond to features
plt.xticks(ticks=np.arange(len(iris.feature_names)), labels=iris.feature_names, rotation=45)
plt.yticks(ticks=np.arange(len(iris.feature_names)), labels=iris.feature_names)

#label axis and title heatmap
plt.title('Covariance Matrix Heatmap for Iris Dataset')
plt.xlabel('Features')
plt.ylabel('Features')

plt.show()


# In[1324]:


#show heatmap of iris data covariance matrix using seaborn
sns.heatmap(cov_iris, cmap='YlOrBr', annot=True, fmt=".2f", cbar_kws={'label': 'Relative Covariance Between Features'})

#label axis and title heatmap
plt.title('Covariance Matrix Heatmap for Iris Dataset')
#use same tickmark scheme for clarity
plt.xticks(ticks=np.arange(len(iris.feature_names)), labels=iris.feature_names, rotation=45)
plt.yticks(ticks=np.arange(len(iris.feature_names)), labels=iris.feature_names, rotation=0)


# I would say the two graphs are the same. Given that the colormaps match, the overlaid heatmaps are nearly identical except for minor formatting details.
