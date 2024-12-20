#!/usr/bin/env python
# coding: utf-8

# **Design and Simulation of a Trial to Evaluate Vaccination Efficacy in a Model Organism**
# 
# To test vaccination efficacy, it is useful to design trials using model organisms to test impacts. This program simulates one such trial using entirely undefined values for probabilities. The probability of becoming infected following an exposure is represented by the variable x for untreated individuals of the model organism. The developer of the vaccination believes that this probability is reduced to y following vaccination. The control and treatment groups have the same number of individuals, represented by N, and the statistical significance is evaluated using a permutation test. This program would allow the user to input the x and y probabilities and simulate trials to find the N to ensure approximately 80% power to detect a significant reduction in infection probability due to the vaccination at the α = 0.05 level.

# In[535]:


#import packages
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()
import scipy.stats as stats

#create function that simulates outcomes of control, treatment
def create_samples(outcomes, control_p, treated_p, repeats, N):
    #repeats = number of simulations, sample_size = number of samples per simulation

    #create samples
    control_sample = np.random.choice(outcomes, size=(repeats, N), p=control_p)
    treated_sample = np.random.choice(outcomes, size=(repeats, N), p=treated_p)

    return control_sample, treated_sample


# In[537]:


control_sample, treated_sample = create_samples(outcomes = [0, 1], control_p = [x, 1-x], treated_p = [y, 1-y], repeats = 1000, N = 22)
control_sample, treated_sample
print("Treated sample mean:", treated_sample.mean())
print("Control sample mean:", control_sample.mean())


# In[539]:


def fisher_permutation(control_sample, treated_sample, split_value, N):
    #combine datasets, create copy for calculations
    combined_data = np.concatenate((treated_sample, control_sample), axis=1)  # Shape: (repeats, N)
    
    #flatten the combined data to ensure shuffling
    combined_data = combined_data.flatten()
    difference = np.zeros(split_value)
    
    for i in range(split_value):
        #shuffle the entire combined dataset
        np.random.shuffle(combined_data)
        
        #simulate treated and control groups for this iteration
        sim_treated = combined_data[:N]  #first 'N' values as treated group
        sim_control = combined_data[N:]  #remainder 'N' as control group
        
        #calculate the difference between simulated means
        difference[i] = sim_treated.mean() - sim_control.mean()

    return difference


# In[541]:


difference = fisher_permutation(control_sample, treated_sample, split_value = 1000, N = 22)


# In[542]:


#use right-sided test to determine p_value

#calculate observed difference
observed_difference = treated_sample.mean() - control_sample.mean()
#check that the relative sample means are the same as before to ensure observed_difference calculation correct
print("Treated sample mean", treated_sample.mean())
print("Control sample mean", control_sample.mean())

#calculate p-value
p_value = np.sum(difference >= observed_difference)/difference.size

#plot p-value
fig = plt.hist(difference, bins=100, align='left')
plt.xlabel('differences')
plt.ylabel('PDF')
plt.axvline(x=observed_difference, color='red');


# In[557]:


#find the best 'N' for 80% threshold that the test passes the statistical significance threshold of alpha = 0.05

#calculate relative power
def calculate_power(N, repeats, split_value, alpha=0.05):
    
    #simulate the control and treatment groups using create_samples
    control_sample, treated_sample = create_samples(outcomes=[0,1], control_p=[x, 1-x], treated_p=[y, 1-y], repeats=repeats, N=N)
    
    #perform permutation test using fisher_permutation
    difference = fisher_permutation(control_sample, treated_sample, split_value=split_value, N=N)
    
    #calculate observed difference
    observed_difference = treated_sample.mean() - control_sample.mean()
    
    #calculate p-value
    p_value = np.sum(difference >= observed_difference) / difference.size
    
    #return if p-value is less than alpha (0.05)
    return p_value < alpha

#calculate power for multiple values of 'N' to choose best 'N'
def power_analysis(N_values, repeats=1000, split_value=1000, alpha=0.05, threshold=0.8):

    #initialize power to store proportions of significant tests
    power = []
    #initialize threshold_N to store when power crosses threshold (80%)
    threshold_N = None 

    #
    for N in N_values:
        #initialize counter to store if p-value < alpha
        significant_count = 0
        
        for _ in range(repeats):
            #use returned equality evals to count significant iterations, count with significant_count
            if calculate_power(N, repeats=1, split_value=split_value, alpha=alpha):
                significant_count += 1

        #calculate proportion of iterations where test rejects null hypothesis (is significant)
        power_value = significant_count / repeats

        #append to 'power' list
        power.append(power_value)
        
        #use threshold_N to find first 'N' past threshold (80%)
        if power_value >= threshold and threshold_N is None:
            threshold_N = N  #capture the first N where power crosses the threshold
    
    return power, threshold_N

#define the range of N values (initial range determined by trying potential 'N' values above)
N_values = range(19, 32, 1)

#calculate the power for each value of N, get the threshold_N
power, threshold_N = power_analysis(N_values, repeats=100, split_value=1000, alpha=0.05)

#plot the power curve
plt.figure(figsize=(8, 6))
plt.plot(N_values, power, marker='o')
plt.xlabel('Sample of Model Organism (N)')
plt.axhline(y=0.8, color='red', linestyle='--', label="Power = 80%")
plt.ylabel('Probability of Significant Result (Power)')
plt.title('Analysis for Vaccination Trial')
plt.grid(True)
plt.show()

#print the power values for each N
for N, p in zip(N_values, power):
    print(f"N={N}, Power={p:.3f}")

# Print the N where power crosses 80%
if threshold_N is not None:
    print(f"The first N where power exceeds 80% is: N={threshold_N}")
else:
    print("No N value exceeds the 80% power threshold.")

Having run my simulation several times, I would choose an initial N value around 24 pigs.