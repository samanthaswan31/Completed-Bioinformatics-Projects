# Neuronal Classifier Using iEEG Data to Predict Choice in a Gambling Task

The project processes a dataset consisting of high-gamma power spectral data for 20 subjects. Each subject's data is analyzed to reduce the dimensionality of the feature space using PCA, followed by trial classification based on Euclidean Distance to the mean of 'safe' and 'not safe' trials. The accuracy of the classification is computed by comparing the predicted classifications with the true labels for each trial. The accuracy scores for different numbers of principal components are plotted to assess the model's performance.

## Code Description

Data Loading: The script begins by iterating through each subject (1 to 20), loading their corresponding dataset from .mat files stored in the Data/ directory. For each subject, the high-gamma power spectral data (powspctrm) and gamble labels (gambles) are extracted and converted into numpy arrays.

Preprocessing: The dataset is transposed to rearrange the dimensions, ensuring the data is in the correct format for PCA. The labels, indicating whether a trial is "safe" (0) or "not safe" (1), are also extracted from the dataset for classification purposes.

PCA Transformation: PCA is applied to the transposed dataset for each subject, varying the number of components from 1 to 10. The dataset is reshaped to a 2D format (time x trials) to apply PCA. The transformed dataset is then reshaped back to a 3D structure (components x time x trials) for further processing.

Classification: For each PCA-transformed dataset, the trials are classified as either "safe" or "not safe" based on their Euclidean distance to the mean of the 'safe' and 'not safe' trials in PCA space. Specifically, each trial is compared to the means of both classes, and the trial is classified based on the smaller Euclidean distance. The classification is compared with the original labels to compute the accuracy.

Accuracy Calculation: The accuracy for each subject is calculated by comparing the predicted trial classifications with the true labels. This process is repeated for each number of PCA components, from 1 to 10, to evaluate the impact of dimensionality reduction on classification performance.

Plotting Results: The accuracy scores for each subject are plotted as a function of the number of principal components. This allows for a visual analysis of how the number of PCA components affects classification accuracy.


This neuronal classifier was developed in collaboration with the (Moxon Lab)[https://moxonlab.bme.ucdavis.edu/] at the University of California Davis. Thank you to Dr. Karen Moxon and Logan Peters for their support with this project.
