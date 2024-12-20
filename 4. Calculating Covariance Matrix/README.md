# Developing Custom Functions for Mean and Covariance Matrix Calculation and Visualization Using the Iris Dataset

This project involves developing custom Python functions for calculating the mean and covariance matrix of a dataset, and visualizing the results with the Iris dataset from scikit-learn. These functions are built from scratch, avoiding reliance on built-in functions for computation, to enhance understanding of the underlying mathematical concepts.

After implementing these functions, the project compares the covariance matrix derived from the custom cov function with the one computed using pandas, validating their equivalence. Additionally, the covariance matrix is visualized through two heatmaps: one created using Matplotlib's imshow and another using Seaborn's heatmap utility. These visualizations help illustrate the feature relationships in the Iris dataset and offer insights into the dataset's internal structure. 

## Functions

mean(input_matrix):

Calculates and returns the mean values of the features in a data matrix.
Handles matrices with uneven row lengths by filling in missing values with placeholders before computation.
Returns the input directly for single-row matrices.

cov(input_matrix):

Computes the covariance matrix for a dataset.
Accounts for missing data by dynamically adjusting calculations.
Internally uses the custom mean function for mean computation.

