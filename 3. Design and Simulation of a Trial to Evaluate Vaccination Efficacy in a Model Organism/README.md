# Design and Simulation of a Veterinary Trial to Evaluate a Vaccine Efficacy in a Model Organism Based Study

This project simulates and analyzes a hypothetical veterinary trial to evaluate the efficacy of a vaccination in reducing the probability of infection in a population of model organisms. The simulation uses statistical methods, including permutation testing and power analysis, to determine the minimum sample size required to detect a significant difference in infection probabilities between a vaccinated and control group.

The analysis is designed to answer the following key question:
What is the smallest number of model organisms 'N' required to achieve at least 80% power to detect a significant difference at an alpha level (α) of 0.05?

The simulation framework is entirely parameterized by:

x: Probability of infection in the control group

y: Probability of infection in the treatment (vaccinated) group

## Features

Sample Generation:

Randomly generates infection outcomes for control and treatment groups based on specified probabilities x and y.


Permutation Testing:

Uses a Fisher-style permutation test to evaluate the statistical significance of observed differences in mean infection rates between groups.


Power Analysis:

Iteratively calculates the power of the test for varying sample sizes 'N' to identify the smallest 'N' that meets the power threshold (80%).


Visualization:

Produces a power curve that visualizes the relationship between sample size and statistical power.


Parameter Flexibility:

Fully customizable inputs for control (x) and treatment (y) probabilities, allowing users to model a wide range of scenarios.


This project was developed with guidance from the Boston University BF550 Fall 2024 teaching team. Thank you to Dr. Ilijia Dukovsky, Yichi Zhang, and Junming Hu for your support and feedback.
