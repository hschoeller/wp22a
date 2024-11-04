#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 12:25:23 2024

@author: schoelleh96
"""

from sklearn.metrics import pairwise_distances

def plot_distance_histogram(dataFrame, bins=50):
    """
    Plots a histogram of all pairwise Euclidean distances between points
    in the input DataFrame.

    Args:
    dataFrame (pd.DataFrame): Input DataFrame with shape (n_samples, n_features).
    bins (int): Number of bins for the histogram.

    Returns:
    None
    """
    # Convert DataFrame to NumPy array
    data = dataFrame.to_numpy()

    # Compute all pairwise Euclidean distances
    distances = pairwise_distances(data, metric='euclidean')

    # Extract the upper triangular part of the distance matrix, excluding the diagonal
    upper_triangle_distances = distances[np.triu_indices(distances.shape[0], k=1)]

    # Plot the histogram of distances
    plt.hist(upper_triangle_distances, bins=bins, edgecolor='black', alpha=0.7)
    plt.title("Histogram of Pairwise Euclidean Distances")
    plt.xlabel("Distance")
    plt.ylabel("Frequency")
    plt.show()