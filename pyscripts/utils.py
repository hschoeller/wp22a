#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 12:25:23 2024

@author: schoelleh96
"""

import pickle
import os
from sklearn.metrics import pairwise_distances
from sklearn.neighbors import BallTree
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from dask import delayed, compute
import scipy.sparse as sps
import scipy.linalg as spl
from scipy.linalg import LinAlgError
import pandas as pd
import seaborn as sb


def plot_distance_histogram(data, bins=50):
    """
    Plot a histogram of all pairwise Euclidean distances.

    Parameters
    ----------
    data (numpy.ndarray): Input data (n_samples, n_features).
    bins (int): Number of bins for the histogram.

    Returns
    -------
    None
    """
    # Compute all pairwise Euclidean distances
    distances = pairwise_distances(data, metric='euclidean')

    # Extract the upper triangular part of the distance matrix,
    # excluding the diagonal
    upper_triangle_distances = distances[np.triu_indices(distances.shape[0],
                                                         k=1)]

    # Plot the histogram of distances
    plt.hist(upper_triangle_distances, bins=bins, edgecolor='black', alpha=0.7)
    plt.title("Histogram of Pairwise Euclidean Distances")
    plt.xlabel("Distance")
    plt.ylabel("Frequency")
    plt.show()


def calculate_distances(data, radius):
    """
    Calculate Euclidean distances for the input dataframe.

    Parameters
    ----------
    data (numpy.ndarray): Input data with shape (n_samples, n_features).
    radius (float): The maximum radius for neighbor queries.

    Returns
    -------
    row_indices, col_indices, dist_values: Three lists containing row, col
        indices and the values according to the sparse distance matrix.
    """
    # Create a BallTree for efficient neighbor queries
    tree = BallTree(data, metric='euclidean')

    @delayed
    def query_point(i):
        """Query radius for a single point and return indices and distances."""
        ind, dist = tree.query_radius([data[i]], r=radius,
                                      return_distance=True)
        return i, ind[0], dist[0]

    # Create a list of delayed query tasks
    delayed_results = [query_point(i) for i in range(data.shape[0])]

    # Execute the graph in parallel
    results = compute(*delayed_results)
    row_indices = []
    col_indices = []
    dist_values = []

    for i, ind, dist in results:
        row_indices.extend([i] * len(ind))  # Row index
        col_indices.extend(ind)  # Column indices are the neighbors
        dist_values.extend(dist)  # Distances for each neighbor

    return row_indices, col_indices, dist_values


def calc_Sums(data, radius, epsilon_opt, d_path, force_calc=False):
    """
    Calculate sum of diffusion distances.

    Parameters
    ----------
    data (numpy.ndarray): Input data with shape (n_samples, n_features).
    radius (float): The maximum radius for neighbor queries.
    epsilon_opt (np.ndarray): array containing the epsilon options
    d_path (str): path to store data

    Returns
    -------
    Sums (pd.DataFrame): Contains diffusion distances across varying eps.

    """
    Sum_list = list()

    if radius == 1e6:
        d_path = d_path + "full"
    if os.path.exists(d_path + "dist.pkl") and not force_calc:
        with open(d_path + "dist.pkl", "rb") as f:
            print("loading " + d_path + "dist.pkl")
            dist = pickle.load(f)
    else:
        row, col, dist = calculate_distances(data, radius)
        print("No. of dists: ")
        print(len(dist))
        print("Avg. dist:")
        print(np.asarray(dist).sum()/(len(dist)-data.shape[0]))
        with open(d_path + "rows.pkl", "wb") as f:
            print(d_path + "rows.pkl")
            pickle.dump(row, f)
        with open(d_path + "cols.pkl", "wb") as f:
            print(d_path + "cols.pkl")
            pickle.dump(col, f)
        with open(d_path + "dist.pkl", "wb") as f:
            print(d_path + "dist.pkl")
            pickle.dump(dist, f)
    i = 0
    for e in epsilon_opt:
        print(e)
        vals = np.exp(-np.square(dist)/e)
        print(vals)
        df = pd.DataFrame({'eps': e,
                           'sums': (vals.sum() / (data.shape[0]**2)),
                           'm': data.shape[0],
                           'dist': np.asarray(dist).mean()}, index=[i])
        i += 1
        Sum_list.append(df)

    Sums = pd.concat(Sum_list, ignore_index=True)

    return Sums


def calc_dd(data):
    """
    Calculate some derived values using diffusion distances.

    Parameters
    ----------
    data (pd.DataFrame): The diffusion distances.

    Returns
    -------
    pd.Series: Contains the derived quantaties.

    """
    leps = np.log(data["eps"]).reset_index(drop=True)
    lS = np.log(data["sums"] * np.square(data["m"])).reset_index(drop=True)

    # This is the derivative of log(S(eps) with respect to log(eps):
    dlS = lS.diff().dropna() / leps.diff().dropna()
    # We assume that the steps in the logarithm of eps are uniform

    e_max = np.argmax(dlS)
    d = 2 * dlS[e_max]

    density = ((np.log(np.pi) + leps[e_max])/2 +
               (np.log(data["m"].reset_index(drop=True)[0]) - lS[e_max])/d)
    rho = density * d

    meandist = np.mean(np.sqrt(-np.log(data["sums"]) * data["eps"]))
    dist = np.mean(data["dist"])
    return pd.Series({"dimension": d, "density": density, "rho": rho,
                      "e_max": np.exp(leps[e_max]), "L": np.exp(density),
                      "meandist": meandist, "dist": dist})


def plot_epsloglog(Sums, e_max, dimension):
    """
    Plot double logarithmic curve for sum of diffusion distances over eps.

    Parameters
    ----------
    Sums (pd.DataFrame): Diffusion distances sums.

    Returns
    -------
    None.

    """
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)

    sb.lineplot(data=Sums, x="eps", y="sums",  # marker=".",
                dashes=False,
                linewidth=.5)
    evaluated_x_values = Sums['eps'].unique()
    sb.rugplot(x=evaluated_x_values, height=0.025, color='green')
    idx = np.argwhere(np.abs(Sums["eps"] - e_max) < 1e-6)
    ax.axline(xy1=(e_max, Sums.iloc[idx[0, 0]]["sums"]), slope=dimension/2,
              color="red", linewidth=.5)
    ax.axline(xy1=(e_max, Sums.iloc[idx[0, 0]]["sums"]),
              slope=np.floor(dimension/2),
              color="orange", linewidth=.5)
    ax.axline(xy1=(e_max, Sums.iloc[idx[0, 0]]["sums"]),
              slope=np.floor(dimension/2)-1,
              color="orange", linewidth=.5)
    ax.axline(xy1=(e_max, Sums.iloc[idx[0, 0]]["sums"]),
              slope=np.ceil(dimension/2),
              color="orange", linewidth=.5)
    ax.axline(xy1=(e_max, Sums.iloc[idx[0, 0]]["sums"]),
              slope=np.ceil(dimension/2)+3,
              color="orange", linewidth=.5)
    ax.set_xlabel("$\\epsilon$ [km\\textsuperscript{2}]")
    ax.set_ylabel("$S_{t}(\\epsilon) m^{-2}$ [--]")

    plt.xscale('log')
    plt.yscale('log')
    # Create custom legend entries
    legend_elements = [
        Line2D([0], [0], color='red', lw=0.5, label=f'Slope = {dimension/2}'),
        Line2D([0], [0], color='orange', lw=0.5, label=f'Slope = {np.floor(
            dimension/2)-1}, {np.floor(dimension/2)}, {np.ceil(dimension/2)},{np.ceil(dimension/2)+1}')
    ]

    # Add the legend to the plot
    ax.legend(handles=legend_elements, loc='best')


def calc_diff_map(eps, N_v, dist_mat):
    """
    Calculate diffusion maps.

    Diffusion maps: eigenvectors of the diffusion transition matrix
                        along with its eigenvalues.

    Parameters
    ----------
    eps (float): diffusion bandwidth.
    N_v (int): how many eigenvalues and -vectors to compute.
    dist_mat (sps.csr_matrix): distance matrix

    Returns
    -------
    vals (np.ndarray): Eigenvalues.
    vecs (dp.ndarray): Eigenvectors (the diffusion maps).

    """
    dist_mat.data = np.exp(-dist_mat.data**2/eps)
    pre_norm_mat = dist_mat / \
        np.outer(dist_mat.sum(axis=1), dist_mat.sum(axis=1))
    norm_mat = pre_norm_mat.multiply(1/pre_norm_mat.sum(axis=1))

    try:
        vals, vecs = sps.linalg.eigs(norm_mat, k=N_v)
    except LinAlgError as e:
        print(f"Eigs (sparse) failed with error {e}, using eig")
        vals, vecs = spl.eig(norm_mat.toarray())
    vals = np.flip(np.sort(np.real(vals)))[:N_v]
    return (vals, np.real(vecs))
