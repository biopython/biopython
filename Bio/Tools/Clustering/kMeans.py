# Copyright 2001 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""This provides code for doing k-Means clustering of data.

k-Means is an algorithm for unsupervised clustering of data.

Glossary:
clusters - A group of closely related data.
centroids - A vector "in the middle" of a cluster.

Functions:
cluster         Cluster a list of data points.

    Distance Functions:
euclidean_dist  The euclidean distance between two points.

"""
import random

try:
    from Numeric import *
except ImportError, x:
    raise ImportError, "This module requires NumPy"

from Bio.Tools import listfns

# Import the euclidean distance function in the kNN class.  This
# should probably be pulled out in a separate module.
from Bio.Tools.Classification.kNN import euclidean_dist, euclidean_dist_py

def random_centroids(data, k):
    """random_centroids(data, k) -> list of centroids

    Return a list of data points to serve as the initial centroids.
    This is k randomly chosen data points.

    """
    centroids = []
    indexes = range(len(data))
    random.shuffle(indexes)
    return take(data, indexes[:k])

def first_k_points_as_centroids(data, k):
    """first_k_points_as_centroids(data, k) -> list of centroids

    Picks the first K points as the initial centroids.  This isn't a
    good method (unless the data is randomized), but does provide
    determinism that's useful for debugging.

    """
    return take(data, range(k))

def _find_closest_centroid(vector, centroids, distance_fn):
    """_find_closest_centroid(vector, centroids, distance_fn) ->
    index of closest centroid

    """
    closest_index = 0
    closest_dist = distance_fn(vector, centroids[0])
    for i in range(1, len(centroids)):
        dist = distance_fn(vector, centroids[i])
        if dist < closest_dist:
            closest_dist = dist
            closest_index = i
    return closest_index

def cluster(data, k, distance_fn=euclidean_dist,
            init_centroids_fn=random_centroids,
            max_iterations=1000, update_fn=None):
    """cluster(data, k[, distance_fn][, max_iterations][, update_fn]) ->
    (centroids, clusters) or None

    Organize data into k clusters.  Return a list of cluster
    assignments between 0-(k-1), where the items in the list
    corresponds to the list of data points.  If the algorithm does not
    converge by max_iterations (default is 1000), returns None.  data
    is a list of data points, which are vectors of numbers.
    distance_fn is a callback function that calculates the distance
    between two vectors.  By default, the Euclidean distance wwill be
    used.  If update_fn is specified, it is called at the beginning of
    every iteration and passed the iteration number, cluster
    centroids, and current cluster assignments.

    """
    # Do some checking to make sure the inputs are reasonable.
    if k < 1:
        raise ValueError, "Please specify a positive number of clusters."
    if not data:
        raise ValueError, "Please pass in some data."
    if len(data) <= k:
        raise ValueError, "Please specify more data points than clusters."
    if max_iterations < 1:
        raise ValueError, "You should have at least one iteraction."
    ndims = len(data[0])
    for i in range(1, len(data)):
        if len(data[i]) != ndims:
            raise ValueError, "All data should have the same dimensionality."

    # Convert the data array into a Numeric array, for speed.
    data = asarray(data, Float)

    # Initialize the clusters without any assignments, and pick the
    # initial centroids.
    clusters = [None] * len(data)
    centroids = init_centroids_fn(data, k)

    for i in range(max_iterations):
        # Call update_fn, if specified.
        if update_fn is not None:
            update_fn(i, centroids, clusters)

        # Assign the clusters.  Each data point is assigned to the
        # closest centroid.
        old_clusters = clusters
        clusters = [None] * len(data)
        for j in range(len(data)):
            clusters[j] = _find_closest_centroid(
                data[j], centroids, distance_fn)

        # Stop if the clusters were the same as the previous iteration.
        if clusters == old_clusters:
            break

        # Calculate the new centroids.  The centroid of a cluster is
        # the average of all its data points.  Calculate the average
        # by summing all the data points in a cluster and dividing by
        # the number of points.  This will result in a divide by zero
        # error and fail if a cluster is empty.  However, this should
        # not happen if the initial centroids were chosen carefully.
        centroids = zeros((k, ndims), Float)
        for j in range(len(data)):
            cluster = clusters[j]
            centroids[cluster] = centroids[cluster] + data[j]
        num_in_cluster = listfns.count(clusters)
        for j in range(k):
            centroids[j] = centroids[j] / num_in_cluster[j]
    else:
        # The loop iterated without converging.
        return None
    return centroids, clusters
