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

"""
import random

try:
    from Numeric import *
except ImportError, x:
    raise ImportError, "This module requires NumPy"

from Bio import listfns
from Bio import distance

def random_centroids(data, k):
    """random_centroids(data, k) -> list of centroids

    Return a list of data points to serve as the initial centroids.
    This is k randomly chosen data points.  Tries to avoid having
    repeated centroies, if possible.

    """
    if k > len(data):
        raise ValueError, "k is larger than the number of data points"
    # Randomize the centroids.
    indexes = range(len(data))
    random.shuffle(indexes)
    # Now get a list of the first k unique data points.
    centroid_indexes = []
    seen = {}
    i = 0
    while len(centroid_indexes) < k and i < len(indexes):
        key = ','.join(map(str, data[i]))   # make the array hashable
        if not seen.has_key(key):
            centroid_indexes.append(i)
            seen[key] = 1
        i += 1
    # If there aren't k unique data points, then just pick the first
    # k.
    if len(centroid_indexes) < k:
        centroid_indexes = indexes[:k]
    return take(data, centroid_indexes)

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

def cluster(data, k, distance_fn=distance.euclidean,
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
        clusters = [_find_closest_centroid(x, centroids, distance_fn) \
                    for x in data]

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
            c, d = clusters[j], data[j]
            centroids[c] = centroids[c] + d
        num_in_cluster = listfns.count(clusters)
        for j in range(k):
            if num_in_cluster.has_key(j):
                centroids[j] = centroids[j] / num_in_cluster[j]
    else:
        # The loop iterated without converging.
        return None
    return centroids, clusters

import ckMeans
_find_closest_centroid = ckMeans._find_closest_centroid
