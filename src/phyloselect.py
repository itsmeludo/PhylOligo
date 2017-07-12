#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" evaluate Phyloligo results

K-medoids by :
    clustering Timo Erkkilä <timo.erkkila@gmail.com>
    Antti Lehmussola <antti.lehmussola@gmail.com>
    License: BSD 3 clause

"""
import matplotlib
matplotlib.use("Agg")

import os, sys, argparse
import tempfile, time

from Bio import SeqIO

import hdbscan
import numpy as np
import h5py
import sklearn.cluster as cluster
from sklearn.manifold import TSNE
from sklearn.base import BaseEstimator, ClusterMixin, TransformerMixin
from sklearn.metrics.pairwise import PAIRWISE_DISTANCE_FUNCTIONS
from sklearn.utils import check_array, check_random_state
from sklearn.utils.validation import check_is_fitted


import matplotlib.pyplot as plt

plot_kwds = {'alpha' : 0.3, 's' : 30, 'linewidths':0}

import warnings


class KMedoids(BaseEstimator, ClusterMixin, TransformerMixin):
    """
    k-medoids class.

    Parameters
    ----------
    n_clusters : int, optional, default: 8
        How many medoids. Must be positive.

    distance_metric : string, optional, default: 'euclidean'
        What distance metric to use.

    clustering : {'pam'}, optional, default: 'pam'
        What clustering mode to use.

    init : {'random', 'heuristic'}, optional, default: 'heuristic'
        Specify medoid initialization.

    max_iter : int, optional, default : 300
        Specify the maximum number of iterations when fitting.

    random_state : int, optional, default: None
        Specify random state for the random number generator.
    """

    # Supported clustering methods
    CLUSTERING_METHODS = ['pam']

    # Supported initialization methods
    INIT_METHODS = ['random', 'heuristic']

    def __init__(self, n_clusters=8, distance_metric='euclidean',
                 clustering_method='pam', init='heuristic',
                 max_iter=300, random_state=None):

        self.n_clusters = n_clusters

        self.distance_metric = distance_metric

        self.init = init

        self.max_iter = max_iter

        self.clustering_method = clustering_method

        self.random_state = random_state

    def _check_init_args(self):

        # Check n_clusters
        if self.n_clusters is None or self.n_clusters <= 0 or \
                not isinstance(self.n_clusters, int):
            raise ValueError("n_clusters has to be nonnegative integer")

        # Check distance_metric
        if callable(self.distance_metric):
            self.distance_func = self.distance_metric
        elif self.distance_metric in PAIRWISE_DISTANCE_FUNCTIONS:
            self.distance_func = \
                PAIRWISE_DISTANCE_FUNCTIONS[self.distance_metric]
        else:
            raise ValueError("distance_metric needs to be " +
                             "callable or one of the " +
                             "following strings: " +
                             "{}".format(PAIRWISE_DISTANCE_FUNCTIONS.keys()) +
                             ". Instead, '{}' ".format(self.distance_metric) +
                             "was given.")

        # Check clustering_method
        if self.clustering_method not in self.CLUSTERING_METHODS:
            raise ValueError("clustering must be one of the following: " +
                             "{}".format(self.CLUSTERING_METHODS))

        # Check init
        if self.init not in self.INIT_METHODS:
            raise ValueError("init needs to be one of " +
                             "the following: " +
                             "{}".format(self.INIT_METHODS))

        # Check random state
        self.random_state_ = check_random_state(self.random_state)

    def fit(self, X, y=None):
        """Fit K-Medoids to the provided data.

        Parameters
        ----------
        X : array-like or sparse matrix, shape=(n_samples, n_features)

        Returns
        -------
        self
        """

        self._check_init_args()

        # Check that the array is good and attempt to convert it to
        # Numpy array if possible
        X = self._check_array(X)

        # Apply distance metric to get the distance matrix
        if self.distance_metric != "precomputed":
            D = self.distance_func(X)
        else:
            D = X

        medoid_ics = self._get_initial_medoid_indices(D, self.n_clusters)

        # Old medoids will be stored here for reference
        old_medoid_ics = np.zeros((self.n_clusters,))

        # Continue the algorithm as long as
        # the medoids keep changing and the maximum number
        # of iterations is not exceeded
        self.n_iter_ = 0
        while not np.all(old_medoid_ics == medoid_ics) and \
                self.n_iter_ < self.max_iter:

            self.n_iter_ += 1

            # Keep a copy of the old medoid assignments
            old_medoid_ics = np.copy(medoid_ics)

            # Get cluster indices
            cluster_ics = self._get_cluster_ics(D, medoid_ics)

            # Update medoids with the new cluster indices
            self._update_medoid_ics_in_place(D, cluster_ics, medoid_ics)

        # Expose labels_ which are the assignments of
        # the training data to clusters
        self.labels_ = cluster_ics

        # Expose cluster centers, i.e. medoids
        self.cluster_centers_ = X.take(medoid_ics, axis=0)

        # Return self to enable method chaining
        return self

    def _check_array(self, X):

        X = check_array(X)

        # Check that the number of clusters is less than or equal to
        # the number of samples
        if self.n_clusters > X.shape[0]:
            raise ValueError("The number of medoids " +
                             "({}) ".format(self.n_clusters) +
                             "must be larger than the number " +
                             "of samples ({})".format(X.shape[0]))

        return X

    def _get_cluster_ics(self, D, medoid_ics):
        """Returns cluster indices for D and current medoid indices"""

        # Assign data points to clusters based on
        # which cluster assignment yields
        # the smallest distance
        cluster_ics = np.argmin(D[medoid_ics, :], axis=0)

        return cluster_ics

    def _update_medoid_ics_in_place(self, D, cluster_ics, medoid_ics):
        """In-place update of the medoid indices"""

        # Update the medoids for each cluster
        for cluster_idx in range(self.n_clusters):

            if sum(cluster_ics == cluster_idx) == 0:
                warnings.warn("Cluster {} is empty!".format(cluster_idx))
                continue

            # Find current cost that is associated with cluster_idx.
            # Cost is the sum of the distance from the cluster
            # members to the medoid.
            curr_cost = np.sum(D[medoid_ics[cluster_idx],
                                 cluster_ics == cluster_idx])

            # Extract the distance matrix between the data points
            # inside the cluster_idx
            D_in = D[cluster_ics == cluster_idx, :]
            D_in = D_in[:, cluster_ics == cluster_idx]

            # Calculate all costs there exists between all
            # the data points in the cluster_idx
            all_costs = np.sum(D_in, axis=1)

            # Find the index for the smallest cost in cluster_idx
            min_cost_idx = np.argmin(all_costs)

            # find the value of the minimum cost in cluster_idx
            min_cost = all_costs[min_cost_idx]

            # If the minimum cost is smaller than that
            # exhibited by the currently used medoid,
            # we switch to using the new medoid in cluster_idx
            if min_cost < curr_cost:

                # Find data points that belong to cluster_idx,
                # and assign the newly found medoid as the medoid
                # for cluster c
                medoid_ics[cluster_idx] = \
                    np.where(cluster_ics == cluster_idx)[0][min_cost_idx]

    def transform(self, X):
        """Transforms X to cluster-distance space.

        Parameters
        ----------
        X : array-like or sparse matrix, shape=(n_samples, n_features)
            Data to transform.

        Returns
        -------
        X_new : array, shape=(n_samples, n_clusters)
            X transformed in the new space.
        """

        check_is_fitted(self, "cluster_centers_")

        # Apply distance metric wrt. cluster centers (medoids),
        # and return these distances
        return self.distance_func(X, Y=self.cluster_centers_)

    def predict(self, X):

        check_is_fitted(self, "cluster_centers_")

        # Check that the array is good and attempt to convert it to
        # Numpy array if possible
        X = check_array(X)

        # Apply distance metric wrt. cluster centers (medoids)
        D = self.distance_func(X, Y=self.cluster_centers_)

        # Assign data points to clusters based on
        # which cluster assignment yields
        # the smallest distance
        labels = np.argmin(D, axis=1)

        return labels

    def inertia(self, X):

        # Map the original X to the distance-space
        Xt = self.transform(X)

        # Define inertia as the sum of the sample-distances
        # to closest cluster centers
        inertia = np.sum(np.min(Xt, axis=1))

        return inertia

    def _get_initial_medoid_indices(self, D, n_clusters):

        if self.init == 'random':  # Random initialization

            # Pick random k medoids as the initial ones.
            medoids = self.random_state_.permutation(D.shape[0])[:n_clusters]

        elif self.init == 'heuristic':  # Initialization by heuristic

            # Pick K first data points that have the smallest sum distance
            # to every other point. These are the initial medoids.
            medoids = list(np.argsort(np.sum(D, axis=1))[:n_clusters])

        else:

            raise ValueError("Initialization not implemented for method: " +
                             "'{}'".format(self.init))

        return medoids


def get_cmd():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", action="store", dest="distmat", required=True, 
                        help="The input matrix file")
    parser.add_argument("-t", action="store_true", dest="performtsne", default=False,
                        help="Perform tsne for visualization and pre-clustering")
    parser.add_argument("-p", action="store", dest="perplexity", default=100, type=int, 
                        help="Change the perplexity value")
    parser.add_argument("-m", action="store", dest="method", required=True, choices=["hdbscan", "kmedoids"], 
                        help="Method to use to compute cluster on transformed distance matrix")
    parser.add_argument("--minclustersize", action="store", dest="min_cluster_size", type=int, 
                        help="Set the minimal cluster size of an HDBSCAN cluster")
    parser.add_argument("--minsamples", action="store", dest="min_samples", type=int, 
                        help="Set the minimal sample size of an HDBSCAN cluster")
    parser.add_argument("-k", action="store", dest="nbk", type=int, 
                        help="Number of cluster")
    parser.add_argument("-f", action="store", dest="fastafile", 
                        help="Path of the original fasta file used for the computation of the distance matrix")
    parser.add_argument("--interactive", action="store_true", dest="interactive", default=False,
                        help="Allow the user to run the script in an interactive mode and change "
                        "clustering parameter on the fly (require -t)")
    parser.add_argument("--large", action="store", choices=["memmap", "h5py"], dest="large", help="Used in combination with "
                        "joblib for large dataset", default=False)
    parser.add_argument("--noX", action="store_true", dest="noX", help="Instead of showing pictures, "
                        "store them in png")
    parser.add_argument("-o", action="store", dest="outputdir", required=True)
    parser.add_argument("-q", "--infreq", action="store", dest="in_freq_file", 
                        help="kmer frequencies input file[default:%(default)s]. If provided, the clustering is performed on the kmer frequency matrix instead of on the contig distance matrix.")
    params = parser.parse_args()
    
    if params.interactive and not params.performtsne:
        print("Error, interactive mode (--interactive) requires tsne (-t)", file=sys.stderr)
        sys.exit(1)
    return params

def read_distmat(path):
    """ read distance matrix
    
    Parameter:
    ----------
    path: string
        path of the input matrix
        
    Return:
    -------
    data: np.ndarray
        the ditance matrix
    """
    return np.loadtxt(path)

def read_freqmat(path):
    """ read distance matrix
    
    Parameter:
    ----------
    path: string
        path of the input matrix
        
    Return:
    -------
    data: np.ndarray
        the ditance matrix
    """
    return np.loadtxt(path)


def transform_matrix_tsne(data, perplexity):
    """ perform tsne on the initial distance matrix for visualization
    
    Parameter:
    ----------
    data: np.ndarray
        the ditance matrix, n by n
    perplexity: int
        the perplexity parameter of t-sne, change dimensionallity reduction (suggested between 30 and 100)
        
    Return:
    -------
    X: np.ndarray
        a n by 2 matrix, each point corresponding to an entry of the distance matrix
    """
    tsne_model = TSNE(n_components=2, random_state=0, perplexity=perplexity, metric="precomputed")
    X = tsne_model.fit_transform(data)
    return X

def find_clusters(data, method, kwargs):
    """ find cluster in the dataset
    
    Parameters:
    -----------
    data: np.ndarray
        bi-dimensional matrix (distance matrix transformed by t-sne)
    method: string
        name of the clustering method to use
    kwargst: dict
        optional arguments of the clustering mthod
    
    Return:
    -------
    labels: np.array
        predicted classes of each matrix entry
    
    """
    if method == "hdbscan":
        clusterer = hdbscan.HDBSCAN(**kwargs)
        clusterer.fit(data.astype(np.float))
        labels = clusterer.labels_
    elif method == "kmedoids":
        clusterer = KMedoids(**kwargs)
        clusterer.fit(data)
        labels = clusterer.labels_
    else:
        print("Error, unknown method {}".format(method), file=sys.stderr)
    return labels
    
def plot_labels(data, labels, algorithm, output):
    """ display resulting labels
    
    Parameters:
    -----------
    data: np.ndarray
        bi-dimensional matrix (distance matrix transformed by t-sne)
    labels: np.array
        predicted classes of each matrix entry
    algorithm: string
        name of the clustering method used
    output: string
        where to save the output plot
    """
    #palette = sns.color_palette('deep', np.unique(labels).max() + 1
    norm = matplotlib.colors.Normalize(vmin=0, vmax=labels.max())
    palette = plt.get_cmap("gist_ncar")
    #palette = plt.get_cmap("viridis")
    colors = [palette(norm(x)) if x >= 0 else (0.0, 0.0, 0.0) for x in labels]
    fig, ax = plt.subplots()
    ax.scatter(data.T[0], data.T[1], c=colors, **plot_kwds)
    frame = plt.gca()
    frame.axes.get_xaxis().set_visible(False)
    frame.axes.get_yaxis().set_visible(False)
    ax.set_title('Clusters found by {}'.format(str(algorithm)), fontsize=24)
    plt.savefig(output, dpi=300)
    plt.close()

def show_labels(data, labels, algorithm, noX, prefix="", dirout=None, verbose=0):
    """ display resulting labels
    
    Parameters:
    -----------
    data: np.ndarray
        bi-dimensional matrix (distance matrix transformed by t-sne)
    labels: np.array
        predicted classes of each matrix entry
    algorithm: string
        name of the clustering method used
    output: string
        where to save the output plot
    """
    #palette = sns.color_palette('deep', np.unique(labels).max() + 1)
    norm = matplotlib.colors.Normalize(vmin=0, vmax=labels.max())
    palette = plt.get_cmap("gist_ncar")
    #palette = plt.get_cmap("viridis")
    colors = [palette(norm(x)) if x >= 0 else (0.0, 0.0, 0.0) for x in labels]
    fig, ax = plt.subplots()
    ax.scatter(data.T[0], data.T[1], c=colors, **plot_kwds)
    frame = plt.gca()
    frame.axes.get_xaxis().set_visible(False)
    frame.axes.get_yaxis().set_visible(False)
    ax.set_title('Clusters found by {}'.format(str(algorithm)), fontsize=24)
    if noX:
        date = time.strftime("%Y%m%d")
        curtime = time.strftime("%H%M%S")
        if dirout != None:
            pathout = tempfile.mktemp(
                prefix="phyloselect_{}_{}_{}_".format(date, curtime, prefix), 
                suffix=".png", dir=dirout)
        else:
            pathout = tempfile.mktemp(
                prefix="phyloselect_{}_{}_{}_".format(date, curtime, prefix), 
                suffix=".png")
        if verbose :
            print("saving file at {}".format(pathout))
        plt.savefig(pathout)
    else:
        if verbose:
            print("Displaing current clustering with algorithm {}".format(algorithm))
        plt.show()
    plt.close()

    
def get_fasta_record(indexes, fastafile):
    """ get specific fasta record from the initial file
    
    Parameters
    ----------
    indexes: np.array
        read index of a specific group
    fastafile: string
        path to the initial fasta file
        
    Return
    ------
    records: list
        selected fasta record corresponding to indexes
    """
    cnt = 0
    records = list()
    for record in SeqIO.parse(fastafile, "fasta"):
        if cnt in indexes:
            records.append(record)
        cnt += 1
    return records
    
def write_fastafile(labels_pred, fastafile, outputdir):
    """ write each fasta sequence of a specific cluster into a fasta file
    
    Parameters
    ----------
    labels_pred: np.array
        group label of each entry of the matrix
    fastafile: string
        path to the initial fasta file
    outputdir: string
        path of the output directory
    """
    all_classes = np.unique(labels_pred)
    for cl in all_classes:
        if cl == -1:
            pathout = os.path.join(outputdir, "data_fasta_unclust.fa".format(cl))
            indexes = np.where(labels_pred == cl)[0]
            records = get_fasta_record(set(list(indexes)), fastafile)
            with open(pathout, "w") as outf:
                SeqIO.write(records, outf, "fasta")
        else:
            pathout = os.path.join(outputdir, "data_fasta_cl{}.fa".format(cl))
            indexes = np.where(labels_pred == cl)[0]
            records = get_fasta_record(set(list(indexes)), fastafile)
            with open(pathout, "w") as outf:
                SeqIO.write(records, outf, "fasta")

def clusterize(data, method, min_cluster_size=None, min_samples=None, nbk=None):
    # cluster points
    kwargs = dict()
    
    if method == "hdbscan":
        if min_cluster_size != None:
            kwargs["min_cluster_size"] = min_cluster_size
        if min_samples != None:
            kwargs["min_samples"] = min_samples
        kwargs["metric"] = "precomputed"
            
    if method == "kmedoids":
        if nbk != None:
            kwargs["n_clusters"] = nbk
        kwargs["distance_metric"] = "precomputed"
        
    labels_pred = find_clusters(data, method, kwargs)
    return labels_pred

def main():
    params = get_cmd()
    
    if not os.path.isdir(params.outputdir):
        os.makedirs(params.outputdir)
    
    # get matrix and transform the data to a bidimensional set of point with tsne
    print("Read matrix")
    if params.large != False:
        if params.large == "memmap":
            # if matrix was created using --large option read it as a memmap matrix
            matrix = np.memmap(params.distmat, dtype=np.float32, mode="r")
            s = matrix.shape[0]
            n = np.sqrt(s)
            if str(n).split(".")[1] != "0":
                print("Error, weird shape for matrix {}".format(params.distmat), file=sys.stderr)
                sys.exit(1)
            matrix = matrix.reshape((int(n), int(n)))
        elif params.large == "h5py":
            # read matrix as a h5py file
            with h5py.File(params.distmat, "r") as hf:
                matrix = hf.get("distances")
                matrix = matrix.value[:]
    else:
        # numpy savetxt matrix
        matrix = read_distmat(params.distmat)
        
    if params.performtsne:
        print("Transform matrix")
        data = transform_matrix_tsne(matrix, params.perplexity)
    else:
        data = matrix
    
    # display original tsne:
    if params.interactive:
        fig, ax = plt.subplots()
        ax.scatter(data[:,0], data[:,1], **plot_kwds)
        if params.noX:
            date = time.strftime("%Y%m%d")
            curtime = time.strftime("%H%M%S")
            pathout = tempfile.mktemp(
                prefix="phyloselect_init_{}_{}".format(date, curtime), 
                suffix=".png", dir=params.outputdir)
            print("Saving dimensionallity reduction to {}".format(pathout))
            plt.savefig(pathout)
        else:
            plt.show()
        plt.close()
    
    print("Clusterize")
    labels_pred = clusterize(matrix, params.method, min_cluster_size=params.min_cluster_size, min_samples=params.min_samples, nbk=params.nbk)
    
    # plot the different classes if reduction of dimentionality
    if params.performtsne and not params.interactive:
        pathout = os.path.join(params.outputdir, "data_tsne_reduc.png")
        print("Save tsne clustering projection in {}".format(pathout))
        plot_labels(data, labels_pred, params.method, pathout)
    elif params.interactive:
        # loop until user is satisfied
        msg = "y"
        method = params.method
        min_cl_size = params.min_cluster_size
        min_samples = params.min_samples
        nbk = params.nbk
        cnt = 0
        while msg != "n":
            labels_pred = clusterize(matrix, method, 
                                     min_cluster_size=min_cl_size, min_samples=min_samples, nbk=nbk)
            show_labels(data, labels_pred, method, params.noX, prefix="clustering_{}".format(cnt), 
                        dirout=params.outputdir, verbose=1)
            cnt += 1
            print("perform an other run? [y/n]")
            msg = input("--> ")
            if msg == "y":
                new_method = ""
                while new_method not in ["hdbscan", "kmedoids", "n"]:
                    print("change method? [hdbscan/kmedoids/n]")
                    new_method = input("--> ")
                if new_method != "n":
                    method = new_method
                if method == "hdbscan":
                    new_min_cl_size = None
                    while new_min_cl_size == None:
                        print("[hdbscan] change min_cluster_size, old value = {}?".format(min_cl_size))
                        new_min_cl_size = input("--> ")
                        if new_min_cl_size == "" and min_cl_size != None :
                            new_min_cl_size = min_cl_size
                        else:
                            try:
                                new_min_cl_size = int(new_min_cl_size)
                            except ValueError:
                                print("Please provide an integer value")
                                new_min_cl_size = -1
                            if new_min_cl_size < 1:
                                new_min_cl_size = None
                            else:
                                min_cl_size = new_min_cl_size
                    new_min_samples = None
                    while new_min_samples == None:
                        print("[hdbscan] change min_samples, old value = {}?".format(min_samples))
                        new_min_samples = input("--> ")
                        if new_min_samples == "" and min_samples != None:
                            new_min_samples = min_samples
                        else:
                            try:
                                new_min_samples = int(new_min_samples)
                            except ValueError:
                                print("Please provide an integer value")
                                new_min_samples = -1
                            if new_min_samples < 1:
                                new_min_samples = None
                            else:
                                min_samples = new_min_samples
                    print("[hdbscan] running hdbscan with parameters min_cluster_size={} min_samples={}".format(min_cl_size, min_samples))
                elif method == "kmedoids":
                    new_nbk = None
                    while new_nbk == None:
                        print("[kmedoids] change nbk, old value = {}?".format(nbk))
                        new_nbk = input("--> ")
                        if new_nbk == "" and nbk != None:
                            new_nbk = nbk
                        else:
                            try:
                                new_nbk = int(new_nbk)
                            except ValueError:
                                print("Please provide an integer value")
                                new_nbk = -1
                            if new_nbk < 1:
                                new_nbk = None
                            else:
                                nbk = new_nbk
                    print("[kmedoids] running kmedoids with parameters nbk={}".format(nbk))
                else:
                    print("Unknown method name {}".format(method))
                    sys.exit(1)
          
    # write cluster indexes
    pathout = os.path.join(params.outputdir, "data_cluster_indexes.dat")
    print("Store cluster indexes in {}".format(pathout))
    all_classes = np.unique(labels_pred)
    with open(pathout, "w") as outf:
        for cl in all_classes:
            indexes = np.where(labels_pred == cl)[0]
            for idx in indexes:
                outf.write("{} {}\n".format(cl, idx))
                
    if params.fastafile:
        print("Write fasta per classes in {}/data_fasta_*.fa".format(params.outputdir))
        write_fastafile(labels_pred, params.fastafile, params.outputdir)
        
    return 0

if __name__ == "__main__":
    ret = main()
    sys.exit(ret)
    
    
