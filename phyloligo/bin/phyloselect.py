#!/usr/bin/env python3
""" evaluate Phyloligo results
"""
import matplotlib
matplotlib.use("Agg")

import os, sys, argparse
import tempfile, time

from Bio import SeqIO

from phyloligo.core.kmedoids import KMedoids
import hdbscan
import numpy as np
import h5py
import sklearn.cluster as cluster
from sklearn.manifold import TSNE

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

sns.set_context('poster')
sns.set_color_codes()
plot_kwds = {'alpha' : 0.3, 's' : 30, 'linewidths':0}


def get_cmd():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", action="store", dest="distmat", required=True, 
                        help="the input matrix file")
    parser.add_argument("-t", action="store_true", dest="performtsne", default=False,
                        help="perform tsne for visualization and pre-clustering")
    parser.add_argument("-p", action="store", dest="perplexity", default=100, type=int, 
                        help="change the perplexity value")
    parser.add_argument("-m", action="store", dest="method", required=True, choices=["hdbscan", "kmedoids"], 
                        help="method to use to compute cluster on transformed distance matrix")
    parser.add_argument("--minclustersize", action="store", dest="min_cluster_size", type=int, 
                        help="set the minimal cluster size of an HDBSCAN cluster")
    parser.add_argument("--minsamples", action="store", dest="min_samples", type=int, 
                        help="set the minimal sample size of an HDBSCAN cluster")
    parser.add_argument("-k", action="store", dest="nbk", type=int, 
                        help="number of cluster")
    parser.add_argument("-f", action="store", dest="fastafile", 
                        help="path of the original fasta file used for the computation of the distance matrix")
    parser.add_argument("--interactive", action="store_true", dest="interactive", default=False,
                        help="allow the user to run the script in an interactive mode and change "
                        "clustering parameter on the fly (require -t)")
    parser.add_argument("--large", action="store", choices=["memmap", "h5py"], dest="large", help="used in combination with "
                        "joblib for large dataset", default=False)
    parser.add_argument("--noX", action="store_true", dest="noX", help="instead of showing pictures, "
                        "store them in pdf")
    parser.add_argument("-o", action="store", dest="outputdir", required=True)
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
    palette = sns.color_palette('deep', np.unique(labels).max() + 1)
    colors = [palette[x] if x >= 0 else (0.0, 0.0, 0.0) for x in labels]
    fig, ax = plt.subplots()
    ax.scatter(data.T[0], data.T[1], c=colors, **plot_kwds)
    frame = plt.gca()
    frame.axes.get_xaxis().set_visible(False)
    frame.axes.get_yaxis().set_visible(False)
    ax.set_title('Clusters found by {}'.format(str(algorithm)), fontsize=24)
    plt.savefig(output, dpi=300)
    plt.close()

def show_labels(data, labels, algorithm, noX, prefix="", dirout="/tmp/", verbose=0):
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
    palette = sns.color_palette('deep', np.unique(labels).max() + 1)
    colors = [palette[x] if x >= 0 else (0.0, 0.0, 0.0) for x in labels]
    fig, ax = plt.subplots()
    ax.scatter(data.T[0], data.T[1], c=colors, **plot_kwds)
    frame = plt.gca()
    frame.axes.get_xaxis().set_visible(False)
    frame.axes.get_yaxis().set_visible(False)
    ax.set_title('Clusters found by {}'.format(str(algorithm)), fontsize=24)
    if noX:
        date = time.strftime("%Y%m%d")
        curtime = time.strftime("%H%M%S")
        pathout = tempfile.mktemp(
            prefix="phyloselect_{}_{}_{}_".format(date, curtime, prefix), 
            suffix=".pdf", dir=dirout)
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
            pathout = os.path.join(outputdir, "data_fasta_noclass.fa")
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
    if params.large != []:
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
                suffix=".pdf", dir="/tmp/")
            #pathout = tempfile.mktemp(suffix=".pdf", prefix="initial_clustering_", dir="/tmp/")
            print("Saving dimensionallity reduction to {}".format(pathout))
            plt.savefig(pathout)
        else:
            plt.show()
        plt.close()
    
    print("Clusterize")
    labels_pred = clusterize(matrix, params.method, min_cluster_size=params.min_cluster_size, min_samples=params.min_samples, nbk=params.nbk)
    
    # plot the different classes if reduction of dimentionality
    if params.performtsne and not params.interactive:
        pathout = os.path.join(params.outputdir, "data_tsne_reduc.pdf")
        print("Save tsne clustering projection in ".format(pathout))
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
                        dirout="/tmp/", verbose=1)
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
                
        
        ²
          
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
    
    
