#!/usr/bin/env python
""" evaluate Phyloligo results
"""

import os, sys, argparse

import hdbscan
import numpy as np
import sklearn.cluster as cluster
import sklearn.metrics as metrics
from sklearn.manifold import TSNE

import matplotlib.pyplot as plt
import seaborn as sns

sns.set_context('poster')
sns.set_color_codes()
plot_kwds = {'alpha' : 0.25, 's' : 80, 'linewidths':0}

def get_cmd():
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", action="store", dest="distmat", required=True)
    parser.add_argument("-l", action="store", dest="labels", required=True)
    parser.add_argument("-e", action="store", dest="evaluation", choices=["tsne+hdbscan", "tsne+kmeans", "tsne+kmeans", "dbscan", "hdbscan"], required=True)
    parser.add_argument("-k", action="store", dest="kwargs", nargs="+", default=[])
    
    params = parser.parse_args()
    return params

def process_kwargs(arguments):
    kwargs = dict()
    for argval in arguments:
        arg, val = argval.split("=")
        kwargs[arg] = val
    return kwargs

def read_distmat(path):
    """ read distance matrix
    """
    return np.loadtxt(path)

def read_labels(path):
    """ read labels for the corresponding matrix
    """
    labels = list()
    group2class = dict()
    with open(path) as inf:
        for i, line in enumerate(inf):
            label = line.strip()
            if label not in group2class:
                cl = len(group2class)
                group2class[label] = cl
            else:
                cl = group2class[label]
            labels.append(cl)
    return group2class, np.array(labels)
            
def clusterize_matrix(data, method, args, kwds):
    X = None
    if method == "hdbscan":
        clusterer = hdbscan.HDBSCAN(*args, metric="precomputed")
        clusterer.fit(data)
        labels = clusterer.labels_
    elif method == "dbscan":
        clusterer = cluster.DBSCAN(*args, metric="precomputed")
        clusterer.fit(data)
        labels = clusterer.labels_
    elif method == "tsne+hdbscan":
        tsne_model = TSNE(n_components=2, random_state=0, perplexity=50, metric="precomputed")
        X = tsne_model.fit_transform(data)
        plt.show()
        clusterer = hdbscan.HDBSCAN(*args, **kwds)
        clusterer.fit(X)
        labels = clusterer.labels_
    elif method == "tsne+kmeans":
        tsne_model = TSNE(n_components=2, random_state=0, perplexity=50, metric="precomputed")
        X = tsne_model.fit_transform(data)
        clusterer = hdbscan.KMEANS(*args, **kwds)
        clusterer.fit(X)
        labels = clusterer.labels_
    else:
        print("Error, unknown method {}".format(method), file=sys.stderr)
    return labels, X

def plot_labels(data, labels, algorithm):
    """ display resulting labels
    """
    palette = sns.color_palette('deep', np.unique(labels).max() + 1)
    colors = [palette[x] if x >= 0 else (0.0, 0.0, 0.0) for x in labels]
    plt.scatter(data.T[0], data.T[1], c=colors, **plot_kwds)
    frame = plt.gca()
    frame.axes.get_xaxis().set_visible(False)
    frame.axes.get_yaxis().set_visible(False)
    plt.title('Clusters found by {}'.format(str(algorithm)), fontsize=24)
    plt.show()

def main():
    params = get_cmd()
    
    clust_args = process_kwargs(params.kwargs)
    
    matrix = read_distmat(params.distmat)
    classes, labels_true = read_labels(params.labels)
    
    labels_pred, X = clusterize_matrix(matrix, params.evaluation, (), clust_args)
    
    # plot the different classes if reduction of dimentionality
    if X!=None and X.shape[1] == 2:
        plot_labels(X, labels_pred, params.evaluation)
    
    # clutering evaluation
    randscore = metrics.adjusted_rand_score(labels_true, labels_pred)
    names = [classes[str(i)] for i in range(len(classes))]
    report = metrics.classification_report(labels_true, labels_pred)
    
    print("Rand Score: {}".format(randscore))
    
    print(report)
          
    return 0

if __name__ == "__main__":
    ret = main()
    sys.exit(ret)
    
    