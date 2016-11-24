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
from matplotlib.backends.backend_pdf import PdfPages

import seaborn as sns

sns.set_context('poster')
sns.set_color_codes()
plot_kwds = {'alpha' : 0.25, 's' : 20, 'linewidths':0}

def get_cmd():
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", action="store", dest="distmat", required=True)
    parser.add_argument("-l", action="store", dest="labels", required=True)
    parser.add_argument("-e", action="store", dest="evaluation", choices=["tsne+hdbscan", "tsne+kmeans", "tsne+kmeans", "dbscan", "hdbscan"], required=True)
    parser.add_argument("-k", action="store", dest="kwargs", nargs="+", default=[])
    parser.add_argument("-o", action="store", dest="outprefix")
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
    d = np.load(path)
    return d["arr_0"]

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
            
def clusterize_matrix(data, method, kwds):
    X = np.array([])
    if method == "hdbscan":
        clusterer = hdbscan.HDBSCAN(metric="precomputed")
        clusterer.fit(data)
        labels = clusterer.labels_
    elif method == "dbscan":
        clusterer = cluster.DBSCAN(metric="precomputed")
        clusterer.fit(data)
        labels = clusterer.labels_
    elif method == "tsne+hdbscan":
        tsne_model = TSNE(n_components=2, random_state=0, perplexity=100, metric="precomputed")
        X = tsne_model.fit_transform(data)
        plt.show()
        clusterer = hdbscan.HDBSCAN(**kwds)
        clusterer.fit(X)
        labels = clusterer.labels_
    elif method == "tsne+kmeans":
        tsne_model = TSNE(n_components=2, random_state=0, perplexity=100, metric="precomputed")
        X = tsne_model.fit_transform(data)
        clusterer = hdbscan.KMEANS(**kwds)
        clusterer.fit(X)
        labels = clusterer.labels_
    else:
        print("Error, unknown method {}".format(method), file=sys.stderr)
    return labels, X

def plot_labels(data, labels, algorithm, output):
    """ display resulting labels
    """
    palette = sns.color_palette('deep', np.unique(labels).max() + 1)
    colors = [palette[x] if x >= 0 else (0.0, 0.0, 0.0) for x in labels]
    fig, ax = plt.subplots()
    ax.scatter(data.T[0], data.T[1], c=colors, **plot_kwds)
    frame = plt.gca()
    frame.axes.get_xaxis().set_visible(False)
    frame.axes.get_yaxis().set_visible(False)
    ax.set_title('Clusters found by {}'.format(str(algorithm)), fontsize=24)
    #plt.show()
    plt.savefig(output, dpi=300)
    plt.close()

def plot_distances(matrix, labels, outprefix):
    """ plot distances
    """
    values = matrix.flatten()
    # plot raw distances
    fig, ax = plt.subplots()
    ax.hist(values, normed=True, bins=20)
    ax.set_xlabel("distance")
    ax.set_ylabel("probability density")
    plt.savefig(outprefix+"_rawdistances.pdf", dpi=600)
    #plt.show()
    plt.close()
    # plot distances inter and intra classes
    classes = list(set(labels))
    with PdfPages(outprefix+"_grpdistances.pdf") as pdf:
        for i, cl1 in enumerate(classes):
            idx1 = np.where(labels==cl1)
            for j in range(i, len(classes)):
                cl2 = classes[j]
                idx2 = np.where(labels==cl2)
                classes_dist = matrix[idx1, idx2].flatten()
                fig, ax = plt.subplots()
                ax.hist(classes_dist, normed=True, bins=20)
                ax.set_xlabel("distance")
                ax.set_ylabel("probability density")
                ax.set_title("{} vs {}".format(cl1, cl2))
                pdf.savefig(dpi=600)
                plt.close()

def evaluate_clustering(true_labels, pred_labels):
    """ evaluation of clustering results
    """
    
    expected_classes = list(set(true_labels))
    
    # for each expected classes count the number of different classes
    evaluation = dict()
    for cl in expected_classes:
        idx = np.where(true_labels == cl)
        pred_classes = pred_labels[idx]
        evaluation[cl] = dict()
        for pred_cl in pred_classes:
            if pred_cl != -1:
                evaluation[cl][pred_cl] = evaluation[cl].get(pred_cl, 0) + 1
                
    # now attribute pred classes to expected classes
    predicted = dict()
    for cl in evaluation:
        for pred_cl in evaluation[cl]:
            predicted.setdefault(pred_cl, list())
            predicted[pred_cl].append((evaluation[cl][pred_cl], cl))
    
    # associated classes
    expected2pred = dict()
    for pred_cl in predicted:
        predicted[pred_cl].sort(reverse=True)
        for val, cl in predicted[pred_cl]:
            if cl not in expected2pred:
                expected2pred[cl] = pred_cl
                break
    
    # now compute the TP, TN etc ...
    for cl in expected_classes:
        pred_cl = expected2pred[cl]
        idx = np.where(true_labels == cl)
        pred_classes = pred_labels[idx]
        tot_pred_cl = np.sum(pred_labels==pred_cl)
        tot_other_pred_cl = len(np.where(pred_labels != pred_cl and pred_labels != -1))
        cl_predcl = np.sum(pred_classes==pred_cl)
        cl_othercl = len(np.where(pred_classes != pred_cl and pred_classes != -1))
        print(cl, pred_cl, cl_predcl, tot_pred_cl, cl_othercl, tot_other_pred_cl)

def main():
    params = get_cmd()
    
    clust_args = process_kwargs(params.kwargs)
    clust_args["min_cluster_size"] = 150
    
    matrix = read_distmat(params.distmat)
    classes, labels_true = read_labels(params.labels)
    
    # plot matrix
    plot_distances(matrix, labels_true, params.outprefix)
    
    labels_pred, X = clusterize_matrix(matrix, params.evaluation, clust_args)
    
    # plot the different classes if reduction of dimentionality
    if len(X) > 0 and X.shape[1] == 2:
        plot_labels(X, labels_pred, params.evaluation, params.outprefix+"_predlabels.pdf")
        plot_labels(X, labels_true, params.evaluation, params.outprefix+"_truelabels.pdf")
        
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
    
    