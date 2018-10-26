import numpy as np
import Bio.Cluster
from sklearn.metrics.pairwise import pairwise_distances
import sys
""" distances must be "imported" to avoid multiprocessing Pool bug
https://stackoverflow.com/questions/41385708/multiprocessing-example-giving-attributeerror
"""

np.seterr(divide='ignore', invalid='ignore')


def posdef_check_value(d):
    d[np.isnan(d)]=0    
    d[np.isinf(d)]=0
    
## Distance functions

def KL(a, b):
    """ compute the KL distance
    """
    if a.ndim == 1 and b.ndim == 1:
        d = a * np.log(a/b)
        posdef_check_value(d)
        res = np.sum(d)
    elif a.ndim == 2 and b.ndim == 2:
        X, Y = check_pairwise_arrays(a, b)
        X = X[:,np.newaxis]
        d = X * np.log(X/Y)
        posdef_check_value(d)
        res = np.sum(d, axis=2).T
    else:
        print("Dimension erro in KL, a={}, b={}".format(a.ndim, b.ndim), file=sys.stderr)
        sys.exit(1)
    return res

def Eucl(a, b):
    """ compute Euclidean distance
    """
    d = pow(a-b,2)
    posdef_check_value(d)
    return np.sqrt(np.sum(d))

def JSD(a, b):
    """ Compute JSD distance
    """
    if a.ndim == 1 and b.ndim == 1:
        h = 0.5 * (a + b)
        d = 0.5 * (KL(a, h) + KL(b, h))
    elif a.ndim==2 and b.ndim == 1:
        h = 0.5 * (a[np.newaxis,:] + b)
        d1 = a[np.newaxis,:] * np.log(a[np.newaxis,:]/h)
        posdef_check_value(d1)
        d1 = np.sum(d1, axis=2)
        d2 = b * np.log(b/h)
        posdef_check_value(d2)
        d2 = np.sum(d2, axis=2)
        d = 0.5 * (d1 + d2)
    else:
        h = 0.5 * (a[np.newaxis,:] + b[:, np.newaxis])
        d1 = a[np.newaxis,:] * np.log(a[np.newaxis,:]/h)
        posdef_check_value(d1)
        d1 = np.sum(d1, axis=2)
        d2 = b[:, np.newaxis] * np.log(b[:, np.newaxis]/h)
        posdef_check_value(d2)
        d2 = np.sum(d2, axis=2)
        d = 0.5 * (d1 + d2)
        #d = d.T
    return d


def KT(a, b):
    """ Compute Kendall_Tau distance
    """
    return(1 - Bio.Cluster.distancematrix((a,b), dist="k")[1][0])

def BC(a, b):
    """ Compute Bray-Curtis distance
    """
    return(pairwise_distances(a,b, metric='braycurtis',n_jobs=1)  )


def SC(a, b):
    """ Compute Spearman correlation distance
    """
    return(1 - spearmanr(a,b).correlation)

def weighted_rank(a, b):
    return 0
