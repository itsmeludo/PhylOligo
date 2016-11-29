#!/usr/bin/env python3

from scoop import futures
  
import os, sys
import subprocess, shlex, shutil
import tempfile
import h5py

from phyloligo.core.phyloutils import remove_folder

## clustering
import numpy as np
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.metrics.pairwise import check_pairwise_arrays
from sklearn.utils import gen_even_slices
from sklearn.externals.joblib import Parallel, delayed
from sklearn.externals.joblib import dump, load

def posdef_check_value(d):
    d[np.isnan(d)]=0    
    d[np.isinf(d)]=0
    
#### Distance functions

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

def euclidean_distances_loc(output, X, s):
    distances = euclidean_distances(X, X[s])
    output[s] = distances.T
    
def JSD_loc(output, X, s):
    X, Y = check_pairwise_arrays(X, X[s])
    d = JSD(X, Y)
    output[s] = d

def compute_JSD_unpack(params):
    """ unpack parameters to compute distances
    """
    i, j, freqi, freqj = params
    d = JSD(freqi, freqj)
    return i, j, d

def compute_Eucl_unpack(params):
    """ unpack parameters to compute distances
    """
    i, j, freqi, freqj = params
    d = Eucl(freqi, freqj)
    return i, j, d

def euclidean_distances_h5py(output_dir, input, s):
    with h5py.File(input, "r") as hf:
        X = hf.get("frequencies")
        dist = euclidean_distances(X, X[s]).T
    
    output = os.path.join(output_dir, "distance_{}_{}".format(s.start, s.stop))
    with h5py.File(output, "w") as hf:
        distances = hf.create_dataset("distances", dist.shape, dtype="float32")
        distances[...] = dist[:]
    

def JSD_h5py(output_dir, input, s):
    with h5py.File(input, "r") as hf:
        X = hf.get("frequencies")
        X, Y = check_pairwise_arrays(X, X[s])
        d = JSD(X, Y)
    
    output = os.path.join(output_dir, "distance_{}_{}".format(s.start, s.stop))
    with h5py.File(output, "w") as hf:
        distances = hf.create_dataset("distances", dist.shape, dtype="float32")
        distances[...] = d[:]
    #hf = h5py.File(output, "r+")
    #distances = hf.get("distances")
    #distances[s] = d
    #hf.close()
    

#### different flavor of parallelisms for distance computation
    
def compute_distances_pickle(frequencies, chunksize, metric="Eucl", n_jobs=1, workdir="/tmp/"):
    """ compute pairwises distances
    
    Parameters
    ----------
    frequencies: np.array
        a list of frequencies, each row corresponding to a sample each column to a "word"
    metric: string
        distance method to use ('JSD', 'Eucl', 'KL')
    n_jobs: int
        number of parallel job to execute
    
    Return:
    -------
    distances: np.array
        (n_samples, n_samples) distance matrix
    """
    #scoop.logger.info("Starting distance computation")
    if metric == "Eucl":
        pathin = tempfile.mktemp(dir=workdir)
        pathout = tempfile.mktemp(dir=workdir)
    
        distances = np.zeros((len(frequencies), len(frequencies)), dtype=float)
        for freqchunk in make_freqchunk(frequencies, chunksize):
            with open(pathin, "wb") as outf:
                pickle.dump(freqchunk, outf)
            # subprocess the computation
            cmd = "python3 -m scoop -n {} phylo_batchdist.py {} {} {}".format(n_jobs, pathin, pathout, metric)
            cmd = shlex.split(cmd)
            try:
                res = subprocess.check_call(cmd)
            except:
                print("Error running phylo_batchdist.py on {} {} /{}".format(pathin, pathout, metric), file=sys.stderr)
                sys.exit(1)
                
            with open(pathout, "rb") as inf:
                res = pickle.load(inf)
            #res = futures.map(compute_Eucl_unpack, freqchunk)
            for i, j, d in res:
                distances[i, j] = distances[j, i] = d
        os.remove(pathin)
        os.remove(pathout)
    elif metric == "JSD":
        pathin = tempfile.mktemp(dir=workdir)
        pathout = tempfile.mktemp(dir=workdir)
    
        distances = np.zeros((len(frequencies), len(frequencies)), dtype=float)
        for freqchunk in make_freqchunk(frequencies, chunksize):
            with open(pathin, "wb") as outf:
                pickle.dump(freqchunk, outf)
            # subprocess the computation
            cmd = "python3 -m scoop -n {} phylo_batchdist.py {} {} {}".format(n_jobs, pathin, pathout, metric)
            cmd = shlex.split(cmd)
            try:
                res = subprocess.check_call(cmd)
            except:
                print("Error running phylo_batchdist.py on {} {} /{}".format(pathin, pathout, metric), file=sys.stderr)
                sys.exit(1)    
            with open(pathout, "rb") as inf:
                res = pickle.load(inf)
            #res = futures.map(compute_JSD_unpack, freqchunk)
            for i, j, d in res:
                distances[i, j] = distances[j, i] = d
        os.remove(pathin)
        os.remove(pathout)
    else:
        print("Error, unknown method {}".format(metric), file=sys.stderr)
        sys.exit(1)
    return distances
    
def compute_distances_scoop(frequencies, chunksize, metric="Eucl", localrun=False, n_jobs=1):
    """ compute pairwises distances
    
    Parameters
    ----------
    frequencies: np.array
        a list of frequencies, each row corresponding to a sample each column to a "word"
    metric: string
        distance method to use ('JSD', 'Eucl', 'KL')
    n_jobs: int
        number of parallel job to execute
    
    Return:
    -------
    distances: np.array
        (n_samples, n_samples) distance matrix
    """
    #scoop.logger.info("Starting distance computation")
    if metric == "Eucl":
        if localrun:
            distances = pairwise_distances(frequencies, metric=Eucl, n_jobs=n_jobs)
        else:
            distances = np.zeros((len(frequencies), len(frequencies)), dtype=float)
            for freqchunk in make_freqchunk(frequencies, chunksize):
                res = futures.map(compute_Eucl_unpack, freqchunk)
                for i, j, d in res:
                    distances[i, j] = distances[j, i] = d
    elif metric == "JSD":
        if localrun:
            distances = pairwise_distances(frequencies, metric=JSD, n_jobs=n_jobs)
        else:
            distances = np.zeros((len(frequencies), len(frequencies)), dtype=float)
            for freqchunk in make_freqchunk(frequencies, chunksize):
                res = futures.map(compute_JSD_unpack, freqchunk)
                for i, j, d in res:
                    distances[i, j] = distances[j, i] = d
    
    else:
        print("Error, unknown method {}".format(metric), file=sys.stderr)
        sys.exit(1)
    return distances


def compute_distances_memmap(frequencies, freq_name, output, metric="Eucl", n_jobs=1):
    """ compute pairwises distances
    
    Parameters
    ----------
    freq_name: string
        path of the local frequency array
    metric: string
        distance method to use ('JSD', 'Eucl', 'KL')
    n_jobs: int
        number of parallel job to execute
    
    Return:
    -------
    distances: np.array
        (n_samples, n_samples) distance matrix
    """
    #scoop.logger.info("Starting distance computation")
    #folder = tempfile.mkdtemp()
    #dist_name = os.path.join(folder, output)

    # Pre-allocate a writeable shared memory map as a container for the
    # results of the parallel computation
    distances = np.memmap(output, dtype=frequencies.dtype, shape=(frequencies.shape[0], frequencies.shape[0]), mode='w+')
    
    # close and reopen it for reference
    #dump(distances, distname)
    #distances = load(dist_name, mmap_mode = "r+")
    del distances
    distances = np.memmap(output, dtype=frequencies.dtype, shape=(frequencies.shape[0], frequencies.shape[0]), mode='r+')

    if metric == "Eucl":
        # execute parallel computation of euclidean distance
        fd = delayed(euclidean_distances_loc)
        Parallel(n_jobs=n_jobs, verbose=0)(fd(distances, frequencies, s) for s in gen_even_slices(frequencies.shape[0], n_jobs))
            
        folder = os.path.dirname(freq_name)
        remove_folder(folder)
        print(freq_name)
        
    elif metric == "JSD":
        fd = delayed(JSD_loc)
        Parallel(n_jobs=n_jobs, verbose=0)(fd(distances, frequencies, s) for s in gen_even_slices(frequencies.shape[0], n_jobs))
        
        folder = os.path.dirname(freq_name)
        remove_folder(folder)
            
    else:
        print("Error, unknown method {}".format(metric), file=sys.stderr)
        sys.exit(1)

def join_distance_results(dist_folder, dist_name, size, n_jobs):
    """ join the different distances computation into a single matrix
    
    Parameters
    ----------
    dist_folder: string
        path to the distance folder storing all distances computation
    dist_name: string
        path to the final matrix
    size: int
        size of the final matrix
    n_jobs: int
        number of jobs used, the names of the files storing distances 
        depend on the number of jobs and on the size of the matrix
    """
    with h5py.File(dist_name, "w") as hf:
        distances = hf.create_dataset("distances", (size, size), dtype="float32")
        for s in gen_even_slices(size, n_jobs):
            pathin = os.path.join(dist_folder, "distance_{}_{}".format(s.start, s.stop))
            with h5py.File(pathin, "r") as inf:
                distances[s] = inf.get("distances").value[:]


def compute_distances_h5py(freq_name, dist_name, metric="Eucl", n_jobs=1):
    """ compute pairwises distances
    
    Parameters
    ----------
    freq_name: string
        path of the local frequency array
    metric: string
        distance method to use ('JSD', 'Eucl', 'KL')
    n_jobs: int
        number of parallel job to execute
    
    Return:
    -------
    distances: np.array
        (n_samples, n_samples) distance matrix
    """
    #scoop.logger.info("Starting distance computation")
    #folder = tempfile.mkdtemp()
    #dist_name = os.path.join(folder, output)
    
    folder = os.path.dirname(freq_name)
    
    with h5py.File(freq_name, "r") as hf:
        frequencies = hf.get("frequencies")
        size = frequencies.shape[0]
    
    dist_folder = os.path.join(folder, "distances")
    if not os.path.isdir(dist_folder):
        os.makedirs(dist_folder)
    
    if metric == "Eucl":
        # execute parallel computation of euclidean distance
        fd = delayed(euclidean_distances_h5py)
        Parallel(n_jobs=n_jobs, verbose=0)(fd(dist_folder, freq_name, s) for s in gen_even_slices(size, n_jobs))
        #for s in gen_even_slices(size, n_jobs):
            #euclidean_distances_h5py(dist_folder, freq_name, s)
            
        #print(freq_name)
        remove_folder(folder)
        
    elif metric == "JSD":
        fd = delayed(JSD_h5py)
        Parallel(n_jobs=n_jobs, verbose=0)(fd(dist_folder, freq_name, s) for s in gen_even_slices(size, n_jobs))
        
        remove_folder(folder)
            
    else:
        print("Error, unknown method {}".format(metric), file=sys.stderr)
        sys.exit(1)
        
    # join distance results
    join_distance_results(dist_folder, dist_name, size, n_jobs)
    

def compute_distances(mthdrun, large, frequencies, freq_name, out_file, dist, threads_max, freqchunksize, workdir):
    """ choose which function to call to compute distances
    """
    if mthdrun == "joblib":
        if large == "memmap":
            res = compute_distances_memmap(frequencies, freq_name, out_file, metric=dist, n_jobs=threads_max)
        if large == "h5py":
            res = compute_distances_h5py(freq_name, out_file, metric=dist, n_jobs=threads_max)
        else:
            res = compute_distances_scoop(frequencies, freqchunksize, metric=dist, localrun=True, n_jobs=threads_max)
    elif mthdrun == "scoop1":
        res = compute_distances(frequencies, freqchunksize, metric=dist, localrun=False, n_jobs=threads_max)
    elif mthdrun == "scoop2":
        res = compute_distances_pickle(frequencies, freqchunksize, metric=dist, n_jobs=threads_max, workdir=workdir)
    else:
        print("Error, method {} is not implemented for pairwise distances computation".format(mthdrun), file=sys.stderr)
        
    return res

