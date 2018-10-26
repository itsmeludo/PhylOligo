#!/usr/bin/env python3

import os, sys, argparse
import numpy as np
import h5py

def read_numpy(path):
    return np.loadtxt(path)

def read_h5py(path):
    with h5py.File(path) as hf:
        mat = hf.get("distances")
        matval = mat.value[:]
    return matval

def read_memmap(path):
    matrix = np.memmap(path, dtype=np.float32, mode="r")
    s = matrix.shape[0]
    n = np.sqrt(s)
    if str(n).split(".")[1] != "0":
        print("Error, weird shape for matrix {}".format(path), file=sys.stderr)
        sys.exit(1)
    matrix = matrix.reshape((int(n), int(n)))
    return matrix

format2fn = {
    "numpy": read_numpy,
    "h5py": read_h5py,
    "memmap": read_memmap
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mat1", action="store", dest="matrix1")
    parser.add_argument("--format1", action="store", dest="format1", choices=["numpy", "memmap", "h5py"])
    parser.add_argument("--mat2", action="store", dest="matrix2")
    parser.add_argument("--format2", action="store", dest="format2", choices=["numpy", "memmap", "h5py"])
    params = parser.parse_args()
    
    mat1 = format2fn[params.format1](params.matrix1)
    mat2 = format2fn[params.format2](params.matrix2)
    print("matrix {}, shape: {}".format(params.matrix1, mat1.shape))
    print("matrix {}, shape: {}".format(params.matrix2, mat2.shape))
    print("Identical matrices?:", np.allclose(mat1, mat2, atol=1e-3))
    print()
    print(mat1)
    print()
    print(mat2)
    sys.exit(0)
    
if __name__ == "__main__":
    main()