#!/usr/bin/env python3

import os, sys
from setuptools import setup, find_packages


def readme():
    with open('README.md') as f:
        return f.read()

script_dir = os.path.join("phylopackage", "bin")
all_scripts = list()
for binf in ["phyloligo.py", "phyloligo_comparemat.py", "phyloselect.py",
             "Kount.py", "phylopreprocess.py", "phyloselect.R",
             "contalocate.R"]: # add script names here
    all_scripts.append(os.path.join(script_dir, binf))

setup(name='phyloligo',
    version='v1.0',
    description='clustering of sequence reads based on oligonucleotide composition',
    long_description=readme(),

    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3',
        'Topic :: Text Processing :: Linguistic',
        ],
    keywords='oligonucleotide reads clustering bioinformatics',
    url='https://github.com/itsmeludo/PhylOligo/',
    author='Ludovic Mallet, Tristan Bitard-Feildel',
    author_email='ludovic.mallet@anses.fr, tristan.bitard-feildel@impmc.upmc.fr',
    license='GPLv3',
    scripts=all_scripts,
    packages=find_packages(exclude=[script_dir,]),
    include_package_data=True,
    install_requires=['biopython>=1.68', 
                      'scikit-learn>=0.19.1',
                      'numpy>=1.11.2',
                      'cython>=0.25.1',
                      'hdbscan>=0.8.2',
                      'matplotlib',
                      'scoop',
                      'h5py',
                      ],
     data_files=[
                ('data/references', [
                        'data/references/M.oryzae_TH12_Eucl_ref.mat.xz',
                        'data/references/M.oryzae_TH12_JSD_ref.mat.xz',
                          ]),
                ('doc', ['doc/PhylOligo_doc.pdf']),
                ('benchmark', ['doc/Benchmark.tgz']),
		],

    zip_safe=False,
)
