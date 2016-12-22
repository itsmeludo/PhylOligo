#!/usr/bin/env python3

          
from setuptools import setup, find_packages

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='phyloligo',
    version='0.1',
    description='clustering of sequence reads based on oligonucleotide composition',
    long_description=readme(),

    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Topic :: Text Processing :: Linguistic',
        ],
    keywords='oligonucleotide reads clustering bioinformatics',
    url='https://github.com/itsmeludo/PhylOligo/',
    author='Ludovic Mallet, Tristan Bitard-Feildel',
    author_email='ludovic.mallet@inra.fr, tristan.bitard-feildel@impmc.upmc.fr',
    license='MIT',
    scripts=['phyloligo.py', 'phyloselect.py',
             'phyloligo_comparemat.py', 
             ],
    #packages=find_packages(exclude=["phyloligo/bin/"]),
    include_package_data=True,
    install_requires=['biopython>=1.68', 
                      'scikit-learn>=0.18.1',
                      'numpy>=1.11.2',
                      'cython>=0.25.1',
                      'hdbscan>=0.8.2',
                      'matplotlib',
                      ],
    zip_safe=False,
)
