# Hybrid_Assembly
This tool is made as a part of my Master's thesis (2023)

## Instalation
The hybdrid assembly tool requires: Numpy, sqlite3, joblib, os, psutil, cython, sys, and setuptools
To download the pipeline please run
```
git clone https://github.com/Esmeetbdb/Hybrid_Assembly/tree/master
```
Then before running the pipeline enter the code directory and run:
```
python setup.py build_ext --inplace
```
to cythonize the code and improve performance.

## Run
To run the pipeline use the following command:
```
python main.py <Fasta> <Cmap> <Enzyme_site> <prefix> [options]
```
Where fasta is the path to the FASTA file containing the long-read fasta information, Cmap is the path to the file containing optical maps, Enzyme site is the enzymme recognition site used to generate optical maps, and prefix is the desired prefix for the output files.

Optional parameters are:
Optional argument | Function | Default
---|---|---
--n-threads | The number of threads that can be used for running the pipeline | 16
--k_mer | K-mer size used for overlapping fasta and optical mapping sites. Larger genomes will require larger k-mer sizes. Larger k-mer sizes lower sensitivity. | 5
--deviationlist | The max deviation applied to find overlaps for each round of alignment. Enter as a comma-separated string. Must be as long as overlap_len | 500,1000,2000
--overlap_len | The number of k-mers that have to overlap at least for an alignment to be considered correct. Enter as a comma-separated string. Must be as long as deviationlist | 6,8,10

For more information about running the tool use:
```
python main.py --help
```
