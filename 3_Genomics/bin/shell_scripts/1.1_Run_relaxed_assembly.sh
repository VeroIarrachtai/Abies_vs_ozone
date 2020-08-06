#! /bin/sh

#SBATCH -w nodo3
#SBATCH -n 4

# Verónica Reyes
# Paper:
# 

# To use ipyrad in my cluster
export LD_LIBRARY_PATH=/opt/miniconda2/lib
export HDF5_USE_FILE_LOCKING=FALSE

# Run steps 1-7 in ipyrad

ipyrad -p ../../data/TMVB_5SNPrad.vcf -s 1234567 -f
