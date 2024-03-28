"""
This file is to be called for shortening the memory size of the data file Results.txt.
It can compress the file size from around 1 GB to around 150 MB.
Run as: python3 data_process_npz.py ../data/foldername/ filename.txt True (folder location)
"""
import numpy as np
from os.path import join as pjoin
import os.path
import sys

data_dir = sys.argv[1]
file_name = "result.txt"

if len(sys.argv)>1:
    all_data = sys.argv[2]
else:
    all_data = False

file_path = pjoin(data_dir, file_name)

if not os.path.exists(file_path):
    print(f"File {file_path} not found.")
    sys.exit(1)

# Load the single file
try:
    #x, nde, ndi, ndn, ndb, phi, EF = np.loadtxt(file_path, unpack=True)
    x, nde, ndi,phi, EF = np.loadtxt(file_path, unpack=True)
except ValueError as e:
    print(f"Error loading data from {file_path}: {e}")
    sys.exit(1)

if all_data:
    pro_data = 'processed_results_all.npz'
    #np.savez_compressed(pjoin(data_dir, pro_data), x=x, nde=nde, ndi=ndi, ndn=ndn, ndb = ndb, phi=phi, EF=EF)
    np.savez_compressed(pjoin(data_dir, pro_data), x=x, nde=nde, ndi=ndi, phi=phi, EF=EF)
#else:
    #pro_data = 'processed_results_E.npz'
    #np.savez_compressed(pjoin(data_dir, pro_data), x=x, EF=EF)
#print(f'Processed data saved to {pjoin(data_dir, pro_data)}')
