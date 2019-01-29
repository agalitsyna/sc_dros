"""
Python script that calculates the TADs  with parameters optimizing the expected TADs size.

Usage:
python 04_tads_calling.py <cooler input file> <prefix for output with TADs> <method name: modularity, armatus> <retrieve matrix: balanced, non-balanced> <expected TADs size> <max size of interTAD>

Example usage:
python 04_tads_calling.py ../DATA/COOL/A6.10.cool ../DATA/TAD/tmp modularity non-balanced 12 3

Output files:
<prefix>.bed bed file with end coordinate -- 0.5*resolution to produce visible results in juicebox with this file
<prefix>.2Dannot -- 2Dannotation file for juicebox
<prefix>.opt_gamma.pickle -- pickle file with dictionary with optimal gammas
<prefix>.opt_segmentation.pickle-- pickle file with dictionary with optimal segmentation

"""

from basic_utils import *
import numpy as np
import pickle

from sys import argv

input_cooler      = argv[1]
output_prefix     = argv[2]
calling_algorithm = argv[3] # armatus of modularity
reading_mode      = argv[4] # balanced or not
tad_size          = int(argv[5]) # median expected TADs size in bins
intertad_size     = int(argv[6]) # max length for segmentation unit to be considered as interTAD
max_tad_size      = 10000 #int(argv[7]) # min length for segmentation unit to be considered as interTAD (too large TAD)

output_bed        = output_prefix + '.bed'
output_2Dannot    = output_prefix + '.2Dannot'
output_segmentation = output_prefix + '.segmentation_all.pickle'

c = cooler.Cooler(input_cooler)
if reading_mode=='balanced':
    balance = True
else:
    balance = False
    
if calling_algorithm == 'modularity':
    step1 = 0.1
    mx1 = 375

elif calling_algorithm == 'armatus':
    step1 = 0.1
    mx1 = 10

gammas = np.arange(step1, mx1, step1)
segmentations_all = {}

for gamma in gammas:
    segmentations = {}
    segmentations_all[gamma] = {}
    for chrom in c.chromnames:

      if chrom=='chrM':
          continue

      mtx = c.matrix(balance=balance).fetch('{0}'.format(chrom)).astype(float)
      mtx[np.isnan(mtx)] = 0
      np.fill_diagonal(mtx, 0)

      if calling_algorithm == 'armatus':
          
          mn = np.percentile(mtx[mtx>0], 1)
          mx = np.percentile(mtx[mtx>0], 99)

          mtx[mtx<=mn] = mn
          mtx[mtx>=mx] = mx
          
          mtx = np.log(mtx)
          mtx = mtx-np.min(mtx)

      segments = produce_segmentation(mtx, gamma, method=calling_algorithm, 
                                  max_intertad_size=intertad_size, max_tad_size=max_tad_size)
      segmentations[chrom] = segments.copy()
      segmentations_all[gamma][chrom] = segments.copy()
    
    #segmentations_to_bed(    segmentations,     output_bed     + '.{:.3f}'.format(gamma),     c.binsize)
    #segmentations_to_2Djuice(segmentations,     output_2Dannot + '.{:.3f}'.format(gamma),     c.binsize)

pickle.dump(segmentations_all, open(output_segmentation, 'wb'))
print("Processed:", input_cooler)
