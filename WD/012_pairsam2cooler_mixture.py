"""
Converts pset of pairsam files for one cell into cooler files. 
TODO: optimize experiments handling.
Example run in bash: 
pref1="sc19"; pref2="sc24"; python 02_pairsam2cooler_mixture.py "../DATA/PAIR/${pref1}_*.pairsam" "$pref1" "../DATA/PAIR/${pref2}_*.pairsam" "$pref2"
"""

# importing modules
from sys import argv
import matplotlib as mpl
mpl.use('agg')
import glob
from basic_utils import *

# Reading command line parameters

mask1 = argv[1]
cell1 = argv[2]
mask2 = argv[3]
cell2 = argv[4]

print(mask1, cell1, mask2, cell2)

# Reading files and filtering
filelist1 = glob.glob(mask1)
exp_list1 = [x.split('/')[-1].split('.')[0] for x in filelist1]
df1 = read_pairsams(filelist1, exp_list1, cell1)
df_filtered1, stats = filter_pair_df(df1)

filelist2 = glob.glob(mask2)
exp_list2 = [x.split('/')[-1].split('.')[0] for x in filelist2]
df2 = read_pairsams(filelist2, exp_list2, cell2)
df_filtered2, stats = filter_pair_df(df2)

# Reading chromosomes
chr_sizes   = {x.split()[0]: int(x.split()[1]) for x in open("../DATA/GENOME/dm3.reduced.chrom.sizes", 'r').readlines()}
chromosomes = [x.split()[0] for x in open("../DATA/GENOME/dm3.reduced.chrom.sizes", 'r').readlines()]

# Creating mixture dataframe
df_filtered = pd.concat([df_filtered1.copy(), df_filtered2.copy()])
total_length = len(df_filtered)
total_length1 = len(df_filtered1)
total_length2 = len(df_filtered2)

# Subsampling
for n_sample in range(0, 10):

    for cell, cell_admixed, keep_num in [(cell1, cell2, total_length1), (cell2, cell1, total_length2)]:

      selected = np.random.choice(np.arange(total_length), keep_num, replace=False)
      df_output = df_filtered.copy().reset_index(drop=True).loc[selected,:].copy()
      
      resolutions = [20, 10, 1]
      create_cooler(df_output, "../DATA_SUBSAMPLE/PAIRIX_MIX/{}.admixed.{}.{}.pairsam".format(cell, cell_admixed, n_sample),
                    "../DATA_SUBSAMPLE/COOL_MIX/{}.admixed.{}.{}.{{}}.cool".format(cell, cell_admixed, n_sample),
                    "../DATA/GENOME/dm3.reduced.chrom.sizes",
                    resolutions_list=resolutions)

      # Writing text files
      res = 10
      cool = "../DATA_SUBSAMPLE/COOL_MIX/{}.admixed.{}.{}.{}.cool".format(cell, cell_admixed, n_sample, res)
      cmap = {0:'#241F20', 1:'#5C006D', 2:'#A30080', 3:'#FB5500', 4:'#FFB300', 5:'#F3F773'}
      for chrom in chromosomes:
        outfile = '../DATA_SUBSAMPLE/IMG/{}.admixed.{}.{}.{}.{}.png'.format(cell, cell_admixed, n_sample, res, chrom)
        cooler2png_chr(cool, outfile, cmap=cmap, chromosome=chrom, remove_diagonal=False, balance=False, scale='linear')
        outfile = '../DATA_SUBSAMPLE/TXT/MTX/{}.{}.admixed.{}.{}.{}.{}.txt'.format('mtx', cell, cell_admixed, n_sample, res, chrom)
        cooler2txt_chr(cool, outfile, fmt='mtx', chromosome=chrom)
        outfile = '../DATA_SUBSAMPLE/TXT/SPARSE/{}.{}.admixed.{}.{}.{}.{}.txt'.format('sparse_bins', cell, cell_admixed, n_sample, res, chrom)
        cooler2txt_chr(cool, outfile, fmt='sparse_bins', chromosome=chrom)
