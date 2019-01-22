"""
Converts a set of pairsam files for a particular cell to cooler files. 
Example run in bash: 
 pref="A6"; python pairsam2cooler.py "../DATA/PAIR/${pref}*.pairsam" "$pref" 

Parameters: 
argv[1] - mask of input pairsam files
argv[2] - cell name (or prefix)
"""

# importing modules
from sys import argv
import glob
from basic_utils import *

# Reading command line parameters

mask          = argv[1]
cell          = argv[2]
out_pairix    = "../DATA/PAIRIX/{}.pairix".format(cell)
out_cool_mask = "../DATA/COOL/{}.{{}}.cool".format(cell)
out_stats     = "../DATA/STATS/{}.filter_stats".format(cell)

chromnames = ['chr4', 'chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R']

logging.debug(mask, cell, out_pairix, out_cool_mask)

filelist = glob.glob(mask)
exp_list = [x.split('/')[-1].split('.')[0] for x in filelist]

logging.debug(list(zip(filelist, exp_list)))

# Reading pairsam files into single dataframe
df = read_pairsams(filelist, exp_list, cell)

# Filtering dataframe by proper ligation junctions
df_filtered, stats = filter_pair_df(df)
stats['cell'] = cell
stats['exp'] = cell
#
## Writing filtering statistics
#with open(out_stats, 'w') as outf:
#    for k in sorted(stats.keys()):
#        outf.write('{}\t{}\n'.format(k, stats[k]))
#
## Writing output to cooler file
#resolutions = [100, 20, 10, 1]
#create_cooler(df_filtered, out_pairix, out_cool_mask,
#              "../DATA/GENOME/dm3.reduced.chrom.sizes", resolutions_list=resolutions)
#
## Writing coolers to images and txt files:
#cmap = {0:'#241F20', 1:'#5C006D', 2:'#A30080', 3:'#FB5500', 4:'#FFB300', 5:'#F3F773'}
#for res in [100, 20, 10]:
#  cool = "../DATA/COOL/{}.{}.cool".format(cell, res)
#  for chrom in chromnames:
#    outfile = '../DATA/IMG/{}.{}.{}.png'.format(cell, res, chrom)
#    cooler2png_chr(cool, outfile, cmap=cmap, chromosome=chrom, remove_diagonal=False, balance=False, scale='linear')
#    outfile = '../DATA/TXT/MTX/{}.{}.{}.{}.txt'.format('mtx', cell, res, chrom)
#    cooler2txt_chr(cool, outfile, fmt='mtx', chromosome=chrom)
#    outfile = '../DATA/TXT/SPARSE/{}.{}.{}.{}.txt'.format('sparse_bins', cell, res, chrom)
#    cooler2txt_chr(cool, outfile, fmt='sparse_bins', chromosome=chrom)



# Randomization of contacts

chr_sizes = {x.split()[0]: int(x.split()[1]) for x in
             open("../DATA/GENOME/dm3.reduced.chrom.sizes", 'r').readlines()}

df_input = df_filtered.query("chrom1==chrom2").reset_index(drop=True).copy()
positions1 = df_input.loc[:, "pos1"].values
positions2 = df_input.loc[:, "pos2"].values
chromosomes = df_input.loc[:, "chrom1"].values

# Incorporation of distance dependence
for mode, index_of_interest in [#('longdist', df_input.query("abs(pos1-pos2)>200000").index), 
                                ['rand', df_input.index]]:

  for n_sample in range(0, 10):
      print(n_sample, exp_list)
      df_output = df_input.copy()
  
      for i in index_of_interest:
          distance = positions2[i] - positions1[i]
          bgn = np.random.randint(0, chr_sizes[chromosomes[i]] - distance)
          end = bgn + distance
  
          df_output.loc[i, "pos1"] = bgn
          df_output.loc[i, "pos2"] = end
  
      resolutions = [20, 10, 1]
      create_cooler(df_output, "../DATA_SHUF/PAIRIX/{}.{}.{}.pairsam".format(cell, mode, n_sample),
                    "../DATA_SHUF/COOL/{}.{}.{}.{{}}.cool".format(cell, mode, n_sample),
                    "../DATA/GENOME/dm3.reduced.chrom.sizes",
                    resolutions_list=resolutions)
      del df_output
  
      # Writing text files
      #res = 10
      #cool = "../DATA_SHUF/COOL/{}.{}.{}.{}.cool".format(cell, mode, n_sample, res)
      #cmap = {0:'#241F20', 1:'#5C006D', 2:'#A30080', 3:'#FB5500', 4:'#FFB300', 5:'#F3F773'}
      #for chrom in chromnames:
      #  outfile = '../DATA_SHUF/IMG/{}.{}.{}.{}.{}.png'.format(cell, mode, n_sample, res, chrom)
      #  cooler2png_chr(cool, outfile, cmap=cmap, chromosome=chrom, remove_diagonal=False, balance=False, scale='linear')
      #  outfile = '../DATA_SHUF/TXT/MTX/{}.{}.{}.{}.{}.{}.txt'.format('mtx', cell, mode, n_sample, res, chrom)
      #  cooler2txt_chr(cool, outfile, fmt='mtx', chromosome=chrom)
      #  outfile = '../DATA_SHUF/TXT/SPARSE/{}.{}.{}.{}.{}.{}.txt'.format('sparse_bins', cell, mode, n_sample, res, chrom)
      #  cooler2txt_chr(cool, outfile, fmt='sparse_bins', chromosome=chrom)

## subsample

#for n_sample in range(0, 10):
#
#    df_output = df_filtered.copy()
#    total_length = len(df_output)
#
#    for keep_prc in np.arange(5, 105, 5):
#
#        keep_num = total_length*keep_prc//100
#
#        selected = np.random.choice(np.arange(total_length), keep_num, replace=False)
#        
#        df_output = df_filtered.copy().reset_index(drop=True)
#        df_output = df_output.loc[selected,:].copy()
#        
#        resolutions = [20, 10, 1]
#        create_cooler(df_output, "../DATA_SUBSAMPLE/PAIRIX/{}.subsample.{}.{}.pairsam".format(cell, keep_prc, n_sample),
#                      "../DATA_SUBSAMPLE/COOL/{}.{{}}.cool.subsample.{}.{}".format(cell, keep_prc, n_sample),
#                      "../DATA/GENOME/dm3.reduced.chrom.sizes",
#                      resolutions_list=resolutions)
#
#        # Writing text files
#        res = 10
#        cool = "../DATA_SUBSAMPLE/COOL/{}.{}.cool.subsample.{}.{}".format(cell, res, keep_prc, n_sample)
#        cmap = {0:'#241F20', 1:'#5C006D', 2:'#A30080', 3:'#FB5500', 4:'#FFB300', 5:'#F3F773'}
#        for chrom in chromnames:
#          outfile = '../DATA_SUBSAMPLE/IMG/{}.subsample.{}.{}.{}.png'.format(cell, keep_prc, n_sample, res, chrom)
#          cooler2png_chr(cool, outfile, cmap=cmap, chromosome=chrom, remove_diagonal=False, balance=False, scale='linear')
#          outfile = '../DATA_SUBSAMPLE/TXT/MTX/{}.subsample.{}.{}.{}.{}.{}.txt'.format('mtx', cell, keep_prc, n_sample, res, chrom)
#          cooler2txt_chr(cool, outfile, fmt='mtx', chromosome=chrom)
#          outfile = '../DATA_SUBSAMPLE/TXT/SPARSE/{}.subsample.{}.{}.{}.{}.{}.txt'.format('sparse_bins', cell, keep_prc, n_sample, res, chrom)
#          cooler2txt_chr(cool, outfile, fmt='sparse_bins', chromosome=chrom)
