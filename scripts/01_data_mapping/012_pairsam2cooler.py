"""
Converts a set of pairsam files for a particular cell into cooler file.
Example run in bash:
 pref="A6"; python pairsam2cooler.py "../DATA/PAIR/${pref}*.pairsam" "$pref" "../data/cool/"

Parameters:
argv[1] - mask of input pairsam files
argv[2] - cell name (or prefix)
"""

from utils import * # from ../../lib/ folder
from sys import argv

mask        = argv[1]
pref = cell = argv[2]
output_path = argv[3]

output_cool = f"{output_path}/{pref}.{{}}.cool" # A mask for writing cool files
output_pairix = f"{output_path}/{pref}.pairix"
output_stats = f"{output_path}/{pref}.stats"

chromnames = ['chr4', 'chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R']

logging.debug(f"Running pairsam2cool for: {mask} {pref}\n Writing into {output_path}/{pref} files")

filelist = glob.glob(mask)
exp_list = [x.split('/')[-1].split('.')[0] for x in filelist]

print(list(zip(filelist, exp_list)))

logging.debug("Reading dataframe...")

# Reading pairsam files into single dataframe
df = read_pairsams(filelist, exp_list, cell)

# Filtering dataframe by proper ligation junctions
df_filtered, stats = filter_pair_df(df)
stats['cell'] = cell
stats['exp'] = cell

# Writing filtering statistics
with open(output_stats, 'w') as outf:
   for k in sorted(stats.keys()):
       outf.write('{}\t{}\n'.format(k, stats[k]))

# Writing output to cooler file
create_cooler(df_filtered, output_pairix, output_cool,
              "../../data/GENOME/dm3.reduced.chrom.sizes",
              resolutions_list=[10])
