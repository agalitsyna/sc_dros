import glob
from basic_utils import *

import numpy as np
import pandas as pd

experiment_ids = 'A1 A2 A3 A4 A5 A6 A7 A8 A9 B1 B10 B11 B12 B13 B14 B15 B16 B17 B18 B19 B2 B20 B21 B22 B23 B25 B26 B27 B29 B3 B30 B31 B32 B33 B34 B35 B36 B37 B38 B39 B4 B40 B41 B42 B43 B44 B45 B47 B48 B49 B5 B50 B51 B53 B54 B6 B7 B8 B9 sc1 sc10 sc11 sc12 sc13 sc14 sc15 sc16 sc17 sc18 sc19 sc2 sc20 sc21 sc22 sc23 sc24 sc25 sc26 sc27 sc28 sc29 sc32 sc34 sc35 sc36 sc4 sc7 sc9'.split()


# Preliminary work: 
#for exp in experiment_ids:
#    fasta_list = glob.glob('../DATA/FASTQ/{}_*R1*.fastq'.format(exp))
#    output = "../DATA/IDS/{}.ids.txt".format(exp)
#
#    command = 'cat'
#    for f in fasta_list:
#        command += " <(awk 'NR % 4 == 1' {})".format(f)
#
#    command += " | gawk '{{match($0, \"@(.+) \", a)}} {{print a[1]}}' > {}".format(output)
#
#    run_command(command)


# Subsampling steps:

step = 100000
niter = 1

df_robustness = {'exp':[], 'iter':[], 'n_contacts_remained':[], 
                 'n_deleted_fasta':[], 'idx_todelete_overlap':[], 
                 'idx_total_pairsam':[], 'idx_total_fasta':[]}

for exp in experiment_ids:

    output = "../DATA/IDS/{}.ids.txt".format(exp)
    idx_fasta = np.loadtxt(output, dtype='S64')

    filelist = glob.glob('../../examples_single_cell_2018/DATA/PAIR/{}_*.pairsam.JJ'.format(exp))
    exp_list = [x.split('/')[-1].split('.')[0] for x in filelist]

    df = read_pairsams(filelist, exp_list, exp)
    df_filtered, stats = filter_pair_df(df)
    idx_pairsam = df_filtered.readID.values.astype('S64')
    
    
    total_length_pairsam = len(idx_pairsam)
    total_length_fasta = len(idx_fasta)
    total_steps = total_length_fasta//step
    
    print(exp, total_steps, len(idx_fasta))

    for it in range(niter):
        n_todelete_count = 0
        idx_todelete_new = []
        idx_remained = idx_fasta.copy()
        for n in range(total_steps):
            idx_remained = np.setdiff1d(idx_remained, idx_todelete_new)
            if len(idx_remained)<=step:
                break
            idx_todelete = np.random.choice(idx_remained, step, replace=False)
            idx_todelete_new = np.concatenate([idx_todelete, idx_todelete_new])

            n_deleted = len(np.intersect1d(idx_pairsam, idx_todelete_new))
            df_robustness['exp'] += [exp]
            df_robustness['iter'] += [it]
            df_robustness['n_contacts_remained'] += [total_length_pairsam-n_deleted]
            df_robustness['n_deleted_fasta'] += [n_todelete_count]
            df_robustness['idx_todelete_overlap'] += [n_deleted]
            df_robustness['idx_total_pairsam'] += [total_length_pairsam]
            df_robustness['idx_total_fasta'] += [total_length_fasta]
            n_todelete_count += step
            
            
resulting_df = pd.DataFrame(df_robustness)

resulting_df.to_csv("../DATA/TMP/robustness_anaylsys2.py")
