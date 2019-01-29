import glob
from basic_utils import *

import numpy as np
import pandas as pd

experiment_ids = 'A2 A6 A1 A3 A4 A5 A6 A7 A8 A9 B1 B10 B11 B12 B13 B14 B15 B16 B17 B18 B19 B2 B20 B21 B22 B23 B25 B26 B27 B29 B3 B30 B31 B32 B33 B34 B35 B36 B37 B38 B39 B4 B40 B41 B42 B43 B44 B45 B47 B48 B49 B5 B50 B51 B53 B54 B6 B7 B8 B9 sc1 sc10 sc11 sc12 sc13 sc14 sc15 sc16 sc17 sc18 sc19 sc2 sc20 sc21 sc22 sc23 sc24 sc25 sc26 sc27 sc28 sc29 sc32 sc34 sc35 sc36 sc4 sc7 sc9'.split()


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


def create_selections(len_initial, start, end, step, add=False, replace=False):
    """
    Creation of random selections of indices
    """
    idx_full = np.arange(len_initial)
    toreturn = []

    assert (end-start)/step>=1

    if add:
        for i in range(start, end + step, step):
            lst = np.random.choice(idx_full, i, replace=replace)
            idx_full = np.setdiff1d(idx_full, lst)

    else:
        for i in range(start, end+step, step):
            lst = np.random.choice(idx_full, i, replace=replace)
            toreturn.append(lst)

    return(toreturn)


nthreads = 10

from multiprocessing import Pool

if __name__ == '__main__':

    niter = 2
    idxfa_dct = {}
    pairsam_dct = {}
    for exp in experiment_ids:
        file_fa = "../DATA/IDS/{}.ids.txt".format(exp)
        idxfa_dct[exp] = np.loadtxt(file_fa, dtype='S64')

        files_pairsam = glob.glob('../DATA/PAIR/{}_*.pairsam.JJ'.format(exp))
        exp_list = [x.split('/')[-1].split('.')[0] for x in files_pairsam]
        pairsam_dct[exp] = read_pairsams(files_pairsam, exp_list, exp)

        len_idxfa = len(idxfa_dct[exp])
        if len_idxfa > 5000000:
            step = 1000000
        else:
            step = 100000

        lst = create_selections(len_idxfa, step, step*(end//step), step, add=False, replace=False)


    def run_filter(exp, idx_faidx, idx_selected):
        """
        Run selection of fasta indices specified by idx_selected list,
        querying of df_pairsam and filtering of unique contacts.
        :param df_pairsam: pandas dataframe after read_pairsams
        :param idx_faidx: numpy list (S64) of all indices of fa
        :param list_selected: numpy list (int) of selected numbers of indexes in fa
        :return: number of unique contacts and dataframe with statistics
        """
        assert len(idx_faidx) == len(idx_selected)

        df = pairsam_dct[exp].copy().reset_index(drop=True)
        df.index = df.readID.values.astype('S64')
        idx_pairsam = df.readID.values.astype('S64')
        df_filtered, stats = filter_pair_df(df.iloc[])

        return n, stats


    lst = [*x for x in create_selections() ]


    p = Pool(nthreads)

    print(p.map(f, [1, 2, 3]))

# Subsampling steps:

step = 100000
niter = 1

df_robustness = {'exp':[], 'iter':[], 'n_contacts_remained':[], 
                 'n_deleted_fasta':[], 'idx_todelete_overlap':[], 
                 'idx_total_pairsam':[], 'idx_total_fasta':[]}

for exp in experiment_ids:

    output = "../DATA/IDS/{}.ids.txt".format(exp)
    idx_fasta = np.loadtxt(output, dtype='S64')
    
    if len(idx_fasta)>5000000:
        step = 1000000
    else:
        step = 100000

    filelist = glob.glob('../DATA/PAIR/{}_*.pairsam.JJ'.format(exp))
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
            print(it, n, len(idx_remained))
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
        pd.DataFrame(df_robustness).to_csv("../DATA/TMP/robustness_anaylsys3.py")
        del idx_todelete_new
        del idx_remained
            
            
resulting_df = pd.DataFrame(df_robustness)

resulting_df.to_csv("../DATA/TMP/robustness_anaylsys3.py")
