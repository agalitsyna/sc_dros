import glob
from basic_utils import *

import numpy as np
import pandas as pd
experiment_ids = 'Cell1 Cell2 Cell3 Cell4 Cell5 Cell6 Cell7 Cell8 Cell9 Cell10 Cell11 Cell12 Cell13 Cell14 Cell15 Cell16 Cell17 Cell18 Cell19 Cell20 Cell21 Cell22 Cell23 Cell24 Cell25 Cell26 Cell27 Cell28 Cell29 Cell30 Cell31 Cell32 Cell33 Cell34 Cell35 Cell36 Cell37 Cell38 Cell39 Cell40 Cell41 Cell42 Cell43 Cell44 Cell45 Cell46 Cell47 Cell48 Cell49 Cell50 Cell51 Cell52 Cell53 Cell54 Cell55 Cell56 Cell57 Cell58 Cell59 Cell60 Cell61 Cell62 Cell63 Cell64 Cell65 Cell66 Cell67 Cell68 Cell69 Cell70 Cell71 Cell72 Cell73 Cell74 Cell75 Cell76 Cell77 Cell78 Cell79 Cell80'.split()

# Preliminary work, selection of ids: 
for exp in experiment_ids:
    fasta_list = glob.glob('../data/FASTQ/{}_*R1*.fastq'.format(exp))
    output = "../data/IDS/{}.ids.txt".format(exp)

    command = 'cat'
    for f in fasta_list:
        command += " <(awk 'NR % 4 == 1' {})".format(f)

    command += " | gawk '{{match($0, \"@(.+) \", a)}} {{print a[1]}}' > {}".format(output)

    run_command(command)

def create_selections(len_initial, start, end, step, add=False, replace=False):
    """
    Creation of random selections of indices
    """
    idx_full = np.arange(len_initial)
    toreturn = []

    assert (end-start)/step>=1

    if add:
        for i in range(start, end + step, step)[::-1]:
            lst = np.random.choice(idx_full, i, replace=replace)
            idx_full = lst.copy() #np.setdiff1d(idx_full, lst)
            toreturn.append(lst)

    else:
        for i in range(start, end+step, step):
            lst = np.random.choice(idx_full, i, replace=replace)
            toreturn.append(lst)

    toreturn.append(np.arange(len_initial))
    return(toreturn)

from multiprocessing import Pool
import os

nthreads = 10
niter = 10

for exp in experiment_ids:
    print(exp)

    lst_total_add = []
    lst_total_noadd = []
    
    pairsam_dct = {}
    idxfa_dct = {}
    
    file_fa = "../data/IDS/{}.ids.txt".format(exp)
    idxfa_dct[exp] = np.loadtxt(file_fa, dtype='S64')

    files_pairsam = glob.glob('../data/PAIR/{}_*.pairsam.JJ'.format(exp))
    exp_list = [x.split('/')[-1].split('.')[0] for x in files_pairsam]
    pairsam_dct[exp] = read_pairsams(files_pairsam, exp_list, exp)

    len_idxfa = len(idxfa_dct[exp])
    if len_idxfa > 5000000:
        step = 1000000
    else:
        step = 100000
        
    lst_add = [(exp, x, n) 
          for n in range(niter)
          for x in create_selections(len_idxfa, step, step*(len_idxfa//step), step, add=True, replace=False)
          ]
    lst_total_add += lst_add
    
    lst_noadd = [(exp, x, n) 
      for n in range(niter)
      for x in create_selections(len_idxfa, step, step*(len_idxfa//step), step, add=False, replace=False)
      ]
    lst_total_noadd += lst_noadd
    
    def run_filter(args):
        """
        Run selection of fasta indices specified by idx_selected list,
        querying of df_pairsam and filtering of unique contacts.
        
        Global parameters used (query by exp): idxfa_dct and pairsam_dct
        
        :param exp: experiment label
        :param idx_selected: numpy list (int) of selected numbers of indexes in fa
        :param niter: metainfo
        :return: stats
        """
        exp, idx_selected, niter = args
        idxfa = idxfa_dct[exp]
    
        df = pairsam_dct[exp].copy()
        idxpairsam = df.readID.values.astype('S64')
        df.index = idxpairsam
    
        idxfa_selected = idxfa[idx_selected].astype('S64')
        idxpairsam_selected = np.intersect1d(idxfa_selected, idxpairsam)
        df_selected = df.loc[idxpairsam_selected,:].reset_index(drop=True)
        df_filtered, stats = filter_pair_df(df_selected)
        
        stats['exp'] = exp
        stats['niter'] = niter
        stats['len_idxfa'] = len(idxfa)
        stats['len_idxfa_selected'] = len(idxfa_selected)
        stats['len_idxpairsam'] = len(idxpairsam)
        stats['len_idxpairsam_selected'] = len(idxpairsam_selected)
        
        return stats
    
    print(len(idxfa_dct[exp]))
    if len(idxfa_dct[exp])>30*1e6:
      res = []
      for i in lst_total_add:
        res.append(run_filter(i))
      res = pd.DataFrame(res)
    else:
      # Parallel analogue of run_filter([exp, idx_selected, 0])
      p = Pool(nthreads)
      res = pd.DataFrame(p.map(run_filter, lst_total_add))

    if not os.path.isfile('../data/TABLES/robustness_analysis_add.txt'):
        res.to_csv('../data/TABLES/robustness_analysis_add.txt')
    else:
        res.to_csv('../data/TABLES/robustness_analysis_add.txt', mode='a', header=False)
    
    if len(idxfa_dct[exp])>30*1e6:
      res = []
      for i in lst_total_noadd:
        res.append(run_filter(i))
      res = pd.DataFrame(res)
    else:
      p = Pool(nthreads)
      res = pd.DataFrame(p.map(run_filter, lst_total_noadd))
    
    if not os.path.isfile('../data/TABLES/robustness_analysis_noadd.txt'):
        res.to_csv('../data/TABLES/robustness_analysis_noadd.txt')
    else:
        res.to_csv('../data/TABLES/robustness_analysis_noadd.txt', mode='a', header=False)

    del lst_total_add
    del lst_total_noadd
    del pairsam_dct
    del idxfa_dct
