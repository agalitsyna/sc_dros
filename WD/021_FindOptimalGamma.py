"""
The script to call optimal TAD segmentation based on pickle with all segmentations.
Inherited from 10_modularity_analysis_universal.ipynb

Input:
- experiment id (cell name)
- mode (just the keyword to recall the setup of settings)
- input pickle name
- folders for image and bed/juicer

Output: image, bed and juicer with optimal segmentation.
"""

from sys import argv
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from basic_utils import *

from sys import argv
infile = argv[1]
if len(argv)>2:
    exp = argv[2]
    mode = argv[3]
    OUTIMG = argv[4]
    OUTDIR = argv[5]

resolution = 10
res_bp = resolution*1000

# Reading input pickle

df = {'bgn': [], 'end': [], 'cell':[], 'state':[], 'length':[], 'ch':[], 'gamma':[], 'mode':[]}
segmentations_tmp_all = pickle.load(open(infile, 'rb'))

for gamma in segmentations_tmp_all.keys():
    segmentations_tmp = segmentations_tmp_all[gamma]
    chrms = list(segmentations_tmp.keys())
    for ch in chrms:
        segments = segmentations_tmp[ch]
        df['bgn'] += list(segments[:, 0])
        df['end'] += list(segments[:, 1])
        df['length'] += list(segments[:, 1] - segments[:, 0])
        df['cell'] += [exp for x in range(len(segments))]
        df['state'] += ['real' for x in range(len(segments))]
        df['ch'] += [ch for x in range(len(segments))]
        df['gamma'] += [gamma for x in range(len(segments))]
        df['mode'] += [mode for x in range(len(segments))]

df = pd.DataFrame(df)

# Determining resolution of pickle

all_gammas = np.sort(np.unique(df.gamma.values))
min_step = abs(all_gammas[1]-all_gammas[0])
epsilon = min_step/2
print(epsilon)

# Finding the best options

df_optimums = {x: [] for x in ['gamma', 'ch', 'exp', 'mode',
                               'optimum_kind', 'median',
                               'mean', 'num', 'cov']}
for ch in chrms:
    if ch == 'chrY':
        continue
    ### Selection of proper chromosome
    query = 'ch=="{}"'.format(ch)
    df_tmp = df.query(query)

    ### Grouping
    grp = df_tmp.groupby(['cell', 'state', 'gamma'])
    d_tmp = {
        'mean': grp.mean()['length'].reset_index(),
        'median': grp.median()['length'].reset_index(),
        'num': grp.count()['length'].reset_index(),
        'cov': grp.sum()['length'].reset_index()
    }

    ### Determing optimal params
    idx = np.argmax(d_tmp['num']['length'])
    opt_max_peak = {
        k: d_tmp[k].loc[idx, 'length'] for k in d_tmp.keys()
    }
    opt_max_peak['gamma'] = d_tmp['num'].loc[idx, 'gamma']

    idx = np.argmin(d_tmp['num'].query('gamma<{}'.format(opt_max_peak['gamma']))['length'])
    opt_min_peak = {
        k: d_tmp[k].loc[idx, 'length'] for k in d_tmp.keys()
    }
    opt_min_peak['gamma'] = d_tmp['num'].loc[idx, 'gamma']

    diff = np.abs(d_tmp['median']['length'].values - 12)
    idxs = np.where(np.abs(diff - min(diff)) < 0.0001)[0]
    idx = idxs[np.argmax(d_tmp['num']['length'].values[idxs])]
    opt_120 = {
        k: d_tmp[k].loc[idx, 'length'] for k in d_tmp.keys()
    }
    opt_120['gamma'] = d_tmp['median'].loc[idx, 'gamma']

    ### Populating the dictionary of optimals
    df_optimums['exp'].append(exp)
    df_optimums['ch'].append(ch)
    df_optimums['mode'].append(mode)
    df_optimums['optimum_kind'].append('optimumMaxPeak')
    for k in list(d_tmp.keys()) + ['gamma']:
        df_optimums[k].append(opt_max_peak[k])

    df_optimums['exp'].append(exp)
    df_optimums['ch'].append(ch)
    df_optimums['mode'].append(mode)
    df_optimums['optimum_kind'].append('optimumMinPeak')
    for k in list(d_tmp.keys()) + ['gamma']:
        df_optimums[k].append(opt_min_peak[k])

    #TODO: write optimum for the center of everything
    df_optimums['exp'].append(exp)
    df_optimums['ch'].append(ch)
    df_optimums['mode'].append(mode)
    df_optimums['optimum_kind'].append('optimum120')
    for k in list(d_tmp.keys()) + ['gamma']:
        df_optimums[k].append(opt_120[k])

    ### Plotting
    plt.figure(figsize=[15, 4])
    for k in d_tmp.keys():
        plt.plot(d_tmp[k]['gamma'],
                 d_tmp[k]['length'] / max(d_tmp[k]['length']),
                 label='{}, max: {:.1f}'.format(k, max(d_tmp[k]['length'])))

    selected_gammas = [
        opt_min_peak['gamma'] + i * (opt_max_peak['gamma'] - opt_min_peak['gamma']) / 6 for i in
        [0, 1, 2, 3, 4, 5, 6]
    ]
    for a in selected_gammas:
        plt.plot([a, a], [0, 1])

    plt.legend()
    plt.xlim([0, 350])

    plt.title("Cell: {} Mode: {} Chromosome: {}".format(exp, mode, ch))

    plt.savefig(
        OUTIMG+"5_gamma-selection_{}_{}_{}.pdf".format(exp, ch, mode))

    ### Saving segmentations

    writing_mode = "w" if ch==chrms[0] else "a"
    print(exp, ch, selected_gammas)

    for i, a in enumerate(selected_gammas):
       df_towrite = df_tmp.query("abs(gamma-{:.1f})<{:.2f}".format(a, epsilon))
       print(len(df_towrite), i, a, exp, mode, ch)
       segmentations_to_bed( df_towrite,
                        OUTDIR+"segmentation_{}_{}_{:.2f}.bed".format(exp, mode, i/6),
                        res_bp, "df", writing_mode)
       segmentations_to_2Djuice( df_towrite,
                        OUTDIR+"segmentation_{}_{}_{:.2f}.2Dannot".format(exp, mode, i/6),
                        res_bp, "df", writing_mode)

df_optimums = pd.DataFrame(df_optimums)

df_optimums.to_csv(OUTDIR+"optimal_gammas_{}_{}.csv".format(exp, mode))
