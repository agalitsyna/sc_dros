import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')

import glob
from utils import * # from ../../lib/ folder

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression

def get_sums(mtx):
    """
    Calculate scaling and marginal distribution for 2D matrix
    """
    v_scaling = np.array([np.diagonal(mtx, i).sum() for i in range(len(mtx))])
    v_marg = np.sum(mtx, axis=0) + np.sum(mtx, axis=1)
    return(v_scaling, v_marg)

def df2mtx(pos_list, l):
    """
    Return 2D matrix from sparse positions list
    """
    mtx = np.zeros([l,l])
    for i, j in pos_list:
        mtx[i, j]+=1
    return mtx

def reconstruct_sparse(mtx, l_bp, binsize=10000):
    """
    Return sparse positions list from 2D matrix
        l_bp is the size of chromosome, needed for proper simulation of the last bin.
    """
    ret = []
    lim = l_bp % binsize
    for i in range(len(mtx)):
        for j in range(i+1, len(mtx)):
            for k in range(int(mtx[i,j])):
                randint1 = np.random.randint(15, binsize-15 if j<len(mtx)-1 else lim-15, 1)[0]
                randint2 = np.random.randint(15, binsize-15 if j<len(mtx)-1 else lim-15, 1)[0]
                ret.append([i*binsize+randint1, j*binsize+randint2])
    return (ret)

from sys import argv

mask        = argv[1]
pref = cell = argv[2]
output_path = argv[3]

out_pairix    = f"{output_path}/{cell}.pairix"
out_cool_mask = f"{output_path}/{cell}.{{}}.cool"
out_stats     = f"{output_path}/{cell}.filter_stats"
out_figure    = f"{output_path}/save_margi_{cell}_{{}}.pdf"


chromnames = ['chr4', 'chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R']

logging.debug("{} {} {} {}".format(mask, cell, out_pairix, out_cool_mask))

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


# Simulation preparation
binsize = 10000
df_filtered.loc[:, 'bin1'] = df_filtered.pos1//binsize
df_filtered.loc[:, 'bin2'] = df_filtered.pos2//binsize
df_filtered = df_filtered.query('chrom1==chrom2')

chrom_sizes = pd.read_csv('../data/GENOME/chrom.sizes_t', sep='\t', header=None, index_col=0)
chromnames = ['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R']

logging.debug("Starting preparation...")

# Smooth scaling for all chromosomes
scaling = np.bincount((df_filtered.bin2-df_filtered.bin1).values.astype(int), 
          minlength=np.max((chrom_sizes//binsize)[1].values.astype(int))+1 )

bins = np.arange(len(scaling))

x, y = np.log(bins*binsize), np.log(scaling)
mask = np.isfinite(x) & np.isfinite(y)
x, y = x[mask], y[mask]
X = x[:,np.newaxis]

polynomial_features = PolynomialFeatures(degree=10,
                                         include_bias=False)
linear_regression = LinearRegression()
pipeline = Pipeline([("polynomial_features", polynomial_features),
                     ("linear_regression", linear_regression)])
pipeline.fit(X, y)

X_new = np.log(bins*binsize)[1:][:,np.newaxis]
y_new = pipeline.predict(X_new)


# Figure 1
plt.figure(figsize=[7,7])
plt.plot(np.exp(X[:,0]), np.exp(y))
plt.plot(np.exp(X_new[:,0]), np.exp(y_new))

plt.xscale('log')
plt.yscale('log')
plt.savefig(out_figure.format("approx_scaling"))
###

scaling_approximated = np.concatenate( [[0], np.exp(y_new)/np.sum(np.exp(y_new))] )


# Figure 2
plt.figure(figsize=[7,7])
for ch in chromnames:
    l = ((chrom_sizes//binsize)[1].astype(int))[ch] + 1
    df_tmp = df_filtered.query('chrom1=="{}"'.format(ch))
    scaling_true = np.bincount((df_tmp.bin2-df_tmp.bin1).values.astype(int), minlength=l  )
    scaling_true = scaling_true/np.sum(scaling_true)

    plt.plot(scaling_approximated[0:l][1:], scaling_true[1:], lw=2, label=ch)
    
plt.plot([0,0.5], [0,0.5], '--', lw=2, color='grey')
plt.legend()
plt.savefig(out_figure.format("per_chr_scaling"))
###

logging.debug("Starting simulation...")
for n_sample in np.arange(0, 10):
    logging.debug("Simulation: {}".format(n_sample))

    df_final = pd.DataFrame()

    for ch in chromnames:

        l = ((chrom_sizes//binsize)[1].astype(int))[ch] + 1

        df_tmp = df_filtered.query('chrom1=="{}"'.format(ch)).query("bin1!=bin2").reset_index(drop=True)

        initial_marginal = np.bincount(df_tmp.bin1, minlength=l) + np.bincount(df_tmp.bin2, minlength=l)

        # Preserving approximated scaling 
        # Probabilities equal to marginals and scaling

        flag = 1
        it = 0
        while flag:
            flag = 0
            print('Iteration: {}, chromosome: {}'.format(it, ch))
            marginal = initial_marginal.copy()
            scaling  = scaling_approximated[0:l]
            mtx_simulated = np.zeros([l, l])

            for bin1 in np.arange(l):
                current_marginal = marginal[bin1:]
                current_scaling = scaling[:l-bin1]*(current_marginal>0).astype(int)
                prob = current_scaling/np.sum(current_scaling)

                counts = marginal[bin1]
                if np.sum(current_scaling>0)==0 and counts>0:
                    print('error {} {}', counts, bin1)
                    flag = 1
                    it += 1
                    break

                if counts>0:
                    if flag or (np.sum((prob>0)*marginal[bin1:])<counts):
                        flag = 1
                        it += 1
                        break
                    it_small=0
                    while True:
                        if it_small<1000:
                            idx_selected = np.random.choice(np.where(prob>0)[0], int(counts), p=prob[prob>0])+bin1
                        elif it_small<2000:
                            idx_selected = np.random.choice(np.where(prob>0)[0], int(counts))+bin1
                        else:
                            flag = 1
                            break
                            
                        upd_marginal = np.zeros(l)
                        for idx in idx_selected:
                            upd_marginal[idx] += 1
                        upd_marginal[bin1] += counts

                        if np.all(marginal-upd_marginal>=0):
                            break
                        if it_small>0:
                            print(it_small)
                        it_small +=1

                    marginal = marginal-upd_marginal
                    for bin2 in idx_selected:
                        mtx_simulated[bin1, bin2] += 1
        

        s1, m1 = get_sums(mtx_simulated)

        if np.sum( m1 != initial_marginal )>0:
            print("Sum is not equal to 0: {}".format(np.sum( m1 != initial_marginal )))


        # Figure 3
        plt.figure(figsize=[7,7])
        plt.plot(np.arange(0, l), s1/np.sum(s1), label="Simulated", alpha=0.5)
        plt.plot(np.arange(0, l), scaling/np.sum(scaling), label="Approximated", alpha=0.5)
        scaling_true = np.bincount((df_tmp.bin2-df_tmp.bin1).values.astype(int), minlength=l  )
        scaling_true = scaling_true/np.sum(scaling_true)
        plt.plot(np.arange(0, l), scaling_true/np.sum(scaling_true), label="Observed {}".format(ch), alpha=0.5)
        plt.xscale('log')
        plt.yscale('log')
        plt.legend()
        plt.savefig(out_figure.format("scaling_{}".format(ch)))


        # Figure 4 
        plt.figure(figsize=[7,5])
        plt.scatter(m1, initial_marginal, s=4)
        plt.xlabel("Simulated")
        plt.xlabel("Initial taget")
        plt.title("Marginals before and after simulation")
        plt.savefig(out_figure.format("compare_marinals_{}".format(ch)))

        # Figure 5
        plt.figure(figsize=[7,5])
        plt.scatter(s1, scaling, s=4)
        plt.xlabel("Simulated")
        plt.xlabel("Approximated taget")
        plt.title("Scalings before and after simulation")
        plt.savefig(out_figure.format("compare_scalings_{}".format(ch)))

        ###
        mtx = df2mtx(df_tmp.loc[:, ['bin1', 'bin2']].values, l)
        s2, m2 = get_sums(mtx)

        # Figure 6
        plt.figure(figsize=[10,10])
        sns.heatmap(mtx[0:200,0:200])
        plt.savefig(out_figure.format("map_{}_orig".format(ch)))
        plt.figure(figsize=[10,10])
        sns.heatmap(mtx_simulated[0:200,0:200])
        plt.savefig(out_figure.format("map_{}_sim".format(ch)))

        # Figure 7
        plt.figure(figsize=[5,5])
        sns.distplot(mtx_simulated[mtx_simulated>0], kde=False)
        sns.distplot(mtx[mtx>0], kde=False)
        plt.savefig(out_figure.format("counts_distr_{}".format(ch)))

        ###
        df_new = pd.DataFrame( reconstruct_sparse(mtx_simulated, ((chrom_sizes)[1].astype(int))[ch] ), columns=['pos1', 'pos2'])
        df_new.loc[:,'chrom1'] = ch
        df_new.loc[:,'chrom2'] = ch
        df_new.loc[:,'strand1'] = df_tmp.loc[:, 'strand1']
        df_new.loc[:,'strand2'] = df_tmp.loc[:, 'strand2']
        df_new.loc[:,'readID'] = df_tmp.loc[:, 'readID']

        df_final = pd.concat([df_final, df_new])

    resolutions = [20, 10, 1]
    mode = "withmarg"
    create_cooler(df_final, "../data_shuf/PAIRIX/{}.{}.{}.pairsam".format(cell, mode, n_sample),
                "../data_shuf/COOL/{}.{}.{}.{{}}.cool".format(cell, mode, n_sample),
                "../data/GENOME/dm3.reduced.chrom.sizes",
                resolutions_list=resolutions)

    del df_final

    res = 10
    cool = "../data_shuf/COOL/{}.{}.{}.{}.cool".format(cell, mode, n_sample, res)
    cmap = {0:'#241F20', 1:'#5C006D', 2:'#A30080', 3:'#FB5500', 4:'#FFB300', 5:'#F3F773'}
    for chrom in chromnames:
      outfile = '../data_shuf/IMG/{}.{}.{}.{}.{}.png'.format(cell, mode, n_sample, res, chrom)
      cooler2png_chr(cool, outfile, cmap=cmap, chromosome=chrom, remove_diagonal=False, balance=False, scale='linear')
      outfile = '../data_shuf/TXT/MTX/{}.{}.{}.{}.{}.{}.txt'.format('mtx', cell, mode, n_sample, res, chrom)
      cooler2txt_chr(cool, outfile, fmt='mtx', chromosome=chrom)
      outfile = '../data_shuf/TXT/SPARSE/{}.{}.{}.{}.{}.{}.txt'.format('sparse_bins', cell, mode, n_sample, res, chrom)
      cooler2txt_chr(cool, outfile, fmt='sparse_bins', chromosome=chrom)
"""
for n_sample in np.arange(0, 10):
    mode = "withmarg"
    res = 10
    cool = "../data_shuf/COOL/{}.{}.{}.{}.cool".format(cell, mode, n_sample, res)
    cmap = {0:'#FFFFFF', 1:'#5C006D', 2:'#A30080', 3:'#FB5500', 4:'#FFB300', 5:'#F3F773'}
    for chrom in chromnames:
      outfile = '../IMG/HEATMAPS/{}.{}.{}.{}.{}.png'.format(cell, mode, n_sample, res, chrom)
      cooler2png_chr(cool, outfile, cmap=cmap, chromosome=chrom, remove_diagonal=False, balance=False, scale='linear')
      #outfile = '../data_shuf/TXT/MTX/{}.{}.{}.{}.{}.{}.txt'.format('mtx', cell, mode, n_sample, res, chrom)
      #cooler2txt_chr(cool, outfile, fmt='mtx', chromosome=chrom)
      #outfile = '../data_shuf/TXT/SPARSE/{}.{}.{}.{}.{}.{}.txt'.format('sparse_bins', cell, mode, n_sample, res, chrom)
      #cooler2txt_chr(cool, outfile, fmt='sparse_bins', chromosome=chrom)
"""
