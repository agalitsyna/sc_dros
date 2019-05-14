import pandas as pd
import glob
files = glob.glob("../DATA/TAD_modularity/*.csv")
df = pd.concat([pd.read_csv(f) for f in files])
df_tmp = df.query('optimum_kind!="optimum120"').drop(["mode", "Unnamed: 0"], axis=1).reset_index(drop=True)
df_tmp.columns = ['gamma', 'chromosome', 'experiment', 'mode', 'medianTADsize', 'meanTADsize', '#TADs', 'coverageTADs']
df_tmp.to_csv("../TABLES/gamma_values_for_TAD_calling.csv")
