import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

class PlotRegion:
    @staticmethod
    def hist_region_len(df, specie):
        sdf = df[(df['specie']==specie) & (df['region_name'].str.contains('IMGT'))]
        sdf = sdf.reset_index()
        sdf['region_name'] = sdf['region_name'].map(lambda x: x.replace('-IMGT', ''))

        fig, ax = plt.subplots(1,3, figsize=(10,3), layout='tight')
        fig.suptitle("Fragment Length of V-regions")
        fig.supxlabel('Length of Amino Acids')
        fig.supylabel('Region')

        i=0
        hdf = sdf[sdf['chain_type']=='H']
        sns.histplot(hdf, x='seq_len', y='region_name', bins=50, binwidth=2, ax=ax[i])
        ax[i].set_title(f"Heavy Chain, {len(hdf)}")
        ax[i].set_xlabel(None)
        ax[i].set_ylabel(None)

        i=1
        kdf = sdf[sdf['chain_type']=='K']
        sns.histplot(kdf, x='seq_len', y='region_name', bins=50, binwidth=2, ax=ax[i])
        ax[i].set_title(f"Kappa Chain,  {len(kdf)}")
        ax[i].set_xlabel(None)
        ax[i].set_ylabel(None)

        i=2
        ldf = sdf[sdf['chain_type']=='L']
        sns.histplot(ldf, x='seq_len', y='region_name', bins=50, binwidth=2, ax=ax[i])
        ax[i].set_title(f"Lambda Chain, {len(ldf)}")
        ax[i].set_xlabel(None)
        ax[i].set_ylabel(None)

        plt.show()