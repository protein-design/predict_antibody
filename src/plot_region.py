import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

class PlotRegion:
    def __init__(self, df:pd.DataFrame, verbose:bool=True):
        self.df = df
        self.verbose = verbose

    def hist_region_len(self, ax, params:dict):
        specie = params['specie']
        chain_type = params['chain_type']
        sdf = self.df[(self.df['specie']==specie) & (self.df['region_name'].str.contains('IMGT'))]
        sdf = sdf.reset_index()
        sdf['region_name'] = sdf['region_name'].map(lambda x: x.replace('-IMGT', ''))
        hdf = sdf[sdf['chain_type']==chain_type]

        # plots
        sns.histplot(hdf, x='seq_len', y='region_name', ax=ax, color='black',
            bins=params.get('bins', 50), binwidth=2)
        ax.set_title(params.get('title', chain_type))
        ax.set_xlabel('Length of fragments')
        ax.set_ylabel(None)
        return ax
