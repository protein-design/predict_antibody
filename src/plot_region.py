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
        frag = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3',]
        sdf = self.df[(self.df['specie']==specie) & (self.df['region_name'].isin(frag))]
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

    def dssp_percent(self, specie:str):
        sgdf = self.df[self.df['specie']==specie]
        if self.verbose:
            print(specie, sgdf.shape)

        g = sns.FacetGrid(sgdf, col='chain_type', row='structure')
        g.map(sns.barplot, 'region_name', 'percent',
            order=['FR1', 'FR2', 'FR3', 'CDR1','CDR2','CDR3'],
            color='black', width=.6
        )
        g.set_titles(col_template="{col_name}", row_template="{row_name}",
            size=8, fontweight='bold')
        g.set(xlabel=None, ylabel=None)
        g.fig.text(0.5, -.03, 'Fragments of V-region',
            ha='center', rotation='horizontal',
            fontsize=10, fontweight='bold'
        )
        g.fig.text(0, 0.5, 'Percentage of amino acids, %',
            va='center', rotation='vertical',
            fontsize=10, fontweight='bold'
        )
        for ax in g.axes.flat:
            ax.tick_params(axis='both', labelsize=8)
            ax.axhline(50, linestyle='--', color='lightgrey')
    
        g.fig.subplots_adjust(wspace=.1, hspace=.6)
        g.fig.set_size_inches(17/2.54, 20/2.54)
        return g