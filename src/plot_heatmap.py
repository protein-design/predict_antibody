import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

class PlotHeatmap:
    def __init__(self, df:pd.DataFrame, verbose:bool=True):
        self.df = df
        self.verbose = verbose
    
    def region_seq_abundance(self, ax, params):
        region_name = params['region_name']
        col_name = params['col']
        row_name = params['row']
        cdf = self.df[self.df['region_name']==region_name]
        cdf = cdf[[col_name, row_name]].groupby([col_name, row_name]).agg(len)
        cdf = cdf.reset_index().rename(columns={0:'count'})
        cdf = cdf.pivot(index=row_name, columns=col_name, values='count')
        cdf = cdf.fillna(0)

        # filter tops in rows: default is seq
        num_row = params.get('num_row', 10)
        row_counts = cdf.apply(sum, axis=1)
        row_counts = row_counts.sort_values(ascending=False)
        top_row_counts = row_counts[:num_row]
        if self.verbose:
            print(top_row_counts.index)
        # filter tops in column
        num_col = params.get('top_col', 10)
        col_counts = cdf.apply(sum, axis=0)
        col_counts = col_counts.sort_values(ascending=False)
        top_col_counts = col_counts[:num_col]
        if self.verbose:
            print(top_col_counts.index)
        
        # build metrics
        fcdf = cdf[top_col_counts.index]
        fcdf = fcdf.loc[top_row_counts.index]
        if self.verbose:
            print(fcdf.shape)
            print(fcdf.head())
            
        # plot
        colors_list = ['gainsboro', 'blue', 'green', 'red']
        cmap = mpl.colors.LinearSegmentedColormap.from_list("custom_red_green_white", colors_list)
        cbar = {'orientation': 'horizontal','location':'top','pad':.04}
        sns.heatmap(fcdf, cmap=cmap, ax=ax, cbar_kws=cbar, linewidth=.5,
            yticklabels=True, xticklabels=True)
        # ax.tick_params(labelsize=8, labelfontfamily='Serif')
        ax.set_xticklabels(ax.get_xticklabels(), size=8, weight='normal')
        ax.set_yticklabels(ax.get_yticklabels(), size=8, weight='normal')
        xlabel = params['xlabel'] if 'xlabel' in params else params['col']
        ax.set_xlabel(xlabel, size=10, fontweight='bold')
        ylabel = params['ylabel'] if 'ylabel' in params else params['row']
        ax.set_ylabel(ylabel, size=10, fontweight='bold')
        ax.text(0, -.1, 1, size=8)
        ax.text(num_col-1, -.1, num_col, size=8)
        ax.text(num_col, 1, 1, size=8)
        ax.text(num_col, num_row, num_row, size=8)
        return ax