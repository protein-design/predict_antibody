import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

class PlotInterface:

    def __init__(self, df, verbose:bool=True):
        self.df = df
        self.verbose = verbose


    def hist_overlap(self, ax, region_name):
        overlap_percent = self.df['overlap']*100
        sns.histplot(overlap_percent, stat='percent',
            bins=50, color='black', ax=ax)
        ax.set_xlabel(f'{region_name}: Overlap rate, %')
        ax.set_ylabel('Percentage, %')
        return ax

    def hist_len(self, ax, region_name, xlim):
        sns.histplot(self.df, x='seq_len', stat='percent',
            bins=50, color='grey', ax=ax)
        ax.set_xlabel(f'Length of {region_name}')
        ax.set_xlim(xlim)
        ax.set_ylabel('Percentage, %')
        return ax