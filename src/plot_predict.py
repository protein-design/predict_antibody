import numpy as np
from numpy.polynomial.polynomial import polyfit
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

class PlotPredict:

    def __init__(self, df, verbose:bool=True):
        self.df = df
        self.verbose = verbose

    def violin_plddt(self, ax):
        pdf = self.df[(self.df['predictor']=='alphafold2') & (self.df['ranking']==1)]

        sns.violinplot(pdf, y='chain_type', x='avg_plddt', ax=ax,
            fill=False, inner='quart', split=True)
        ax.set_xlabel('Average pLDDT')
        ax.set_ylabel('Chain ')
        return ax
    
    def violin_pdb_rmsd(self, ax):
        self.df['rmsd_adj'] = self.df['rmsd'] + 1e-5

        sns.violinplot(self.df, x='chain_status', y='rmsd_adj', ax=ax,
            log_scale=True, color='black', fill=False)
        ax.set_xlim(-.6, 2.5)
        ax.set_xlabel('PDB Post-processing')
        ax.set_ylabel('Log RMSD')
        x_offset= -.55
        n=5
        ax.axhline(np.log(n), linestyle='--', color='grey')
        ax.text(x_offset, np.log(n)-1, n)
        n=50
        ax.axhline(np.log(n), linestyle='--', color='grey')
        ax.text(x_offset, np.log(n)+1.2, 50)
        return ax
    

    def violin_pdb_tm(self, ax):
        sns.violinplot(self.df, x='chain_status', y='tm1', ax=ax,
            log_scale=True, color='black', fill=False)
        ax.set_xlim(-.6, 2.5)
        ax.set_ylim(-.1, 1.4)
        ax.set_xlabel('PDB Post-processing')
        ax.set_ylabel('TM-align')
        ax.axhline(1, linestyle='--', color='grey')
        return ax
    


    def dot_plddt_rmsd(self, ax, chain_type=None, q=.95):
        sdf = self.df[self.df['chain_type']==chain_type] \
            if chain_type else self.df
        
        col1 = 'avg_plddt'
        col2 = 'rmsd'
        sns.scatterplot(sdf, x=col1, y=col2, ax=ax,
            alpha=.1, color='black', s=10)
        ax.set_xlabel('Average pLDDT')
        ax.set_ylabel('RMSD')
        ax.set_xlim(60, 100)
        ax.set_ylim(-1, 30)
        ax.invert_yaxis()

        # fit
        b, m = polyfit(sdf[col1], sdf[col2], 1)
        print(f"intercept={b}, coef={m}")
        ax.plot(sdf[col1], b + m * sdf[col1],
            '-', color='blue')
        coef, p_value = pearsonr(sdf[col1], -sdf[col2])
        ax.set_title(f"$R^2$ = {coef:.2f}")

        # quantile
        q1 = np.quantile(sdf[col1], 1-q)
        print(q1)
        ax.axvline(q1, linestyle='--', color='grey')
        q2 = np.quantile(sdf[col2], q)
        print(q2)
        ax.axhline(q2, linestyle='--', color='grey')
        return ax

    def dot_plddt_tm(self, ax, chain_type=None, q=.95):
        sdf = self.df[self.df['chain_type']==chain_type] \
            if chain_type else self.df
        
        sns.scatterplot(sdf, x='avg_plddt', y='tm1', ax=ax,
            alpha=.1, color='black', s=10)
        ax.set_xlabel('Average pLDDT')
        ax.set_ylabel('TM-align')
        ax.set_xlim(60, 100)

        # fit
        b, m = polyfit(sdf['avg_plddt'], sdf['tm1'], 1)
        print(f"intercept={b}, coef={m}")
        ax.plot(sdf['avg_plddt'], b + m * sdf['avg_plddt'],
            '-', color='blue')
        coef, p_value = pearsonr(sdf['avg_plddt'], sdf['tm1'])
        ax.set_title(f"$R^2$ = {coef:.2f}")

        # quantile
        q_plddt = np.quantile(sdf['avg_plddt'], 1-q)
        print(q_plddt)
        ax.axvline(q_plddt, linestyle='--', color='grey')
        q_tm = np.quantile(sdf['tm1'], q)
        print(q_tm)
        ax.axhline(q_tm, linestyle='--', color='grey')
        return ax
    
    def dot_plddt_rmsd2(self, ax, chain_type=None, q=.95):
        sdf = self.df[self.df['chain_type']==chain_type] \
            if chain_type else self.df
        
        col1 = 'avg_plddt'
        col2 = 'rmsd'
        sns.scatterplot(sdf, x=col1, y=col2, ax=ax,
            alpha=.1, color='black', s=10)
        ax.set_xlim(60, 100)
        ax.set_xlabel('average pLDDT')
        ax.set_ylabel('RMSD')

        # fit
        b, m = polyfit(sdf[col1], sdf[col2], 1)
        print(f"intercept={b}, coef={m}")
        ax.plot(sdf[col1], b + m * sdf[col1],
            '-', color='blue')
        coef, p_value = pearsonr(sdf[col1], sdf[col2])
        ax.text(60, b+m*60, round(coef,2))
        ax.invert_yaxis()

        # quantile
        q_plddt = np.quantile(sdf[col1], 1-q)
        print(q_plddt)
        ax.axvline(q_plddt, linestyle='--', color='grey')
        q2 = np.quantile(sdf[col2], q)
        print(q2)
        ax.axhline(q2, linestyle='--', color='grey')
        return ax
