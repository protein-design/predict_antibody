import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import py3Dmol


class PlotPdb:
    
    @staticmethod
    def bar_count_pdb(df, ymax=None):
        counts = df['structure_method'].value_counts()
        num_pdb = len(df['pdb_id'].unique())

        fig, ax = plt.subplots(1, figsize=(6,3), layout='tight')
        bars = ax.bar(counts.index, counts)
        ax.set_title(f'{num_pdb} PDB')
        ax.set_ylabel('Number of PDB')
        if ymax:
            ax.set_ylim(0, ymax)
        for j, bar in enumerate(bars):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2, height,
                str(counts.iloc[j]), ha='center', va='bottom')
        return fig, ax

    @staticmethod
    def violin_resolution(df):
        fig, ax = plt.subplots(1, figsize=(6,4), layout='tight')
    
        rdf = df[df['resolution'].notna()]
        sns.violinplot(rdf, x='structure_method', y='resolution', log_scale=True,
            color='grey', fill=False, inner='stick', ax=ax)
        ax.set_title('Distribution of Resolution')
        ax.set_xlim(-1, 3)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=20, ha='right')

        # add lines
        r = 5
        ax.axhline(np.log(r), color='green', linestyle='--')
        ax.text(3, np.log(r+1), str(r)+ r' $\AA$')
        r = 20
        ax.axhline(np.log(r), color='orange', linestyle='--')
        ax.text(3, np.log(r+2), str(r) + r' $\AA$')
        ax.set_xlabel(None)
        ax.set_ylabel(r'Resolution, log-$\AA$')
        return fig, ax

    @staticmethod
    def hist_resolution(df, height=200):
        '''
        select best resolution
        '''
        fig, ax = plt.subplots(1, figsize=(8,4), layout='tight')

        rdf = df[df['structure_method'].isin(['x-ray diffraction', 'electron microscopy'])]
        sns.histplot(rdf, x='resolution', hue='structure_method', alpha=.5, ax=ax)
        num_pdb = len(rdf['pdb_id'].unique())
        ax.set_title(f'Resolution Selection from {num_pdb} PDB data')
        ax.set_xlim(0, 10)

        # text
        resolution = rdf['resolution'][rdf['resolution'].notna()]
        q = .99
        qr = round(np.quantile(resolution, q), 1)
        ax.axvline(qr, linestyle='--')
        ax.text(qr+.1, height, f"quantile={q}\nresolution={qr}"+ ' $\AA$')
        q = .95
        qr = round(np.quantile(resolution, q), 1)
        ax.axvline(qr, linestyle='--')
        ax.text(qr+.1, height, f"quantile={q}\nresolution={qr}"+ ' $\AA$')
        ax.set_xlabel(r'Resolution, $\AA$')
        ax.set_ylabel('Number of PDB')
        return fig, ax

    @staticmethod
    def pie_complex(df):
        num_pdb = len(df['pdb_id'].unique())
        g = df.groupby('pdb_id').agg({'chain_no':'nunique'})
        g = g.value_counts()
        n = 6
        counts = g[:n]
        counts.index = ['monomer','dimer','trimer','tetramer', 'pentamer', 'heptamer']
        counts['other'] = sum(g[n:])
        counts = counts.reset_index()
        print(counts)

        fig, ax = plt.subplots(1, figsize=(6,4), layout='tight')
        # light color
        colors = plt.cm.Pastel1(np.arange(len(counts)))
        explode = [0, 0, 0, 0.1, 0.2, 0.3, 0.1]
        ax.pie(counts['count'], labels=counts['index'], autopct='%.1f%%', 
            colors=colors, explode=explode, startangle=90)
        ax.set_title(f'Protein Complex in {num_pdb} PDB data')
        return fig, ax

    @staticmethod
    def dot_bfactor(df, title_label:str=None, quantile=.9, xlim_max:int=None):
        fig, ax = plt.subplots(1, figsize=(8,4), layout='tight')

        sns.scatterplot(df, y='resolution', x='avg_bfactor', alpha=.1, s=10, ax=ax)
        ax.set_xlabel('Average B-factor')
        ax.set_ylabel(r'Resolution, $\AA$')
        title = f"PDB: {len(df)}"
        if title_label:
            title += title_label + f', quantile={quantile}'
        ax.set_title(title)
        if xlim_max:
            ax.set_xlim(0,xlim_max)
        # quantile
        bq = np.quantile(df['avg_bfactor'], quantile)
        print(f'Quantile of {quantile}, B-factor = {bq}')
        ax.axvline(bq, linestyle='--', color='lightgreen')
        rq = np.quantile(df['resolution'], quantile)
        print(f'Quantile of {quantile}, resolution = {rq}')
        ax.axhline(rq, linestyle='--', color='lightgreen')
        rect = patches.Rectangle((-5, -.1), bq+5, rq+.1, facecolor='lightgreen', alpha=.3)
        ax.add_patch(rect)

        # plt.yscale('log')
        # invert axis
        plt.gca().invert_xaxis()
        plt.gca().invert_yaxis()
        return fig, ax

'''
    @staticmethod
    def (df):
        fig, ax = plt.subplots(1, figsize=(10,4), layout='tight')
        return fig, ax
'''