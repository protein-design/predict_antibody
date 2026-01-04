import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

class PlotSeq:

    def __init__(self, df:pd.DataFrame, verbose:bool=True):
        self.df = df
        self.verbose = verbose

    def pie_specie_counts(self, ax, params):
        # counts
        cdf = self.df[self.df['specie'].notna()]
        values = cdf['specie'].value_counts()
        values = values.reset_index()
        topn = params['n']
        counts = pd.DataFrame(values.iloc[:topn,:])
        other = sum(values['count'][topn:])
        if other > 0:
            topn += 1
            counts.loc[topn] = ['Others', other]

        # draw pie
        wedges, texts, autotexts = ax.pie(
            counts['count'],
            labels=counts['specie'],
            autopct='%1.1f%%',
            pctdistance=params.get('pctdistance', .5),
            explode=params['explode'],
            labeldistance=params.get('labeldistance', 1.0),
            colors=['lightgrey',] * topn,
        )
        # set labels
        for text in texts:
            text.set_rotation(0)
            text.set_horizontalalignment('left')
            text.set_fontsize(8)
            text.set_fontweight('bold')
        # fontstyle except 'others'
        for text in texts[:-1]:
            text.set_fontstyle('italic')
        # position label of 'human' and mouse
        texts[0].set_position(((-.3, .6)))
        texts[1].set_position(((-.9, -.2)))
        return ax

    def table_species(self, ax, params):
        # counts
        cdf = self.df[self.df['specie'].notna()]
        values = cdf['specie'].value_counts()
        values = values.reset_index()
        topn = params.get('n', 0)

        # other counts
        counts = pd.DataFrame()
        if n:
            counts = pd.DataFrame(values.iloc[:topn,:])
            # other counts
            other = sum(values['count'][topn:])
            if other > 0:
                topn += 1
                counts.loc[topn] = ['Others', other]
        else:
            counts = values.copy()

        # total counts
        topn += 1
        counts.loc[topn] = ['Total', sum(values['count'])]
        
        # draw table
        ax.axis('off')
        table = pd.plotting.table(
            ax,
            counts,
            loc="center",
            cellLoc="left",
            colWidths=[0.5,] * counts.shape[1]
        )
        return ax

    def plot_summary(self, n):
        top_names = list(self.df['specie'].value_counts()[:n].index)
        if self.verbose:
            print(top_names)
        sdf = self.df[self.df['specie'].isin(top_names)]

        # count numbers
        counts = sdf.groupby(['specie', 'chain_type']).agg({'chain_seq': len})
        counts = counts.reset_index()
        counts = counts.pivot(index='chain_type', columns='specie', values='chain_seq')
        
        # plot
        fig, ax = plt.subplots(1, figsize=(8,5), layout='tight')
        fig.suptitle(f"V-Region Sequences")
        sns.boxplot(sdf, x='specie', y='pro_len', hue='chain_type', notch=True, ax=ax)
        ax.set_ylabel('Length of AA')
        ax.set_ylim(0, 600)
        pd.plotting.table(ax, counts, loc="top", cellLoc="center", colWidths=[0.2,]*counts.shape[1])
        plt.legend(loc='lower right')
        plt.show()
        return counts
   
    def bar_gene_family(self, ax, params):
        specie = params['specie']
        sdf = self.df[self.df['specie']==specie]
        if self.verbose:
            print(sdf.shape)
        g = sdf.groupby('chain_type').agg({'gene_name':'nunique'})
        g = g.reset_index()

        # barplot
        sns.barplot(g, x='chain_type', y ='gene_name', ax=ax, color='grey')
        ax.set_xlabel(None)
        xtick_labels = ax.get_xticklabels()
        ax.set_xticks(range(len(xtick_labels)))
        ax.set_xticklabels(xtick_labels, rotation=30, ha='right')
        ax.tick_params(axis='x', pad=0)
        ax.set_ylabel('Number of\ngene families')
        return ax

    def bar_allele_name(self, ax, params):
        specie = params['specie']
        sdf = self.df[self.df['specie']==specie]
        if self.verbose:
            print(sdf.shape)
        g = sdf.groupby('chain_type').agg({'allele_name':'nunique'})
        g = g.reset_index()

        # barplot
        sns.barplot(g, x='chain_type', y ='allele_name', ax=ax, color='grey')
        ax.set_xlabel(None)
        xtick_labels = ax.get_xticklabels()
        ax.set_xticks(range(len(xtick_labels)))
        ax.set_xticklabels(xtick_labels, rotation=30, ha='right')
        ax.tick_params(axis='x', pad=0)
        ax.set_ylabel('Number of\nallele names')
        return ax
    
    def bar_chain_seq(self, ax, params):
        specie = params['specie']
        sdf = self.df[self.df['specie']==specie]
        if self.verbose:
            print(sdf.shape)
        g = sdf.groupby('chain_type').agg({'chain_seq':'nunique'})
        g = g.reset_index()

        # barplot
        sns.barplot(g, x='chain_type', y ='chain_seq', ax=ax, color='grey')
        ax.set_xlabel(None)
        xtick_labels = ax.get_xticklabels()
        ax.set_xticks(range(len(xtick_labels)))
        ax.set_xticklabels(xtick_labels, rotation=30, ha='right')
        ax.tick_params(axis='x', pad=0)
        ax.set_ylabel('Number of\nunique sequences')
        return ax

    def box_chain_len(self, ax, params):
        specie = params['specie']
        sdf = self.df[self.df['specie']==specie]
        if self.verbose:
            print(sdf.shape)

        sns.violinplot(sdf, x='chain_type', y ='pro_len', ax=ax, color='lightgrey')
        ax.set_xlabel('Chain of antibody')
        ax.set_ylabel('Length of AA')
        return ax
    
    def hist_chain_len(self, ax, params):
        specie = params['specie']
        sdf = self.df[self.df['specie']==specie]
        if 'min_len' in params:
            sdf = sdf[sdf['pro_len']>=params['min_len']]
        if 'max_len' in params:
            sdf = sdf[sdf['pro_len']<params['max_len']]
        if self.verbose:
            print(sdf.shape)

        bins=params.get('bins', 50)
        sns.histplot(sdf, x ='pro_len', hue='chain_type', ax=ax, bins=bins, alpha=.5)
        ax.set_xlabel('Length of AA')
        return ax
       
    def pie_specie_chain(self, ax, specie, chain_type:str, key:str, topn:int, params):
        cdf = self.df[(self.df['specie']==specie)&(self.df['chain_type']==chain_type)]
        values = cdf[key].value_counts()
        values = values.reset_index()
        counts = pd.DataFrame(values.iloc[:topn,:])
        other = sum(values['count'][topn:])
        if other > 0:
            topn += 1
            counts.loc[topn] = ['Others', other]
        if self.verbose:
            print(f"{specie}, number of pies: {len(counts)}")

        # plot
        explode = params.get('explode', [.1,] * topn)
        wedges, texts, autotexts = ax.pie(
            counts['count'],
            labels=counts[key],
            autopct='%1.1f%%',
            textprops={'fontsize': 8},
            pctdistance=params.get('pctdistance', .7),
            explode=explode,
            labeldistance=params.get('labeldistance', 1.1),
            colors = ['lightgrey'] * topn,
        )
        ax.set_xlabel(chain_type + ' chain', labelpad=0)
        # set labels
        for text in texts:
            text.set_rotation(0)
            text.set_fontsize(8)
            text.set_horizontalalignment(params.get('ha', 'left'))
        # adjust label positions
        for i, xy in params.get('label_pos', {}).items():
            texts[i].set_position(xy)
        return ax

    def table_specie_chain(self, ax, specie, chain_type:str, key:str, topn:int):
        cdf = self.df[(self.df['specie']==specie)&(self.df['chain_type']==chain_type)]
        values = cdf[key].value_counts()
        values = values.reset_index()
        counts = pd.DataFrame(values.iloc[:topn,:])
        other = sum(values['count'][topn:])
        if other > 0:
            topn += 1
            counts.loc[topn] = ['Others', other]

        # add table
        topn += 1
        counts.loc[topn] = ['Total', sum(values['count'])]
        ax.axis('off')
        table = pd.plotting.table(ax, counts, loc="center", cellLoc="left", \
            colWidths=[0.5,]*counts.shape[1])
        return ax

    def bar_specie_chain(self, ax, specie, chain_type:str, key:str, topn:int, params):
        cdf = self.df[(self.df['specie']==specie)&(self.df['chain_type']==chain_type)]
        values = cdf[key].value_counts()
        values = values.reset_index()
        counts = pd.DataFrame(values.iloc[:topn,:])
        counts['percent'] = counts['count'] * 100 /sum(counts['count'])

        width = params.get('bar_width', None)
        sns.barplot(counts, x='percent', y=key, color='grey', width=width, ax=ax)
        ax.set_title(chain_type + ' chain')
        ax.set_xlabel('Percentage, %', fontsize=8)
        ax.set_ylabel(params.get('ylabel', key), fontsize=10)
        ax.tick_params(axis='x', labelsize=8)
        ax.tick_params(axis='y', labelsize=8)
        return ax
