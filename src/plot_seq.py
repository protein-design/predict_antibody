import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

class PlotSeq:

    @staticmethod
    def pie_species(df:pd.DataFrame, n:int, explode:list, args):
        # counts
        cdf = df[df['specie'].notna()]
        values = cdf['specie'].value_counts()
        values = values.reset_index()
        counts = pd.DataFrame(values.iloc[:n,:])
        other = sum(values['count'][n:])
        counts.loc[n+1] = ['Others', other]

        fig, ax = plt.subplots(1,2, figsize=(9,6), layout='tight')
        fig.suptitle(f'Aligned IMGT Genes')

        # draw pie
        wedges, texts, autotexts = ax[0].pie(
            counts['count'],
            labels=counts['specie'],
            autopct='%1.1f%%',
            pctdistance=args.get('pctdistance', .5),    # Move percentages closer to the center
            explode=explode,
            labeldistance=args.get('labeldistance', 1.0)
        )
        for text in texts[:-1]:
            text.set_rotation(args.get('angle', 10))
            text.set_horizontalalignment(args.get('ha', 'center'))
            text.set_fontsize(8)
            text.set_fontstyle('italic')
            text.set_fontweight('bold')

        # add table
        ax[1].axis('off')
        counts.loc[n+2] = ['Total', sum(values['count'])]
        table = pd.plotting.table(ax[1], counts, loc="center", cellLoc="left", \
                    colWidths=[0.5,]*counts.shape[1])
        plt.show()

    @staticmethod
    def plot_summary(df:pd.DataFrame):
        # combine data
        # data = self.df_combined()

        # count numbers
        counts = df.groupby(['specie', 'chain_type']).agg({'chain_seq': len})
        counts = counts.reset_index()
        counts = counts.pivot(index='chain_type', columns='specie', values='chain_seq')
        
        # plot
        fig, ax = plt.subplots(1, figsize=(8,5), layout='tight')
        fig.suptitle(f"V-Region Sequences")
        sns.boxplot(df, x='specie', y='pro_len', hue='chain_type', notch=True, ax=ax)
        ax.set_ylabel('Length of AA')
        ax.set_ylim(0, 600)
        pd.plotting.table(ax, counts, loc="top", cellLoc="center", colWidths=[0.2,]*counts.shape[1])
        plt.legend(loc='lower right')
        plt.show()
        return df, counts
    
    @staticmethod    
    def plot_summary_specie(df:pd.DataFrame, specie):
        sdf = df[df['specie']==specie]
        print(sdf.shape)
        g = sdf.groupby('chain_type').agg({'gene_name':'nunique', 'allele_name':'nunique', 'chain_seq': 'nunique'})
        g = g.reset_index()

        fig, ax = plt.subplots(2,2, figsize=(8,5), layout='tight')
        fig.suptitle(f"Chain Sequences of Antibodies, {specie}")

        sns.barplot(g, x='chain_type', y ='gene_name', ax=ax[0][0])
        ax[0][0].set_title('Number of IMGT Gene Family')
        ax[0][0].set_ylabel('Unique Number')
        sns.barplot(g, x='chain_type', y ='allele_name', ax=ax[0][1])
        ax[0][1].set_title('Number of IMGT Allele Name')
        ax[0][1].set_ylabel('Unique Number')
        sns.barplot(g, x='chain_type', y ='chain_seq', ax=ax[1][0])
        ax[1][0].set_title('Number of Chains')
        ax[1][0].set_ylabel('Number')

        sns.boxplot(sdf, x='chain_type', y ='pro_len', notch=True, ax=ax[1][1])
        ax[1][1].set_title('Length of Sequences')
        ax[1][1].set_ylabel('Length of AA')

   
    @staticmethod
    def pie_gene(df:pd.DataFrame, specie, chain_type:str, key:str, n:int, args):
        cdf = df[(df['specie']==specie)&(df['chain_type']==chain_type)]
        values = cdf[key].value_counts()
        values = values.reset_index()
        counts = pd.DataFrame(values.iloc[:n,:])
        other = sum(values['count'][n:])
        explode = [args.get('explode_unit', .1),] * n
        if other > 0:
            counts.loc[n+1] = ['Others', other]
            explode += [0,]

        fig, ax = plt.subplots(1,2, figsize=(9,6), layout='tight')
        fig.suptitle(f'Aligned IMGT Genes, {key}, {specie}, chain: {chain_type}')
        wedges, texts, autotexts = ax[0].pie(
            counts['count'],
            labels=counts[key],
            autopct='%1.1f%%',
            pctdistance=args.get('pctdistance', .5),    # Move percentages closer to the center
            explode=explode,
            labeldistance=args.get('labeldistance', 1.2)
        )
        for text in texts[:-1]:
            text.set_rotation(args.get('angle', 20))
            text.set_horizontalalignment(args.get('ha', 'center'))
        # add table
        ax[1].axis('off')
        counts.loc[n+2] = ['Total', sum(values['count'])]
        table = pd.plotting.table(ax[1], counts, loc="center", cellLoc="left", \
                    colWidths=[0.5,]*counts.shape[1])
        plt.show()
