import numpy as np
from numpy.polynomial.polynomial import polyfit
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from sklearn.metrics import precision_recall_curve, roc_curve, auc

class PlotBinding:

    def __init__(self, data, verbose:bool=True):
        self.data = data
        self.verbose = verbose

    def hist_contacts(self, ax, cutoff):
        df = self.data[cutoff]
        sns.histplot(df, x='num_ca', stat='percent', ax=ax, bins=100, color='black')
        ax.set_title(f'threshold = {cutoff}')
        ax.set_xlabel('Number of Ca contacts')
        ax.set_ylabel('Percentage, %')
    
    def hist_binding_affinity(self, ax):
        df = self.data
        sns.histplot(df, x='binding_affinity', stat='percent', ax=ax, bins=100, color='grey')
        ax.set_ylim(0,12)
        ax.set_xlabel('Binding affinity, kcal/mol')
        ax.set_ylabel('Percentage, %')
        q = np.quantile(df['binding_affinity'], .95)
        ax.axvline(q, linestyle='--', color='black')
        ax.text(-18, 11, f"95% = {q:.1f}", fontsize=8)

    def hist_kd(self, ax):
        df = self.data
        sns.histplot(df, x='log-kd', stat='percent', ax=ax, bins=100, color='grey')
        ax.set_ylim(0, 12)
        ax.set_xlabel('Dissociation constant, -logM')
        ax.set_ylabel('Percentage, %')
        q = np.quantile(df['dissociation_constant'], .95)
        ax.axvline(np.log(q), linestyle='--', color='black')
        ax.text(-10, 11, f"95% = {q}", fontsize=8, ha='right')

    def precision_recall(self, ax, col):
        droc = self.data[col]
        precisions, recalls, thresholds = precision_recall_curve(
            droc['y'], -droc[col]
        )
        ax.plot(thresholds, recalls[:-1], label='Recall',
            color='black', linestyle='-')
        ax.plot(thresholds, precisions[:-1], label='Precision',
            color='black', linestyle='--')
        ax.axvline(-10, linestyle='--', color='grey')
        ax.axvline(-20, linestyle='--', color='grey')
        ax.set_xlim(-100, 0)
        ax.set_xlabel('Minus minimum distance of Ca, -$\AA$')
        ax.legend(loc='center left', fontsize=8)

    def roc(self, ax, col):
        droc = self.data[col]
        fpr, tpr, thresholds = roc_curve(droc['y'], -droc[col])
        roc_auc = auc(fpr, tpr)
        ax.plot(fpr, tpr, color='black')
        ax.set_xlabel('False positive rate')
        ax.set_ylabel('True positive rate')
        ax.set_title(f"AUC = {roc_auc:.2f}", fontsize=8)

    def box_contacts_distance(self, ax, col):
        df = self.data[self.data[col]<=30].reset_index()
        df[col] = df[col].round(0).astype(int)
        sns.boxplot(df, x=col, y='binding_affinity', ax=ax, 
            notch=True, fill=False, color='black')
        ax.set_xlabel(f"Minimum distance of the top {col} ranked Ca, $\AA$")
        ax.set_ylabel('Binding affinity, kcal/mol')
        ax.axhline(-5.5, color='grey', linestyle='--')