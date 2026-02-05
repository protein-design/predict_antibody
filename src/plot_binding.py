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
        return ax

    def roc(self, ax, col):
        droc = self.data[col]
        fpr, tpr, thresholds = roc_curve(droc['y'], -droc[col])
        roc_auc = auc(fpr, tpr)
        ax.plot(fpr, tpr, color='grey')
        ax.set_xlabel('False positive rate')
        ax.set_ylabel('True positive rate')
        ax.set_title(f"AUC = {roc_auc:.2f}", fontsize=8)
        return ax