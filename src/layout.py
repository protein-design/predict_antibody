'''
figure 6.5x8 inch (16.5x20 cm)
'''
import numpy as np
import pandas as pd
import seaborn as sns
import string

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'sans-serif']

class Layout:
    # unit: centimeters
    # one column, 1.5 column, full width
    col_width = (6, 8.5, 11.4, 17.4)
    panel_label = string.ascii_uppercase

    def __init__(self, args:dict={}):
        self.args = args
        width_level = self.args.get('width_level', 1)
        self.width = self.col_width[width_level]
        self.height = self.args.get('height', self.width * .9)
        # unit: inch
        self.figsize = (self.width / 2.54, self.height / 2.54)
        # print(self.figsize)
        print(f"figure size: {round(self.width, 1)} x {round(self.height, 1)} cm")

    def one(self):
        fig, ax = plt.subplots(1, figsize=self.figsize, layout='tight')
        return fig, ax
    
    def row(self, label_x:tuple, label_y:float):
        '''
        example label_x = (-5,-5), label_y = 10
        '''
        fig = plt.figure(figsize=self.figsize, constrained_layout=True)
        n = len(label_x)
        wrt = self.args.get('width_ratios', None)
        gs = gridspec.GridSpec(1, n, figure=fig, width_ratios=wrt)
        gs.update(wspace=self.args.get('space', .1))
        axes = []
        for i in range(n):
            ax = fig.add_subplot(gs[0, i])
            xytext = (label_x[i], label_y)
            ax.annotate(self.panel_label[i], (0,1), xytext=xytext,
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='top', fontsize=10, fontweight='bold')
            axes.append(ax)
        return fig, tuple(axes)

    def column(self, label_x:float, label_y:tuple):
        '''
        example label_x = (-5,-5), label_y = 10
        '''
        fig = plt.figure(figsize=self.figsize, constrained_layout=True)
        n = len(label_y)
        gs = gridspec.GridSpec(n, 1, figure=fig)
        axes = []
        for i in range(n):
            ax = fig.add_subplot(gs[i, 0])
            xytext = (label_x, label_y[i])
            ax.annotate(self.panel_label[i], (0,1), xytext=xytext,
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='top', fontsize=10, fontweight='bold')
            axes.append(ax)
        return fig, tuple(axes)

    def grid(self, nrow:int, ncol:int, panel_xytext:list):
        fig = plt.figure(figsize=self.figsize, constrained_layout=True)
        gs = gridspec.GridSpec(ncols=ncol, nrows=nrow, figure=fig)
        # spaces between panels
        gs.update(
            wspace=self.args.get('wspace', .1),
            hspace=self.args.get('hspace', .1)
        )
        # add panels
        axes = []
        for i in range(nrow):
            for j in range(ncol):
                ax = fig.add_subplot(gs[i, j])
                axes.append(ax)

        # add panel labels
        n = nrow * ncol
        for i in range(n):
            ax = axes[i]
            ax.annotate(
                self.panel_label[i],
                (0,1), 
                xytext=panel_xytext[i],
                xycoords='axes fraction',
                textcoords='offset points',
                ha='left',
                va='top',
                fontsize=10,
                fontweight='bold'
            )
        return fig, tuple(axes)
        
    def row1_row2(self, panel_xytext:list):
        '''
        1st row: one figure, 2nd row: two figures
        '''
        fig = plt.figure(figsize=self.figsize)
        hrt=self.args.get('height_ratios', [1, 1])
        wrt=self.args.get('width_ratios', [1, 1])
        gs = gridspec.GridSpec(nrows=2, ncols=2, height_ratios=hrt, width_ratios=wrt)
        # spaces between panels
        gs.update(
            wspace=self.args.get('wspace', 0.4),
            hspace=self.args.get('hspace', 0.4)
        )

        # Create axes objects that span specific grid locations
        ax1 = fig.add_subplot(gs[0, :])
        ax2 = fig.add_subplot(gs[1, 0])
        ax3 = fig.add_subplot(gs[1, 1])
        axes = [ax1, ax2, ax3]

        # add panel labels
        for i in range(3):
            ax = axes[i]
            ax.annotate(self.panel_label[i], (0,1), xytext=panel_xytext[i],
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='top', fontsize=10, fontweight='bold')
        return fig, tuple(axes)

    def row2_row3(self, panel_xytext:list):
        '''
        1st row: 2 figures, 2nd row: 3 figures
        '''
        fig = plt.figure(figsize=self.figsize)
        gs = gridspec.GridSpec(nrows=2, ncols=6, 
            height_ratios=self.args.get('height_ratios', [1,1,]),
            width_ratios=self.args.get('width_ratios', [1,]*6),
            # spaces between panels
            wspace=self.args.get('wspace', .2),
            hspace=self.args.get('hspace', .2),
        )
        # Create axes objects that span specific grid locations
        ax1 = fig.add_subplot(gs[0, :3])
        ax2 = fig.add_subplot(gs[0, 3:])
        ax3 = fig.add_subplot(gs[1, 0:2])
        ax4 = fig.add_subplot(gs[1, 2:4])
        ax5 = fig.add_subplot(gs[1, 4:6])
        axes = [ax1, ax2, ax3, ax4, ax5]

        # add panel labels
        for i in range(5):
            ax = axes[i]
            ax.annotate(self.panel_label[i], (0,1), xytext=panel_xytext[i],
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='top', fontsize=10, fontweight='bold')
        return fig, tuple(axes)
    
    def col1_col2(self, panel_xytext:list):
        '''
        1st col: one figure, 2nd col: two figures in two rows
        '''
        fig = plt.figure(figsize=self.figsize)
        gs = gridspec.GridSpec(
            nrows=2,
            ncols=2,
            height_ratios=self.args.get('height_ratios', [1, 1]),
            width_ratios=self.args.get('width_ratios', [1, 1])
        )
        # spaces between panels
        gs.update(
            wspace=self.args.get('wspace', 0.4),
            hspace=self.args.get('hspace', 0.4)
        )

        # Create axes objects that span specific grid locations
        ax1 = fig.add_subplot(gs[:, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[1, 1])
        axes = [ax1, ax2, ax3]

        # add panel labels
        for i in range(3):
            ax = axes[i]
            ax.annotate(self.panel_label[i], (0,1), xytext=panel_xytext[i],
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='top', fontsize=10, fontweight='bold')
        return fig, tuple(axes)