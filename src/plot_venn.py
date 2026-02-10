
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from venny4py.venny4py import venny4py
from matplotlib.patches import Circle

class PlotVenn:

    @staticmethod
    def venn_chain(sets, ax, labels_xy):
        colors = ['grey',] * len(sets)
        venny4py(sets=sets, size=3, colors=colors, edge_color='black',
            column_spacing=2, line_width=.5, dpi=300, asax=ax)
        # set labels
        legend = ax.legend()
        legend.set_visible(False)
        for i in range(len(sets)):
            x, y = labels_xy[i]
            label = list(sets)[i]
            ax.text(x, y, label, fontsize=10)
        return ax