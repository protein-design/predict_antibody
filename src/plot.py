import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


class Plot:

    @sstaticmethod:
    def pie(ax, values, labels, params):
        explode = params['explode'] if 'explode' in params \
            else [0.01,] * len(labels)
        colors = params['colors'] if 'colors' in params \
            else ['lightgrey',] * params.get('topn', len(labels))
        wedges, texts, autotexts = ax.pie(
            values,
            labels=labels,
            autopct='%1.1f%%',
            pctdistance=params.get('pctdistance', .5),
            labeldistance=params.get('labeldistance', 1.0),
            explode=explode,
            colors=colors,
        )
        for text in texts:
            text.set_rotation(params.get('angle', 0))
            text.set_horizontalalignment(params.get('ha', 'center'))
            text.set_fontsize(params.get('fontsize', 0))
            text.set_fontstyle(params.get('fontstyle', 'normal'))
            text.set_fontweight(params.get('fontweight', 'bold'))
        return ax