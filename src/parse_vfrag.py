import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

class ParseVfrag:
    data_dir = '/home/yuan/output/pdb'

    def __init__(self, vregion, dist, verbose:bool=False):
        self.vregion = vregion
        self.dist = dist
        self.verbose = verbose
        self.data = None
    
    def process(self):
        vg = self.vregion.groupby(['pdb_id', 'chain_no', 'region_name'])
        vg = dict(tuple(vg))

