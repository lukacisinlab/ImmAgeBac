import scanpy as sc
import muon as mu
import time
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit, fsolve, least_squares
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors
import math
import statsmodels.api as sm
from statsmodels.formula.api import ols
import sys
from contextlib import contextmanager
from constantsForData import Constants

from cellCounts import CellCounts

class CellVisualization(CellCounts):
    def __init__(self):
        super().__init__()

    def get_counts(self, donor = None):
        counts_dict = {}
        data = self.adata
        if donor != None:
            data = self.adata[self.full_data.obs.donor == donor]
        markers = self.cell_into_markers()
        self.intersection = self.getIntersections()

        for cell, marker in markers.items():
            try:
                res = self.count_rows_with_multiple_conditions(data, marker, "CD27")
                counts_dict[cell] = res
            except:
                print("Failed at cell type:", cell, "\nWith markers:", marker)
            print("\n")

        return counts_dict


    def plot_donor_cells_histogram(self, n_rows, n_cols, save = False):
        fig = plt.figure(figsize=(12, 30))
        max_y = 0
        for i, donor in enumerate(self.donors):
            print("----", donor, "-----")
            ax = fig.add_subplot(n_rows, n_cols, i + 1)
            counts_dict = self.get_counts(donor)
            total_donor_count = self.adata[self.full_data.obs.donor == "P1"].n_obs

            ax.bar(counts_dict.keys(), [x / total_donor_count for x in counts_dict.values()], color='g')
            ax.set_title(donor)

            ax.set_xticks(range(len(counts_dict)))
            ax.set_xticklabels(counts_dict.keys(), rotation=90, ha='right')

            ax.set_ylim(0, 1)

        plt.tight_layout(pad=6.0)
        if save:
            plt.savefig("cell_bar_plots.png")
        else:
            plt.show()


    def swarm_plot(self, save = False):
        counts_per_donor = {donor: {} for donor in self.donors}

        for donor in self.donors:
            counts_dict = self.get_counts(donor)
            total_donor_count = self.adata[self.full_data.obs.donor == donor].n_obs

            for cell_type, count in counts_dict.items():
                counts_per_donor[donor][cell_type] = count / total_donor_count

        df = pd.DataFrame.from_dict(counts_per_donor, orient='index')
        df_melted = df.reset_index().melt(id_vars='index', var_name='Cell_Type', value_name='Count')
        df_melted.rename(columns={'index': 'Donor'}, inplace=True)

        plt.figure(figsize=(12, 8))
        sns.swarmplot(data=df_melted, x='Cell_Type', y='Count', hue='Donor', palette='Set2')
        plt.title('Swarm Plot of Immune Cell Counts per Donor')
        plt.ylabel('Count freqency')
        plt.legend(title='Donor')
        plt.xticks(rotation=45, ha='right')
        if save:
            plt.savefig("swaarm_plot.png")
        else:
            plt.show()
        return df_melted

#test = CellVisualization()
#test.plot_donor_cells_histogram(4,2,True)