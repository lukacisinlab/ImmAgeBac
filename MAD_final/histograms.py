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

class Histograms:
    def __init__(self):
        self.full_data = mu.read("Satija.h5mu")
        self.adata = self.full_data['ADT']

        self.expected_dictionary = Constants.EXPECTED_DICTIONARY
        self.manual_histograms = Constants.MANUAL_HISTOGRAMS
        self.missing_data = Constants.MISSING_DATA
        self.manual_intersection = Constants.MANUAL_INTERSECTION
        self.cell_type_dict = Constants.CELL_TYPE_DICT
        self.CELLS = Constants.CELLS

        self.metadata_df = pd.read_excel("mmc3.xlsx")
        self.metadata_df.set_index("Donor", inplace=True)
        self.donors = sorted(list({x for x in self.full_data.obs.donor}))


    def draw_histograms_manual(self):
        variables = self.manual_histograms
        n_rows = math.ceil(len(variables) / 2)
        n_cols = 2
        fig = plt.figure(figsize = (12,12))
        for i, var_name in enumerate(variables):
            ax = fig.add_subplot(n_rows, n_cols, i + 1)
            data = self.adata.X[: ,self.adata.var.index == var_name]
            y,x,_ = ax.hist(data, bins = 100, alpha = .3)
            ax.set_title(var_name)

        fig.tight_layout()
        plt.show()


    def gauss(self,x, mu, sigma, A):
        return A * np.exp(-(x - mu)**2 / 2 / sigma**2)

    def bimodal(self, x, mu1, sigma1, A1, mu2, sigma2, A2):
        return self.gauss(x, mu1, sigma1, A1) + self.gauss(x, mu2, sigma2, A2)

    def find_intersections(self, gauss1_params, gauss2_params, x_range):
        mu1, sigma1, A1 = gauss1_params
        mu2, sigma2, A2 = gauss2_params

        def difference(x):
            return self.gauss(x, mu1, sigma1, A1) - self.gauss(x, mu2, sigma2, A2)

        x_initial_guesses = np.linspace(x_range[0], x_range[1], 500)
        intersections = []

        for x0 in x_initial_guesses:
            root, info, ier, msg = fsolve(difference, x0, full_output=True)
            if ier == 1 and x_range[0] <= root <= x_range[1]:
                intersections.append(root[0])

        intersections = np.unique(intersections)
        return intersections


    def draw_histograms(self, save = False):
        variables = self.expected_dictionary.keys()
        n_rows = math.ceil(len(variables) / 2)
        n_cols = 2

        intersection_dict = {}
        fig = plt.figure(figsize=(12,12))
        for i, var_name in enumerate(variables):
            ax = fig.add_subplot(n_rows, n_cols, i + 1)
            data = self.adata.X[:, self.adata.var.index == var_name]
            if data.size == 0:
                print(var_name + " is empty\n")
                continue
            y, x, _ = ax.hist(data, bins=100, alpha=.3)
            x = (x[1:] + x[:-1]) / 2

            expected = self.expected_dictionary[var_name]
            try:
                params, cov = curve_fit(self.bimodal, x, y, expected)
            except RuntimeError:
                print(var_name)
                continue
            sigma = np.sqrt(np.diag(cov))
            x_fit = np.linspace(x.min(), x.max(), 500)

            ax.plot(x_fit, self.bimodal(x_fit, *params), color='red', lw=0)

            ax.plot(x_fit, self.gauss(x_fit, *params[:3]), color='red', lw=1, ls="--")
            ax.plot(x_fit, self.gauss(x_fit, *params[3:]), color='green', lw=1, ls=":")

            intersections = self.find_intersections(params[:3], params[3:], (x.min(), x.max()))
            for x_int in intersections:
                ax.axvline(x=x_int, color='blue', linestyle='--', lw=1)

            intersection_dict[var_name] = intersections

            ax.set_title(var_name)


        fig.tight_layout()
        if save == True:
            plt.savefig("positivityThreshold.png")
        else:
            plt.show()
        return intersection_dict


##test = Histograms()
##test.draw_histograms(True)
