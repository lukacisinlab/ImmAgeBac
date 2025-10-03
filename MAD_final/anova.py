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

from cellVisualizations import CellVisualization

class ANOVA(CellVisualization):
    def __init__(self):
        super().__init__()

    def create_melted_df(self):
        counts_per_donor = {donor: {} for donor in self.donors}

        for donor in self.donors:
            counts_dict = self.get_counts(donor)
            total_donor_count = self.adata[self.full_data.obs.donor == donor].n_obs

            for cell_type, count in counts_dict.items():
                counts_per_donor[donor][cell_type] = count / total_donor_count

        df = pd.DataFrame.from_dict(counts_per_donor, orient='index')
        df_melted = df.reset_index().melt(id_vars='index', var_name='Cell_Type', value_name='Count')
        df_melted.rename(columns={'index': 'Donor'}, inplace=True)

        df_melted["CMV_positive"] = df_melted["Donor"].apply(
            lambda x: True if self.metadata_df.loc[x].CMV == "positive" else False)
        df_melted["Gender"] = df_melted["Donor"].apply(
            lambda x: "Male" if self.metadata_df.loc[x].Gender == "M" else "Female")
        df_melted["Treatment"] = df_melted["Donor"].apply(lambda x: self.metadata_df.loc[x].Treatment)

        return df_melted


    def anova_get_p_value(self):
        counts_long_df = self.create_melted_df()
        cell_types = set(counts_long_df.iloc[:, 1])
        anova_var = {'Gender', 'Treatment', 'CMV_positive'}
        anova_to_p_vals_dict = {}
        p_vals = {}

        for column in anova_var:
            print("Now running anova for:", column, "\n")

            for cell_type in cell_types:
                curr_data = counts_long_df[counts_long_df['Cell_Type'] == cell_type]
                formula = f'Count ~ C({column})'
                model = ols(formula, data=curr_data).fit()
                anova_results = sm.stats.anova_lm(model, typ=2)

                p_value = anova_results["PR(>F)"][f"C({column})"]
                if p_value <= 0.05:
                    print("Null hypothesis rejected with cell type:", cell_type, "and p value:", str(round(p_value, 5)))
                p_vals[cell_type] = p_value

            print("\nFinished running anova for:", column, "\n")

            anova_to_p_vals_dict[column] = p_vals
            p_vals = {}
        return anova_to_p_vals_dict

    def save_anova_res(self):
        @contextmanager
        def redirect_stdout(to):
            original_stdout = sys.stdout
            sys.stdout = to
            try:
                yield
            finally:
                sys.stdout = original_stdout

        with open('anova_res.txt', 'w') as f:
            with redirect_stdout(f):
                anova_vals = self.anova_get_p_value()

##test = ANOVA()
##test.save_anova_res()