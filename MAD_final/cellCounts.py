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

from histograms import Histograms

class CellCounts(Histograms):
    def __init__(self):
        super().__init__()

    def cell_into_markers(self):
        markers = {}
        for cell in self.CELLS:
            res = ""
            args = cell.split()
            for arg in args:
                if arg in self.cell_type_dict:
                    res += " ".join(self.cell_type_dict[arg]) + " "
                else:
                    res += arg
                    if arg != args[-1]:
                        res += " "
            markers[cell] = res
        return markers

    def print_markers(self):
        markers = self.cell_into_markers()
        for cell, marker in markers.items():
            total_width = 79
            cell_len = len(cell)
            buffer_size = total_width - cell_len

            left_buffer_size = buffer_size // 2
            right_buffer_size = buffer_size - left_buffer_size

            print("-" * left_buffer_size + cell + "-" * right_buffer_size)
            print(marker, "\n")

    def getIntersections(self):
        variables = self.expected_dictionary.keys()
        n_rows = math.ceil(len(variables) / 2)
        n_cols = 2

        intersection_dict = {}
        for i, var_name in enumerate(variables):
            data = self.adata.X[:, self.adata.var.index == var_name]
            if data.size == 0:
                print(var_name + " is empty\n")
                continue
            y, x = np.histogram(data, bins=100)
            x = (x[1:] + x[:-1]) / 2

            expected = self.expected_dictionary[var_name]
            try:
                params, cov = curve_fit(self.bimodal, x, y, expected)
            except RuntimeError:
                print(var_name)
                continue
            sigma = np.sqrt(np.diag(cov))
            x_fit = np.linspace(x.min(), x.max(), 500)
            intersections = self.find_intersections(params[:3], params[3:], (x.min(), x.max()))

            intersection_dict[var_name] = intersections


        return intersection_dict

    def count_rows_with_multiple_conditions(self, adata, conditions, CCR7_replacement = "CCR7"):
        intersection = self.intersection
        mask = np.ones(adata.shape[0], dtype=bool)
        args = conditions.split()

        for i in range(len(args)):
            column_name = args[i]

            if column_name in ["+", "-"]:
                continue
            if i == len(args) - 1:
                break
            if column_name == "CCR7" and CCR7_replacement != "CCR7":
                column_name = CCR7_replacement

            if column_name not in adata.var_names:
                print("MARKER NOT FOUND:", column_name)
                raise ValueError(f"Column '{column_name}' not found in adata.var_names")

            operator = args[i + 1]
            column_index = adata.var_names.get_loc(column_name)
            try:
                value = intersection[column_name][0]
                if column_name == "CD25":
                    value = intersection[column_name][1]

            except KeyError:
                value = self.manual_intersection[column_name]

            column_data = adata.X[:, column_index]

            if operator == "+":
                mask &= (column_data > value)
            elif operator == "-":
                mask &= (column_data < value)

        count = np.sum(mask)

        return count

    def print_counts(self, donor = None):
        data = self.adata
        if donor != None:
            data = self.adata[self.full_data.obs.donor == donor]
        markers = self.cell_into_markers()
        self.intersection = self.getIntersections()
        for cell, marker in markers.items():
            try:
                res = self.count_rows_with_multiple_conditions(data, marker, "CD27")
                print("Success at cell type:", cell, "\nWith markers:", marker, "\nCount:", res)
            except:
                print("Failed at cell type:", cell, "\nWith markers:", marker)
            print("\n")

#test = CellCounts()
#test.print_markers()

#test = CellCounts()
#test.print_counts("P1")