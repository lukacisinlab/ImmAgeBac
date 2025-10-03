
class Constants:
    EXPECTED_DICTIONARY = {"CD57": (2, .1, 250, 3, .1, 125),
                                "CD14": (0.5, .1, 250, 1.7, .1, 125),
                                "CD4-1": (2, .1, 250, 5, .1, 125),
                                "CD8": (0.5, .1, 250, 3.5, .1, 125),
                                "CD45RA": (0.5, .1, 250, 1.8, .1, 125),
                                "CD3-1": (0.7, .1, 250, 2.7, .1, 125),
                                "CD3-2": (0.7, .1, 250, 3, .1, 125),
                                "CD127": (0.6, .1, 250, 2, .1, 125),
                                "CD16": (1, .1, 250, 3, .1, 125),
                                "CD161": (0, .1, 250, 1.2, .1, 125),
                                "CD27": (0.5, .1, 250, 3, .1, 125),
                                "CD28": (1, .1, 250, 2.1, .1, 125),
                                }
    MANUAL_HISTOGRAMS = ["CD4-2", "CD20", "CD38-1", "CD38-2", "CD185", "HLA-DR", "CD19",
                                                  "CD28", "CD25", "CD196", "CD185"]
    MISSING_DATA = ["CD33", "CD183"]

    MANUAL_INTERSECTION = {
        "HLA-DR": 1.7,
        "CD185": 0.6,
        "CD38-2": 2.1,
        "CD38-1": 2.6,
        "CD20": 0.8,
        "CD4-2": 1.1,
        "CD19": 3,
        "CD28": 1.6,
        "CD25": 1.0,
        "CD196": 1.7,
        "CD185": 0.8
    }
    CELL_TYPE_DICT = {
        "T": ["CD3-1 +", "CD3-2 +", "CD14 -"],
        "Effector": ["CD45RA +", "CCR7 -"],
        "Treg_cells": ["CD3-1 +", "CD3-2 +", "CD4-1 +", "CD4-2 +", "CD25 +", "CD127 +"],
        "Monocytes": ["CD14 +", "CD33 +"],
        "NK": ["CD3-1 -", "CD3-2 -", "CD16 +"],
        "Central_memory": ["CD45RA -", "CCR7 +"],
        "Effector_memory": ["CD45RA -", "CCR7 -"],
        "Lymphocytes": ["CD14 -", "CD33 -"],
        "Plasmablast": ["CD3-1 -", "CD3-2 -", "CD20 -", "CD27 +", "CD38-1 +", "CD38-2 -"],
        "Naive": ["CD45RA +", "CCR7 +"],
        "B": ["CD20 +", "CD19 +"]
    }

    CELLS = [
        "CD57 + CD8 + T",
        "Effector CD8 + T",
        "CD28 - CD8 + T",
        "Effector_memory CD8 + T",
        "Treg_cells",
        "Effector_memory CD4-1 + CD4-2 + T",
        "Monocytes",
        "NK",
        "CD57 + NK",
        "Central_memory CD4-1 + CD4-2 + T",
        "HLA-DR - CD38-1 + CD38-2 + CD4-1 + CD4-2 + T",
        "T",
        "Lymphocytes",
        "CD161 + NK",
        "CD8 + T",
        "Plasmablast",
        "Naive CD4-1 + CD4-2 + T",
        "B",
        "CD27 + CD8 + T",
        "CD161 - CD45RA + CD4-1 + CD4-2 + Treg_cells",
        "CD185 + CD4-1 + CD4-2 + T",
        "Naive CD8 + T"
    ]