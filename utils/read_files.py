import os
import sys
import settings
import pandas as pd


def process_smina_atom_terms_block(block):
    lines = block.strip().split('\n')[1:-1]  # Pomijamy "atomid" i "END"
    data = [line.split() for line in lines]

    atomid = [int(row[0]) for row in data]
    element = [row[1] for row in data]
    coordinates = [row[2].replace('<','').replace('>','').split(',') for row in data]
    scoring_components = [list(map(float, row[3:])) for row in data]

    # Zliczanie sum dla każdej kolumny składowej funkcji skorującej
    sum_gauss_1 = sum(row[0] for row in scoring_components)
    sum_gauss_2 = sum(row[1] for row in scoring_components)
    sum_repulsion = sum(row[2] for row in scoring_components)
    sum_hydrophobic = sum(row[3] for row in scoring_components)
    sum_non_dir_h_bond = sum(row[4] for row in scoring_components)

    scoring_components_sum = [sum_gauss_1, sum_gauss_2, sum_repulsion, sum_hydrophobic, sum_non_dir_h_bond]

    return atomid, element, coordinates, scoring_components_sum