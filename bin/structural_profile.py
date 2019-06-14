#!/usr/bin/env python3
"""
Module containing functions to calculate structural environment profiles of AAs
"""
import numpy as np
import Bio.PDB

def make_residue_distance_matrix(residues, ref_atom='CA'):
    """
    Generate a distance matrix from an iterator of Bio.PDB residues. Hetero residues
    are ignored. ref_atom is the atom to measure distances between, with no checking
    for atom presence (i.e. C-beta in Gly will fail)
    """
    residues = [r for r in residues if r.id[0] == ' ']

    dist = np.zeros((len(residues), len(residues)))

    for i, res1 in enumerate(residues):
        for j, res2 in enumerate(residues[i + 1:]):
            dist[i, j + i + 1] = dist[j + i + 1, i] = res1[ref_atom] - res2[ref_atom]

    return dist
