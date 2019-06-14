#!/usr/bin/env python3
"""
Module containing functions to calculate structural environment profiles of AAs
"""

# TODO Tests to make sure these are producing correct protiles

from collections import defaultdict
import numpy as np
import Bio.PDB
from Bio.SeqUtils import seq1, seq3
from Bio.Alphabet.IUPAC import protein as protein_alphabet

AA_INDEX = {aa: i for i, aa in enumerate(protein_alphabet.letters)}

def residue_distance_matrix(residues, ref_atom='CA'):
    """
    Generate a distance matrix from an iterable of Bio.PDB residues.
    There is no checking whether the specified atom makes sense (i.e. using C-beta
    with Gly will fail)

    residues: iterable of Bio.PDB residues to get distances from
    ref_atom: atom to measure distances from
    """

    dist = np.zeros((len(residues), len(residues)))

    for i, res1 in enumerate(residues):
        for j, res2 in enumerate(residues[i + 1:]):
            dist[i, j + i + 1] = dist[j + i + 1, i] = res1[ref_atom] - res2[ref_atom]

    return dist

# TODO exclude hetero atoms properly
# currently doesn't work if any are included with a distance matrix
def k_nearest_residues_profile(res_index, residues, distance_matrix=None, k=10):
    """
    Calculate chemical environment of a residue as the make up of the k nearest AAs.
    Hetero residues are excluded.

    res_index: residue to process
    residues: list of all residues to consider
    distance_matrix: numpy matrix of distances between residues, with rows/columns in
                     that order
    k: count the k nearest residues
    """

    if distance_matrix is None:
        residues = [r for r in residues if r.id[0] == ' ']
        distance_matrix = residue_distance_matrix(residues)

    dists = distance_matrix[res_index,]
    dists[res_index] = dists.max() + 1 # don't want the residue itself

    nearest_k = [residues[i] for i in np.argpartition(dists, k)[:k]]

    counts = defaultdict(lambda: 0)
    for i in nearest_k:
        counts[seq1(i.get_resname())] += 1

    return np.array([counts[aa] for aa in protein_alphabet.letters])

# TODO exclude hetero properly
# TODO good distance to choose as default?
def within_distance_profile(res_index, residues, distance_matrix=None, max_dist=10):
    """
    Calculate chemical environment of a residue as the residues within max_dist
    angstroms of it. Hetero residues are excluded

    res_index: residue to process
    residues: list of all residues to consider
    distance_matrix: numpy matrix of distances between residues, with rows/columns in
                     that order
    max_dist: maximum distance to count within
    """

    if distance_matrix is None:
        residues = [r for r in residues if r.id[0] == ' ']
        distance_matrix = residue_distance_matrix(residues)

    dists = distance_matrix[res_index,]
    dists[res_index] = max_dist + 1 # ignore residue itself

    res_within_dist = [residues[i] for i in np.argwhere(dists < max_dist)[:,0]]

    counts = defaultdict(lambda: 0)
    for i in res_within_dist:
        counts[seq1(i.get_resname())] += 1

    return np.array([counts[aa] for aa in protein_alphabet.letters])