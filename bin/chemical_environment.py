#!/usr/bin/env python3
"""
Module containing functions to calculate structural environment profiles of AAs
"""

# TODO Tests to make sure these are producing correct profiles

from collections import defaultdict
import numpy as np
from Bio.SeqUtils import seq1#, seq3
from Bio.Alphabet.IUPAC import protein as protein_alphabet

AA_INDEX = {aa: i for i, aa in enumerate(protein_alphabet.letters)}

def k_nearest_residues_profile(res_index, residues, distance_matrix=None, k=10):
    """
    Calculate chemical environment of a residue as the make up of the k nearest AAs.
    Hetero atoms are included so must be dropped separately if desired.

    res_index: residue to process
    residues: list of all residues to consider
    distance_matrix: numpy matrix of distances between residues, with rows/columns in
                     that order
    k: count the k nearest residues

    returns: structural profile of the residue (np.array)
    """

    if distance_matrix is None:
        distance_matrix = residue_distance_matrix(residues)

    dists = distance_matrix[res_index,]
    dists[res_index] = dists.max() + 1 # don't want the residue itself

    nearest_k = [residues[i] for i in np.argpartition(dists, k)[:k]]

    counts = defaultdict(lambda: 0)
    for i in nearest_k:
        counts[seq1(i.get_resname())] += 1

    return np.array([counts[aa] for aa in protein_alphabet.letters])

# TODO good distance to choose as default?
def within_distance_profile(res_index, residues, distance_matrix=None, max_dist=10):
    """
    Calculate chemical environment of a residue as the residues within max_dist
    angstroms of it. Hetero atoms are included so must be dropped separately if desired.

    res_index: residue to process
    residues: list of all residues to consider
    distance_matrix: numpy matrix of distances between residues, with rows/columns in
                     that order
    max_dist: maximum distance to count within

    returns: structural profile of the residue (np.array)
    """

    if distance_matrix is None:
        distance_matrix = residue_distance_matrix(residues)

    dists = distance_matrix[res_index,]
    dists[res_index] = max_dist + 1 # ignore residue itself

    res_within_dist = [residues[i] for i in np.argwhere(dists < max_dist)[:, 0]]

    counts = defaultdict(lambda: 0)
    for i in res_within_dist:
        counts[seq1(i.get_resname())] += 1

    return np.array([counts[aa] for aa in protein_alphabet.letters])

def get_profiles(residues, profile_func=k_nearest_residues_profile, drop_hetero=True):
    """
    Calculate chemical environment profiles for a list of residues

    residues: list of Bio.PDB residues to retrieve profiles for
    profile_func: function to calculate profiles with sig f(res_ind, res_list, dist_matrix)
    drop_hetero: filter out hetero atoms (bool)

    returns: profiles (list of np.array)
    """
    if drop_hetero:
        residues = drop_hetero_atoms(residues)

    distance_matrix = residue_distance_matrix(residues)

    return [profile_func(i, residues, distance_matrix) for i in range(len(residues))]

def residue_distance_matrix(residues, ref_atom='CA'):
    """
    Generate a distance matrix from an iterable of Bio.PDB residues.
    There is no checking whether the specified atom makes sense (i.e. using C-beta
    with Gly will fail)

    residues: iterable of Bio.PDB residues to get distances from
    ref_atom: atom to measure distances from

    returns: distance matrix, rows and columns in order of residues (np.array)
    """

    dist = np.zeros((len(residues), len(residues)))

    for i, res1 in enumerate(residues):
        for j, res2 in enumerate(residues[i + 1:]):
            dist[i, j + i + 1] = dist[j + i + 1, i] = res1[ref_atom] - res2[ref_atom]

    return dist

def drop_hetero_atoms(residues):
    """
    Drop hetero atoms.

    residues: iterable of Bio.PDB residues

    returns: residues without hetero atoms (list)
    """
    return [r for r in residues if r.id[0] == ' ']
