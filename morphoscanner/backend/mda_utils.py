#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 30 21:01:00 2023

@author: lillux
"""

import MDAnalysis as mda

# TODO: integrate all the MDAnalysis functionalities
def make_universe(topology:str, trajectory:str=None, in_memory=False):
    '''
    Instantiate the MDAnalysis.Universe, that is used to parse the trajectory data
    
    Parameters
    ----------
    topology : str
        system path of topology file (gro for GROMACS).
    trajectory : str, optional
        system path of trajectory file (xtc or trr for GROMACS).
        can be provided as a single file or
        a list of consecutive trajectory files,
        as [part1.trr, part2.trr, ...].
        Can be left empty to load only the system configuration,
        in case you are loading a single structural file (eg. pdb)
    in_memory : bool, optional
        The default is False.
        Move data to memory for faster (~100x faster)
        frames coordinate retrival.
        It use much more memory (RAM).

    Returns
    -------
    universe : MDAnalysis.Universe()
        
    '''
    if trajectory:
        universe = mda.Universe(topology, trajectory, in_memory=in_memory)
    else:
        universe = mda.Universe(topology, in_memory=in_memory)
    return universe


def getAtomInUniverse(universe) -> dict:
    '''
    Read atoms in MDAnalysis.Universe, and count their occurences

    Parameters
    ----------
    universe : MDAnalysis.Universe
        An instance of MDAnalysis.Universe.

    Returns
    -------
    count_dict : dict
        A dict with key:value pairs:
            atomType:count
        It describe how many occurence of each atomType is present in a frame
        of the Universe.
    '''
    
    count_dict = {}
    for name in universe.atoms.names:
        try:
            count_dict[name] += 1
        except KeyError:
            count_dict[name] = 1
    return count_dict