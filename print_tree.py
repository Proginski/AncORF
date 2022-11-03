#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 12:19:21 2022

@author: paul.roginski
"""

import argparse
import dendropy


# Arguments parsing
parser = argparse.ArgumentParser()
parser.add_argument("-tree", required=True, help="newick file for the input phylogeny tree")
args = parser.parse_args()

  
# newick file for the phylogeny tree
tree_file = args.tree

# tree_file = "/home/me/STORE/Eukaryotes/PHYLO/SCER_NCBI/GENOMES/Saccharomyces_species_corr_doubleoutput.nwk"


tree = dendropy.Tree.get(path=tree_file, schema='newick', preserve_underscores=True)

print(tree.as_ascii_plot())
        