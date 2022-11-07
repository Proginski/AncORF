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
parser.add_argument("-names", required=True, help="list of names to retain", default=False)
parser.add_argument("-tree", required=True, help="newick file for the input phylogeny tree")
parser.add_argument("-out", required=True, help="output newick file")
parser.add_argument("-preserve_underscores", required = False, 
                    help="should underscores be kept when reading the tree, or considered as spaces", 
                    default = False, type=bool)
args = parser.parse_args()


# newick file for the phylogeny tree
tree_file = args.tree
# tree_file = "/home/me/STORE/Eukaryotes/PHYLO/SCER_SGD/GENOMES/Saccharomyces_species.nwk"
tree = dendropy.Tree.get(path=tree_file, schema='newick', preserve_underscores = args.preserve_underscores)

print("Original tree : ", tree)
print(tree.as_ascii_plot())

with open(args.names) as f:
    names = [name.strip() for name in f]
# names = ["Scer_SGD","Spar","Smik","Skud","Sarb","Sbay"]

taxa = [ tree.taxon_namespace[i].label for i in range(0,len(tree.taxon_namespace)) ]
extra_taxa = [ label for label in taxa if label not in names]
if len(extra_taxa) > 0 :
    print("Extra taxa are : {}".format(extra_taxa))
    
new_tree = tree.extract_tree_with_taxa_labels(labels=names)
print("New tree : ", new_tree)
print(new_tree.as_ascii_plot())


new_tree.write(path=args.out, schema="newick")