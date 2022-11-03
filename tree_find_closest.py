#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 12:19:21 2022

@author: paul.roginski

The idea behind this script is that all possible most recent common ancestors 
of : a given species on the first hand, and its neighbors on the other hand,
are simply all the (node-represented) ancestors of the given species.
AND, all the ancestors-nodes of a given species, can be order by counting the 
number of edges between them and the root.
"""


import argparse
import dendropy

# Arguments parsing
parser = argparse.ArgumentParser()
parser.add_argument("-tree", required=True, help="newick file for the input phylogeny tree")
parser.add_argument("-focal", required=True, help="name of the focal species")
parser.add_argument("-names", required=True, help="list of species among which you want to determine which are (topologivally) the closest to the focal")
parser.add_argument("-out", required=False, help="output file name")
args = parser.parse_args()

  
# newick file for the phylogeny tree
tree_file = args.tree

# tree_file = "/home/me/STORE/Eukaryotes/PHYLO/SCER_NCBI/GENOMES/Saccharomyces_species_corr_doubleoutput.nwk"
# tree_file = "/home/me/STORE/Eukaryotes/PHYLO/SCER_NCBI/GENOMES/Saccharomyces_species.nwk"
# tree_file = "/home/me/STORE/Eukaryotes/PHYLO/SCER_NCBI/lolilol.nwk"
# tree_file = "/home/me/STORE/Eukaryotes/PHYLO/PV/GENOMES/cladogram.txt"
# tree_file = "/home/me/STORE/Eukaryotes/PHYLO/PV/GENOMES/Plasmo_tree.txt"

tree = dendropy.Tree.get(path=tree_file,
                         schema='newick',
                         preserve_underscores=True,
                         rooting='force-rooted')



print(tree.as_ascii_plot())

with open(args.names) as f:
    taxon_labels = [name.strip() for name in f]


# Retrieve in the taxon_namespace attribute of the tree object, the element 
# corresponding with the focal species.
focal_name = args.focal
# focal_name = "Scer_NCBI"
# focal_name = "Plasmodium_vivax"
print("focal name : {}".format(focal_name))
print("names : {}".format(taxon_labels))


index = [ i for i in range(0,len(tree.taxon_namespace)) if tree.taxon_namespace[i].label == focal_name ][0]
focal = tree.taxon_namespace[index]


# Set every edge length to 1. By doing so, we are sure the calc_node_root_distances can be used
# Besides, we do not care about the actual value of each edge lentgh in this case.
for node in tree.preorder_node_iter():
    node.edge.length = 1

# Adds attribute “root_distance” to each node, with value set to the sum of 
# edge lengths from the node to the root. Returns list of distances. 
tree.calc_node_root_distances(return_leaf_distances_only=False) 

# For each taxon provided in names, get its most recent common ancestor with 
# the focal species (= a node), and get its distance to the tree.
root_distance = { name:tree.mrca(taxon_labels=[focal_name, name]).root_distance for name in taxon_labels if name != focal_name }
print(root_distance)
max_value = max(root_distance.values())

print("max number of node-represented-ancestors in common with {} = {}".format(focal_name,max_value))

closest_ones = [ name for name in root_distance.keys() if root_distance[name] == max_value ] 

print(closest_ones)
        
if args.out : 
    with open(args.out, 'w') as fp:
        fp.write("%s\n" % focal_name)
        for name in closest_ones:
            fp.write("%s\n" % name)
    print('{} generated.'.format(args.out))