#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 16:40:15 2022

@author: paul.roginski
"""


import argparse
import dendropy
import pandas as pd


# Arguments parsing
parser = argparse.ArgumentParser()
parser.add_argument("-names", required=True, help="csv file matching names (col1) and phylip_names (col2)", default=False)
parser.add_argument("-tree", required=True, help="newick file for the input phylogeny tree")
parser.add_argument("-out", required=True, help="output newick file")
args = parser.parse_args()
    
# newick file for the phylogeny tree
tree_file = args.tree

names = pd.read_csv(args.names, names = ["name", "phylip_name"])

# names = ["Saccharomyces cerevisiae","Saccharomyces paradoxus","Saccharomyces bayanus"]
# tree_file = "/run/user/4766/gvfs/smb-share:server=store.intra.i2bc.paris-saclay.fr,share=equipes/BIM/MEMBERS/paul.roginski/Eukaryotes/PHYLO/ORFdate/liste_carvunis.nwk"
# names = ["Scer_NCBI","Spar","Sbay"]
# names = pd.read_csv("/run/user/4766/gvfs/smb-share:server=store.intra.i2bc.paris-saclay.fr,share=equipes/BIM/MEMBERS/paul.roginski/Eukaryotes/PHYLO/SCER_NCBI/test", names = ["name", "phylip_name"])
# tree_file = "/run/user/4766/gvfs/smb-share:server=store.intra.i2bc.paris-saclay.fr,share=equipes/BIM/MEMBERS/paul.roginski/Eukaryotes/PHYLO/SCER_NCBI/GENOMES/Saccharomyces_species.nwk"


tree = dendropy.Tree.get(path=tree_file, schema='newick', preserve_underscores=True)

for name in names["name"] :
    
    node_to_change = tree.find_node_with_taxon_label(name)
    node_to_change.taxon.label = names.loc[names['name'] == name, 'phylip_name'].item()


# print(tree)
# tree.write(path="/run/user/4766/gvfs/smb-share:server=store.intra.i2bc.paris-saclay.fr,share=equipes/BIM/MEMBERS/paul.roginski/Eukaryotes/PHYLO/SCER_NCBI/lol", schema="newick")
tree.write(path=args.out, schema="newick")