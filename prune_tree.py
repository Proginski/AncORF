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
parser.add_argument("-verbose", required=False, help="Whether or not to print a lot of information")
args = parser.parse_args()

# argsnames = "/run/user/4766/gvfs/smb-share:server=store.intra.i2bc.paris-saclay.fr,share=equipes/BIM/MEMBERS/paul.roginski/Eukaryotes/PHYLO/SCER_NCBI/test"
with open(args.names) as f:
    names = [name.strip() for name in f]
    
# newick file for the phylogeny tree
tree_file = args.tree


# names = ["Saccharomyces cerevisiae","Saccharomyces paradoxus","Saccharomyces bayanus"]
# tree_file = "/run/user/4766/gvfs/smb-share:server=store.intra.i2bc.paris-saclay.fr,share=equipes/BIM/MEMBERS/paul.roginski/Eukaryotes/PHYLO/ORFdate/liste_carvunis.nwk"
# names = ["Scer_NCBI","Spar","Sbay"]
# tree_file = "/home/me/STORE/Eukaryotes/PHYLO/SCER_NCBI/GENOMES/Saccharomyces_species_corr_doubleoutput.nwk"


tree = dendropy.Tree.get(path=tree_file, schema='newick', preserve_underscores=True)

taxa = [ tree.taxon_namespace[i].label for i in range(0,len(tree.taxon_namespace)) ]
extra_taxa = [ label for label in taxa if label not in names]

if args.verbose : print("Extra taxa are : {}".format(extra_taxa))

if len(extra_taxa) > 0 :
    
    kept_taxa = [ taxon for taxon in taxa if taxon not in extra_taxa ]
    
    if len(kept_taxa) > 0 :
        
        tree.prune_taxa_with_labels(extra_taxa)
        # Also need to correct the "taxon_namespace" attribute
        tree.taxon_namespace = [ tree.taxon_namespace[i] for i in range(0,len(tree.taxon_namespace)) if tree.taxon_namespace[i].label in names]
        
        if args.verbose : 
            print("Corrected tree")
            print(tree.as_ascii_plot())
            
    else : print("\nThere is no taxon to keep.\nExiting."); exit
        



tree.write(path=args.out, schema="newick")