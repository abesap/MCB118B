#!/usr/bin/env python

"""
Author      : Yi-Chieh Wu
Studens     : Abel Sapirstein and Emily Patroni
Class       : HMC MCB 118b
Date        : 2018 Mar 22
Description : Tree Representations

You do not need to submit this file. It is an optional extension for the exercises.
"""

import sys
from treelib import *

#========================================
# representation functions

def canonical(tree):
    """Returns a canonical tree.
    
    Parameters
    --------------------
        tree  -- RLR tree
    """
    
    ### ========== TODO : START ========== ###
    # add your code here
    # you can modify 'tree' as much as you like
    
    ### ========== TODO : END ========== ###
    
    return tree

#========================================
# main

def main():
    if len(sys.argv) == 2:
        fn = sys.argv[1]
        
        tree = read_tree(fn, use_list=True)
        print RLR_to_newick(tree)
        
        ctree = canonical(tree)
        print RLR_to_newick(ctree)
    
    else:
        raise Exception("%s: error: invalid command" % sys.argv[0])

if __name__ == '__main__':
    sys.exit(main())

