#!/usr/bin/env python

"""
Author      : Yi-Chieh Wu
Class       : HMC MCB 118b
Date        : 2018 Mar 22
Description : Nearest Neighbor Interchange
"""

import sys
import treelib

#========================================
# NNI functions

def get_nnis(tree):
    """Given a tree ((A,B),(C,D), find all NNI neighbors
    along the edge connecting (A,B) and (C,D).
    
    Parameters
    --------------------
        tree      -- RLR tree
    
    Return
    --------------------
        nni_trees -- list of neighboring RLR trees
    """
    
    nni_trees = []
    
    root, left, right = tree
    lt_name, A, B = left
    rt_name, C, D = right
    
    # nothing to swap -- (A,(B,C)) or ((A,B),C)
    if len(A) == 0 or len(C) == 0:
        return nni_trees
    
    # swap B and C -> ((A,C),(B,D))
    nni_tree = (root,
                      (lt_name, A, C),
                      (rt_name, B, D)
               )
    nni_trees.append(nni_tree)
    
    # swap B and D -> ((A,D),(C,B))
    nni_tree = (root,
                     (lt_name, A, D),
                     (rt_name, C, B)
               )
    nni_trees.append(nni_tree)
    
    return nni_trees


def get_all_rerootings(tree):
    """Given a tree, find all rerootings.
    
    Parameters
    --------------------
        tree  -- RLR tree
    
    Return
    --------------------
        trees -- list of rerooted RLR trees
    """
    
    def reroot_left(t):
        """Helper function to move root to left subtree.
        ((A,B),right) => (A,(B,right)) and (B,(A,right))"""
        root, left, right = t
        lt_name, A, B = left
        if len(A) == 0: # cannot move root
            return
        
        # create rerooted trees
        t1 = (root,
                    A,
                    (lt_name, B, right)
             )
        t2 = (root,
                    B,
                    (lt_name, A, right)
             )
        
        # add to rerooted trees
        trees.append(t1)
        trees.append(t2)
        
        # recur to keep moving down left subtree
        reroot_left(t1)
        reroot_left(t2)
    
    def reroot_right(t):
        """Helper function to move root to right subtree.
        (left,(C,D)) => ((left,D),C) and ((left,C),D)"""
        root, left, right = t
        rt_name, C, D = right
        if len(C) == 0: # cannot move root
            return
        
        # create rerooted trees
        t1 = (root,
                   (rt_name, left, D),
                   C
             )
        t2 = (root,
                   (rt_name, left, C),
                   D
             )
        # add to rerooted trees
        trees.append(t1)
        trees.append(t2)
        
        # recur to keep moving down right subtree
        reroot_right(t1)
        reroot_right(t2)
    
    trees = [tree]
    reroot_left(tree)
    reroot_right(tree)
    return trees


def get_all_nnis(tree):
    """Given a tree ((A,B),(C,D), find all NNI neighbors
    along all edges.
    
    Parameters
    --------------------
        tree      -- RLR tree
    
    Return
    --------------------
        nni_trees -- list of neighboring RLR trees
    """
    
    nni_trees = []
    trees = get_all_rerootings(tree)
    for current_tree in trees:
        current_nni_trees = get_nnis(current_tree)
        nni_trees.extend(current_nni_trees)
    return nni_trees


#========================================
# main

def main():
    if len(sys.argv) == 2:
        fn = sys.argv[1]
        
        tree = treelib.read_tree(fn)
        nni_trees = get_all_nnis(tree)
        for nni_tree in nni_trees:
            print treelib.RLR_to_newick(nni_tree)
    
    else:
        raise Exception("%s: error: invalid command" % sys.argv[0])

#if __name__ == '__main__':
#    sys.exit(main())

