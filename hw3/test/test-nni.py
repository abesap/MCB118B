#! /usr/bin/env python

"""
Author      : Yi-Chieh Wu
Students    : Abel Sapirstein and Emily Petroni
Class       : HMC MCB 118b
Date        : 2017 Mar 31
Description : Tests for NNIs
"""

import os, sys

# add current working directory to path
sys.path.insert(0, os.getcwd())

import treelib
import nni

#========================================
# test code
# do not worry if you do not understand this code

import unittest

class NNITests(unittest.TestCase):
    
    def test_get_nnis(self):
        def helper(newick, values):
            """Check newick string 'newick' against NNI list 'values'."""
            tree = treelib.newick_to_RLR(newick)
            RLR_trees = nni.get_nnis(tree)
            
            # same number of nni trees?
            m = len(values)
            n = len(RLR_trees)
            self.assertEqual(m, n,
                             "tree %s --> found %d nni trees (correct %d)" % \
                             (newick, n, m))
            
            # same nni trees?
            newick_trees = [treelib.RLR_to_newick(RLR_tree) for RLR_tree in RLR_trees]
            for value in values:
                flag = value in newick_trees
                self.assertTrue(flag,
                                "tree %s --> missing nni tree %s" % \
                                (newick, value))
                newick_trees.remove(value)
        
        ### ========== TODO : START ========== ###
        # determine the output for the given tree
        # we have filled out one case for you
        
        # case 0 -- let newick = "(A,(B,C))"
        #helper("(A,(B,C))", [])
        # case 1 -- let newick = "((A,B),(C,D))"
        #helper("((A,B),(C,D))",
            #   ["((A,C),(B,D))", "((A,D),(C,B))"])
        # case 2 -- let newick = "((A,B),(C,(D,E)))"
        helper("((A,B),(C,(D,E)))", 
            ["((A,C),(B,(D,E)))", "((A,(D,E)),(C,B))"])
        # case 3 -- let newick = "(((A,B),(C,D)),((E,F),(G,H)))"
        helper("(((A,B),(C,D)),((E,F),(G,H)))", 
            ["(((A,B),(E,F)),((C,D),(G,H)))", "(((A,B),(G,H)),((E,F),(C,D)))"])
        ### ========== TODO : END ========== ###


#========================================
# main

def main():
    return unittest.main()

if __name__ == '__main__':
    sys.exit(main())
