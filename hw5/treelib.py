#!/usr/bin/env python

"""
Author      : Yi-Chieh Wu
Class       : HMC MCB 118b
Date        : 2018 Mar 28
Description : Tree Library
"""

import sys
import re, ast

#========================================
# input / output

def RLR_to_newick(rlr_tree):
    """Convert RLR tree to newick string.
    
    Parameters
    --------------------
        rlr_tree -- RLR tree
    
    Return
    --------------------
        tree     -- newick string
    """
    
    root, left, right = rlr_tree
    
    if len(left) == 0:
         # leaf node
        return root
    else:
        # internal node
        return "(%s,%s)" % (RLR_to_newick(left), RLR_to_newick(right))


def newick_to_RLR(string, use_list=False):
    """Convert newick string to RLR tree.
    
    Parameters
    --------------------
        string   -- newick string
        use_list -- True to encode node as list rather than tuple
    
    Return
    --------------------
        tree     -- RLR tree
    """
    
    def helper(newick_tree, ancs=None):
        """Convert newick tree to RLR tree.
        
        Parameters
        --------------------
          newick_tree -- newick tree
          ancs        -- list of ancestor numbers
        
        Return
        --------------------
          tree        -- RLR tree
        """
        
        if use_list:
            make_node = lambda toks: list(toks)
        else:
            make_node = lambda toks: toks
        
        if ancs is None:
            ancs = [-1]
        
        if type(newick_tree) != tuple:
            return (newick_tree, (), ())
        else:
            anc = max(ancs) + 1
            ancs.append(anc)
            return make_node(("Anc" + str(anc),
                              helper(newick_tree[0], ancs),
                              helper(newick_tree[1], ancs)))
    
    # wrap names in single quotes 
    string2 = re.sub(r"(\w+)", r"'\1'", string)
    
    # convert to tuple of tuples
    string3 = ast.literal_eval(string2)
    
    # convert to RLR tree
    tree = helper(string3)
    
    return tree

parse_newick = newick_to_RLR


def read_tree(fn, use_list=False):
    """Read tree from file.
    
    Parameters
    --------------------
        fn       -- filename
        use_list -- True to encode node as list rather than tuple
    
    Return
    --------------------
        tree     -- RLR tree
    """
    
    with open(fn, "r") as f:
        line = f.read()
        line = line.replace('\n', '')   # strip new lines
        line = re.sub(r'\s+', '', line) # strip white space
        if line[-1] == ";":
            line = line[:-1]            # strip semicolon
        return parse_newick(line, use_list=use_list)


#========================================
# tree functions

def get_name(tree):
    """Returns the name of the root node of 'tree'."""
    return tree[0]


def set_name(tree, label):
    """Sets the name of the root node of 'tree' to 'label'."""
    tree[0] = label


def is_leaf(tree):
    """Returns True if the root node of 'tree' is a leaf node and False otherwise."""
    return len(tree[1]) == 0


def left_subtree(tree):
    """Returns the left subtree of 'tree'."""
    return tree[1]


def right_subtree(tree):
    """Returns the right subtree of 'tree'."""
    return tree[2]


def swap_children(tree):
    """Swaps the left and right subtree of 'tree'."""
    tree[1], tree[2] = tree[2], tree[1]

