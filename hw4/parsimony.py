import sys
import random
import fasta
import treelib
import nni
INF = float('inf')
CHARS = ("A", "T", "C", "G")
def best_cost_character(tree, character, pos, mapping, memo): 
    if (tree,character,pos) in memo.keys():
        return memo[(tree,character,pos)]
    root, left, right = tree
# base case
    if len(left) == 0: # leaf node
        #print(mapping[root][pos] ) d
        if mapping[root][pos] == character:
            return 0
        return INF
        # recursion
    else:
        costs=[]
        for left_char in CHARS:
            left_cost = best_cost_character(left, left_char, pos, mapping, memo) 
            for right_char in CHARS:
                right_cost = best_cost_character(right, right_char, pos, mapping, memo)
                # add cost along immediate branches
                cost=left_cost + right_cost + int(character != left_char)+ int(character != right_char)
                memo[(tree,character,pos)] = cost
                costs.append(cost)
    return min(costs)
def best_cost_tree(tree, leaf_mapping):
    anc = tree
    while anc[0] not in leaf_mapping.keys():
        anc=anc[1]
    cost = 0 
    memo = {}
    length = len(leaf_mapping[anc[0]])
    for pos in range(length):
        temp = []
        for char in CHARS:
            
            temp.append(best_cost_character(tree, char, pos, leaf_mapping, memo))
        cost+= min(temp)
        print(pos,cost)
    return cost