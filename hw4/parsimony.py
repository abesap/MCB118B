import sys
import random
import fasta
import treelib
import nni
import numpy as np

INF = float('inf')
CHARS = ("A", "T", "C", "G")
def best_cost_character(tree, character, pos, mapping, memo): 
    if (tree,character,pos) in memo:
        #print("hitting a memoized one")
        return memo[(tree,character,pos)]
    root, left, right = tree
# base case
    if len(left) == 0: # leaf node
        #print(mapping[root][pos] ) d
        if mapping[root][pos] == character:
            memo[(tree,character,pos)] = 0
            return 0
        memo[(tree,character,pos)] = INF
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
                
                costs.append(cost)
        memo[(tree,character,pos)] = min(costs)
        return min(costs)

def best_cost_tree(tree, leaf_mapping):
    #print(tree)
    cost = 0 
    memo = {}
    length = len(list(leaf_mapping.values())[0])
    for pos in range(length):
        temp = []
        for char in CHARS:
            temp.append(best_cost_character(tree, char, pos, leaf_mapping, memo))
        cost+= min(temp)
    return cost

def helper(leaf_list):
    if len(leaf_list) == 2:
        return leaf_list
    else:
        leaf_list[-2] = (leaf_list[-2],leaf_list[-1])
        leaf_list = leaf_list[:-1]
        return helper(leaf_list)
def make_one_tree(leaf_list):
    leafs = leaf_list
    tree = helper(leafs)
    toreturn =""
    for x in str(tuple(tree)):
        if x != '\'':
            toreturn = toreturn + x 
    return treelib.newick_to_RLR(toreturn)
def max_parsimony(leaf_mapping, sample_size):
    best_tree = make_one_tree(list(leaf_mapping.keys()))
    best = best_cost_tree(best_tree, leaf_mapping)
    nnis= nni.get_all_nnis(best_tree)
    while len(nnis)>0:
        if len(nnis) > sample_size:
            select_nnis = random.sample(nnis, sample_size)
        else:
            select_nnis = nnis
        score_list = [best_cost_tree(x, leaf_mapping) for x in select_nnis]
        if best <= min(score_list):
            break
        best = min(score_list)
        best_tree = select_nnis[np.argmin(score_list)]
        nnis= nni.get_all_nnis(best_tree)
    return best_tree,best
def main():
    #print(len(sys.argv))
    if len(sys.argv) == 3:
        print(sys.argv)
        sample_size = int(sys.argv[2])
        leaf_mapping = dict(fasta.load(str(sys.argv[1])))
        best_tree,best = max_parsimony(leaf_mapping, sample_size)
        print(treelib.RLR_to_newick(best_tree))
        print(best)
	return best_tree,best
    
    else:
        raise Exception("%s: error: invalid command" % sys.argv[0])
if __name__ == '__main__':
    sys.exit(main())
 """"-cores 32 --trimmomatic /data/mcfadden/aquattrini/PROGRAMS/Trimmomatic-0.35/trimmomatic-0.35.jar
17:16 [asapirstein@purves:~/Desktop/118b/hw4]
25 % (('Anc0', ('Anc5', ('Anc4', ('Puma', (), ()), ('Anc6', ('Leopard', (), ()), ('BlackBear', (), ()))), ('Anc7', ('Homotherium', (), ()), ('AmericanCat', (), ()))), ('Anc8', ('Anc3', ('AfricanWildCat', (), ()), ('Anc2', ('ChineseDesertCat', (), ()), ('Anc1', ('Sabertooth', (), ()), ('SpottedHyena', (), ())))), ('Anc9', ('Anc10', ('Anc11', ('Lion', (), ()), ('Tiger', (), ())), ('CaveBear', (), ())), ('Anc12', ('Anc17', ('StripedHyena', (), ()), ('Anc16', ('BrownBear', (), ()), ('Anc15', ('Dog', (), ()), ('Wolf', (), ())))), ('Anc14', ('Cheetah', (), ()), ('Anc13', ('FelisCat',(), ()), ('WildCat', (), ()))))))), 455)

[3]    exit 1     python parsimony.py data/animals.fa 10 > 10.txt
17:16 [asapirstein@purves:~/Desktop/118b/hw4]
25 % save 5.txt
zsh: correct 'save' to 'osage' [nyae]? n
zsh: command not found: save
17:20 [asapirstein@purves:~/Desktop/118b/hw4]
26 %
17:20 [asapirstein@purves:~/Desktop/118b/hw4]
26 % (('Anc0', ('Anc9', ('Anc11', ('Lion', (), ()), ('Anc10', ('CaveBear', (), ()), ('Tiger', (), ()))), ('Anc6', ('Anc7', ('Homotherium', (), ()), ('Leopard', (), ())), ('Anc5', ('Anc4', ('AfricanWildCat', (), ()), ('Anc2', ('Puma',(), ()), ('ChineseDesertCat', (), ()))), ('Anc1', ('Sabertooth', (), ()), ('Anc3', ('SpottedHyena', (), ()), ('BlackBear', (), ())))))), ('Anc8', ('AmericanCat', (), ()), ('Anc13', ('Anc17', ('Anc14', ('Anc16', ('BrownBear', (), ()),('Anc15', ('Dog', (), ()), ('Wolf', (), ()))), ('StripedHyena', (), ())), ('Cheetah', (), ())), ('Anc12', ('WildCat', (), ()), ('FelisCat', (), ()))))), 457)

[5]  - exit 1     python parsimony.py data/animals.fa 75 > 75.txt
17:21 [asapirstein@purves:~/Desktop/118b/hw4]
26 % (('Anc0', ('Anc9', ('Anc11', ('Lion', (), ()), ('Anc10', ('CaveBear', (), ()), ('Tiger', (), ()))), ('Anc6', ('Anc7', ('Homotherium', (), ()), ('Leopard', (), ())), ('Anc5', ('Anc4', ('AfricanWildCat', (), ()), ('Anc2', ('Puma',(), ()), ('ChineseDesertCat', (), ()))), ('Anc1', ('Sabertooth', (), ()), ('Anc3', ('SpottedHyena', (), ()), ('BlackBear', (), ())))))), ('Anc8', ('AmericanCat', (), ()), ('Anc13', ('Anc17', ('Anc14', ('Anc16', ('BrownBear', (), ()),('Anc15', ('Dog', (), ()), ('Wolf', (), ()))), ('StripedHyena', (), ())), ('Cheetah', (), ())), ('Anc12', ('WildCat', (), ()), ('FelisCat', (), ()))))), 457)

[4]  - exit 1     python parsimony.py data/animals.fa 50 > 50.txt
17:21 [asapirstein@purves:~/Desktop/118b/hw4]
26 % (('Anc0', ('Anc9', ('Anc11', ('Lion', (), ()), ('Anc10', ('CaveBear', (), ()), ('Tiger', (), ()))), ('Anc6', ('Anc7', ('Homotherium', (), ()), ('Leopard', (), ())), ('Anc5', ('Anc4', ('AfricanWildCat', (), ()), ('Anc2', ('Puma',(), ()), ('ChineseDesertCat', (), ()))), ('Anc1', ('Sabertooth', (), ()), ('Anc3', ('SpottedHyena', (), ()), ('BlackBear', (), ())))))), ('Anc8', ('AmericanCat', (), ()), ('Anc13', ('Anc17', ('Anc14', ('Anc16', ('BrownBear', (), ()),('Anc15', ('Dog', (), ()), ('Wolf', (), ()))), ('StripedHyena', (), ())), ('Cheetah', (), ())), ('Anc12', ('WildCat', (), ()), ('FelisCat', (), ()))))), 457)

[6]  + exit 1     python parsimony.py data/animals.fa 100 > 100.txt
17:21 [asapirstein@purves:~/Desktop/118b/hw4]
26 % (('Anc0', ('Anc3', ('Anc6', ('Anc7', ('Anc1', ('Sabertooth', (), ()), ('Homotherium', (), ())), ('Anc5', ('BlackBear', (), ()), ('SpottedHyena', (), ()))), ('Leopard', (), ())), ('Anc10', ('Anc8', ('Anc2', ('Puma', (), ()), ('AmericanCat', (), ())), ('ChineseDesertCat', (), ())), ('Anc13', ('Anc15', ('Anc17', ('StripedHyena', (), ()), ('Anc14', ('BrownBear', (), ()), ('Anc16', ('Wolf', (), ()), ('Dog', (), ())))), ('Cheetah', (), ())), ('Anc9', ('CaveBear', (), ()), ('Anc11', ('Lion', (), ()), ('Tiger', (), ())))))), ('Anc4', ('Anc12', ('FelisCat', (), ()), ('WildCat', (),())), ('AfricanWildCat', (), ()))), 425)"""