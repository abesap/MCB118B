import treelib as tl
def get_nnis(tree):
    """INPUT
    tree- a RLR tree
    RETURNS
        - A list of NNI trees
    """
    root,left,right=tree
    if left[1]==():
        return []
    else:
        one = (root, (left[0],left[1], right[1]), (right[0],left[2],right[2])) 
        two = (root, (left[0],left[1], right[2]), (right[0], right[1], left[2]) )
        #print(one)
        return [one,two]

def get_all_rerootings(tree):
    """INPUT
    tree- a RLR tree
    RETURNS
        - A list of NNI trees
    """
    if tree == (): #empty tree
        return 1
    elif tl.is_leaf(tree[1]): #left is leaf
        return 1
    elif tl.is_leaf(tree[2]): #right is leaf
        return 1
    else:
        left =left_to_root(tree)
        right = right_to_root (tree)
        ll = get_all_rerootings(left)
        rr = get_all_rerootings(right)
        return [tree, left, right, ll, rr]
def left_to_root(tree):
    """" Inputs a Tree
    sets the nearest left node to root"""
    root, left, right=tree
    loot , lleft, light = left
    new_tree = (loot, lleft, (root,light,right))
    return new_tree
def right_to_root(tree):
    """" Inputs a Tree
    sets the nearest left node to root"""
    root, left, right=tree
    right_root, right_left, right_right = right
    new_tree = (right_root, (root, left, right_left) ,right_right)
    return new_tree