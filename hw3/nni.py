import treelib as tl
def get_nnis(tree):
    """INPUT
    tree- a RLR tree
    RETURNS
        - A list of NNI trees
    """
    root,left,right=tree
    if tl.is_leaf(left) or tl.is_leaf(right):
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
    else:
        output=landr_to_root(tree)
        return output
def landr_to_root(tree):
    """ Helper function to append to loop"""
    output=[tree]
    def left_to_root(tree):
        """" Inputs a Tree
        sets the nearest left node to root
        recurses on right and left subtrees"""
        if len(tree) == 0: #if tree is empty
            return 
        if tl.is_leaf(tree):
            return
        if tl.is_leaf(tree[1]): #if left is leaf
            return 
        root, left, right=tree
        loot , lleft, light = left
        new_tree1 = (loot, lleft, (root,light,right))
        new_tree2 = (loot, light, (root,lleft,right))
        [output.append(x)for x in [new_tree1,new_tree2, left_to_root(new_tree1), left_to_root(new_tree2)]]#, right_to_root(new_tree)]
    def right_to_root(tree):
        """" Inputs a Tree
        sets the nearest right node to root
        recurses on right and left subtrees"""
        if len(tree) == 0: #empty
            return    
        if tl.is_leaf(tree):
            return
        if tl.is_leaf(tree[2]): #right is leaf
            return  
        root, left, right=tree
        if len(right) == 0:
            return []
        right_root, right_left, right_right = right
        new_tree1 = (right_root, (root, left, right_left) ,right_right)
        new_tree2 = (right_root, (root, left, right_right) ,right_left)
        [output.append(x) for  x in [new_tree1,new_tree2, right_to_root(new_tree1), right_to_root(new_tree2)]]
    left_to_root(tree)
    right_to_root(tree)
    return list(filter(None, output))

def get_all_nnis(tree):
    """" INPUTS 
        tree    an RLR format tree to find NNI's of
        OUTPUTS
        output  a list of all NNI trees, in RLR format
        """"
    output=[]
    reroots=get_all_rerootings(tree)
    for tree in reroots:
        addee = get_nnis(tree)
        if addee not in output:
            output.extend(addee)
    return output