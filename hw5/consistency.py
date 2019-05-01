import numpy as np
import matplotlib.pyplot as plt
import treelib
import parsimony
import hash_tree
import sys
CHARS = np.array(("A", "T", "C", "G"))
def random_seq(n):
    return ''.join(list(np.random.choice(CHARS,n)))

def evolve_seq(seq,lambd):
    out=list(seq)
    for n in range(len(out)):
        changes=np.random.poisson(lambd,1)
        for y in range(changes):
            chars = [x for x in CHARS if x != out[n]]
            out[n]=str(np.random.choice(chars,1)[0])
    return ''.join(out)

def seq_diff(seq1,seq2):
    if (len(seq1) != len(seq2)):
        return "Sequences are not of the same length"
    s1=np.array(list(seq1))
    s2=np.array(list(seq2))
    return sum(s1!=s2)/float(len(seq1))

def jukes_cantor_model(lambd):
    return -float(3)/4*(np.exp(lambd*-(float(4)/3))-1)
def plot_seq_changes(max_changes, step_changes, seq_length):
    x = np.arange(0,max_changes, step_changes)
    seq = random_seq(seq_length)
    jc = map(jukes_cantor_model,x)
    o = [evolve_seq(seq, y) for y in x]
    obs = [seq_diff(y,seq)for y in o]
    plt.close()
    plt.plot(jc,'ro')
    plt.plot(obs)
    plt.xlabel("Number of Maximum Changes")
    plt.ylabel("Predicted number of Observable Changes per Site")
    plt.savefig("Seq"+str(seq_length)+"max"+str(max_changes)+".png")

def calc_parsimony_consistency(tree, mu,
                           seq_lengths, sample_size, reps):
    return sample_size

def correct_evolved_seq(time,mu,seq):
    lambd = time*mu 
    cor_lamb = jukes_cantor_model(lambd)
    return evolve_seq(seq, cor_lamb)
def prog_seq(tree,mu,leaf_dic, anc_seq):
    if not treelib.is_leaf(tree):
        root, left, right = tree
        name,time = root
        left_name = left[0][0]
        right_name = right[0][0]
        ls= correct_evolved_seq(left[0][1],mu,anc_seq)
        rs = correct_evolved_seq(right[0][1], mu, anc_seq)
        prog_seq(left, mu,leaf_dic, ls)
        prog_seq(right,mu,leaf_dic, rs)
    else:
        root, left, right = tree
        name,time = root
        leaf_dic[name] = anc_seq

def clean(tree):
    if treelib.is_leaf(tree):
        return (tree[0][0],tree[1],tree[2])
    else:
        return (tree[0][0],clean(tree[1]),clean(tree[2]))
def calc_parsimony_consistency(tree, mu, seq_lengths, sample_size, reps):
    out=[]
    for n in seq_lengths:
        count=0
        for y in range(reps):
            anc_seq=random_seq(n)
            leaf_dic={}
            prog_seq(tree,mu,leaf_dic, anc_seq)
            proposed = parsimony.max_parsimony(leaf_dic,sample_size)[0]
            real = clean(tree)
            if hash_tree.canonical(proposed, True)==hash_tree.canonical(real, True):
                count+=1
        out.append(count)
    return out
def plot_parsimony_consistency(tree, mu, seq_lengths, sample_size, reps):
    x = calc_parsimony_consistency(tree, mu, seq_lengths, sample_size, reps)
    plt.plot(seq_lengths,x)
    plt.xscale('log')
    plt.xlabel("Sequence Length")
    plt.ylabel("Number of Correct Occurences of "+ str(reps)+"with mu of"+ str(mu))
    plt.ylim(0, reps*1.1)
    plt.title("Sample Size vs Number of Correct Trees, of "+str(reps)+"Reps")
    plt.savefig(str(seq_lengths)+str(mu)+".png")
    plt.close()
    return x
def main():
    seq_lengths = [10,30,50,100,500,1000,5000]
    LONGBRANCH = (("",0),(("",0.01), (("A",0.01),(),()), (("B",0.4),(),())),(("",0.01), (("C",0.01),(),()), (("D",0.4),(),())))
    mu = [.00001,0.0005,0.0001,0.005,0.01,0.04,.1]
    increase_seq = False 
    for x in mu: 
        num_reps = 100
        vals = plot_parsimony_consistency(LONGBRANCH, x, seq_lengths, 10, num_reps)
        if vals[-1] is not max(vals) and vals[-1] != num_reps:#If inconsistent
            break
    
    else:
        raise Exception("%s: error: invalid command" % sys.argv[0])
if __name__ == '__main__':
    sys.exit(main())





