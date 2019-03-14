"""
Author      : Matina Donaldson-Matasci, Yi-Chieh Wu, Eliot Bush
Class       : HMC MCB 118b
Date        : 2019 Mar 12
Description : Fasta module
"""

import sys

#========================================
# functions

def load(filename):
    """Load fasta or multifasta file.
    
    Parameters
    --------------------
      fn -- name of file containing sequences in fasta format
    
    Return
    --------------------
      fa -- list of tuples (name, sequence)
    """
    
    ### ========== TODO : START ========== ###
    fa = []
    
    # YOUR CODE HERE
    with open(filename,"r") as f:
        for line in f:
            pass
    
    return fa
    ### ========== TODO : END ========== ###


#========================================
# main

def main():
    if len(sys.argv) == 2:
        # process input
        fn = sys.argv[1]
        
        # load file
        fa = load(fn)
        
        # print back out what was read
        for key, seq in fa:
            print ">%s" % key
            print seq
    
    else:
        raise Exception("%s: error: invalid command" % sys.argv[0])

if __name__ == '__main__':
    sys.exit(main())
