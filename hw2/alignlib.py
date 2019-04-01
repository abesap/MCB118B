#! /usr/bin/env python

"""
Author      : Yi-Chieh Wu and Matina Donaldson-Matasci (adapted from Eliot Bush and Brian Tjaden)
Class       : HMC MCB 118b
Date        : 2018 Mar 21
Description : Alignment Library
"""

import sys, fasta

#========================================
# globals

STOP = -1
UP_LEFT = 0
UP = 1
LEFT = 2


#========================================    
# helper functions

### ========== TODO : START ========== ###
# part b.1: add any helper functions here


### =========== TODO : END =========== ###


#========================================
# main alignment function

def align(v, w, match=1, mismatch=-2, indel=-3):
    """Compute the optimal pairwise alignment of two genomic sequences.
    
    Parameters
    --------------------
        v, w                 -- strings to align
        match                -- score for match (default=1)
        mismatch             -- score for mismatch (default=-2)
        indel                -- score for a gap/indel (default=-3)
    
    Return
    --------------------
        score_array          -- DP table (2D array)
        aligned_v, aligned_w -- aligned strings
        max_score            -- alignment score
    """
    
    ### ========== TODO : START ========== ###
    # you can modify anything in this TODO block
    #
    # setup: copy the body your global_align function from last week here
    #        modify it to also return the maximum alignment score
    #
    # part a.1: support local alignment
    # part b.2: support matrix scoring
    
    # initialize DP tables
    score_array = [[]]
    
    
    # complete DP table
    
    
    # find alignment score
    max_score = 0
    
    
    # traceback to create aligned strings
    aligned_v = ""
    aligned_w = ""
    
    
    # return
    return score_array, aligned_v, aligned_w, max_score
    ### =========== TODO : END =========== ###


#========================================
# main

def main():
    if len(sys.argv) == 6:
        
        # global or local alignment?
        ### ========== TODO : START ========== ###
        # part a.2: parse sys.argv[1]
        
        ### =========== TODO : END =========== ###
        
        # fixed or matrix scoring?
        match = None
        mismatch = None
        matrixfile = None
        if ',' in sys.argv[2]: # fixed
            # read the match,mismatch parameters as a single argument and parse
            toks = sys.argv[2].split(',')
            if len(toks) != 2:
                raise Exception("%s: error: third argument must be of form 'match,mismatch'")
            match = int(toks[0])
            mismatch = int(toks[1])
        else:                  # matrix
            ### ========== TODO : START ========== ###
            # part b.3: parse sys.argv[2]
            pass
            ### =========== TODO : END =========== ###
        
        # linear gap
        indel = int(sys.argv[3])
        
        # filenames
        vFileName = sys.argv[4]
        wFileName = sys.argv[5]
        v = fasta.load(vFileName)[0][1]
        w = fasta.load(wFileName)[0][1]
        
        # align
        ### ========== TODO : START ========== ###
        # part a.2: support global or local alignment
        # part b.3: support matrix scoring
        
        score_array, aligned_v, aligned_w, max_score = \
            align(v, w, match, mismatch, indel)
        ### =========== TODO : END =========== ###
        
        # print aligned strings
        print aligned_v
        print aligned_w
        
        # print alignment score
        print max_score
        
        # print DP table, use if numpy is not available
        # for L in score_array:
        #    print(L)
        
        # for prettier printing of arrays, convert to numpy array, then print
        import numpy as np
        print(np.array(score_array))
    
    else:
        raise Exception("%s: error: invalid command" % sys.argv[0])

if __name__ == '__main__':
    sys.exit(main())
