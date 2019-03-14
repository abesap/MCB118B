#! /usr/bin/env python

"""
Authors     : Matina Donaldson-Matasci, Yi-Chieh Wu, Eliot Bush
Class       : HMC MCB 118b
Date        : 2019 Mar 12
Description : Global Alignment
"""

import sys, fasta

#========================================
# globals

STOP = -1
UP_LEFT = 0
UP = 1
LEFT = 2

#========================================
# functions

def global_align(v, w, match, mismatch, indel):
    """Given two strings and alignment parameters,
    return the global alignment DP table.
    
    Parameters
    --------------------
      v, w                   -- strings to align
      match, mismatch, indel -- alignment parameters
    
    Return
    --------------------
      score_array            -- DP table (2D array)
      aligned_v              -- aligned string
      aligned_w              -- aligned string
    """
    
    ### ========== TODO : START ========== ###
    # you can modify anything in this TODO block
    
    # initialize DP tables
    score_array = [[]]
    
    
    # complete DP table
    
    
    # traceback to create aligned strings
    aligned_v = ""
    aligned_w = ""
    
    
    # return
    return score_array, aligned_v, aligned_w
    ### ========== TODO : END ========== ###


#========================================
# main

def run(match, mismatch, indel, vFileName, wFileName):
    # load each sequence in fasta files
    v = fasta.load(vFileName)[0][1]
    w = fasta.load(wFileName)[0][1]
    
    # compute DP table
    scoreArray = global_align(v, w, match, mismatch, indel)
    
    return scoreArray


def main():
    if len(sys.argv) == 6:
        # process input
        match = int(sys.argv[1])
        mismatch = int(sys.argv[2])
        indel = int(sys.argv[3])
        vFileName = sys.argv[4]
        wFileName = sys.argv[5]
        
        # compute DP table
        scoreArray, aligned_v, aligned_w = run(match, mismatch, indel, vFileName, wFileName)
        
        # print DP table, use if numpy is not available
        #for L in scoreArray:
        #    print(L)
        
        # for prettier printing of arrays, convert to numpy array, then print
        import numpy as np
        print(np.array(scoreArray))
        
        # print aligned strings
        print aligned_v
        print aligned_w
    
    else:
        raise Exception("%s: error: invalid command" % sys.argv[0])

if __name__ == '__main__':
    sys.exit(main())
