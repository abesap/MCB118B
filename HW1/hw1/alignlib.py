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
    nrow=len(v)+1
    ncol=len(w)+1
    score_array = [[0 for n in xrange(ncol)] for n in xrange(nrow)]
    # complete DP table
    for row in xrange(nrow):
        for col in xrange(ncol):
            score = 0
            if col > 0 and row > 0:
                diag = score_array[row-1][col-1]
                up = score_array[row][col-1]+indel
                left = score_array[row-1][col]+indel
                if v[row-1]==w[col-1]:
                    diag = diag+match
                else:
                    diag= diag+mismatch
                score = max(diag,up,left)
                ### enter conditionals
            else:
                score= indel*(row+col)
            score_array[row][col]= score
    
    # traceback to create aligned strings
    aligned_v = ""
    aligned_w = ""
    last_row = nrow-1
    last_col = ncol-1
    v_pos = nrow-1
    w_pos = ncol-1
    while last_row>0 or last_col>0:
        current=score_array[last_row][last_col]
        diag = score_array[last_row-1][last_col-1]
        if v[last_row-1]==w[last_col-1]:
            diag = diag+match
        else:
            diag= diag+mismatch
        up = score_array[last_row-1][last_col]+indel
        left = score_array[last_row][last_col-1]+indel
        val= max(diag,up,left)
        if val == diag:
            aligned_v = v[v_pos-1]+aligned_v
            aligned_w = w[w_pos-1]+aligned_w
            v_pos = v_pos-1
            w_pos = w_pos-1
            last_row = last_row-1
            last_col = last_col-1
        elif val == up:
            aligned_w = "-"+aligned_w
            aligned_v = v[v_pos-1]+aligned_v
            v_pos-=1
            last_row = last_row-1

        else:
            aligned_v = "-"+aligned_v
            aligned_w = w[w_pos-1]+aligned_w
            w_pos-=1
            last_col = last_col-1
    
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

#if __name__ == '__main__':
#   sys.exit(main())
