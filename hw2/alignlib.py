#! /usr/bin/env python

"""
Author      : Yi-Chieh Wu and Matina Donaldson-Matasci (adapted from Eliot Bush and Brian Tjaden)
Students    : Emily Petroni and Abel Sapirstein, Assignment 2
Class       : HMC MCB 118b
Date        : 2018 Mar 21
Description : Alignment Library
"""

import sys, fasta
import numpy as np
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

def align(v, w, match=1, mismatch=-2, indel=-3 ,aligntype="global"):
    """Compute the optimal pairwise alignment of two genomic sequences.
    
    Parameters
    --------------------
        v, w                 -- strings to align
        match                -- score for match (default=1)
        mismatch             -- score for mismatch (default=-2)
        indel                -- score for a gap/indel (default=-3)
        aligntype            -- aligment type, can be local or global
    
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
    # initalize DP tble 
    nrow=len(v)+1
    ncol=len(w)+1
    score_array = [[0 for n in xrange(ncol)] for n in xrange(nrow)]
    parent_array = [["" for n in xrange(ncol)] for n in xrange(nrow)]
    if aligntype=="global":

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
                    if score == diag:
                        p = "d"
                    elif score == up:
                        p = "u"
                    elif score == left:
                        p = "l"
                    ### enter conditionals
                else:
                    score= indel*(row+col)
                score_array[row][col]= score
                parent_array[row][col] = p
    if aligntype=="local":
        # complete DP table
        for row in xrange(nrow):
            for col in xrange(ncol):
                score = 0
                p = ""
                if col > 0 and row > 0:
                    diag = score_array[row-1][col-1]
                    up = score_array[row][col-1]+indel
                    left = score_array[row-1][col]+indel
                    if v[row-1]==w[col-1]:
                        diag = diag+match
                    else:
                        diag= diag+mismatch
                    score = max(diag,up,left,0)
                    if score == diag:
                        p = "d"
                    elif score == left:
                        p = "l"
                    elif score == up:
                        p = "u"
                    
                    ### enter conditionals
                score_array[row][col]= score
                parent_array[row][col] = p
    print(np.array(parent_array))
    if aligntype == 'global':
        a,b,c = globaltraceback(v,w,match,mismatch,indel,score_array)
        return score_array,a,b,c
    if aligntype == 'local':
        a,b,c=localtraceback(v,w,match,mismatch,indel,score_array,parent_array)
        return score_array,a,b,c
        
    
    
            
        
        ### ========== TODO : END ========== ###

def localtraceback(v,w,match, mismatch, indel, scoreArray,parent_array):
    ""
    maxrow = 0
    maxcol = 0 
    maxval = 0 
    for row in range(len(scoreArray)):
        for col in range(len(scoreArray[row])):
            if scoreArray[row][col] > maxval:
                maxval = scoreArray[row][col]
                maxrow = row
                maxcol = col
    val = maxval
    current_row = maxrow
    current_col = maxcol
    aligned_w = ""
    aligned_v = ""
    v_pos = current_row
    w_pos = current_col
    while val > 0 and current_col>0 and current_row>0:
        print(val)
        p = parent_array[current_row][current_col]
        if p == 'd':
            print("Diag")
            aligned_v = v[v_pos-1]+aligned_v
            aligned_w = w[w_pos-1]+aligned_w
            v_pos = v_pos-1
            w_pos = w_pos-1
            current_row = current_row-1
            current_col = current_col-1
        elif p == 'u':
            print("Up")
            aligned_w = "-"+aligned_w
            aligned_v = v[v_pos-1]+aligned_v
            v_pos-=1
            current_row = current_row-1
            

        elif p == 'l':
            print("left")
            aligned_v = "-"+aligned_v
            aligned_w = w[w_pos-1]+aligned_w
            w_pos-=1
            current_col = current_col-1
        elif val ==0: 
            print("Breaking")
            break 
        val = scoreArray[current_row][current_col]
    return aligned_v, aligned_w, maxval

def globaltraceback(v,x,match,mismatch, indel, scoreArray):
    aligned_v = ""
    aligned_w = ""
    last_row = nrow-1
    last_col = ncol-1
    v_pos = nrow-1
    w_pos = ncol-1
    maxval = scoreArray[last_row][last_col]
    while last_row>0 or last_col>0:
        current=score_array[last_row][last_col]
        diag = score_array[last_row-1][last_col-1]
        if v[last_row-1]==w[last_col-1]:
            diag = val+match
        else:
            diag= val+mismatch
        up = scoreArray[last_row-1][last_col]+indel
        left = scoreArray[last_row][last_col-1]+indel
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
    return  aligned_v, aligned_w, maxval

#========================================
# main

def run(match, mismatch, indel, vFileName, wFileName):
    # load each sequence in fasta files
    v = fasta.load(vFileName)[0][1]
    w = fasta.load(wFileName)[0][1]
    
    # compute DP table
    scoreArray = global_align(v, w, match, mismatch, indel)
    
    
    
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
        aligntype = sys.argv[1]
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
        indel = int(sys.argv[3])
        vFileName = sys.argv[4]
        wFileName = sys.argv[5]
        # part b.3: support matrix scoring
        
        score_array, aligned_v, aligned_w, max_score = \
            align(v, w, match, mismatch, indel, aligntype)
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
        
        print(np.array(score_array))
    
    else:
        raise Exception("%s: error: invalid command" % sys.argv[0])

if __name__ == '__main__':
    sys.exit(main())
