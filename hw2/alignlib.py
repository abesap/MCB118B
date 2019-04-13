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
import pandas as pd
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
def read_matrix_file(matrixfile):
    """Given a matrix file name return a dictionary of dictionaries for
        score lookup. 
    Parameters
    --------------------
        matrixfile   -- an fname to be read in to generate
                        the matrix
    
    Return
    --------------------
        scoring_matrix  -- DOD's representing scores

    Approach: Create a dictionary of all characters represented
    for each key store another dictionary of possible values, 
    keys being the characters (eg A,C,T,G) and values being the mis/match
    score for that particular potion. 
    """
    matrix=[]
    d = {}
    dod = {}
    with open(matrixfile,"r") as f: # Read in Data
        for line in f:
            row=line.split()
            matrix.append(row)      
    names = matrix[0]
    for x in range(len(names)): # Loop to create DoDs 
        row = matrix[x+1]
        d= {}
        for y in range(1,len(names)+1):
            d[names[y-1]]=row[y]
        dod[row [0]] = d

    return dod
### =========== TODO : END =========== ###


#========================================
# main alignment function

def align(v, w, match=1, mismatch=-2, indel=-3 ,aligntype="global", matrixfile = None):
    """Compute the optimal pairwise alignment of two genomic sequences.
    
    Parameters
    --------------------
        v, w                 -- strings to align
        match                -- score for match (default=1)
        mismatch             -- score for mismatch (default=-2)
        indel                -- score for a gap/indel (default=-3)
        aligntype            -- aligment type, can be local or global
        matrixfile           -- file containing score data
    
    Return
    --------------------
        score_array          -- DP table (2D array)
        aligned_v, aligned_w -- aligned strings
        max_score            -- alignment score
    """
    
    ### ========== TODO : START ========== ###
    if matrixfile: #If using a table scoring schema
        scoring_matrix = read_matrix_file(matrixfile)
    nrow=len(v)+1
    ncol=len(w)+1
    #Generate blank arrays, chars for parent, 0's for score
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
                p = ""
                #if not on the edge 
                if col > 0 and row > 0:
                    #Look at potential parent cells
                    diag = score_array[row-1][col-1]
                    up = score_array[row][col-1]+indel
                    left = score_array[row-1][col]+indel
                    #If regular scoring algorithim 
                    if matrixfile == None:
                        if v[row-1]==w[col-1]:
                            diag = diag+match
                        else:
                            diag= diag+mismatch
                    #If using matrix to score
                    else:
                        #print(scoring_matrix[v[row-1]][w[col-1]])
                        diag = diag+int(scoring_matrix[v[row-1]][w[col-1]])
                    score = max(diag,up,left)
                    #Set parent arrays based on max score
                    if score == diag:
                        p = "d"
                    elif score == up:
                        p = "u"
                    elif score == left:
                        p = "l"
                # if on the edge
                else:
                    score= indel*(row+col)
                score_array[row][col]= score
                parent_array[row][col] = p
    if aligntype=="local":
        for row in xrange(nrow):
            for col in xrange(ncol):
                score = 0
                p = ""
                if col > 0 and row > 0:
                    diag = score_array[row-1][col-1]
                    up = score_array[row][col-1]+indel
                    left = score_array[row-1][col]+indel
                    if matrixfile == None:
                        if v[row-1]==w[col-1]:
                            diag = diag+match
                        else:
                            diag= diag+mismatch
                    else:
                        diag = diag+int(scoring_matrix[v[row-1]][w[col-1]])
                    #Scoring approach now includes 0
                    score = max(diag,up,left,0)
                    if score == diag:
                        p = "d"
                    elif score == left:
                        p = "l"
                    elif score == up:
                        p = "u"
                score_array[row][col]= score
                parent_array[row][col] = p
    if aligntype == 'global':
        a,b,c = globaltraceback(v,w,match,mismatch,indel,score_array, parent_array)
        return score_array,a,b,c
    if aligntype == 'local':
        a,b,c=localtraceback(v,w,match,mismatch,indel,score_array,parent_array)
        return score_array,a,b,c
        
    
    
            
        
        ### ========== TODO : END ========== ###

def localtraceback(v,w,match, mismatch, indel, scoreArray,parent_array):
    """ Parameters
    --------------------
        v, w                 -- strings to align
        match                -- score for match (default=1)
        mismatch             -- score for mismatch (default=-2)
        indel                -- score for a gap/indel (default=-3)
        scoreArray           -- Array of scores 
        parent_array         -- Array of pointers back to parent cell 
    
    Return
    --------------------
        score_array          -- DP table (2D array)
        aligned_v, aligned_w -- aligned strings
        max_score            -- alignment score
    """ 
    maxrow = 0
    maxcol = 0 
    maxval = 0 
    for row in range(len(scoreArray)):
        for col in range(len(scoreArray[row])):
            if scoreArray[row][col] > maxval: #finds maximum score to start at
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
        p = parent_array[current_row][current_col] #starts at maximum value
        if p == 'd': #diagonal parent
            aligned_v = v[v_pos-1]+aligned_v #adds corresponding character from strings (match or mismatch)
            aligned_w = w[w_pos-1]+aligned_w
            v_pos = v_pos-1
            w_pos = w_pos-1
            current_row = current_row-1
            current_col = current_col-1
        elif p == 'l': #left parent
            aligned_w = "-"+aligned_w
            aligned_v = v[v_pos-1]+aligned_v
            v_pos-=1
            current_row = current_row-1
            

        elif p == 'u': #upper parent
            aligned_v = "-"+aligned_v
            aligned_w = w[w_pos-1]+aligned_w
            w_pos-=1
            current_col = current_col-1
        elif val ==0: #when alignment score is zero, stop
            break 
        val = scoreArray[current_row][current_col]
    return aligned_v, aligned_w, maxval

def globaltraceback(v,w,match,mismatch, indel, score_array,parent_array):
    """ Parameters
    --------------------
        v, w                 -- strings to align
        match                -- score for match (default=1)
        mismatch             -- score for mismatch (default=-2)
        indel                -- score for a gap/indel (default=-3)
        scoreArray           -- Array of scores 
        parent_array         -- Array of pointers back to parent cell 
    
    Return
    --------------------
        score_array          -- DP table (2D array)
        aligned_v, aligned_w -- aligned strings
        max_score            -- alignment score
    """ 
    aligned_v = ""
    aligned_w = ""
    nrow = len(v)
    ncol = len(w)
    current_row = nrow-1
    current_col = ncol-1
    v_pos = nrow
    w_pos = ncol
    maxval = score_array[current_row][current_col]
    while current_col>0 and current_row>0: #while in table; stop at top left corner
        p = parent_array[current_row][current_col]
	if p == 'd': #diagonal parent, adds corresponding letters from strings (match or mismatch)
            aligned_v = v[v_pos-1]+aligned_v
            aligned_w = w[w_pos-1]+aligned_w
            v_pos = v_pos-1
            w_pos = w_pos-1
            current_row = current_row-1
            current_col = current_col-1
        elif p == 'l': #left parent, adds gap
            aligned_w = "-"+aligned_w
            aligned_v = v[v_pos-1]+aligned_v
            v_pos-=1
            current_row = current_row-1
        else: #upper parent, adds gap
            aligned_v = "-"+aligned_v
            aligned_w = w[w_pos-1]+aligned_w
            w_pos-=1
            current_col = current_col-1

    
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
            matrixfile=sys.argv[2]
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
            align(v, w, match, mismatch, indel, aligntype,matrixfile)
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
