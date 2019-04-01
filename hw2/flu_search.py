#! /usr/bin/env python

"""
Author      : Yi-Chieh Wu (adapted from Eliot Bush)
Class       : HMC MCB 118b
Date        : 2017 Mar 13
Description : Flu Search
"""

import sys, random
import fasta, search
from Bio.SubsMat.MatrixInfo import pam250

#========================================
# functions

def run_search(query_filename, db_filename,
               matrix, gap_open, gap_length,
               num_reps):
    """Load sequences from file and run a search printing results to standard out.
    
    Parameters
    --------------------
      query_filename  -- query filename, string
      db_filename     -- database filename, string
      matrix          -- scoring matrix, BioPython matrix
      gap_open        -- gap open penalty, int
      gap_length      -- gap length penalty, int
      num_reps        -- number of repetitions, int
    """
    
    # load stuff
    query_key, query_seq = fasta.load(query_filename)[0]
    db = fasta.load(db_filename)
    
    # find matches
    matches = search.search_one(query_seq, db,
                                matrix, gap_open, gap_length,
                                num_reps)
    
    # print results
    
    # ...the search seq first
    print '>' + query_key
    print query_seq
    
    # ... the matches
    for key, seq in matches:
        print '>' + key
        print seq


#========================================
# main

if __name__ == "__main__":
    # fixed alignment and search parameters
    matrix = pam250
    gap_open = -10
    gap_length = -2
    num_reps = 20
    
    # user inputs
    query_filename = sys.argv[1]
    db_filename = sys.argv[2]
    
    # search
    run_search(query_filename, db_filename, 
               matrix, gap_open, gap_length, num_reps)
