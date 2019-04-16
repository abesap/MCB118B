"""
Author      : Yi-Chieh Wu (adapted from Eliot Bush)
Class       : HMC MCB 118b
Date        : 2018 Feb 15
Description : Fasta module
"""

import sys

#========================================
# functions

def load(filename):
    """Load fasta or multifasta.
    
    Parameters
    --------------------
      fn -- filename of BLAST hits
    
    Return
    --------------------
      fa -- list of tuples (key, seq)
    """
    
    fa = []
    
    with open(filename,"r") as f:
        header = None
        for line in f:
            # strip whitespace from right (e.g. strip newline)
            line = line.rstrip()
            
            # skip empty lines
            if line == "":
                continue
            
            if line[0] == ">": # header
                if header is not None:
                    # put together previous (key,seq)
                    seq = "".join(temp_seq)
                    fa.append((header, seq))
                
                # reset
                header = line[1:]
                temp_seq = []
            else :              # seq
                temp_seq.append(line)
        
        # put together last (key,seq)
        seq = "".join(temp_seq)
        fa.append((header, seq))
    
    return fa


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
