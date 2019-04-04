from Bio import pairwise2
import random

def search_one(query, db, matrix, gap,gap_length, num_reps):
    """parameters
            query 	-- query sequence, string
            db 		-- database of sequences, list of (key, sequence) tuples
            matrix 	-- scoring matrix, BioPython matrix
            gap 		-- gap open penalty for affine gap, int
            gap_length 	-- gap length penalty for affine gap, int
            num_reps 	-- number of random shuffles, int
        return values
            matches 	-- significant matches, list of (key,seq) tuples"""
    print("Starting search")
    scores=[pairwise2.align.localds(query,x[1],matrix,gap,gap_length, 
            penalize_extend_when_opening = True, score_only = True
            ) for x in db]
    print("Done Querying")
    scoreslist=[]
    for n in range(num_reps):
        temp_s=list(query)
        random.shuffle(temp_s)
        new_s=''.join(temp_s)
        [scoreslist.append(pairwise2.align.localds(new_s,x[1],matrix,gap,gap_length, 
            penalize_extend_when_opening = True, score_only = True
            )) for x in db]
    scores= [1 if x> max(scoreslist) else 0 for x in scores ]
    toreturn = []
    for x in range(len(scores)):
        if scores[x]:
            toreturn.append(db[x])
    return toreturn
