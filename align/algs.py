def score_seqs(blosum, seq_m, seq_n, gapO, gapE):
    """
    This function fills in a scoring matrix for two sequences to align, as well
    as a matrix that stores the previous state of that score (whether it was
    derived from up, left, or diagonal state).

    Adapted from pseudocode here: https://github.com/elizabethfong/SparkSmithWaterman/wiki/Smith-Waterman-Algorithm

    Inputs:
        blosum: the substitution matrix being used
        seq_n, seq_m: sequences to align
    Output:
        score_mat: matrix containing maximum alignment scores for each cell
        state_mat: matrix containing direction from which maximum was derived
    """

    # Initialize empty matrices
    score_mat = [[None]*len(seq_n) for _ in range(0,len(seq_m))]
    state_mat = [[None]*(len(seq_n)+1) for _ in range(0,len(seq_m)+1)]

    # Set the first row of score_mat equal to zero
    score_mat.insert(0,[0] * len(seq_n))
    # Set the first column equal to zero
    for row in score_mat:
        row.insert(0,0)

    # Start scoring
    for i in range(1,len(seq_m)+1):
        for j in range(1, len(seq_n)+1):

            # Initialize values
            maxS = 0
            state = ""

            # Calculate deletion score (west cell)
            # (Depends on whether gap is opening or extending)
            if state_mat[i][j-1] == "del":
                score = score_mat[i][j-1] + gapE
            else:
                score = score_mat[i][j-1] + gapO

            # Accept new score if it's greater than 0
            if score > maxS:
                maxS = score
                state = "del"

            # Calculate insertion score (north cell)
            if state_mat[i-1][j] == "ins":
                score = score_mat[i-1][j] + gapE
            else:
                score = score_mat[i-1][j] + gapO

            if score > maxS:
                maxS = score
                state = "ins"

            # Calculate alignment score (northwest cell)
            aa_i = seq_m[i-1]
            aa_j = seq_n[j-1]
            alignmentScore = blosum.loc[aa_i,aa_j]
            score = score_mat[i-1][j-1] + alignmentScore
            if score > maxS:
                maxS = score
                state = "align"

            # Assign new value to maximum of above
            score_mat[i][j] = maxS
            state_mat[i][j] = state

    # Return the score matrix
    return score_mat, state_mat

def traceback(score_mat, state_mat, max_seen, max_list, seq_m, seq_n):
    """
    This function accepts two m+1 by n+1 matrices. It locates the coordinates
    of the maximum alignment score (given by the score matrix) and traces back
    (using the state matrix) to return the alignment of the two sequences.

    Inputs:
        score_mat: scoring matrix
        state_mat: matrix indicating the source of each cell's maximum score
        ("align" if it's from an alignment, "ins" if it's from insertion i.e. from
        north cell, or "del" if it's from deletion i.e. from west cell)

    Output: consensus alignment of the two sequences

    """

    # Find optimal alignments for each max cell
    alignments = []
    score = max_seen

    for cell in max_list:
        residues_m = [] #keep track of residues in alignment
        residues_n = []
        i, j = cell

        while score > 0:

            # Alignment score
            if state_mat[i][j] == "align":
                residues_m.append(seq_m[i-1])
                residues_n.append(seq_n[j-1])
                i = i-1
                j = j-1

            # Insertion score (came from north cell)
            elif state_mat[i][j] == "ins":
                residues_m.append("-")
                residues_n.append(seq_n[j-1])
                i = i-1

            # Deletion score (came from west cell)
            elif state_mat[i][j] == "del":
                residues_m.append(seq_m[i-1])
                residues_n.append("-")
                j = j-1

            # Update score of focal cell
            score = score_mat[i][j]

    return list(reversed(residues_m)), list(reversed(residues_n))
