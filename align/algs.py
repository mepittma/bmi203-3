import pandas as pd

def read_blosum(filepath):
    """
    This function accepts the filepath of a substitution matrix that dictates the
    penalty for swapping two distinct amino acids and returns a pandas dataframe to
    refer to during alignment.

    Input: filepath to whitespace-delimited substitution matrix
    Output: pandas dataframe with row/column names of amino acids
    """

    # Read in the matrix
    data = pd.read_csv(filepath, comment="#", delim_whitespace=True)

    # Rename the rows to AA abbreviations
    idx = pd.Series(list(data))
    data = pd.concat([data, idx], axis=1)
    data.set_index(0, inplace=True)
    data = data.rename_axis(None)

    return data

def read_prot(filepath):
    """
    This function accepts the filepath of a protein to align, ignores the first line
    (proceeded by '>' char), strips newlines, and returns the protein as a single string.
    """
    seq = ""
    with open(filepath) as f:
        for line in f:
            if not line.startswith(">"):
                seq += (line.rstrip())

    return seq

def score_seqs(blosum, seq_m, seq_n):
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

    # Set the gap penalty etc.
    gapO = -5 #opening a gap
    gapE = -2 #extending a gap

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
            print("Score in row ", i, ", column ",j,": ",maxS)
            state_mat[i][j] = state

    # Return the score matrix
    return score_mat, state_mat
