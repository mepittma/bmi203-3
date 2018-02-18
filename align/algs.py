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

def scoring_mat(blosum, seq_m, seq_n):
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
    state_mat = [[None]*len(seq_n) for _ in range(0,len(seq_m))]
    print(len(score_mat))
    print
