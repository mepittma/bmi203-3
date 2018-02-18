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


def scoring_mat(blosum, seq_n, seq_m, score_mat, state_mat):
    """
    This function fills in a scoring matrix for two sequences to align, as well
    as a matrix that stores the "up", "left", or "right" state of each
    """
