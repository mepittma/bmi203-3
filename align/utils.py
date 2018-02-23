import pandas as pd
from .algs import score_seqs, traceback

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


def max_score(score_mat):
    """
    Input: scoring matrix
    Output: tuple containing (maximum score value, [list of indexes for this score])
    """
    # Loop through the list of lists, saving the max value and locations
    max_seen = 0
    max_list = []
    for i in range(1,len(score_mat)):
        for j in range(1, len(score_mat[0])):
            if score_mat[i][j] > max_seen:
                max_seen = score_mat[i][j]
                max_list = [(i,j)]
            elif score_mat[i][j] == max_seen:
                max_list.append((i,j))
    return max_seen, max_list[0] #only return index of first maximum.


def get_cutoff(base_dir, gapO, gapE, threshold, blosum):
    """
    Inputs:
        filepath: path to list of positive alignment hits
        gapO: penalty for opening a gap
        gapE: penalty for extending a gap
        threshold: desired TP rate
        blosum: substitution matrix to use

    Output: the score threshold required to accept the desired number of TPs
    """

    # Read in the list of true positives, creating a list for alignment scores
    TP_scores = []
    with open(base_dir + "/HW3_due_02_23/Pospairs.txt") as f:

        for line in f:
            prot_1, prot_2 = line.split()

            # Read in these two sequences
            seq_m = read_prot(base_dir + "/HW3_due_02_23/" + prot_1)
            seq_n = read_prot(base_dir + "/HW3_due_02_23/" + prot_2)

            # Create the alignment score matrix
            score_mat, state_mat = score_seqs(blosum, seq_m.upper(), seq_n.upper(), -gapO, -gapE)
            max_seen, max_list = max_score(score_mat)
            TP_scores.append(max_seen)

    # Find the 70% cutoff for acceptance of true positives
    n_to_keep = int(round(0.7 * len(TP_scores)))
    TP_scores.sort(reverse = True)
    keep = TP_scores[0:n_to_keep]

    return keep[-1]

    #print("All max scores: ", TP_scores)
    #print("Score cutoff for opening penalty of ", gapO, " and extension penalty of ",gapE,": ",t)
    # Confirm that this is a TP rate of 0.7
    #TP = sum(i > t for i in TP_scores)
    #print(TP/pos_count)

def score_all_prots(subst, filepath, gapO, gapE):
    """
    Inputs: sustitution matrix and the filepath of the protein pair file to read,
    plust gap penalty scheme for opening/extension
    Output: list of max alignment scores for each of these pairs
    """

    base_dir = "/Users/student/Documents/BMI206/bmi203-3"

    with open(filepath) as f:

        scores = []

        for line in f:
            prot_1, prot_2 = line.split()

            # Read in these two sequences
            seq_m = read_prot(base_dir + "/HW3_due_02_23/" + prot_1)
            seq_n = read_prot(base_dir + "/HW3_due_02_23/" + prot_2)

            # Create the alignment score matrix
            score_mat, state_mat = score_seqs(subst, seq_m.upper(), seq_n.upper(), gapO, gapE)
            max_seen, max_list = max_score(score_mat)
            scores.append(max_seen)

    return scores

def normal_score(subst, filepath, gapO, gapE):
    """
    Inputs: sustitution matrix and the filepath of the protein pair file to read,
    plust gap penalty scheme for opening/extension
    Output: list of max alignment scores for each of these pairs, normalized by length of shortest sequence
    """

    base_dir = "/Users/student/Documents/BMI206/bmi203-3"

    with open(filepath) as f:

        scores = []

        for line in f:
            prot_1, prot_2 = line.split()

            # Read in these two sequences
            seq_m = read_prot(base_dir + "/HW3_due_02_23/" + prot_1)
            seq_n = read_prot(base_dir + "/HW3_due_02_23/" + prot_2)

            # Find the minimum sequence length
            norm = min(len(seq_m),len(seq_n))

            # Create the alignment score matrix
            score_mat, state_mat = score_seqs(subst, seq_m.upper(), seq_n.upper(), gapO, gapE)
            max_seen, max_list = max_score(score_mat)
            scores.append(max_seen/norm)

    return scores

def align_and_write(filepath,outpath,subst,gapO,gapE):
    """

    """
    base_dir = "/Users/student/Documents/BMI206/bmi203-3"
    with open(outpath,'w') as w:
        with open(filepath) as f:

            for line in f:
                prot_1, prot_2 = line.split()

                w.write("Alignment of {} and {}:\n".format(prot_1,prot_2))

                # Read in these two sequences
                seq_m = read_prot(base_dir + "/HW3_due_02_23/" + prot_1)
                seq_n = read_prot(base_dir + "/HW3_due_02_23/" + prot_2)

                # Create the alignment score matrix
                score_mat, state_mat = score_seqs(subst, seq_m.upper(), seq_n.upper(), gapO, gapE)
                max_seen, max_list = max_score(score_mat)
                alignment = traceback(score_mat, state_mat, max_seen, max_list, seq_m, seq_n)

                w.write("{}\n{}\n\n".format(''.join(alignment[0]),''.join(alignment[1])))

def calc_ROC(pos_scores,neg_scores):
    """
    """

    # Select FPR cut-off points to plot
    TPR = []
    FPR = []
    for c in np.arange(0.05, 1.01, 0.01):

        # Find the number of negatives to keep for a FPR of i
        neg_scores.sort(reverse=True)
        n_to_keep = int(round(c * len(neg_scores)))
        keep = neg_scores[0:n_to_keep]
        cutoff = keep[-1]

        # Find TPR at this cutoff
        TPR.append(sum(i > cutoff for i in pos_scores)/len(pos_scores))
        # Find FRP at this cutoff (should be close to i, but doubles and rounding do happen)
        FPR.append(sum(i > cutoff for i in neg_scores)/len(neg_scores))

    # return
    return TPR,FPR
