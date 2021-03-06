from align import algs
from align import utils

# Testing the alignment of the following two sequences:
seq_m = "ARND"
seq_n = "ARNCD"

# Read in substitution matrix
filepath = "/Users/student/Documents/BMI206/bmi203-3/HW3_due_02_23/BLOSUM50"
blosum = utils.read_blosum(filepath)
gapO = -5
gapE = -2

def test_scoring():
    # ASSUMES BLOSUM50 matrix

    score_mat, state_mat = algs.score_seqs(blosum, seq_m, seq_n, gapO, gapE)
    assert score_mat == [[0, 0, 0, 0, 0, 0], [0, 5, 0, 0, 0, 0], [0, 0, 12, 7, 5, 3],
    [0, 0, 7, 19, 14, 12], [0, 0, 5, 14, 15, 22]]
    assert state_mat == [[None, None, None, None, None, None],
    [None, 'align', '', '', '', ''], [None, '', 'align', 'del', 'del', 'del'],
    [None, '', 'ins', 'align', 'del', 'del'],
    [None, '', 'ins', 'ins', 'align', 'align']]

def test_max_finding():
    # Make sure that the function returns 22 as the max score, and returns (4,5) as its index.
    max_seen, max_list = utils.max_score(algs.score_seqs(blosum, seq_m, seq_n, gapO, gapE)[0])
    assert max_seen == 22
    assert max_list == [(4,5)]

def test_traceback():
    # Make sure that the function returns the expected alignment:
    #['A', 'R', 'N', 'N', 'D']
    #['A', 'R', 'N', '-', 'D']
    score_mat, state_mat = algs.score_seqs(blosum, seq_m, seq_n, gapO, gapE)
    max_seen, max_list = utils.max_score(score_mat)
    assert algs.traceback(
        score_mat, state_mat, max_seen, max_list, seq_m, seq_n) == (['A', 'R', 'N', 'N', 'D'], ['A', 'R', 'N', '-', 'D'])

def test_thresholding():
    # Make sure my thresholding method is doing what it should
    # Find the 70% cutoff for acceptance of true positives
    TP_scores = [0,1,2,3,4,5,6,7,8,9]
    n_to_keep = int(round(0.7 * len(TP_scores)))
    TP_scores.sort(reverse = True)
    keep = TP_scores[0:n_to_keep]

    assert keep == [9,8,7,6,5,4,3]
    assert keep[-1] == 3
