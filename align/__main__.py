# this file is only called when the package is called from the command
# line
from .algs import score_seqs, traceback
from .utils import read_blosum, read_prot, max_score, get_cutoff

# # # # # # # # # # # # # #   PART I   # # # # # # # # # # # # # #

# # # # # # # Question 1 # # # # # # #
# Consider the false positive rate (proportion of negative pairs with scores that
# exceed a score threshold) when the true positive rate (proportion of positive
# pairs with scores above the threshold) is 0.7. What's the best false positive
# rate that you can achieve with varying both gap opening (from 1 to 20) and
# extension penalties (from 1 to 5) with the BLOSUM50 matrix? What is the best gap
# penalty combination?

# Constants: set filepath, read in the BLOSUM50 matrix
base_dir = "/Users/student/Documents/BMI206/bmi203-3"
blosum = read_blosum(base_dir + "/HW3_due_02_23/BLOSUM50")

# Test different gap opening/extension penalties on the positives
params = []

for gapO in range(1,21):
    for gapE in range(1,6):

        t = get_cutoff(base_dir, gapO, gapE, 0.7, blosum)

        # Find how many negatives are hits
        negs = []

        filepath = base_dir + "/HW3_due_02_23/Negpairs.txt"
        with open(filepath) as f:
            neg_count = 0
            for line in f:
                prot_1, prot_2 = line.split()

                # Read in these two sequences
                seq_m = read_prot(base_dir + "/HW3_due_02_23/" + prot_1)
                seq_n = read_prot(base_dir + "/HW3_due_02_23/" + prot_2)

                # Create the alignment score matrix
                score_mat, state_mat = score_seqs(blosum, seq_m.upper(), seq_n.upper(), -gapO, -gapE)
                max_seen, max_list = max_score(score_mat)
                negs.append(max_seen)

                # Record how many total negative comparisons there are
                neg_count += 1

        # Calculate false positive rate
        print("All scores for negative pairs: ", negs)
        FP = sum(i > t for i in negs)
        print("Number of FP for opening penalty of ", gapO, " and extension penalty of ",gapE,": ",FP)
        FPR = FP/neg_count
        print("FPR: ", FPR,"\n")

        params.append((gapO,gapE,FPR))

print("Parameter results: ", params)


# # # # # # # Question 2 # # # # # # #
# Using the gap penalties you determined from question 1, which of the provided
# scoring matrices performs the best, in terms of false positive rate (at a true
# positive rate of 0.7)? What are the performance rates of each of the matrices?
# Create a Receiver Operator Curve (ROC) graph which shows the fraction of true
# positives on the Y axis and the fraction of false positives on the X axis.
# Include on the graph data for each of the provided matrices. Please take care
# to make your ROC graphs square,with both X and Y axes limited to the range [0:1].

# # # # # # # Question 3 # # # # # # #
# How does the performance change if you normalize the Smith-Waterman scores by
# the length of the shortest sequence in a pair (i.e. divide the raw score by the
# min length)?  Show the ROC curves for your best matrix and for the same matrix
# with normalized scores. Are the false positive rates better or worse? Why do
# you think this is so?
