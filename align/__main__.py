# this file is only called when the package is called from the command
# line
from .algs import score_seqs, traceback
from .utils import read_blosum, read_prot, max_score, get_cutoff, score_all_prots, normal_score
from .optimize import obj_score
import matplotlib.pyplot as plt
import numpy as np

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

"""
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
"""

# # # # # # # Question 2 # # # # # # #
# Using the gap penalties you determined from question 1, which of the provided
# scoring matrices performs the best, in terms of false positive rate (at a true
# positive rate of 0.7)? What are the performance rates of each of the matrices?
# Create a Receiver Operator Curve (ROC) graph which shows the fraction of true
# positives on the Y axis and the fraction of false positives on the X axis.
# Include on the graph data for each of the provided matrices. Please take care
# to make your ROC graphs square,with both X and Y axes limited to the range [0:1].
"""
gapO = -6
gapE = -5

# For each substitution matrix
roc_items = []
mat_list = ["BLOSUM50","BLOSUM62","MATIO","PAM100","PAM250"]
for mat in mat_list:

    # Read in the substitution matrix
    subst = read_blosum(base_dir + "/HW3_due_02_23/" + mat)
    print("Read in substitution matrix ", mat,"\n")

    # Get back scores for positive pairs for this matrix
    filepath = base_dir + "/HW3_due_02_23/Pospairs.txt"
    pos_scores = score_all_prots(subst, filepath, gapO, gapE)
    print("Found positive scores: ", pos_scores,"\n")

    # Get back scores for negative pairs for this matrix
    filepath = base_dir + "/HW3_due_02_23/Negpairs.txt"
    neg_scores = score_all_prots(subst, filepath, gapO, gapE)
    print("Found negative scores: ", neg_scores,"\n")

    # Write these out to a file in case something goes horribly wrong
    outpath = base_dir + "/HW3_due_02_23/" + mat + "_maxAlignmentScores.txt"
    print("Printing out to ", outpath)
    with open(outpath,'w') as f:
        f.write("positive {}\nnegative {}".format(pos_scores,neg_scores))

    # Select TPR cut-off points to plot
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

    # Save to ROC list
    print("True Positive Rates: ",TPR,"\n")
    print("False Positive Rates: ", FPR,"\n")
    roc_items.append((FPR,TPR))

    # Save out
    with open(outpath, 'a') as f:
        f.write("\nTPR: {} \nFPR: {}".format(TPR, FPR))

    # Test ROC curve for this substitution matrix
    #plt.plot(FPR,TPR)
    #plt.show()

for i in range(0,len(mat_list)):
    plt.plot(roc_items[i][0],roc_items[i][1])
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.gca().set_aspect('equal', adjustable='box')
plt.legend(mat_list, loc='lower right')
plt.show()
"""
# # # # # # # Question 3 # # # # # # #
# How does the performance change if you normalize the Smith-Waterman scores by
# the length of the shortest sequence in a pair (i.e. divide the raw score by the
# min length)?  Show the ROC curves for your best matrix and for the same matrix
# with normalized scores. Are the false positive rates better or worse? Why do
# you think this is so?
"""
gapO = -6
gapE = -5
mat = "BLOSUM50"
roc_items = []

# Read in the substitution matrix, get positive and negative scores
subst = read_blosum(base_dir + "/HW3_due_02_23/" + mat)
filepath = base_dir + "/HW3_due_02_23/Pospairs.txt"
pos_scores = score_all_prots(subst, filepath, gapO, gapE)
filepath = base_dir + "/HW3_due_02_23/Negpairs.txt"
neg_scores = score_all_prots(subst, filepath, gapO, gapE)

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

# Save to ROC list
print("True Positive Rates: ",TPR,"\n")
print("False Positive Rates: ", FPR,"\n")
roc_items.append((FPR,TPR))

# # # # # # # NORMALIZED # # # # # # #
# Read in the substitution matrix, get positive and negative scores
subst = read_blosum(base_dir + "/HW3_due_02_23/" + mat)
filepath = base_dir + "/HW3_due_02_23/Pospairs.txt"
pos_scores = normal_score(subst, filepath, gapO, gapE)
filepath = base_dir + "/HW3_due_02_23/Negpairs.txt"
neg_scores = normal_score(subst, filepath, gapO, gapE)

# Repeat with normalized values
nTPR = []
nFPR = []
for c in np.arange(0.05, 1.01, 0.01):

    # Find the number of negatives to keep for a FPR of i
    neg_scores.sort(reverse=True)
    n_to_keep = int(round(c * len(neg_scores)))
    keep = neg_scores[0:n_to_keep]
    cutoff = keep[-1]

    # Find TPR at this cutoff
    nTPR.append(sum(i > cutoff for i in pos_scores)/len(pos_scores))
    # Find FRP at this cutoff (should be close to i, but doubles and rounding do happen)
    nFPR.append(sum(i > cutoff for i in neg_scores)/len(neg_scores))

# Save to ROC list
print("True Positive Rates, normal: ",TPR,"\n")
print("False Positive Rates, normal: ", FPR,"\n")
roc_items.append((nFPR,nTPR))

for i in range(0,2):
    plt.plot(roc_items[i][0],roc_items[i][1])
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.title("ROC curve: BLOSUM50 raw vs normalized")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.gca().set_aspect('equal', adjustable='box')
plt.legend(["Raw scores","Normalized scores"], loc='lower right')
plt.show()
"""
# # # # # # # # # # # # # #   PART I   # # # # # # # # # # # # # #

# Write out the alignments of positive and negative pairs.

# # # # # # # Question 1 # # # # # # #
# Devise an optimization algorithm to modify the values in a starting score matrix
# such as to maximize the following objective function: sum of TP rates for FP
# rates of 0.0, 0.1, 0.2, and 0.3. The maximum value for the objective function
# is 4.0 (where you are getting perfect separation of positive and negative pairs
# even at the lowest false positive rate). You should use the gap and extension
# penalties derived from Part 1. Remember, you must maintain symmetry in your matrix.
# You can make use of real-valued scores in the matrices if desired (this is
# probably a good idea).



# # # # # # # Question 2 # # # # # # #
# Beginning from the best matrix from above (that which produced the alignments),
# run your optimization algorithm to maximize the fitness of the new matrix.
# How much improvement do you see in the fitness?  Show the full ROC curves for
# the original matrix and the optimized matrix. What happens when you now realign
# the sequences using the new matrix and rescore? Show the new ROC curve following
# realignment on the same graph as above. Qualitatively discuss how your matrix
# and your alignments change following optimization.

# # # # # # # Question 3 # # # # # # #
# Beginning from the MATIO matrix, but using the same initial sequence alignments,
# re-run the optimization. Show the same ROC plots as for (2). Discuss the
# relationship between the results you see here and the results you saw for (2).
