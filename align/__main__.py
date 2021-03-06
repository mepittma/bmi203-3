# this file is only called when the package is called from the command
# line
from .algs import score_seqs, traceback
from .utils import read_blosum, read_prot, max_score, get_cutoff, score_all_prots, normal_score, align_and_write, calc_ROC
from .optimize import obj_score, random_permut
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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
        FP = sum(i >= t for i in negs)
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
        TPR.append(sum(i >= cutoff for i in pos_scores)/len(pos_scores))
        # Find FRP at this cutoff (should be close to i, but doubles and rounding do happen)
        FPR.append(sum(i >= cutoff for i in neg_scores)/len(neg_scores))

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
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.grid()
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

# # # # # # # Question 1 # # # # # # #
# Devise an optimization algorithm to modify the values in a starting score matrix
# such as to maximize the following objective function: sum of TP rates for FP
# rates of 0.0, 0.1, 0.2, and 0.3. The maximum value for the objective function
# is 4.0 (where you are getting perfect separation of positive and negative pairs
# even at the lowest false positive rate). You should use the gap and extension
# penalties derived from Part 1. Remember, you must maintain symmetry in your matrix.
# You can make use of real-valued scores in the matrices if desired (this is
# probably a good idea).

# See optimize.py for optimization functions

# # # # # # # Question 2 # # # # # # #
# Beginning from the best matrix from above (that which produced the alignments),
# run your optimization algorithm to maximize the fitness of the new matrix.
# How much improvement do you see in the fitness?  Show the full ROC curves for
# the original matrix and the optimized matrix. What happens when you now realign
# the sequences using the new matrix and rescore? Show the new ROC curve following
# realignment on the same graph as above. Qualitatively discuss how your matrix
# and your alignments change following optimization.

# Write out the alignments of positive and negative pairs with original BLOSUM50 matrix.
"""
filepath = base_dir + "/HW3_due_02_23/Pospairs.txt"
outpath = base_dir + "/output/Pos_align.txt"
subst = read_blosum(base_dir + "/HW3_due_02_23/BLOSUM50")
gapO = -6
gapE = -5

align_and_write(filepath,outpath,subst,gapO,gapE)

filepath = base_dir + "/HW3_due_02_23/Negpairs.txt"
outpath = base_dir + "/output/Neg_align.txt"
align_and_write(filepath,outpath,subst,gapO,gapE)


# Find the best-scoring negative pairs and worst-scoring positive pairs, saving out to file
subst = read_blosum(base_dir + "/HW3_due_02_23/BLOSUM50")
gapO = -6
gapE = -5

def get_pair_score(subst,filepath,gapO,gapE):
        base_dir = "/Users/student/Documents/BMI206/bmi203-3"

        with open(filepath) as f:

            score_pairs = []

            for line in f:
                prot_1, prot_2 = line.split()

                # Read in these two sequences
                seq_m = read_prot(base_dir + "/HW3_due_02_23/" + prot_1)
                seq_n = read_prot(base_dir + "/HW3_due_02_23/" + prot_2)

                # Create the alignment score matrix
                score_mat, state_mat = score_seqs(subst, seq_m.upper(), seq_n.upper(), gapO, gapE)
                max_seen, max_list = max_score(score_mat)

                # Pair-score tuple
                score_pairs.append(((prot_1, prot_2), max_seen))

        return score_pairs

pos_f = base_dir + "/HW3_due_02_23/Pospairs.txt"
neg_f = base_dir + "/HW3_due_02_23/Negpairs.txt"
pos_pair_scores = get_pair_score(subst, pos_f, gapO, gapE)
pos_pair_scores.sort(key=lambda x: x[1])
neg_pair_scores = get_pair_score(subst, neg_f, gapO, gapE)
neg_pair_scores.sort(key=lambda x: x[1], reverse=True)

def save_out_max(scores, n, filename):
    scores = scores[:n]
    with open(filename, "w") as f:
        for pair in scores:
            f.write(" ".join([pair[0][0], pair[0][1]]))
            f.write("\n")

save_out_max(pos_pair_scores, 15, base_dir + "/output/Pos_max.txt")
save_out_max(neg_pair_scores, 15, base_dir + "/output/Neg_max.txt")
"""

"""
# Optimize the matrix for true positives at all FPR.
print("Reading in the blosum matrix...")
subst = read_blosum(base_dir + "/HW3_due_02_23/BLOSUM50")
opt_mat = random_permut(start_mat = subst, max_iter = 200)
# Save out as csv so we never have to do that again
opt_mat[0].to_csv(base_dir+"/output/optimized_matrix.csv", sep=' ', header=True)


# Find ROC values for the original and optimized matrix; plot
gapO = -6
gapE = -5
roc_items = []

# Filenames
subst = read_blosum(base_dir + "/HW3_due_02_23/BLOSUM50")
pos_f = base_dir + "/HW3_due_02_23/Pospairs.txt"
neg_f = base_dir + "/HW3_due_02_23/Negpairs.txt"

# Read in the original BLOSUM matrix, get positive and negative scores
pos_scores = score_all_prots(subst, pos_f, gapO, gapE)
neg_scores = score_all_prots(subst, neg_f, gapO, gapE)
roc_items.append(calc_ROC(pos_scores,neg_scores))

# Now read in the optimized BLOSUM matrix, get TPR and FPR
#opt_mat = read_blosum(base_dir + "/output/optimized_matrix")
pos_scores = score_all_prots(opt_mat[0], pos_f, gapO, gapE)
neg_scores = score_all_prots(opt_mat[0], neg_f, gapO, gapE)
roc_items.append(calc_ROC(pos_scores,neg_scores))

# Align and write out
filepath = base_dir + "/HW3_due_02_23/Pospairs.txt"
outpath = base_dir + "/output/BLOSUM_opt_Pos_align.txt"
align_and_write(filepath,outpath,opt_mat[0],gapO,gapE)

filepath = base_dir + "/HW3_due_02_23/Negpairs.txt"
outpath = base_dir + "/output/BLOSUM_opt_Neg_align.txt"
align_and_write(filepath,outpath,opt_mat[0],gapO,gapE)

for i in range(0,2):
    plt.plot(roc_items[i][1],roc_items[i][0])
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.title("ROC curve: BLOSUM50 raw vs optimized")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.gca().set_aspect('equal', adjustable='box')
plt.legend(["Raw scores","Normalized scores"], loc='lower right')
plt.savefig(base_dir+"/output/optimized_ROC.pdf")


print("Done optimizing BLOSUM50. Moving on to MATIO...")
"""
# # # # # # # Question 3 # # # # # # #
# Beginning from the MATIO matrix, but using the same initial sequence alignments,
# re-run the optimization. Show the same ROC plots as for (2). Discuss the
# relationship between the results you see here and the results you saw for (2).


# Write out the alignments using this matrix.
"""
# Optimize the matrix for true positives at all FPR.
print("Reading in the MATIO matrix...")
subst = read_blosum(base_dir + "/HW3_due_02_23/MATIO")
opt_mat = random_permut(start_mat = subst, max_iter = 200)
# Save out as csv so we never have to do that again
opt_mat[0].to_csv(base_dir+"/output/optimized_matrix.csv", sep=' ', header=True)
"""
print("reading opt_mat")
# Read in the matrix
filepath = base_dir + "/output/optimized_matrix.csv"
opt_mat = pd.read_csv(filepath, delim_whitespace=True)
print("opt mat: ",opt_mat)
# Find ROC values for the original and optimized matrix; plot
gapO = -6
gapE = -5
roc_items = []

# Filenames
subst = read_blosum(base_dir + "/HW3_due_02_23/MATIO")
pos_f = base_dir + "/HW3_due_02_23/Pospairs.txt"
neg_f = base_dir + "/HW3_due_02_23/Negpairs.txt"

# Read in the original BLOSUM matrix, get positive and negative scores
print("scoring original MATIO")
pos_scores = score_all_prots(subst, pos_f, gapO, gapE)
neg_scores = score_all_prots(subst, neg_f, gapO, gapE)
roc_items.append(calc_ROC(pos_scores,neg_scores))

# Now read in the optimized BLOSUM matrix, get TPR and FPR
print("scoring optimal MATIO")
pos_scores = score_all_prots(opt_mat, pos_f, gapO, gapE)
neg_scores = score_all_prots(opt_mat, neg_f, gapO, gapE)
roc_items.append(calc_ROC(pos_scores,neg_scores))

print("plotting")
for i in range(0,2):
    plt.plot(roc_items[i][1],roc_items[i][0])
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.title("ROC curve: BLOSUM50 raw vs optimized")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.gca().set_aspect('equal', adjustable='box')
plt.legend(["Raw scores","Optimized scores"], loc='lower right')
plt.savefig(base_dir+"/output/optimized_MATIO_ROC.pdf")


# Align and write out
filepath = base_dir + "/HW3_due_02_23/Pospairs.txt"
outpath = base_dir + "/output/MATIO_opt_Pos_align.txt"
align_and_write(filepath,outpath,opt_mat,gapO,gapE)

filepath = base_dir + "/HW3_due_02_23/Negpairs.txt"
outpath = base_dir + "/output/MATIO_opt_Neg_align.txt"
align_and_write(filepath,outpath,opt_mat,gapO,gapE)
