from .utils import score_all_prots
from copy import copy
import numpy as np
import pandas as pd
import pickle as pkl
import random

def obj_score(subst, mini=False):
    """
    Objective function to optimize a matrix.
    Input: substitution matrix for which to test objective score
    Output: objective score of that matrix
    """
    base_dir = "/Users/student/Documents/BMI206/bmi203-3"
    gapO = -6
    gapE = -5

    # Test matrix on positive and negative scores
    print("\tScoring positives...")
    filepath = base_dir + "/HW3_due_02_23/Pospairs.txt"
    if mini == True:
        filepath = base_dir + "/output/Pos_max.txt"
    pos_scores = score_all_prots(subst, filepath, gapO, gapE)

    print("\tScoring negatives...")
    filepath = base_dir + "/HW3_due_02_23/Negpairs.txt"
    if mini == True:
        filepath = base_dir + "/output/Neg_max.txt"
    neg_scores = score_all_prots(subst, filepath, gapO, gapE)
    neg_scores.sort(reverse=True)

    # Select FPR cut-off points to plot
    print("\tFinding cutoff...")
    score = 0.0
    for c in np.arange(0, 0.31, 0.1):

        # Find the number of negatives to keep for a FPR of i
        n_to_keep = int(round(c * len(neg_scores)))
        keep = neg_scores[0:n_to_keep]

        # If there were any to keep at all, return them
        if len(keep) > 0:
            cutoff = keep[-1]
        else:
            #If we aren't accepting any negative alignments, the cutoff is the highest neg alignment score +1
            cutoff = neg_scores[0] + 1

        # Find TPR at this cutoff, add to score
        score += (sum(i >= cutoff for i in pos_scores)/len(pos_scores))

    return score

def random_permut(start_mat, max_iter):
    """
    """

    base_dir = "/Users/student/Documents/BMI206/bmi203-3"

    # Initialize
    #start_mat is pandas
    #res = [[pandas, score],[pandas,score]]
    print("Calculating the objective score of the original BLOSUM matrix...")
    prev_best = obj_score(start_mat)#, mini=True)
    print("Original score: ", prev_best)

    # Interesting function of symmetry: we try addition of one step and addition of two steps at the same time
    print("Calculating first permutation of the BLOSUM matrix...")

    for c in range(0,max_iter):

        names = list(start_mat.columns.values)
        del names[-1] # we don't care about the *
        step = random.sample([-5,-4,-3,-2,-1,1,2,3,4,5],1)
        x = copy(start_mat)                     # make a shallow copy of the last matrix to be seen
        i,j = random.sample(names, 2)           # choose a random cell to alter step size (excluding * cells)
        print("Looking at row ", i, ", col ",j)
        print("Value at this index: ", x.loc[i,j])
        x.loc[i,j] = x.loc[i,j] + step          # add the chosen step to this parameter
        x.loc[i,j] = x.loc[i,j]                 # keep symmetric!
        score = obj_score(x)#, mini=True)         # score this new matrix
        print("Score of new matrix: ", score)

        if score > prev_best:
            start_mat = x
            prev_best = score
            print("Accepting new matrix.")

        if score >= 4.0:
            print("Maximized function.")
            return(x, score)

        c += 1

    print("Reached max iteration of ", max_iter)
    return(x, score)
