from .utils import score_all_prots
from copy import copy
import numpy as np
import pandas as pd
import pickle as pkl
import random

def move_matrix(coef,mat_1,mat_2,add_mat):
    """
    Do matrix calculations without messing with numpy.

    Inputs:
        coef: coefficient for the current matrix combination
        centroid: the centroid matrix
        point: the point by which to alter the centroid

    Output: a single combined matrix, representing the new centroid
    """
    # Subtract the point from the centroid to get the new point
    new_point = coef * [map(lambda x, y: x-y, ii, jj) for ii, jj in zip(mat_1, mat_2)]
    # Add the old centroid to the new point and return
    return [map(lambda x,y: x+y, ii, jj) for ii, jj in zip(add_mat, new_point)]


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
        filepath = base_dir + "/HW3_due_02_23/worst_pospairs.txt"
    pos_scores = score_all_prots(subst, filepath, gapO, gapE)

    print("\tScoring negatives...")
    filepath = base_dir + "/HW3_due_02_23/Negpairs.txt"
    if mini == True:
        filepath = base_dir + "/HW3_due_02_23/best_negpairs.txt"
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

def nelder_mead(start_mat, step, improv_thr=0.0001, improv_iter=5, max_iter=20, refl=1, expn=2, contr=0.5, shrnk=0.5):
    """
    (adapted from https://github.com/fchollet/nelder-mead/blob/master/nelder_mead.py)

    Inputs:
        start_mat = substitution matrix to optimize for true positives - pandas dataframe
        step = distance to move to try new matrix
        improv_thresh, improv_iter = after improv_iter number of iterations without
            improvement greater than improv_thresh, return results
        max_iter = always break after this number of iterations
        refl, expn, contr, shrnk = reflection, expansion, contraction and shrink coefficients
            #currently defaults specified by Wikipedia, subject to change

    Output: tuple(best substitution matrix, best score)

    """

    base_dir = "/Users/student/Documents/BMI206/bmi203-3"

    # Initialize
    #start_mat is pandas
    #res = [[pandas, score],[pandas,score]]
    print("Calculating the objective score of the original BLOSUM matrix...")
    prev_best = obj_score(start_mat, mini=True)
    res = [[start_mat, prev_best]]    # list of lists to keep track of matrix/score pairs
    no_improv = 0                     # set a count to keep track of how many times this exact matrix was seen

    # Interesting function of symmetry: we try addition of one step and addition of two steps at the same time
    print("Calculating first permutation of the BLOSUM matrix...")

    """
    Ideal, but takes way too long to run
    for i,row in start_mat.iterrows():
        for j in start_mat.columns:
            x = copy(start_mat)                     # make a shallow copy of the last matrix to be seen
            print("Looking at row ", i, ", col ",j)
            print("Value at this index: ", x.loc[i,j])
            x.loc[i,j] = x.loc[i,j] + step          # add the chosen step to this parameter
            x.loc[i,j] = x.loc[i,j]                 # keep symmetric!
            score = obj_score(x)                    # score this new matrix
            res.append([x, score])                  # append this matrix and score to the running list
    """

    names = list(start_mat.columns.values)
    del names[-1]
    print("Names: ", names)
    for i in range(0,2):#range(0,200):
        x = copy(start_mat)                     # make a shallow copy of the last matrix to be seen
        i,j = random.sample(names, 2)           # choose a random cell to alter step size (excluding * cells)
        print("Looking at row ", i, ", col ",j)
        print("Value at this index: ", x.loc[i,j])
        x.loc[i,j] = x.loc[i,j] + step          # add the chosen step to this parameter
        x.loc[i,j] = x.loc[i,j]                 # keep symmetric!
        score = obj_score(x, mini=True)                    # score this new matrix
        res.append([x, score])                  # append this matrix and score to the running list



    # Save out as a pkl incase something goes wrong downstream
    iters = 0
    pkl.dump( res, open( (base_dir + "/output/permutations_{}.p".format(iters)), "wb" ) )

    # Run through this loop until a break condition is met
    while 1:

        # Order by the second value (objective score) in each list in the res running list
        res.sort(key=lambda x: x[1],reverse=True)
        best = res[0][1]

        # If we've hit the maximum iteration threshold, stop here and return the best matrix/value
        if iters >= max_iter:
            print("Max iteration of ",max_iter,"reached. Returning.")
            return res[0]

        # If we've hit the maximum possible optimization, stop here and return the best matrix/value
        if best == 4.0:
            print("Reached maximum value. Returning.")
            return res[0]


        # If we haven't hit the max, we're going through another iteration
        iters += 1
        print("Best score at iteration ", iters, ": ", best)

        # If the best is better than the previous by designated margin,
        # reset improvement stagnation at 0 and update "prev_best"
        if best > prev_best + improv_thr:
            print("Score improved.")
            improv_iter = 0
            prev_best = best

        # Otherwise, add to tally of iterations that yielded no improvement
        else:
            print("Score not improved.")
            improv_iter += 1

        # If the objective score hasn't improved for the past 10 trials, break
        if no_improv >= improv_iter:
            print("no_improv: ", no_improv)
            print("improv_iter: ", improv_iter)
            print("No improvement for", improv_iter, "rounds. Returning.")
            return res[0]

        # Calculate centroid
        # Create a centroid matrix of zeroes (a list of lists of zeros)
        centroid = pd.DataFrame(0, index=start_mat.columns.values, columns=start_mat.columns.values)
        print("Centroid matrix: ", centroid)

        # Iteratively calculate the centroid using every matrix but the last
        # DOUBLE THE NECESSARY TIME/SPACE, should probably fix
        for tup in res[:-1]:

            mat = tup[0]
            for i,row in start_mat.iterrows():
                for j in start_mat.columns:
                    centroid.loc[i,j] += mat.loc[i,j]/(len(res)-1)

        # Attempt reflection - does reflecting the worst point through the centroid improve the objective score?
        #xr = move_matrix(refl, centroid, res[-1][0], centroid)
        xr = centroid + refl*(centroid - res[-1][0])
        #print("Here's what the centroid matrix looks like: ", centroid)
        #print("Here's what the worst matrix looks like: ", res[-1][0])
        #print("Here's what the reflected matrix looks like: ", xr)

        # If the score of the reflected matrix is better, keep it
        rscore = obj_score(xr, mini=True)
        if res[0][1] >= rscore > res[-2][1]:
            # if the score of the best matrix is greater than/equal to the new rscore,
            # which is greater than the second-to-last worst score, replace the worst-scoring matrix
            del res[-1]
            res.append([xr, rscore])
            continue

        # If the reflection point scores the best so far, try moving the list of solutions in that direction
        if rscore > res[0][1]:
            #xe = move_matrix(expn, res[-1][0], centroid, centroid)
            xe = centroid + expn*(res[-1][0] - centroid)
            #print("Here's what the centroid matrix looks like: ", centroid)
            #print("Here's what the worst matrix looks like: ", res[-1][0])
            #print("Here's what the expanded matrix looks like: ", xe)
            escore = obj_score(xe,mini=True)

            # Keep whichever score was better
            if escore > rscore:
                del res[-1]
                res.append([xe, escore])
                continue
            else:
                del res[-1]
                res.append([xr, rscore])
                continue

        # Contraction
        #xc = move_matrix(contr, res[-1][0], centroid, centroid)
        xc = centroid + contr*(res[-1][0] - centroid)
        #print("Here's what the centroid matrix looks like: ", centroid)
        #print("Here's what the worst matrix looks like: ", res[-1][0])
        #print("Here's what the contracted matrix looks like: ", xc)
        cscore = obj_score(xc,mini=True)

        # If the contraction score is better than the worst score, replace
        if cscore > res[-1][1]:
            del res[-1]
            res.append([xc, cscore])
            continue

        # Reduction - replace all points except the best with x_i - x_1
        x1 = res[0][0] # This is the top scoring matrix
        nres = []
        for pair in res:
            redx = x1 + shrnk*(pair[0] - x1)
            #print("Here's what the centroid matrix looks like: ", centroid)
            #print("Here's what the best matrix looks like: ", res[0][0])
            #print("Here's what the reduced matrix looks like: ", redx)

            score = obj_score(redx,mini=True)
            nres.append([redx, score])
        res = nres
