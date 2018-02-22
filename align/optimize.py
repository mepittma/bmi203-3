from .utils import score_all_prots
from copy import copy

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


def obj_score(subst):
    """
    Objective function to optimize a matrix.
    Input: substitution matrix for which to test objective score
    Output: objective score of that matrix
    """

    # Test matrix on positive and negative scores
    filepath = base_dir + "/HW3_due_02_23/Pospairs.txt"
    pos_scores = score_all_prots(subst, filepath, gapO, gapE)

    filepath = base_dir + "/HW3_due_02_23/Negpairs.txt"
    neg_scores = score_all_prots(subst, filepath, gapO, gapE)

    # Select FPR cut-off points to plot
    score = 0.0
    for c in np.arange(0, 0.31, 0.1):

        # Find the number of negatives to keep for a FPR of i
        neg_scores.sort(reverse=True)
        n_to_keep = int(round(c * len(neg_scores)))
        keep = neg_scores[0:n_to_keep]
        cutoff = keep[-1]

        # Find TPR at this cutoff, add to score
        score += (sum(i >= cutoff for i in pos_scores)/len(pos_scores))

    return score

def nelder_mead(start_mat, step, improv_thr=0.0001, improv_iter=10, max_iter=200, refl=1, expn=2, contr=0.5, shrnk=0.5):
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

    # Initialize
    #start_mat is pandas
    #res = [[pandas, score],[pandas,score]]
    prev_best = obj_score(start_mat)
    res = [[start_mat, prev_best]]    # list of lists to keep track of matrix/score pairs
    no_improv = 0                     # set a count to keep track of how many times this exact matrix was seen

    for i, row in start_mat.iterrows():
        for j, col in row.items():
            x = copy.copy(start_mat)                # make a shallow copy of the last matrix to be seen
            x.iloc[i,j] = x.iloc[i,j] + step        # add the chosen step to this parameter
            x.iloc[j,i] = x.iloc[j,i] + step        # keep symmetric!
            score = obj_score(x)                            # score this new matrix
            res.append([x, score])                  # append this matrix and score to the running list

    # Run through this loop until a break condition is met
    iters = 0
    while 1:

        # Order by the second value (objective score) in each list in the res running list
        res.sort(key=lambda x: x[1],reverse=True)
        best = res[0][1]

        # If we've hit the maximum iteration threshold, stop here and return the best matrix/value
        if iters >= max_iter:
            return res[0]

        # If we haven't hit the max, we're going through another iteration
        iters += 1
        print('...best so far:', best)

        # If the best is better than the previous by designated margin,
        # reset improvement stagnation at 0 and update "prev_best"
        if best > prev_best + improv_thr:
            improv_iter = 0
            prev_best = best

        # Otherwise, add to tally of iterations that yielded no improvement
        else:
            improv_iter += 1

        # If the objective score hasn't improved for the past 10 trials, break
        if no_improv >= no_improv_break:
            return res[0]

        # Calculate centroid
        # Create a centroid matrix of zeroes (a list of lists of zeros)
        centroid = pd.DataFrame(0, index=start_mat.columns.values, columns=start_mat.columns.values)
        print("Rows in centroid matrix: ", len(centroid))
        print("Columns in centroid matrix: ", len(centroid[0]))

        # Iteratively calculate the centroid using every matrix but the last
        # DOUBLE THE NECESSARY TIME/SPACE, should probably fix
        for tup in res[:-1]:

            mat = tup[0]
            for i, row in start_mat.iterrows():
                for j, col in row.items():
                    centroid.iloc[i,j] += mat.iloc[i,j]/(len(res)-1)

        # Attempt reflection - does reflecting the worst point through the centroid improve the objective score?
        #xr = move_matrix(refl, centroid, res[-1][0], centroid)
        xr = centroid + alpha*(centroid - res[-1][0])
        print("Here's what the centroid matrix looks like: ", centroid)
        print("Here's what the worst matrix looks like: ", res[-1][0])
        print("Here's what the reflected matrix looks like: ", xr)

        # If the score of the reflected matrix is better, keep it
        rscore = obj_score(xr)
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
            escore = obj_score(xe)

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
        cscore = obj_score(xc)

        # If the contraction score is better than the worst score, replace
        if cscore > res[-1][1]:
            del res[-1]
            res.append([xc, cscore])
            continue

        # Reduction - replace all points except the best with x_i - x_1
        x1 = res[0][0] # This is the top scoring matrix
        nres = []
        for pair in res:
            redx = move_matrix(shrnk, pair[0], x1, x1)
            redx = x1 + shrnk*(pair[0] - x1)
            score = obj_score(redx)
            nres.append([redx, score])
        res = nres
