# this file is only called when the package is called from the command
# line
from .algs import read_blosum, read_prot, score_seqs#, scoring_mat, traceback

# Read in substitution matrix
filepath = "/Users/student/Documents/BMI206/bmi203-3/HW3_due_02_23/BLOSUM50"
blosum = read_blosum(filepath)

# Read in protein sequence
#filepath = "/Users/student/Documents/BMI206/bmi203-3/HW3_due_02_23/sequences/prot-0004.fa"
#seq_m = read_prot(filepath)
seq_m = "ARND"
print(seq_m)

#filepath = "/Users/student/Documents/BMI206/bmi203-3/HW3_due_02_23/sequences/prot-0008.fa"
#seq_n = read_prot(filepath)
seq_n = "ARNCD"
print(seq_n)

score_mat, state_mat = score_seqs(blosum, seq_m, seq_n)
print(score_mat)
print(state_mat)

"""
# Initialize empty matrices
print("Length of m: ", len(seq_m))
print("Length of n: ", len(seq_n))
score_mat = [[None]*len(seq_n) for _ in range(0,len(seq_m))]
state_mat = [[None]*(len(seq_n)+1) for _ in range(0,len(seq_m)+1)]
print("Number of rows: ",len(score_mat))
print("Number of columns: ", len(score_mat[0]),"\n")

# Set the first row of score_mat equal to zero
score_mat.insert(0,[0] * len(seq_n))
# Set the first column equal to zero
for row in score_mat:
    row.insert(0,0)

# Set the gap penalty etc.
gapO = -5 #opening a gap
gapE = -2 #extending a gap

# Start scoring
for i in range(1,len(seq_m)+1):
    for j in range(1, len(seq_n)+1):

        # Initialize values
        maxS = 0
        state = ""

        # Calculate deletion score (west cell)
        # (Depends on whether gap is opening or extending)
        if state_mat[i][j-1] == "del":
            score = score_mat[i][j-1] + gapE
        else:
            score = score_mat[i][j-1] + gapO

        # Accept new score if it's greater than 0
        if score > maxS:
            maxS = score
            state = "del"

        # Calculate insertion score (north cell)
        if state_mat[i-1][j] == "ins":
            score = score_mat[i-1][j] + gapE
        else:
            score = score_mat[i-1][j] + gapO

        if score > maxS:
            maxS = score
            state = "ins"

        # Calculate alignment score (northwest cell)
        aa_i = seq_m[i-1]
        aa_j = seq_n[j-1]
        alignmentScore = blosum.loc[aa_i,aa_j]
        score = score_mat[i-1][j-1] + alignmentScore
        if score > maxS:
            maxS = score
            state = "align"

        # Assign new value to maximum of above
        score_mat[i][j] = maxS
        print("Score in row ", i, ", column ",j,": ",maxS)
        state_mat[i][j] = state

# Return the score matrix
print(score_mat)
"""
