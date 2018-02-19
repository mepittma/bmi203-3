# this file is only called when the package is called from the command
# line
from .algs import read_blosum, read_prot, score_seqs, max_score, traceback

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

# set gap penalties
gapO = -5
gapE = -2

score_mat, state_mat = score_seqs(blosum, seq_m, seq_n, gapO, gapE)
print(score_mat)
print(state_mat)

max_seen, max_list = max_score(score_mat)

print(traceback(score_mat, state_mat, max_seen, max_list, seq_m, seq_n))

"""
# Find optimal alignments for each max cell
alignments = []
score = max_seen

for cell in max_list:
    residues_m = [] #keep track of residues in alignment
    residues_n = []
    i, j = cell

    while score > 0:

        # Alignment score
        if state_mat[i][j] == "align":
            residues_m.append(seq_m[i-1])
            residues_n.append(seq_n[j-1])
            i = i-1
            j = j-1

        # Insertion score (came from north cell)
        elif state_mat[i][j] == "ins":
            residues_m.append("-")
            residues_n.append(seq_n[j-1])
            i = i-1

        # Deletion score (came from west cell)
        elif state_mat[i][j] == "del":
            residues_m.append(seq_m[i-1])
            residues_n.append("-")
            j = j-1

        # Update score of focal cell
        score = score_mat[i][j]

print(list(reversed(residues_m)))
print(list(reversed(residues_n)))
"""
