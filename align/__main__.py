# this file is only called when the package is called from the command
# line
from .algs import read_blosum, read_prot#, scoring_mat, traceback

# Read in substitution matrix
filepath = "/Users/student/Documents/BMI206/bmi203-3/HW3_due_02_23/BLOSUM50"
blosum = read_blosum(filepath)

# Read in protein sequence
filepath = "/Users/student/Documents/BMI206/bmi203-3/HW3_due_02_23/sequences/prot-0004.fa"
seq_m = read_prot(filepath)
print(seq_m)

filepath = "/Users/student/Documents/BMI206/bmi203-3/HW3_due_02_23/sequences/prot-0008.fa"
seq_n = read_prot(filepath)
print(seq_n)

# Initialize empty matrices
#score_mat = [[None]*len(seq_n) for _ in range(0,len(seq_m))]
#state_mat = [[None]*len(seq_n) for _ in range(0,len(seq_m))]
#print(len(score_mat))
#print(len(score_mat[1]))
