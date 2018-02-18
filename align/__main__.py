# this file is only called when the package is called from the command
# line
from .algs import read_blosum#, scoring_mat, traceback

filepath = "/Users/student/Documents/BMI206/bmi203-3/HW3_due_02_23/BLOSUM50"
#basename = os.path.basename(filepath)
#name = os.path.splitext(basename)


df = read_blosum(filepath)
print(df)
