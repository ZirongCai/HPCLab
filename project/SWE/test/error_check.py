#!/usr/bin/python3

import numpy as np
import sys

if (len(sys.argv) == 1):
    file_ref = "out_ref.txt"
    file_test = "out_test.txt"
elif (len(sys.argv) == 3):
    file_ref = str(sys.argv[1])
    file_test = str(sys.argv[2])

fh = open(file_ref, "r")
for line in fh:
    ref=line.rstrip().split()
fh.close()
ref = np.array(ref); 
ref = ref.astype(float)

fh = open(file_test, "r")
for line in fh:
    test=line.rstrip().split()
fh.close()
test = np.array(test); 
test = test.astype(float)

error = np.sum(np.absolute(test-ref))

print("Error: %5.10f\n" % error)
