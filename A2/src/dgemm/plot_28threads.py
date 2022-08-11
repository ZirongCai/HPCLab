import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import pandas

matplotlib.use("Agg")

size=[]
gemm1=[]
file = open("output_28Threads.txt", "r")
file = file.readlines()[1:]
for line in file:
    data=line.rstrip().split()
    size.append(data[0])
    gemm1.append(float(data[3]))


plt.plot(size,gemm1,marker='.')

plt.xlabel('Size')
#plt.xscale('log', base=2)
plt.ylabel('Performance (GFLOPs)')
plt.grid()
plt.legend()
#plt.title('Initial performance on laptop')
#plt.show();
plt.savefig("../../figs/4.28Threads.png",dpi=300,bbox_inches="tight")
#plt.close()

