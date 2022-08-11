import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import pandas

matplotlib.use("Agg")

thread=[]
micro1=[]
gebp1=[]
gemm1=[]
file = open("output_28Threads_nTeams.txt", "r")
file = file.readlines()[1:]
for line in file:
    data=line.rstrip().split()
    thread.append(data[0])
    gemm1.append(float(data[3]))


plt.plot(thread,gemm1,marker='.')

plt.xlabel('Threads per Team')
#plt.xscale('log', base=2)
plt.ylabel('Performance (GFLOPs)')
plt.grid()
plt.legend()
#plt.title('Initial performance on laptop')
#plt.show();
plt.savefig("../../figs/4.ThreadsPerTeam.png",dpi=300,bbox_inches="tight")
#plt.close()

