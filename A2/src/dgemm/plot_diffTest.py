import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import pandas

matplotlib.use("Agg")

thread=[]
micro1=[]
gebp1=[]
gemm1=[]
file = open("output_2048s.txt", "r")
file = file.readlines()[1:]
for line in file:
    data=line.rstrip().split()
    thread.append(int(data[0]))
    micro1.append(float(data[1]))
    gebp1.append(float(data[2]))
    gemm1.append(float(data[3]))


plt.plot(thread,micro1,marker='.', label='Microkernel')
plt.plot(thread,gebp1,marker='.', label='GEBP')
plt.plot(thread,gemm1,marker='.', label='GEMM')

plt.xlabel('nThreads')
plt.xscale('log', base=2)
plt.ylabel('Performance (GFLOPs)')
plt.grid()
plt.legend()
#plt.title('Initial performance on laptop')
#plt.show();
plt.savefig("../../figs/4.1_differenceTests.png",dpi=300,bbox_inches="tight")
#plt.close()

