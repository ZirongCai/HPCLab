import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import pandas

matplotlib.use("Agg")

thread=[]
gemm1=[]
gemm2=[]
file = open("output_2048s.txt", "r")
file = file.readlines()[1:]
for line in file:
    data=line.rstrip().split()
    thread.append(int(data[0]))
    gemm1.append(float(data[3]))


file = open("output_2048s_noft.txt", "r")
file = file.readlines()[1:]
for line in file:
    data=line.rstrip().split()
    gemm2.append(float(data[3]))

tmp1 = gemm1[0]
for i in range(len(gemm1)):
    gemm1[i] = gemm1[i]/tmp1 
    gemm2[i] = gemm2[i]/tmp1


plt.plot(thread,gemm1,marker='.', label='FT Optimisation')
plt.plot(thread,gemm2,marker='.', label='No FT')
plt.xlabel('nThreads')
plt.xscale('log', base=2)
plt.ylabel('Speedup')
plt.grid()
plt.legend()
#plt.title('Initial performance on laptop')
#plt.show();
plt.savefig("../../figs/4.2_CompareNoFT.png",dpi=300,bbox_inches="tight")
#plt.close()

