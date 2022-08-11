import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import pandas

matplotlib.use("Agg")

size=[]
niters=[]
residum=[]
time=[]
file = open("output_gridsizes.txt", "r")
file = file.readlines()[1:]
for line in file:
    data=line.rstrip().split()
    size.append(float(data[0]))
    niters.append(float(data[1]))
    residum.append(float(data[2]))
    time.append(float(data[3]))

my_xticks=[]
for i in range(len(time)):
    my_xticks.append('$\mathregular{2^{-%d}}$'%i)

#--- TIME PLOT -----------------------------------------------------
time = np.array(time)
time = time/time[0]
plt.plot(size, time, marker='.')

plt.xlabel('10 * Gridsize', fontsize=14)
plt.xscale('log', base=2)
plt.yscale('log', base=2)
plt.ylabel('Normalized Time', fontsize=14)
plt.xlim(2**(-4), 2**(-14))
plt.grid()
plt.xticks(size, my_xticks)
plt.savefig("../../figs/5.gridsizes_1.png",dpi=300,bbox_inches="tight")

# create figure and axis objects with subplots()
fig,ax = plt.subplots()
# make a plot
ax.plot(size, niters, color="red", marker="o")
plt.xscale('log', base=2)
plt.xlim(2**(-4), 2**(-14))
plt.xticks(size, my_xticks)
plt.grid()
# set x-axis label
ax.set_xlabel("10 * Gridsize",fontsize=14)
# set y-axis label
ax.set_ylabel("Number of Iterations",color="red",fontsize=14)
ax.tick_params(axis='y', colors='red')

# twin object for two different y-axis on the sample plot
ax2=ax.twinx()
# make a plot with different y-axis using second axis object
ax2.plot(size, residum,color="blue",marker="o")
ax2.set_ylabel("Residum",color="blue",fontsize=14)
ax2.tick_params(axis='y', colors='blue')
plt.xscale('log', base=2)
plt.yscale('log')
#plt.xlim(2**(-4), 2**(-14))
plt.xticks(size, my_xticks)
plt.grid()
# save the plot as a file
plt.savefig("../../figs/5.gridsizes_2.png",dpi=300,bbox_inches="tight")

