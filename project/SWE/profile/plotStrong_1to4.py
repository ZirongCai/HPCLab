import numpy as np
import matplotlib.pyplot as plt
#import pandas as pd
import sys

numThread = np.array([1,2,4])

#--------------------64---------------------------------
with open("./profile/outputs/StrongScale_64.txt","r") as out:
	lines = out.readlines()

cpuTime_64 = []
commTime_64 = []
for line in lines:
	cpuTime_64.append(float(line.split()[4]))
	commTime_64.append(float(line.split()[10]))

out.close()

#----------------------128-------------------------------
with open("./profile/outputs/StrongScale_128.txt","r") as out:
	lines = out.readlines()

cpuTime_128 = []
commTime_128 = []
for line in lines:
	cpuTime_128.append(float(line.split()[4]))
	commTime_128.append(float(line.split()[10]))

out.close()

#--------------------256---------------------------------
with open("./profile/outputs/StrongScale_256.txt","r") as out:
	lines = out.readlines()

cpuTime_256 = []
commTime_256 = []
for line in lines:
	cpuTime_256.append(float(line.split()[4]))
	commTime_256.append(float(line.split()[10]))

out.close()

#--------------------512---------------------------------
with open("./profile/outputs/StrongScale_512.txt","r") as out:
	lines = out.readlines()

cpuTime_512 = []
commTime_512 = []
for line in lines:
	cpuTime_512.append(float(line.split()[4]))
	commTime_512.append(float(line.split()[10]))

out.close()

#--------------------1024---------------------------------
with open("./profile/outputs/StrongScale_1024.txt","r") as out:
	lines = out.readlines()

cpuTime_1024 = []
commTime_1024 = []
for line in lines:
	cpuTime_1024.append(float(line.split()[4]))
	commTime_1024.append(float(line.split()[10]))

out.close()

#--------------------2048---------------------------------
with open("./profile/outputs/StrongScale_2048.txt","r") as out:
	lines = out.readlines()

cpuTime_2048 = []
commTime_2048 = []
for line in lines:
	cpuTime_2048.append(float(line.split()[4]))
	commTime_2048.append(float(line.split()[10]))

out.close()

#--------------------4096---------------------------------
with open("./profile/outputs/StrongScale_4096.txt","r") as out:
	lines = out.readlines()

cpuTime_4096 = []
commTime_4096 = []
for line in lines:
	cpuTime_4096.append(float(line.split()[4]))
	commTime_4096.append(float(line.split()[10]))

out.close()

#--------------------8192---------------------------------
with open("./profile/outputs/StrongScale_8192.txt","r") as out:
	lines = out.readlines()

cpuTime_8192 = []
commTime_8192 = []
for line in lines:
	cpuTime_8192.append(float(line.split()[4]))
	commTime_8192.append(float(line.split()[10]))

out.close()


# SpeedUp_5000000 = np.divide(np.array(serial_time[0]),np.array(exeTime_5000000))
# SpeedUp_10000000 = np.divide(np.array(serial_time[1]),np.array(exeTime_10000000))
# SpeedUp_20000000 = np.divide(np.array(serial_time[2]),np.array(exeTime_20000000))
# SpeedUp_40000000 = np.divide(np.array(serial_time[3]),np.array(exeTime_40000000))


# fig, axe = plt.subplots() # two axes on figure
# fig.suptitle('Strong Scaling with Different DataSize')


# plt.plot(numThread, cpuTime_64, color='Blue', marker='.',label="64")
# plt.plot(numThread, cpuTime_128, color='Red', marker='.',label="128")
# plt.plot(numThread, cpuTime_256, color='Green', marker='.', label="256")
plt.plot(numThread, cpuTime_512, color='Violet', marker='.', label="512")
plt.plot(numThread, cpuTime_1024, marker='.', label="1024")
plt.plot(numThread, cpuTime_2048, marker='.', label="2048")
plt.plot(numThread, cpuTime_4096, marker='.', label="4096")
plt.plot(numThread, cpuTime_8192, marker='.', label="8192")

plt.xticks(np.arange(min(numThread), max(numThread)+1, 1.0))

plt.ylabel('CpuTime',fontsize=14)
plt.xlabel('Number of Nodes',fontsize=14)
plt.grid()

plt.legend()
plt.savefig("./profile/strongscalingCPU.png",dpi=300,box_inches="tight")


plot2 = plt.figure(2)
plt.plot(numThread, commTime_512, color='Violet', marker='.', label="512")
plt.plot(numThread, commTime_1024, marker='.', label="1024")
plt.plot(numThread, commTime_2048, marker='.', label="2048")
plt.plot(numThread, commTime_4096, marker='.', label="4096")
plt.plot(numThread, commTime_8192, marker='.', label="8192")

plt.xticks(np.arange(min(numThread), max(numThread)+1, 1.0))

plt.ylabel('Cpu+Communication Time',fontsize=14)
plt.xlabel('Number of Nodes',fontsize=14)
plt.grid()

plt.legend()

# #plt.show();
plt.savefig("./profile/strongscalingCOMM.png",dpi=300,box_inches="tight")