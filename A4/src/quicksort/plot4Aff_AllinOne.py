import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

numThread = np.array([1,2,4,7,8,14,16,28,32,56])

#--------------------5000000---------------------------------
#--------------------balanced-----------------------------
with open("../../result/quicksort/StrongScale_5000000.txt","r") as out:
	lines = out.readlines()

exeTime_5000000 = []
for line in lines:
	exeTime_5000000.append(float(line.split()[2]))

out.close()



#--------------------10000000---------------------------------
#--------------------balanced-----------------------------
with open("../../result/quicksort/StrongScale_10000000.txt","r") as out:
	lines = out.readlines()

exeTime_10000000 = []
for line in lines:
	exeTime_10000000.append(float(line.split()[2]))

out.close()


#--------------------20000000-----------------------------
#--------------------balanced-----------------------------
with open("../../result/quicksort/StrongScale_20000000.txt","r") as out:
	lines = out.readlines()

exeTime_20000000 = []
for line in lines:
	exeTime_20000000.append(float(line.split()[2]))

out.close()


#--------------------40000000-----------------------------
#--------------------balanced-----------------------------
with open("../../result/quicksort/StrongScale_40000000.txt","r") as out:
	lines = out.readlines()

exeTime_40000000 = []
for line in lines:
	exeTime_40000000.append(float(line.split()[2]))

out.close()


#-----------------serial time -----------------------------

with open("../../result/quicksort/StrongScale_Serial.txt","r") as out:
	lines = out.readlines()

serial_time = []
for line in lines:
	serial_time.append(float(line.split()[1]))
out.close()


SpeedUp_5000000 = np.divide(np.array(serial_time[0]),np.array(exeTime_5000000))
SpeedUp_10000000 = np.divide(np.array(serial_time[1]),np.array(exeTime_10000000))
SpeedUp_20000000 = np.divide(np.array(serial_time[2]),np.array(exeTime_20000000))
SpeedUp_40000000 = np.divide(np.array(serial_time[3]),np.array(exeTime_40000000))


fig, axe = plt.subplots() # two axes on figure
fig.suptitle('Strong Scaling with Different DataSize')


axe.plot(numThread, SpeedUp_5000000, color='Blue', marker='.',label="5M")
axe.plot(numThread, SpeedUp_10000000, color='Red', marker='.',label="10M")
axe.plot(numThread, SpeedUp_20000000, color='Green', marker='.', label="20M")
axe.plot(numThread, SpeedUp_40000000, color='Violet', marker='.', label="40M")

axe.set_ylabel('SpeedUp',fontsize=8)
axe.set_xlabel('NumThread',fontsize=8)
axe.legend()

# #plt.show()
# axes[0,1].plot(numThread, speedup_balanced_4000, color='Blue', marker='.',label='balanced')
# axes[0,1].plot(numThread, speedup_compact_4000, color='Red', marker='.',label='compact')
# axes[0,1].plot(numThread, speedup_scatter_4000, color='Green', marker='.', label='scatter')
# axes[0,1].get_xaxis().set_visible(False)
# axes[0,1].legend(loc='upper right',prop={'size': 4})


# axes[1,0].plot(numThread, speedup_balanced_8000, color='Blue', marker='.',)
# axes[1,0].plot(numThread, speedup_compact_8000, color='Red', marker='.',)
# axes[1,0].plot(numThread, speedup_scatter_8000, color='Green', marker='.', )
# axes[1,0].set_ylabel('SpeedUp',fontsize=8)

# axes[1,1].plot(numThread, speedup_balanced_16000, color='Blue', marker='.',)
# axes[1,1].plot(numThread, speedup_compact_16000, color='Red', marker='.',)
# axes[1,1].plot(numThread, speedup_scatter_16000, color='Green', marker='.', )

# axes[2,0].plot(numThread, speedup_balanced_28000, color='Blue', marker='.',)
# axes[2,0].plot(numThread, speedup_compact_28000, color='Red', marker='.',)
# axes[2,0].plot(numThread, speedup_scatter_28000, color='Green', marker='.', )
# axes[2,0].set_ylabel('SpeedUp',fontsize=8)
# axes[2,0].set_xlabel('NumThread',fontsize=8)

# axes[2,1].plot(numThread, speedup_balanced_56000, color='Blue', marker='.',)
# axes[2,1].plot(numThread, speedup_compact_56000, color='Red', marker='.',)
# axes[2,1].plot(numThread, speedup_scatter_56000, color='Green', marker='.', )
# axes[2,1].set_xlabel('NumThread',fontsize=8)

# #plt.show();
fig.savefig("../../figs/quicksort_strongscaling.png",dpi=300,box_inches="tight")