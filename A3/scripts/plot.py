import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys


with open("../result_logs/result_loginnode.txt","r") as out:
	line = out.readlines()

out = "";
i=0;
lineIter = iter(line)
for line in lineIter:
	if(line[0:3]=="  5"):
		out += line + next(lineIter) + next(lineIter) + next(lineIter) + next(lineIter)


output = [x.split(';') for x in out.split('\n')][:-1] 


i = 0; res = []; flops = []
for x in output:
	res.append(int(x[0]))
	flops.append(float(x[2]))
	i+=1

resbatch = res;
flopsbatch = flops;

#----------------------------------
with open("../result_logs/result_slurm.txt","r") as out:
	line = out.readlines()

out = "";
i=0;
lineIter = iter(line)
for line in lineIter:
	if(line[0:3]=="  5"):
		out += line + next(lineIter) + next(lineIter) + next(lineIter) + next(lineIter)


output = [x.split(';') for x in out.split('\n')][:-1] #pd.read_csv(out,sep=';',header=None).values;
#print(output[0])


i = 0; reslogin = []; flopslogin = []
for x in output:
	reslogin.append(int(x[0]))
	flopslogin.append(float(x[2]))
	i+=1
#----------------------------------



with open("../result_logs/flags_description.txt","r") as compilerflags:
	f = compilerflags.read().splitlines()[2:-1]

#print(f,'\n')
flags = []
for flag in f:
	flags.append(flag.split(' ')[6:])

print(flags)


#batch vs login node


color = ['Yellow','Blue','Red','Violet']
plt.figure(dpi=120)
i=50
plt.plot(resbatch[i:i+5], flopsbatch[i:i+5], color='Blue', )
plt.text(resbatch[i], flopsbatch[i], 'Batch -O2 -xhost -fno-alias', fontsize=9)
i=110
plt.plot(resbatch[i:i+5], flopsbatch[i:i+5], color='Violet')
plt.text(resbatch[i], flopsbatch[i], 'Batch O3 -xCORE-AVX512 -fno-alias', fontsize=9)
i=50
plt.plot(reslogin[i:i+5], flopslogin[i:i+5], color='Red')
plt.text(reslogin[i], flopslogin[i], 'Login -O2 -xhost -fno-alias', fontsize=9)
i=110
plt.plot(reslogin[i:i+5], flopslogin[i:i+5], color='Orange')
plt.text(reslogin[i+5], flopslogin[i+5], 'Login -O3 -xCORE-AVX512 -fno-alias', fontsize=9)


#plt.legend(loc='upper right',prop={'size': 8})
plt.xlabel('Resolution',fontsize=15)
plt.ylabel('MFlops',fontsize=15)
plt.ylim(0,6000)
plt.title('Batch Run Performance',fontsize=15)
plt.xticks(fontsize=13, rotation=0)
plt.yticks(fontsize=13, rotation=0)
plt.show()


"""

i=0
j=0
color = ['Yellow','Blue','Red','Violet']
plt.figure(dpi=120)
for flag in flags:
	plt.plot(res[i:i+5], flops[i:i+5], color=color[j], label=flag)
	plt.text(res[i+4], flops[i+4], flag, fontsize=11)
	i+=5	
	if(i%40==0): j+=1

#plt.legend(loc='upper right',prop={'size': 8})
plt.xlabel('Resolution',fontsize=15)
plt.ylabel('MFlops',fontsize=15)
plt.ylim(0,6000)
plt.title('Login Node Performance',fontsize=15)
plt.xticks(fontsize=13, rotation=0)
plt.yticks(fontsize=13, rotation=0)
#plt.show()

"""

#one plot
"""
i=0
j=0
k=0
fig, axs = 	plt.subplots(2,2)
for flag in flags:
	axs[j,k].plot(res[i:i+5], flops[i:i+5], label=flag)
	axs[j,k].legend(loc='upper right',prop={'size': 8})
	axs[j,k].set_xlabel('Resolution')
	axs[j,k].set_ylabel('MFlops')
	axs[j,k].set_ylim(0,6000)
	#axs[j,k].grid(True)
	i+=5
	if(i%40==0):
		k+=1
	if(i%80==0):
		j+=1
		k=0
	

plt.xlabel('Resolution')
plt.ylabel('MFlops')
fig.suptitle('Compiler Flags vs Performance')
plt.show()

"""
#multiplots

"""

i=0
j=0
for flag in flags:
	plt.figure(j,figsize=(20,10))
	plt.plot(res[i:i+5], flops[i:i+5], label=flag)
	plt.legend(loc='upper right',prop={'size': 14})
	plt.xlabel('Resolution')
	plt.ylabel('MFlops')
	plt.ylim(0,6000)
	plt.title('Compiler Flags vs Performance')
	#axs[j,k].grid(True)
	i+=5
	if(i%40==0):
		j+=1
		name = 'slurm1' + str(j) + '.png'
		plt.savefig(name)
	
	

#plt.xlabel('Resolution')
#plt.ylabel('MFlops')
#fig.suptitle('Compiler Flags vs Performance')
#plt.show()

"""
