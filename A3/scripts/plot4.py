import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import pandas

matplotlib.use("Agg")

res=[500,1500,2500,3500,4500]
flops=[467.35,304.00,278.83,249.97,226.66]
plt.plot(res,flops,marker='.')
plt.xlabel('Resolution')
plt.ylabel('MFlops')
plt.grid()
plt.title('Initial performance on laptop')
#plt.show();
plt.savefig("../figures/initial_performance.png",dpi=300,box_inches="tight")
#plt.close()

