import numpy  as np
import matplotlib.pyplot as plt
data = np.loadtxt('potencialExpDirichlet.txt')
data = np.loadtxt('potencialExpFonte.txt')

print(data.shape)
d=data.shape

for i in range(0,d[0]):
    p = data[i,:]
    plt.plot(p)

plt.show()

