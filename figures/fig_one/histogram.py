import matplotlib.pyplot as plt
import numpy as np


plt.figure(figsize=[6, 6])
a = np.random.random((8, 15))
plt.imshow(a, cmap='hot', interpolation='nearest')
plt.axis('off')
plt.gca().set_position([0, 0, 1, 1])
plt.savefig("test.svg")

