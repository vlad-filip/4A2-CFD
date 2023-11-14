import matplotlib.pyplot as plt
import numpy as np

# Recorded data for sfac = 0.25 and CFL = 0.4
size = [15, 25, 53, 75 , 100, 250]
it = [100, 500, 2000, 6000, 10000, 17500]

# Calculate time for one iterations
# Record 20,000 its and devide
t = 48.13/20000

for i in range(np.size(it)):
	it[i] = it[i] * t

# Plot Bump Case
plt.yscale('log')
plt.xscale('log')
plt.plot(size, it, '-o', label = 'Bump Case')

# Set plot
plt.xlabel('Number i-direction cells')
plt.ylabel('Run time/seconds')
plt.xlim([10,750])
plt.grid()
plt.legend()
plt.show()
