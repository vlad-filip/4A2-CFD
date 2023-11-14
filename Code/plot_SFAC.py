import matplotlib.pyplot as plt

# Recorded data for sfac = 0.25 and 53x37
sfac = [0.075, 0.1 , 0.25, 0.5, 0.75, 0.9]
# Note: 0.05 and 0.0 1 diverged with NaNs
dro = [0.476E-04,0.564E-04, 0.114E-03,  0.174E-03, 0.202E-03, 0.216E-03]

# Plot Bump Case
plt.yscale('log')
plt.xscale('log')
plt.plot(sfac, dro, '-o', label = 'Bump Case')

# Set plot
plt.xlabel('sfac')
plt.ylabel('Average residual density error')
plt.grid()
plt.legend()
plt.show()
