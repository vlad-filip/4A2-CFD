import matplotlib.pyplot as plt

# Recorded data for sfac = 0.25 and CFL = 0.4
size = [15, 25, 50, 75 , 100, 250]
dro = [0.519E-03,0.289E-03,0.125E-03,0.611E-04,0.341E-04,0.993E-05]

# Plot order lines
plt.plot([10,1000],[ 14.00E-04, 14.00E-04], 'g', label = '0th order')
plt.plot([10,1000],[ 14.00E-04, 14.00E-06], color = 'orange', label = '1st order')
plt.plot([10,1000],[ 14.00E-04, 14.00E-08], 'r', label = '2nd order')

# Plot Bump Case
plt.yscale('log')
plt.xscale('log')
plt.plot(size, dro, '-o', label = 'Bump Case')

# Set plot
plt.xlabel('Number i-direction cells')
plt.ylabel('Average residual density error')
plt.xlim([10,750])
plt.grid()
plt.legend()
plt.show()
