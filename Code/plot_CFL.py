import matplotlib.pyplot as plt

# Recorded data for sfac = 0.25 and 53x37
CFL = [0.025, 0.05, 0.075 , 0.1, 0.4, 0.5, 0.75, 0.9]
dro = [0.107E-04, 0.192E-04, 0.254E-04, 0.320E-04 , 0.114E-03, 0.143E-03, 0.215E-03, 0.258E-03]

# Plot order lines
plt.plot([1,0.01],[ 3.00E-04, 3.00E-04], 'g', label = '0th order')
plt.plot([1,0.01],[ 3.00E-04, 3.00E-06], color = 'orange', label = '1st order')
plt.plot([1,0.01],[ 3.00E-04, 3.00E-08], 'r', label = '2nd order')

# Plot Bump Case
plt.yscale('log')
plt.xscale('log')
plt.plot(CFL, dro, '-o', label = 'Bump Case')

# Set plot
plt.xlabel('CFL')
plt.ylabel('Average residual density error')
plt.grid()
plt.legend()
plt.show()
