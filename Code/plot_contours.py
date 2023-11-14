#
#   plot_contours
#                               
#   Script to plot a converged flowfield from the 4A2 solver
#
#   Change to the directory you want to execute the script within and execute 
#   with "python path_to_script/plot_contours.py casename"

# Import modules and functions
from routines import *

def main():

    # Construct full filenames to read the run data
    inname = 'input_' + sys.argv[-1] + '.txt'
    outname = 'out_final_' + sys.argv[-1] + '.bin'

    # Read the settings and the case from file
    av = read_settings(inname)
    g = read_case(outname)

    # When presenting results all values should be non-dimensionalised. Two
    # variables of interest might be:
    #    1. Static pressure coefficient, (p - p_ref) / (pstag_ref - p_ref)
    #    2. Mach number, v / (ga * rgas * t)**0.5

    # First complete the "calc_secondary" function within "routines.py" to
    # calculate static pressure and Mach number, and any others you want!
    g = calc_secondary(av,g)    

    # Use the "cut_i", "mass_av" AND "area_av" functions to calculate the
    # reference pressures at the inlet plane and therefore the static pressure
    # coefficient
    # INSERT
    c = cut_i(g,0)
    pstag, m = mass_av(c,'pstag')
    pstat, l = area_av(c,'pstat')
    g['cp'] = (g['pstat']-pstat)/(pstag - pstat)
    
    # Plot mass flux in i direction
    mass_i = np.zeros(np.size(g['pstat'][:,0]))
    tstag = np.zeros(np.size(g['pstat'][:,0]))
    
    for i in range(np.size(g['pstat'][:,0])):
    	c = cut_i(g,i)
    	pstag, m = mass_av(c,'pstag')
    	tstag[i], m = area_av(c,'tstag')
    	mass_i[i] = np.sum(m)
    
    # Normalise
    mass_i = mass_i/mass_i[0]
    tstag = tstag/tstag[0]
    	
    fig = plt.figure(figsize=[9.6,7.2]); ax = plt.axes();
    ax.scatter(range(np.size(g['pstat'][:,0])),tstag)
    ax.plot(range(np.size(g['tstat'][:,0])),tstag)
    ax.set_xlabel('Section index in the ith direction')
    ax.set_ylabel('Stagnation temperature ratio')
    ax.grid()
    
    # Specify the parameters to plot
    fieldnames = ['pstag', 'mach']; 
    colnames = ['Stagnation pressure','Mach number']
    
    # Plot residual density error, match dimensions
    g['dro'] = np.vstack([g['dro'],g['dro'][51,:]])
    g['dro'] = np.column_stack((g['dro'],g['dro'][:,35]))

    # Plot the calculated non-dimensional parameters to show the flow solution
    for n,name in enumerate(fieldnames):

        # Open figure window
        fig = plt.figure(figsize=[9.6,7.2]); ax = plt.axes();
    
        # Set aspect ratio as equal and remove axes labels
        ax.set_aspect('equal',adjustable='box'); ax.axis('off')
 
        # Plot filled contour levels
        hc = ax.pcolormesh(g['x'],g['y'],abs(g[name]),shading='gouraud')

        # Add colorbar with variable name
        colorbar(hc,colnames[n])

        # Add Mach = 1 contours
        if name == 'mach':
            ax.contour(g['x'],g['y'],g['mach'],[1.0],colors='w',
                linewidths=0.5) 

        # Draw the walls of the block
        plot_wall(ax,g)

    # Show all the plots
    plt.show()

    
main()


