import math
import caesar
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from matplotlib.lines import Line2D


#%% 
model = 'm50n512'
size = 50
colors = ['orange','darkmagenta','dodgerblue', 'tomato',  'seagreen' ]
fb_fols  = ['nofb','noagn','nojet','nox','7jk']  
fb_types= [ 'no fb',  'stellar','AGN winds', 'jets','x-ray']  
snaps = [50,78,104,151]

#%% computation functions

# compute mass function
def mass_func(masses_lg,Vol_Mpc, step, bins):
        # create histgram from that octant octants 
        N, bin_edges = np.histogram(masses_lg,bins=bins)

        # compute the stellar mass function based on that 
        SMF = np.log10(N/(Vol_Mpc*step))
        
        return SMF


# DEFINING JACKKNIFE SAMPLING

def jknife_resample(galaxy_pos,galaxy_masses,L,bins,step):
    
    # define octants
    octants = [[0,0,0],[L/2,0,0],[0,L/2,0],[0,0,L/2]
              ,[L/2,L/2,0],[L/2,0,L/2],[0,L/2,L/2],[L/2,L/2,L/2]]
    
    SMFs = []
    
    
    for corner in octants:
        
        new_mask = []
        for pos in galaxy_pos:
            
            #check if the galaxy lies in the octant  
            if corner[0] <= pos[0] <= corner[0]+ L/2 and corner[1] <= pos[1] <= corner[1]+ L/2 and corner[2] <= pos[2] <= corner[2]+ L/2:
                
                #exclude current octant from the catalog
                new_mask.append(True)
                
            else:
                new_mask.append(False)
            
        # apply mask to galaxy masses to only include current octant 
        jk_mass = galaxy_masses[new_mask]
        
        jk_mass_lg = np.log10(jk_mass) 
    
        Vol = ((L/1000)**3)/8
        
        
        SMF = mass_func(jk_mass_lg, Vol,step,bins)               
        SMFs.append(SMF)
    
    # transpose to give a range of bin contents for the different octants 
    SMF_range = np.array(SMFs).T
    SMF_range  = np.ma.masked_invalid(SMF_range)
    
    var  = 1/math.sqrt(7)*np.ma.std(SMF_range, axis=1)     
    
    return var

  
def master_SMF(model,fb_type,snap,size):
    infile = '/Users/luciescharre/Desktop/Uni/MPhys/data_scripts/cats/%s/%s_%03d.hdf5' %(fb_type,model,snap)
    
    sim = caesar.load(infile)
    h = sim.simulation.hubble_constant    
    z = sim.simulation.redshift
    Vol_Mpc = (size/h)**3
    if z<1:
        z_round =int(round(z))
    else:
        z_round =int(round(z)) 
    print(z)
    print(' ')
    print(' ')

    ################################################### PROCESSING AND PLOTTING
    # Collect the galaxy masses 
    galaxy_masses = [i.masses['stellar'] for i in sim.galaxies] 
    galaxy_masses=np.array(galaxy_masses)
    
    # take the logarithm of masses 
    galaxy_masses_lg = np.log10(galaxy_masses) 
    
    # create the bins between the min and max of the galaxy masses
    Nbins = 15
    bins = np.linspace(galaxy_masses_lg.min(),galaxy_masses_lg.max(),Nbins) 
    step = ((galaxy_masses_lg.max()-galaxy_masses_lg.min())/Nbins)
  
    # compute the SMF per volume
    SMF = mass_func(galaxy_masses_lg, Vol_Mpc, step, bins)

    # collect the positions of glaxies in comoving kiloparsec
    galaxy_pos = [np.array(i.pos) for i in sim.galaxies] 
    galaxy_pos = np.array(galaxy_pos)
    
    # define side lengthof box
    L=size*1000/h # kpc cm
    
    # call resampling function to find the cosmic variance
    var = jknife_resample(galaxy_pos, galaxy_masses, L,bins,step)
    
    # find the bincentres 
    bincen = np.zeros(len(bins)-1)
    for i in range(len(bins)-1):
        bincen[i] = 0.5*(bins[i]+bins[i+1])
    
    return bincen, SMF, var, z_round     
        
 
         
#%% Subplots of different fb runs at different redshifts

model = 'm50n512'
snaps = [36,42,50,62,78,104,151]
ccolors = ['r','gold','lime','g','b','darkviolet','magenta']
size = 50

fb_types = ['nofb','noagn','nojet','nox','7jk', 'legend']  
fb_titles = [ 'no fb',  '+ stellar','+ AGN winds', '+ jets','+ x-ray','legend']  

def type_divided_subplots():
    # loop calling the master_SMF function
    fig, ax = plt.subplots(nrows=3, ncols=2)
    plt.subplots_adjust(wspace=0, hspace=0)
    fig.set_size_inches(9,12)

    for j in range(len(fb_types)):
        fb_type = fb_types[j]
        fb_title = fb_titles[j]
        #fig, ax = plt.subplots(1) 
        
        if j !=5:
            for i in range(len(snaps)): 
                snap = snaps[i]
                color = ccolors[i]
                
                bincen, SMF, var, z_round = master_SMF(model,fb_type,snap,size)
                ax.flat[j].plot(bincen,SMF,label = 'z=%s' %(z_round),color=color)
                ax.flat[j].fill_between(bincen, SMF-var, SMF+var,facecolor= 'lightgrey')
                at = AnchoredText('%s' %(fb_title), prop=dict(size=12), frameon=True, loc='upper right')
                at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
                ax.flat[j].add_artist(at)
  
        else:
                colors = ['r','gold','lime','g','b','darkviolet','magenta']
                custom_lines = [Line2D([0], [0], color=colors[0], lw=2),
                    Line2D([0], [0], color=colors[1], lw=2),
                    Line2D([0], [0], color=colors[2], lw=2),
                    Line2D([0], [0], color=colors[3], lw=2),
                    Line2D([0], [0], color=colors[4], lw=2),
                    Line2D([0], [0], color=colors[5], lw=2),
                    Line2D([0], [0], color=colors[6], lw=2)]
    
                custom_labels = ['z=6','z=5','z=4','z=3','z=2','z=1','z=0']
                ax.flat[5].legend(custom_lines, custom_labels,loc = 'center')
                ax.flat[5].set_yticklabels([])
                #ax.flat[5].set_xticklabels([])
                ax.flat[5].plot(SMF,color = 'w')
                ax.flat[5].set_ylim(-5,-0.75)
                ax.flat[5].set_xlim(8.7,13)
        
        ax.flat[j].set_ylim(-5,-0.75)
        ax.flat[j].set_xlim(8.7,13)
        #ax.flat[j].annotate('%s' %(fb_title),xycoords='data', xy=(-2,12),xytext=(0.5, 0.5), textcoords='axes fraction', color ='k')
        ax.flat[j].tick_params(right=True, top=True)   
        ax.flat[j].tick_params(axis="y",direction="inout")
        ax.flat[j].tick_params(axis="x",direction="inout")

        
        if j == 0 or j==1 or j==2 or j==3:
            ax.flat[j].set_xticklabels([], direction='inout')
        
        if j == 1 or j == 3:    
            ax.flat[j].set_yticklabels([], direction='inout')
            
        
        #fig, ax = plt.subplots()

    
    ax.flat[5].tick_params(labelright=False, labeltop=False, labelbottom = True)   
    #ax.flat[5].xaxis.set_tick_coords(.5, .08)
      
    ax.flat[5].set_xlabel("log $M_{\star}$ [$M_{\odot}$]")
    ax.flat[4].set_xlabel("log $M_{\star}$ [$M_{\odot}$]")
    ax.flat[0].set_ylabel("log $\phi$ [$Mpc^{-3}$]")
    ax.flat[2].set_ylabel("log $\phi$ [$Mpc^{-3}$]")
    ax.flat[4].set_ylabel("log $\phi$ [$Mpc^{-3}$]")
    

    #ax.flat[3].xaxis.set_label_coords(.5, .08)
    
    #plt.title("Stellar mass function at z=%s, box size %s/h Mpc" %(z_round,size))
    #plt.title("%s" %(fb_type))
    #plt.savefig("SMF_%s_%3d_%s.png" %(size,snap,fb_type))

    plt.savefig("SMF_subgrid.pdf", bbox_inches='tight')
    plt.show()   
        

#type_divided_subplots()  


#%% 

def redshift_divided_subplots():

    fig, ax = plt.subplots(nrows=2, ncols=2)
    plt.subplots_adjust(wspace=0, hspace=0)
    fig.set_size_inches(11,9)
    colors = ['orange','darkmagenta','dodgerblue', 'tomato',  'seagreen' ]
    fb_fols  = ['nofb','noagn','nojet','nox','7jk']  
    fb_types= [ 'no fb',  'stellar','AGN winds', 'jets','x-ray']  
    snaps = [50,78,104,151]
        
    for i in range(len(snaps)): 
        snap = snaps[i]
        
        for j in range(len(fb_types)):
            fb_type = fb_fols[j]
            fb_title = fb_types[j]
            color = colors[j]
        
            size =50
            model = 'm50n512'
            bincen, SMF, var, z_round = master_SMF(model,fb_type,snap, size)
            plot=ax.flat[i].plot(bincen,SMF,label = '%s' %(fb_title),color=color)
            ax.flat[i].fill_between(bincen, SMF-var, SMF+var,facecolor= color, alpha = 0.2)
            at = AnchoredText('z=%s' %(z_round), prop=dict(size=12), frameon=True, loc='upper right')
            at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
            ax.flat[i].add_artist(at)
            
            """
            # plot the m100n1024 to test convergence
            if fb_type == '7jk':
                model = 'm100n1024'
                size = 100
                bincen, SMF, var, z_round = master_SMF(model,fb_type,snap,size)
                plot=ax.flat[i].plot(bincen,SMF,label = 'm100n1024',color='k')
                ax.flat[i].fill_between(bincen, SMF-var, SMF+var,facecolor= color, alpha = 0.2)
                """
        
        ax.flat[i].set_ylim(-5,-0.75)
        ax.flat[i].set_xlim(8.7,13)
        #ax.flatj].annotate('%s' %(fb_title),xycoords='data', xy=(-2,12),xytext=(0.5, 0.5), textcoords='axes fraction', color ='k')
        ax.flat[i].tick_params(right=True, top=True)   
        ax.flat[i].tick_params(axis="y",direction="inout")
        ax.flat[i].tick_params(axis="x",direction="inout")

        
       
    ax.flat[0].set_xticklabels([], direction='inout')
    ax.flat[1].set_xticklabels([], direction='inout')
    
    ax.flat[1].set_yticklabels([], direction='inout')
    ax.flat[3].set_yticklabels([], direction='inout')

      
    ax.flat[2].set_xlabel("log $M_{\star}$ [$M_{\odot}$]")
    ax.flat[3].set_xlabel("log $M_{\star}$ [$M_{\odot}$]")
    ax.flat[0].set_ylabel("log $\phi$ [$Mpc^{-3}$]")
    ax.flat[2].set_ylabel("log $\phi$ [$Mpc^{-3}$]")
    
    ax.flat[3].legend(loc='lower left')
    

    plt.savefig("SMF_subgrid_z.pdf", bbox_inches='tight')
    plt.show()   
        

#redshift_divided_subplots() 

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
