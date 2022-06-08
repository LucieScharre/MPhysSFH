import numpy as np
import matplotlib.pyplot as plt
import caesar
from astropy.cosmology import WMAP9 as cosmo
import re
from matplotlib.offsetbox import AnchoredText



#%%
fb_fol = 's50'
model = "m50n512"
size = 50
snaps = [151, 142,125, 104,90, 78,71, 62, 50,42, 36,30,26,22]

# on server
fb_types = ['7jk', 'nox','nojet',  'noagn', 'nofb']
fb_fols = ['s50', 's50nox','s50nojet',  's50noagn', 's50nofb']

# my local folders 
fb_types = ['x-ray','jets','AGN winds',  'stellar', 'no fb']
fb_fols = ['7jk', 'nox','nojet',  'noagn', 'nofb']
colors = [ 'seagreen' ,'tomato', 'dodgerblue','darkmagenta','orange']


#%% compute SFRD 

def SFRD(model, fb_fol, snap, size):
    # LOADING CATALOGUES & PRINTING IMPORTANT INFO
    #infile = '/disk04/rad/sim/%s/%s/Groups/%s_%03d.hdf5' % (
    #    model, fb_fol, model, snap)
    infile = '/Users/luciescharre/Desktop/Uni/MPhys/data_scripts/cats/%s/%s_%03d.hdf5' %(fb_fol,model,snap)
    
    print(infile)
    sim = caesar.load(infile)

    h = sim.simulation.hubble_constant
    z = sim.simulation.redshift
    if z > 1.99:
        z_round = round(z)
    elif 2 > z > 0:
        z_round = round(z, 2)
    print(z)
    print(' ')
    print(' ')

    SFR = np.array([i.sfr for i in sim.halos])

    SFR_sum = np.sum(SFR)

    Vol_Mpc = (size/h)**3
    SFR_density = SFR_sum/Vol_Mpc
    age = cosmo.age(z).value
    return SFR_density, z, age




#%% plot SFRD

def SFRD_plot(fb_fols, fb_types, colors, ratio = False):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    
    for i in range(len(fb_fols)):
        fb_type = fb_types[i]
        fb_fol = fb_fols[i]
        color = colors[i]
        SFRDs = []
        redshifts = []
        ages = []
    
        for j in range(len(snaps)):
            snap = snaps[j]

            SFR_density, z, age = SFRD(model, fb_fol, snap, size)
            SFRDs.append(SFR_density)
            redshifts.append(z)
            ages.append(age)
        
        
        ########## Ratio to fiducial run
        # first save the fiducial run for the ratio plot, make sure to put 7jk first in loop
        if ratio == True and fb_fol == '7jk':
            SFRD_fid = SFRDs
            
        elif ratio == True and fb_fol != '7jk':
            # compute ratio to fiducial run for remaining plots
            SFRD_ratio = np.array(SFRDs)/np.array(SFRD_fid)
            ax1.plot(ages, SFRD_ratio, label=fb_type, color=color)
            
            
        # if ratio == false, just plot as normal SFRD
        else: 
            ax1.plot(ages, np.log10(SFRDs), label=fb_type, color=color)

    
    if ratio ==False:
    ####### add  Observations ########
    #filein = open(
    #    "/home/rad/gizmo-analysis/caesar/Observations/md14_data.txt", "r")
        fig.set_size_inches(8,6)
        filein = open(
            "/Users/luciescharre/Desktop/Uni/MPhys/data_scripts/SFRD/md14_data.txt", "r")
        
        z_obs = []
        logSFRD = []
        sighi = []
        siglo = []
        
        for line in filein.readlines():  # iterates through the lines in the file
            # check if the current line
            # starts with "#"
            if line.startswith("#"):
                continue
            else:
                tokens = re.split(' +', line)
        
                # appends data in their respective lists
                z_obs.append(float(tokens[0]))
                logSFRD.append(float(tokens[1]))
        
                sighi.append(float(tokens[2]))  # appends array to position list
                siglo.append(abs(float(tokens[3])))
        filein.close()
        obs_err = np.array(tuple(zip(siglo, sighi)))
        obs_err = [siglo, sighi]
        imf_factor = 1.8
        logSFRD = logSFRD - np.log10(imf_factor)
        
        ax1.errorbar(cosmo.age(z_obs).value, logSFRD, obs_err,
                     label='MD14', fmt='o', color='k')
    
        ax1.set_ylabel('log SFRD ($M_{\odot} y r^{-1} M p c^{-3}$)')
        ax1.set_ylim(-2.65, -0.4)  
        filename = 'SFRD_pure_%s.pdf' %fb_type

    elif ratio == True:
        ax1.set_ylabel('SFRD Ratio')
        #ax1.set_ylim(-2.65, -0.4)  
        filename = 'SFRD_ratio.png'    
    
    
    redshift_ticks = [0, 1, 2, 4, 6]
    age_ticks = cosmo.age(redshift_ticks).value
    ax1.set_xlabel('Cosmic Time [Gyr]')
    
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xticks(age_ticks)
    ax2.set_xlabel('Redshift')
    ax2.set_xticklabels(redshift_ticks)
    
    ax1.legend(loc='lower center')
    #ax1.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.savefig(filename, bbox_inches='tight', dpi=300)


#SFRD_plot(fb_fols, fb_types, colors, ratio = True)


#%% compute SFRD in mass and SFR bins, called in SFRD_binned_plots
    
start = 11
stop = 14.5
step = 0.5
mass_bins = np.arange(start,stop,step)

start = -1
stop = 4
step = 1
SFR_bins = np.arange(start,stop,step)

def SFRD_binned(model, fb_fol, snap, size,mass_bins, SFR_bins):
    # LOADING CATALOGUES & PRINTING IMPORTANT INFO
    infile = '/Users/luciescharre/Desktop/Uni/MPhys/data_scripts/cats/%s/%s_%03d.hdf5' %(fb_fol,model,snap)
    print(infile)
    sim = caesar.load(infile)

    h = sim.simulation.hubble_constant
    z = sim.simulation.redshift
    if z > 1.99:
        z_round = round(z)
    elif 2 > z > 0:
        z_round = round(z, 2)
    print(z)
    print(' ')
    print(' ')
    
    
    mH_tot = np.array([i.halo.masses['total']
                      for i in sim.galaxies if i.central == 1])
    mH_tot_lg = np.log10(mH_tot)

    SFR = np.array([i.halo.sfr for i in sim.galaxies if i.central == 1])
    
    
    SFR_nonzero = SFR[SFR>0]
    SFR_lg = np.log10(SFR_nonzero)
    #print(min(SFR),max(SFR))
    #print(min(SFR_lg),max(SFR_lg))
    
    Vol_Mpc = (size/h)**3
    

    SFRD_massbinned = []
    SFRD_SFRbinned = []
    

   
    # use masks to bin masses and apply to SFR
    
    for bin_i in mass_bins:
            masked_mH = np.ma.masked_inside(mH_tot_lg, bin_i-step, bin_i+step)
            mH_mask = masked_mH.mask
            SFR_mass_masked = SFR[mH_mask]
            SFRD_massbinned.append(np.sum(SFR_mass_masked)/Vol_Mpc)
    
    
    for bin_i in SFR_bins:
            masked_SFR_lg = np.ma.masked_inside(SFR_lg, bin_i-step, bin_i+step)
            SFR_mask = masked_SFR_lg.mask
            SFR_SFR_masked = SFR_nonzero[SFR_mask]
            #print(bin_i,np.sum(SFR_SFR_masked))
            SFRD_SFRbinned.append(np.sum(SFR_SFR_masked)/Vol_Mpc)
   
    SFRD_all = np.sum(SFR)/Vol_Mpc
    age = cosmo.age(z).value
    return SFRD_massbinned, SFRD_SFRbinned, SFRD_all, z, age



#%% plot SFRD in mass and SFR bins in individual plots

def SFRD_binned_plots(model, fb_fols, fb_types,size,mass_bins, SFR_bins):

    # make different plots for different fb types
    
    for i in range(len(fb_fols)):
        #fig = plt.figure()
        #fig.set_size_inches(8,6)
        #ax1 = fig.add_subplot(111)
        #ax2 = ax1.twiny()
        
        fb_type = fb_types[i]
        fb_fol = fb_fols[i]
        #color = colors[i]
        SFRDs = []
        redshifts = []
        ages = []
        SFRDs_massbinned =[]
        SFRDs_SFRbinned =[]
        
        for j in range(len(snaps)):
            snap = snaps[j]

            SFRD_massbinned,SFRD_SFRbinned, SFRD_all, z, age = SFRD_binned(model, fb_fol, snap, size, mass_bins, SFR_bins)
            
            # [[massbinned_SFRD_z1],..., [massbinned_SFRD_zn]]
            
            SFRDs_massbinned.append(SFRD_massbinned)
            SFRDs_SFRbinned.append(SFRD_SFRbinned)
            SFRDs.append(SFRD_all)
            
            redshifts.append(z)
            ages.append(age)
        
        # need to transpose, [[z_SFRD_m1],..., [z_SFRD_mn]], for each mass bin i now have a z evolution of SFRD
        SFRDs_massbinned_z = np.array(SFRDs_massbinned).T
        SFRDs_SFRbinned_z = np.array(SFRDs_SFRbinned).T
        
        
        ####### mass binned plots
        
        fig = plt.figure()
        fig.set_size_inches(8,5.5)
        ax1 = fig.add_subplot(111)
        ax2 = ax1.twiny()
        
        colormap = plt.cm.gist_ncar
        plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.plasma(np.linspace(0, 1, len(mass_bins)))))

        
        for i in range(len(SFRDs_massbinned_z)):
            ax2.plot(ages, np.log10(SFRDs_massbinned_z[i]), label=str(mass_bins[i]))
        
        ax2.plot(ages, np.log10(SFRDs), label='total', color ='k')
        ax1.set_ylabel('log SFRD ($M_{\odot} y r^{-1} M p c^{-3}$)')
        ax1.set_ylim(-4, -0.4)  
        filename = 'SFRD_massbinned_%s.png'
    

        redshift_ticks = [0, 1, 2, 4, 6]
        age_ticks = cosmo.age(redshift_ticks).value
        ax2.set_xlabel('Cosmic Time [Gyr]')
        
        
        ax1.set_xlim(ax2.get_xlim())
        ax1.set_xticks(age_ticks)
        ax1.set_xlabel('Redshift')
        ax1.set_xticklabels(redshift_ticks)
        
        #ax2.legend(loc='upper right',bbox_to_anchor=(1.04,1))
        plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', title = 'log($M_{h}$/$M_{\odot})$')
        plt.tight_layout()
        plt.title('%s' %(fb_type))
        plt.savefig(filename %(fb_type), bbox_inches='tight')
        
        
        
        ####### SFR binned plots
                
        fig = plt.figure()
        fig.set_size_inches(8,5.5)
        ax1 = fig.add_subplot(111)
        
        colormap = plt.cm.gist_ncar
        plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.coolwarm(np.linspace(0, 1, len(SFR_bins)))))

        ax1.plot(ages, np.log10(SFRDs), label='total', color ='k')
        for i in range(len(SFRDs_SFRbinned_z)):
            ax1.plot(ages, np.log10(SFRDs_SFRbinned_z[i]), label=SFR_bins[i])
        
        
        ax1.set_ylabel('log SFRD ($M_{\odot} y r^{-1} M p c^{-3}$)')
        ax1.set_ylim(-4, -0.4)  
        filename = 'SFRD_SFRbinned_%s.png'
    
        
        redshift_ticks = [0, 1, 2, 4, 6]
        age_ticks = cosmo.age(redshift_ticks).value
        ax1.set_xlabel('Cosmic Time [Gyr]')
        ax2 = ax1.twiny()
        ax2.set_xlim(ax1.get_xlim())
        ax2.set_xticks(age_ticks)
        ax2.set_xlabel('Redshift')
        ax2.set_xticklabels(redshift_ticks)
        
        #ax2.legend(loc='upper right',bbox_to_anchor=(1.04,1))
        #plt.title('%s' %(fb_type))
        at = AnchoredText('%s' %(fb_type), prop=dict(size=12), frameon=True, loc='upper right')
        at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
        ax1.add_artist(at)
        ax1.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', title = 'log(SFR/$M_{\odot} yr^{-1})$')
        plt.tight_layout()
        plt.savefig(filename %(fb_type), bbox_inches='tight', dpi=300)
        

#SFRD_binned_plots(model, fb_fols, fb_types,size,mass_bins, SFR_bins)
     
        
#%% plot SFRD in mass and SFR bins in subplots

def SFRD_binned_subplots(model, fb_fols, fb_types,size,mass_bins, SFR_bins):
    
    # loop calling the master_SMF function
    fig, ax = plt.subplots(nrows=3, ncols=2)
    plt.subplots_adjust(wspace=0, hspace=0)
    fig.set_size_inches(9,12)
    
    fig1, ax1 = plt.subplots(nrows=3, ncols=2)
    plt.subplots_adjust(wspace=0, hspace=0)
    fig1.set_size_inches(9,12)

    # make different plots for different fb types
    
    for i in range(len(fb_fols)):
        iylist = [0,2,4]
        ixlist = [4,5]
        if i !=5:
            fb_type = fb_types[i]
            fb_fol = fb_fols[i]
            #color = colors[i]
            SFRDs = []
            redshifts = []
            ages = []
            SFRDs_massbinned =[]
            SFRDs_SFRbinned =[]
            
            
            for j in range(len(snaps)):
                snap = snaps[j]
    
                SFRD_massbinned,SFRD_SFRbinned, SFRD_all, z, age = SFRD_binned(model, fb_fol, snap, size, mass_bins, SFR_bins)
                
                # [[massbinned_SFRD_z1],..., [massbinned_SFRD_zn]]
                
                SFRDs_massbinned.append(SFRD_massbinned)
                SFRDs_SFRbinned.append(SFRD_SFRbinned)
                SFRDs.append(SFRD_all)
                
                redshifts.append(z)
                ages.append(age)
            
            # need to transpose, [[z_SFRD_m1],..., [z_SFRD_mn]], for each mass bin i now have a z evolution of SFRD
            SFRDs_massbinned_z = np.array(SFRDs_massbinned).T
            SFRDs_SFRbinned_z = np.array(SFRDs_SFRbinned).T
            
            colormap = plt.cm.gist_ncar
            ax.flat[i].set_prop_cycle(plt.cycler('color', plt.cm.plasma(np.linspace(0, 1, len(mass_bins)))))
            
            mass_lines = []
            for ii in range(len(SFRDs_massbinned_z)):
                line = ax.flat[i].plot(ages, np.log10(SFRDs_massbinned_z[ii]), label=str(mass_bins[ii]))
                mass_lines.append(line)
                
            ax.flat[i].plot(ages, np.log10(SFRDs), label='total', color ='k')
            ax.flat[i].set_ylim(-4, -0.3)  
            
            at = AnchoredText('%s' %(fb_type), prop=dict(size=12), frameon=True, loc='upper right')
            at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
            ax.flat[i].add_artist(at)
            #
                    
            colormap = plt.cm.gist_ncar
            ax1.flat[i].set_prop_cycle(plt.cycler('color', plt.cm.coolwarm(np.linspace(0, 1, len(SFR_bins)))))
    
            SFR_lines = []
            #labels = []
            for ii in range(len(SFRDs_SFRbinned_z)):
                line = ax1.flat[i].plot(ages, np.log10(SFRDs_SFRbinned_z[ii]), label=SFR_bins[ii])
                SFR_lines.append(line)
                #labels.label()
                
            ax1.flat[i].plot(ages, np.log10(SFRDs), label='total', color ='k')
            #ax1.set_ylabel('log SFRD ($M_{\odot} y r^{-1} M p c^{-3}$)')
            ax1.flat[i].set_ylim(-4, -0.3)  
            at = AnchoredText('%s' %(fb_type), prop=dict(size=12), frameon=True, loc='upper right')
            at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
            ax1.flat[i].add_artist(at)
            
            #filename = 'SFRD_SFRbinned_%s.png'
            ax.flat[5].plot(ages, np.log10(SFRDs), label='total', color ='w')
            ax1.flat[5].plot(ages, np.log10(SFRDs), label='total', color ='w')
            
            ax.flat[i].tick_params(right=True, top=True)   
            ax.flat[i].tick_params(axis="y",direction="inout")
            ax.flat[i].tick_params(axis="x",direction="inout")
            
            ax1.flat[i].tick_params(right=True, top=True)   
            ax1.flat[i].tick_params(axis="y",direction="inout")
            ax1.flat[i].tick_params(axis="x",direction="inout")
        
    
    redshift_ticks = [0, 1, 2, 4, 6]
    age_ticks = cosmo.age(redshift_ticks).value
    
    ax_up0 = ax.flat[0].twiny()
    ax_up0.set_xlim(ax.flat[0].get_xlim())
    ax_up0.set_xticks(age_ticks)
    ax_up0.set_xlabel('Redshift')
    ax_up0.set_xticklabels(redshift_ticks)
    
    ax_up1 = ax.flat[1].twiny()
    ax_up1.set_xlim(ax.flat[1].get_xlim())
    ax_up1.set_xticks(age_ticks)
    ax_up1.set_xlabel('Redshift')
    ax_up1.set_xticklabels(redshift_ticks)
    
    ax_up0 = ax1.flat[0].twiny()
    ax_up0.set_xlim(ax.flat[0].get_xlim())
    ax_up0.set_xticks(age_ticks)
    ax_up0.set_xlabel('Redshift')
    ax_up0.set_xticklabels(redshift_ticks)
    
    ax_up1 = ax1.flat[1].twiny()
    ax_up1.set_xlim(ax.flat[1].get_xlim())
    ax_up1.set_xticks(age_ticks)
    ax_up1.set_xlabel('Redshift')
    ax_up1.set_xticklabels(redshift_ticks)
    
        
    ax.flat[0].set_xticklabels([], direction='inout')
    ax.flat[1].set_xticklabels([], direction='inout')   
    ax.flat[2].set_xticklabels([], direction='inout')
    ax.flat[3].set_xticklabels([], direction='inout')          
    ax.flat[4].set_xlabel('Cosmic Time [Gyr]')
    ax.flat[5].set_xlabel('Cosmic Time [Gyr]')
    ax.flat[5].set_xlim(ax.flat[4].get_xlim())
    
    ax.flat[0].set_ylabel('log SFRD ($M_{\odot} y r^{-1} M p c^{-3}$)')
    ax.flat[1].set_yticklabels([], direction='inout')
    ax.flat[2].set_ylabel('log SFRD ($M_{\odot} y r^{-1} M p c^{-3}$)')
    ax.flat[3].set_yticklabels([], direction='inout')
    ax.flat[4].set_ylabel('log SFRD ($M_{\odot} y r^{-1} M p c^{-3}$)')
    ax.flat[5].set_yticklabels([], direction='inout')
    ax.flat[5].set_ylim(ax.flat[4].get_ylim())
    
    ax1.flat[0].set_xticklabels([], direction='inout')
    ax1.flat[1].set_xticklabels([], direction='inout')   
    ax1.flat[2].set_xticklabels([], direction='inout')
    ax1.flat[3].set_xticklabels([], direction='inout')   
    ax1.flat[4].set_xlabel('Cosmic Time [Gyr]')
    ax1.flat[5].set_xlabel('Cosmic Time [Gyr]')
    ax1.flat[5].set_xlim(ax1.flat[4].get_xlim())
    
    ax1.flat[0].set_ylabel('log SFRD ($M_{\odot} y r^{-1} M p c^{-3}$)')
    ax1.flat[1].set_yticklabels([], direction='inout')
    ax1.flat[2].set_ylabel('log SFRD ($M_{\odot} y r^{-1} M p c^{-3}$)')
    ax1.flat[3].set_yticklabels([], direction='inout')
    ax1.flat[4].set_ylabel('log SFRD ($M_{\odot} y r^{-1} M p c^{-3}$)')
    ax1.flat[5].set_yticklabels([], direction='inout')
    ax1.flat[5].set_ylim(ax1.flat[4].get_ylim())
    
    
    filename_mass = 'SFRD_massbinned.pdf'
    filename_SFR = 'SFRD_SFRbinned.pdf'
    
    fig.legend(mass_lines,labels = mass_bins,title = 'log($M_{h}$/$M_{\odot})$',loc= [0.71,0.12])
    fig1.legend(SFR_lines,labels = SFR_bins,title = 'log(SFR/$M_{\odot} yr^{-1})$',loc= [0.69,0.14])

    fig.savefig(filename_mass, bbox_inches='tight')
    fig1.savefig(filename_SFR, bbox_inches='tight')

#SFRD_binned_subplots(model, fb_fols, fb_types,size,mass_bins, SFR_bins)

#%% only plot SFRD in SFR bins for AGN winds, jets and x-ray runs 
    
def SFRD_SFR_subplots(model,size,SFR_bins):
    
    fig1, ax1 = plt.subplots(nrows=2, ncols=2)
    plt.subplots_adjust(wspace=0, hspace=0)
    fig1.set_size_inches(11,9)

    # make different plots for different fb types
    fb_fols  = ['nojet','nox','7jk']  
    fb_types= ['AGN winds', 'jets','x-ray']  
    
    for i in range(len(fb_fols)):

        if i !=3:
            fb_type = fb_types[i]
            fb_fol = fb_fols[i]
            #color = colors[i]
            SFRDs = []
            redshifts = []
            ages = []
            SFRDs_SFRbinned =[]
            
            
            for j in range(len(snaps)):
                snap = snaps[j]
    
                SFRD_massbinned,SFRD_SFRbinned, SFRD_all, z, age = SFRD_binned(model, fb_fol, snap, size, mass_bins, SFR_bins)
                
                # [[massbinned_SFRD_z1],..., [massbinned_SFRD_zn]]
                
                SFRDs_SFRbinned.append(SFRD_SFRbinned)
                SFRDs.append(SFRD_all)
                
                redshifts.append(z)
                ages.append(age)
            
            # need to transpose, [[z_SFRD_m1],..., [z_SFRD_mn]], for each mass bin i now have a z evolution of SFRD

            SFRDs_SFRbinned_z = np.array(SFRDs_SFRbinned).T
            
            
            colormap = plt.cm.gist_ncar
            ax1.flat[i].set_prop_cycle(plt.cycler('color', plt.cm.coolwarm(np.linspace(0, 1, len(SFR_bins)))))
    
            SFR_lines = []
            #labels = []
            for ii in range(len(SFRDs_SFRbinned_z)):
                line = ax1.flat[i].plot(ages, np.log10(SFRDs_SFRbinned_z[ii]), label=SFR_bins[ii])
                SFR_lines.append(line)
                #labels.label()
                
            ax1.flat[i].plot(ages, np.log10(SFRDs), label='total', color ='k')
            #ax1.set_ylabel('log SFRD ($M_{\odot} y r^{-1} M p c^{-3}$)')
            ax1.flat[i].set_ylim(-4, -0.3)  
            at = AnchoredText('%s' %(fb_type), prop=dict(size=12), frameon=True, loc='upper right')
            at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
            ax1.flat[i].add_artist(at)
            
            #filename = 'SFRD_SFRbinned_%s.png'
            ax1.flat[3].plot(ages, np.log10(SFRDs), label='total', color ='w')
            
            ax1.flat[i].tick_params(right=True, top=True)   
            ax1.flat[i].tick_params(axis="y",direction="inout")
            ax1.flat[i].tick_params(axis="x",direction="inout")
        
    
    redshift_ticks = [0, 1, 2, 4, 6]
    age_ticks = cosmo.age(redshift_ticks).value

    ax_up0 = ax1.flat[0].twiny()
    ax_up0.set_xlim(ax1.flat[0].get_xlim())
    ax_up0.set_xticks(age_ticks)
    ax_up0.set_xlabel('Redshift')
    ax_up0.set_xticklabels(redshift_ticks)
    
    ax_up1 = ax1.flat[1].twiny()
    ax_up1.set_xlim(ax1.flat[1].get_xlim())
    ax_up1.set_xticks(age_ticks)
    ax_up1.set_xlabel('Redshift')
    ax_up1.set_xticklabels(redshift_ticks)
    
    ax1.flat[0].set_xticklabels([], direction='inout')
    ax1.flat[1].set_xticklabels([], direction='inout')   
    ax1.flat[2].set_xlabel('Cosmic Time [Gyr]')
    ax1.flat[3].set_xlabel('Cosmic Time [Gyr]')
    ax1.flat[3].set_xlim(ax1.flat[2].get_xlim())
    
    ax1.flat[0].set_ylabel('log SFRD ($M_{\odot} y r^{-1} M p c^{-3}$)')
    ax1.flat[1].set_yticklabels([], direction='inout')
    ax1.flat[2].set_ylabel('log SFRD ($M_{\odot} y r^{-1} M p c^{-3}$)')
    ax1.flat[3].set_yticklabels([], direction='inout')
    ax1.flat[3].set_ylim(ax1.flat[2].get_ylim())

    filename_SFR = 'SFRD_SFRbinned_reduced.pdf'
    
    fig1.legend(SFR_lines,labels = SFR_bins,title = 'log(SFR/$M_{\odot} yr^{-1})$',loc= [0.7,0.19])

    fig1.savefig(filename_SFR, bbox_inches='tight', dpi=300)
    
#SFRD_SFR_subplots(model,size,SFR_bins)      
    


#%% ratio of contribution from galaxies at high SFR
    
# idea, take ratio of overall nox and 7jk and take ratio of SFR bins over 10^1.5
def SFRD_SFR_ratio(model,size):
    # make different plots for different fb types
    fb_fols  = ['nox','7jk']  
    fb_types= [ 'jets','x-ray']  
    snaps = [151, 142,125, 104,90, 78,71, 62, 50,42]
    snaps =[51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,78,90, 104, 125, 142, 151]
   
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    SFRDs_high_x =[]
    SFRDs_high_jets =[]
    SFRDs_x= []
    SFRDs_jets = []
    for i in range(len(fb_fols)):
        fb_type = fb_types[i]
        fb_fol = fb_fols[i]
        redshifts = []        
        ages = []        
        for j in range(len(snaps)):
            snap = snaps[j]
            infile = '/Users/luciescharre/Desktop/Uni/MPhys/data_scripts/cats/%s/%s_%03d.hdf5' %(fb_fol,model,snap)
            print(infile)
            sim = caesar.load(infile)
        
            h = sim.simulation.hubble_constant
            z = sim.simulation.redshift
            if z > 1.99:
                z_round = round(z)
            elif 2 > z > 0:
                z_round = round(z, 2)
            print(z)
            print(' ')
            print(' ')
            
        
            SFR = np.array([i.halo.sfr for i in sim.galaxies if i.central == 1])
            
            SFR_bin = 1
            SFR_high = SFR[SFR>10**SFR_bin]
            
            Vol_Mpc = (size/h)**3
            
            
            SFR_sum = np.sum(SFR)
            SFR_sum_high = np.sum(SFR_high)
            
            SFRD_all = SFR_sum/Vol_Mpc
            SFRD_high = SFR_sum_high/Vol_Mpc

            if fb_fol == '7jk': 
                SFRDs_high_x.append(SFRD_high)
                SFRDs_x.append(SFRD_all)
                
            elif fb_fol == 'nox': 
                SFRDs_high_jets.append(SFRD_high)
                SFRDs_jets.append(SFRD_all)
                
                
            redshifts.append(z)
            age = cosmo.age(z).value
            ages.append(age)
            

    
    ax1.plot(ages, np.array(SFRDs_x)/np.array(SFRDs_jets), label='all galaxies', color='darkturquoise')
    ax1.plot(ages,  np.array(SFRDs_high_x)/np.array(SFRDs_high_jets), label='SFR > $10^{%s}  M_{\odot} yr^{-1}$' %SFR_bin, color='crimson')

    redshift_ticks = [0, 1, 2, 4]
    age_ticks = cosmo.age(redshift_ticks).value
    ax1.set_xlabel('Cosmic Time [Gyr]')
    
    ax1.set_ylabel('$SFRD_{x-ray}/SFRD_{jets}$')
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xticks(age_ticks)
    ax2.set_xlabel('Redshift')
    ax2.set_xticklabels(redshift_ticks)
    
    ax1.legend()
    plt.savefig('SFR_ratio_x_jets.png', bbox_inches='tight', dpi=300)

#SFRD_SFR_ratio(model,size)

#%% SFRD convergence test
            
def SFRD_runs_plot(models, sizes,snaps):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    
    for i in range(len(models)):
        model = models[i]
        color = colors[i]
        size = sizes[i]
        SFRDs = []
        redshifts = []
        ages = []
    
        for j in range(len(snaps)):
            snap = snaps[j]
            SFR_density, z, age = SFRD(model, fb_fol, snap, size)
            SFRDs.append(SFR_density)
            redshifts.append(z)
            ages.append(age)
            ax1.plot(ages, np.log10(SFRDs), label=str(model), color=color)
            

    ax1.set_ylabel('log SFRD ($M_{\odot} y r^{-1} M p c^{-3}$)')
    ax1.set_ylim(-2.65, -0.4)  
    filename = 'SFRD_pure.pdf'

    redshift_ticks = [0, 1, 2, 4, 6]
    age_ticks = cosmo.age(redshift_ticks).value
    ax1.set_xlabel('Cosmic Time [Gyr]')
    
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xticks(age_ticks)
    ax2.set_xlabel('Redshift')
    ax2.set_xticklabels(redshift_ticks)
    
    ax1.legend()
    plt.savefig(filename, bbox_inches='tight')


# models on server
models = ['m50n512','m25n512','m25n256', 'm100n1024']
sizes = [50,25,25,100]
snaps = [151, 142,125, 104,90, 78,71, 62, 50,42, 36,30,26,22]


#SFRD_runs_plot(models,sizes,snaps)


