
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 12:10:14 2021

@author: luciescharre
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 10:56:25 2021

@author: luciescharre
"""

import numpy as np
import matplotlib.pyplot as plt
import caesar
from astropy.cosmology import WMAP9 as cosmo
from matplotlib.offsetbox import AnchoredText


model = "m50n512"
size = 50

# server 
fb_types = ['+ x-ray','+ jets','+ AGN winds',  '+ stellar', 'no fb']
fb_fols = ['s50',  's50nox', 's50nojet', 's50noagn', 's50nofb'] 

# local files
fb_fols  = ['nofb','noagn','nojet','nox','7jk']  
fb_types= [ 'no fb',  'stellar','AGN winds', 'jets','x-ray']  

snaps = [36, 50, 62, 78, 84, 90, 97, 104, 125, 137, 142, 151]

#%% computing the bins

def SFR_halobins(model, fb_fol, snap, size):
    # LOADING CATALOGUES & PRINTING IMPORTANT INFO
    # infile = '/disk04/rad/sim/%s/%s/Groups/%s_%03d.hdf5' % (
    #    model, fb_fol, model, snap)
    infile = '/Users/luciescharre/Desktop/Uni/MPhys/data_scripts/cats/%s/%s_%03d.hdf5' %(fb_fol,model,snap)
    #infile = '/disk04/rad/sim/%s/%s/Groups/%s_%03d.hdf5' % (
    #    model, fb_fol, model, snap)
    print(infile)
    sim = caesar.load(infile)
    z = sim.simulation.redshift
    print(z)
    print(' ')
    print(' ')

    # select halos with central galaxies
    # find the mass bins 1e11,12,13,14, include .5 to each side
    # for all haloIDs in those mass bins, take the median, upper and lower percentile of the SFR

    mH_tot = np.array([i.halo.masses['total']
                      for i in sim.galaxies if i.central == 1])
    mH_tot_lg = mH_tot_lg = np.log10(mH_tot)
    mH_stellar = np.array([i.halo.masses['stellar']
                          for i in sim.galaxies if i.central == 1])
    hal_SFR = np.array([i.halo.sfr for i in sim.galaxies if i.central == 1])

    mH_11 = []
    mH_12 = []
    mH_13 = []
    mH_14 = []

    # appends to lists of indices belonging to relevant mass bins
    for mi in range(len(mH_tot_lg)):

        if 10.5 < mH_tot_lg[mi] < 11.5:
            mH_11.append(mi)

        elif 11.5 < mH_tot_lg[mi] < 12.5:
            mH_12.append(mi)

        elif 12.5 < mH_tot_lg[mi] < 13.5:
            mH_13.append(mi)

        elif 13.5 < mH_tot_lg[mi] < 14.5:
            mH_14.append(mi)

    mH_bins_ids = [mH_11, mH_12, mH_13, mH_14]
    # print(mH_11)

    SFR_bins_med = [[], [], [], []]
    sSFR_bins_med = [[], [], [], []]
    nSFR_bins_med = [[], [], [], []]

    SFR_bins_16 = [[], [], [], []]
    sSFR_bins_16 = [[], [], [], []]
    nSFR_bins_16 = [[], [], [], []]

    SFR_bins_84 = [[], [], [], []]
    sSFR_bins_84 = [[], [], [], []]
    nSFR_bins_84 = [[], [], [], []]

    # using the indices, I can then find the median and percentiles of the SFR
    for j in range(len(mH_bins_ids)):
        
        mH_bin = mH_bins_ids[j]
        SFR_temp = []
        sSFR_temp = []
        nSFR_temp = []

        # collect the SFR
        for i in mH_bin:
            SFR_temp.append(hal_SFR[i])
            sSFR_temp.append(hal_SFR[i]/mH_stellar[i])
            nSFR_temp.append(hal_SFR[i]/mH_tot[i])

        # compute the percentiles

        if len(SFR_temp) == 0:
            SFR_bins_med[j] = float("nan")
            sSFR_bins_med[j] = float("nan")
            nSFR_bins_med[j] = float("nan")
            SFR_bins_16[j] = float("nan")
            sSFR_bins_16[j] = float("nan")
            nSFR_bins_16[j] = float("nan")
            SFR_bins_84[j] = float("nan")
            sSFR_bins_84[j] = float("nan")
            nSFR_bins_84[j] = float("nan")

        else:
            SFR_bins_med[j] = np.percentile(SFR_temp, 50)
            sSFR_bins_med[j] = np.percentile(sSFR_temp, 50)
            nSFR_bins_med[j] = np.percentile(nSFR_temp, 50)

            SFR_bins_16[j] = np.percentile(SFR_temp, 16)
            sSFR_bins_16[j] = np.percentile(sSFR_temp, 16)
            nSFR_bins_16[j] = np.percentile(nSFR_temp, 16)

            SFR_bins_84[j] = np.percentile(SFR_temp, 84)
            sSFR_bins_84[j] = np.percentile(sSFR_temp, 84)
            nSFR_bins_84[j] = np.percentile(nSFR_temp, 84)

    SFR_data = [SFR_bins_med, SFR_bins_16, SFR_bins_84]
    sSFR_data = [sSFR_bins_med, sSFR_bins_16, sSFR_bins_84]
    nSFR_data = [nSFR_bins_med, nSFR_bins_16, nSFR_bins_84]

    age = cosmo.age(z).value
    return z, age, SFR_data, sSFR_data, nSFR_data


#%% plotting in subplots
def halobins_plot(model, fb_fols, fb_types,snaps,size):
    
    # loop calling the master_SMF function
    fig, ax = plt.subplots(nrows=3, ncols=2)
    plt.subplots_adjust(wspace=0, hspace=0)
    fig.set_size_inches(9,12)

    # make different plots for different fb types
    
    for i in range(len(fb_fols)):
        if i !=5:
            fb_type = fb_types[i]
            fb_fol = fb_fols[i]

            redshifts = []
            ages = []

            nSFR_bins_z = []    
            nSFR_bins_16_z = []
            nSFR_bins_84_z = []
    
            for snap in snaps:
                z, age, SFR_data, sSFR_data, nSFR_data = SFR_halobins(
                    model, fb_fol, snap, size)
        
                redshifts.append(z)
                ages.append(age)
        
                nSFR_bins_z.append(nSFR_data[0])
                nSFR_bins_16_z.append(nSFR_data[1])
                nSFR_bins_84_z.append(nSFR_data[2])
        
            # Transpose arrays to gather values for the bins together, should contain one
            # series of values with changing z in 4 arrays for wach of the mass bins
            nSFR_bins_z = np.array(nSFR_bins_z).T
            nSFR_bins_16_z = np.array(nSFR_bins_16_z).T
            nSFR_bins_84_z = np.array(nSFR_bins_84_z).T
    
            
            
            line1 = ax.flat[i].plot(ages, nSFR_bins_z[3], color='g')
            ax.flat[i].fill_between(
                ages, nSFR_bins_16_z[3], nSFR_bins_84_z[3], facecolor='lightgreen', alpha=0.9)
    
            line2 = ax.flat[i].plot(ages, nSFR_bins_z[2], color='b')
            ax.flat[i].fill_between(
                ages, nSFR_bins_16_z[2], nSFR_bins_84_z[2], facecolor='lightblue', alpha=0.7)
    
            line3 =ax.flat[i].plot(ages, nSFR_bins_z[1], color='r')
            ax.flat[i].fill_between(
                ages, nSFR_bins_16_z[1], nSFR_bins_84_z[1], facecolor='salmon', alpha=0.5)
    
            line4 =ax.flat[i].plot(ages, nSFR_bins_z[0], color='k')
            ax.flat[i].fill_between(
                ages, nSFR_bins_16_z[0], nSFR_bins_84_z[0], facecolor='lightgrey', alpha=0.5)
            
            lines = [line1,line2,line3,line4]
            labels = ['14','13','12','11']
            
            ax.flat[i].set_ylim(1e-14, 2e-10)  # nSFR
            ax.flat[i].set_yscale('log')
            
            ax.flat[i].tick_params(right=True, top=True)   
            ax.flat[i].tick_params(axis="y",direction="inout")
            ax.flat[i].tick_params(axis="x",direction="inout")
            
            at = AnchoredText('%s' %(fb_type), prop=dict(size=12), frameon=True, loc='upper right')
            at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
            ax.flat[i].add_artist(at)
            
    ax.flat[5].set_yscale('log')
    ax.flat[5].plot(ages, nSFR_bins_z[0], color='w')
    ax.flat[5].tick_params(right=True, top=True)   
    ax.flat[5].tick_params(axis="y",direction="inout")
    ax.flat[5].tick_params(axis="x",direction="inout")

        
    ax.flat[0].set_xticklabels([], direction='inout')
    ax.flat[1].set_xticklabels([], direction='inout')   
    ax.flat[2].set_xticklabels([], direction='inout')
    ax.flat[3].set_xticklabels([], direction='inout')          
    ax.flat[4].set_xlabel('Cosmic Time [Gyr]')
    ax.flat[5].set_xlabel('Cosmic Time [Gyr]')
    ax.flat[5].set_xlim(ax.flat[4].get_xlim())
    
    ax.flat[0].set_ylabel('nSFR [$\mathrm{yr}^{-1}$]')
    ax.flat[1].set_yticklabels([], direction='inout')
    ax.flat[2].set_ylabel('nSFR [$\mathrm{yr}^{-1}$]')
    ax.flat[3].set_yticklabels([], direction='inout')
    ax.flat[4].set_ylabel('nSFR [$\mathrm{yr}^{-1}$]')
    ax.flat[5].set_yticklabels([], direction='inout')
    ax.flat[5].set_ylim(ax.flat[4].get_ylim())

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
    
    fig.legend(lines,labels=labels,title = 'log($M_{h}$/$M_{\odot})$',loc= [0.71,0.14])
    plt.savefig('nSFR_halobins_subplots.pdf' , bbox_inches='tight')

     


halobins_plot(model, fb_fols, fb_types,snaps,size) 
#%%    
    
def halobins_plot_ind(model, fb_fols, fb_types,snaps,size):    
    for i in range(len(fb_fols)):
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        fig.set_size_inches(8,6)
        
        if i !=5:
            fb_type = fb_types[i]
            fb_fol = fb_fols[i]

            redshifts = []
            ages = []

            nSFR_bins_z = []    
            nSFR_bins_16_z = []
            nSFR_bins_84_z = []
    
            for snap in snaps:
                z, age, SFR_data, sSFR_data, nSFR_data = SFR_halobins(
                    model, fb_fol, snap, size)
        
                redshifts.append(z)
                ages.append(age)
        
                nSFR_bins_z.append(nSFR_data[0])
                nSFR_bins_16_z.append(nSFR_data[1])
                nSFR_bins_84_z.append(nSFR_data[2])
        
            # Transpose arrays to gather values for the bins together, should contain one
            # series of values with changing z in 4 arrays for wach of the mass bins
            nSFR_bins_z = np.array(nSFR_bins_z).T
            nSFR_bins_16_z = np.array(nSFR_bins_16_z).T
            nSFR_bins_84_z = np.array(nSFR_bins_84_z).T
    
            
            
            line1 = ax1.plot(ages, nSFR_bins_z[3], color='g',label='$\mathrm{M}_{\mathrm{h}}=10^{14} \mathrm{M}_{\odot}$')
            ax1.fill_between(
                ages, nSFR_bins_16_z[3], nSFR_bins_84_z[3], facecolor='lightgreen', alpha=0.9)
    
            line2 = ax1.plot(ages, nSFR_bins_z[2], color='b',label='$\mathrm{M}_{\mathrm{h}}=10^{13} \mathrm{M}_{\odot}$')
            ax1.fill_between(
                ages, nSFR_bins_16_z[2], nSFR_bins_84_z[2], facecolor='lightblue', alpha=0.7)
    
            line3 =ax1.plot(ages, nSFR_bins_z[1], color='r',label='$\mathrm{M}_{\mathrm{h}}=10^{12} \mathrm{M}_{\odot}$')
            ax1.fill_between(
                ages, nSFR_bins_16_z[1], nSFR_bins_84_z[1], facecolor='salmon', alpha=0.5)
    
            line4 =ax1.plot(ages, nSFR_bins_z[0], color='k',label='$\mathrm{M}_{\mathrm{h}}=10^{11} \mathrm{M}_{\odot}$')
            ax1.fill_between(
                ages, nSFR_bins_16_z[0], nSFR_bins_84_z[0], facecolor='lightgrey', alpha=0.5)
            
            #lines = [line1,line2,line3,line4]
            #labels = ['14','13','12','11']
            
            ax1.set_ylim(1e-14, 2e-10)  # nSFR
            ax1.set_yscale('log')
            
            ax1.tick_params(right=True, top=True)   
            ax1.tick_params(axis="y",direction="inout")
            ax1.tick_params(axis="x",direction="inout")
            
            at = AnchoredText('%s' %(fb_type), prop=dict(size=12), frameon=True, loc='upper right')
            at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
            ax1.add_artist(at)
            
    
                
            ax1.set_xlabel('Cosmic Time [Gyr]')
            ax1.set_ylabel('nSFR [$\mathrm{yr}^{-1}$]')
           
            redshift_ticks = [0, 1, 2, 4, 6]
            age_ticks = cosmo.age(redshift_ticks).value
            
            ax_up0 = ax1.twiny()
            ax_up0.set_xlim(ax1.get_xlim())
            ax_up0.set_xticks(age_ticks)
            ax_up0.set_xlabel('Redshift')
            ax_up0.set_xticklabels(redshift_ticks)
           # plt.title('%s' %fb_type)
            
            ax1.legend(loc = 'lower left')
            plt.savefig('nSFR_halobins_%s.png' %fb_type , bbox_inches='tight', dpi=300)


halobins_plot_ind(model, fb_fols, fb_types,snaps,size)




