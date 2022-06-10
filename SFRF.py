#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 15:28:58 2022

@author: luciescharre
"""

import caesar
import numpy as np
import matplotlib.pyplot as plt    

#%% compute SFRF
def SFRFunc(model,fb_type,snap):
    infile = '/Users/luciescharre/Desktop/Uni/MPhys/data_scripts/cats/%s/%s_%03d.hdf5' %(fb_type,model,snap)
    
    sim = caesar.load(infile)

    h = sim.simulation.hubble_constant    
    z = sim.simulation.redshift
    #print(h)
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
    galaxy_sfr = np.array([i.sfr for i in sim.galaxies if i.central==1])
    
    # take the logarithm of masses 
    galaxy_sfr_lg = np.log10(galaxy_sfr[galaxy_sfr > 0]) 
    
    #galaxy_sfr_lg = galaxy_sfr
    # collect SFR
    
    Nbins = 15
    bins = np.linspace(galaxy_sfr_lg.min(),galaxy_sfr_lg.max(),Nbins) 
    
    step = ((galaxy_sfr_lg.max()-galaxy_sfr_lg.min())/Nbins)
    
    # create histgram from that octant octants 
    N, bin_edges = np.histogram(galaxy_sfr_lg,bins=bins)

    # compute the stellar mass function based on that 
    SFRF = np.log10(N/step)
    
    # find the bincentres 
    bincen = np.zeros(len(bins)-1)
    for i in range(len(bins)-1):
        bincen[i] = 0.5*(bins[i]+bins[i+1])
    
    return bincen, SFRF,  z_round     
    
#%%
model = 'm50n512'
size = 50

fb_types = ['7jk','nox']  
fb_titles = ['x-ray','jets']  
linestyles = ['solid', '--']

snaps = [50,62,78]
colors = ['orangered','orange','g']

# plots SFRF for the nox and 7jk run
def plot_SFRF(model, snaps, colors,fb_types, fb_titles, linestyles):    
    for j in range(len(fb_types)):
        fb_type = fb_types[j]
        #fig, ax = plt.subplots(1) 
        linestyle = linestyles[j]
        fb_title = fb_titles[j]
        
        
        for i in range(len(snaps)): 
            snap = snaps[i]
            color = colors[i]
            
            sfrcen, SFRF, z_round = SFRFunc(model,fb_type,snap)
            plt.plot(10**sfrcen, 10**SFRF, label = '%s z=%s' %(fb_title, z_round),color=color, linestyle = linestyle)
            #plt.fill_between(bincen, SMF-var, SMF+var,facecolor= 'lightgrey')

    plt.legend()
    plt.xlabel("SFR [$\mathrm{M}_{\odot} \mathrm{yr}^{-1}$]")
    plt.ylabel("Number of galaxies")
    plt.loglog()
    plt.savefig("SFRF.pdf")
    plt.show()        


# plots  SFRF ratio of nox and 7jk run
def plot_SFRF_ratio(model, snaps, colors):
    
    for i in range(len(snaps)): 
        snap = snaps[i]
        color = colors[i]
        
        fb_type = '7jk'
        sfrcen_fid, SFRF_fid, z_round_fid = SFRFunc(model,fb_type,snap)
        
        fb_type = 'nox'
        sfrcen_nox, SFRF_nox, z_round_nox = SFRFunc(model,fb_type,snap)
        
        ratio = SFRF_fid/SFRF_nox
        
        plt.plot(sfrcen_fid, ratio, label = 'z=%s' %(z_round_fid),color=color)

    plt.legend()
    plt.xlabel("log(SFR /($\mathrm{M}_{\odot} \mathrm{yr}^{-1}$))")
    plt.ylabel("log($N_{nox}/N_{fiducial}/$)")
    plt.savefig("SFRF_ratio.png")
 
#plot_SFRF(model, snaps, colors,fb_types, fb_titles, linestyles)
    