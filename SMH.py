#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 13:20:55 2021

@author: luciescharre
"""

import caesar
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mlp
from scipy import stats
from matplotlib.offsetbox import AnchoredText


#%%
# function to read in batch of quantities from catalogues 

def read_cat(sim, centralOnly = False):
    if centralOnly == True:
        mH_tot = np.array([i.halo.masses['total'] for i in sim.galaxies if i.central==1])
        #i.masses only instead to do only the mass of the central galaxy
        # do i.halo.masses to do total mass of halo that has a central galaxy
        mH_stellar = np.array([i.halo.masses['stellar'] for i in sim.galaxies if i.central==1])
        mH_BH = np.array([i.masses['bh'] for i in sim.galaxies if i.central==1])
        SFR = np.array([i.halo.sfr for i in sim.galaxies if i.central==1])
        SFR_100 = np.array([i.halo.sfr_100 for i in sim.galaxies if i.central==1])
        bhmdot = np.array([i.halo.bhmdot for i in sim.galaxies if i.central==1])
        bh_fedd= np.array([i.halo.bh_fedd for i in sim.galaxies if i.central==1])
        
    else:
            
        mH_tot = np.array([i.masses['total'] for i in sim.halos])     # Collect the halo masses 
        mH_stellar = np.array([i.masses['stellar'] for i in sim.halos] )
        SFR = np.array([i.sfr for i in sim.halos])
        SFR_100 = np.array([i.sfr_100 for i in sim.halos])
        bhmdot = np.array([i.bhmdot for i in sim.halos]) #Central black hole accretion rate in Msun/yr, 
        bh_fedd = np.array([i.bh_fedd for i in sim.halos])  #central BH eddington ratio
        
        
    return  mH_tot, mH_stellar, SFR, SFR_100, bhmdot, bh_fedd, mH_BH
             

#%%
def percentile16(y):
   return(np.percentile(y,16))

def percentile84(y):
   return(np.percentile(y,84))
    

#%% creating the heatmaps   
model = "m50n512"  
size=50   


fb_titles = ['+ stellar','+ AGN winds','+ jets','+ x-ray'] 
fb_types = ['noagn','nojet','nox','7jk'] 

snaps = [151,78]

def SMH_ratio(model, size, fb_types, fb_titles, snaps):
    for i in range(len(snaps)):
        snap = snaps[i]
        count = 9
        fig, ax = plt.subplots(nrows=2, ncols=2)
        plt.subplots_adjust(wspace=0, hspace=0)
        fig.set_size_inches(12,9)
       
        for j in range(len(fb_types)):
            fb_type = fb_types[j]
            fb_title = fb_titles[j]
            
            count -= 1
        #################################### LOADING CATALOGUES & PRINTING IMPORTANT INFO    
            infile = '/Users/luciescharre/Desktop/Uni/MPhys/data_scripts/cats/%s/%s_%03d.hdf5' %(fb_type,model,snap)
            sim = caesar.load(infile)
            
            h = sim.simulation.hubble_constant    
            z = sim.simulation.redshift
            if z > 2:
                z_round =float(int(round(z)))
            elif 2 > z > 0:
                z_round =round(z,1)
            print(z)
            print(' ')
            print(' ')
            
            
            #mH_tot, mH_stellar, SFR, SFR_100, bhmdot, bh_fedd,mH_BH = read_cat(sim,  centralOnly = True)
            
            mH_tot = np.array([i.halo.masses['total'] for i in sim.galaxies if i.central==1])
            #i.masses only instead to do only the mass of the central galaxy
            # do i.halo.masses to do total mass of halo that has a central galaxy
            mH_stellar = np.array([i.masses['stellar'] for i in sim.galaxies if i.central==1])
            SFR = np.array([i.halo.sfr for i in sim.galaxies if i.central==1])
            
            sSFR = SFR[mH_stellar > 0]/mH_stellar[mH_stellar > 0]
            
            sSFRmin = min(sSFR)
            sSFRmax = max(sSFR)
            #print(sSFRmin, sSFRmax)
            
            sSFRmin = 1e-11
            sSFRmax = 1e-8
            
            if z < 1:
                res = 40
                #sSFRmin = 5e-11
            
            elif 1 <= z <= 2: 
                res = 30
                #sSFRmin = 5e-11
                
            elif 2 < z < 3: 
                res = 30
                #sSFRmin = 1e-10
            else: 
                res = 20
            
            res = 40
            
            Nbins1 = res
            Nbins2 = res
            
            
            
            if snap == 151:
            # hist will be an array with size Nbins1 X Nbins2; 
            # each element contains the number of haloes falling in the corresponding 2D bin of Mhalo and Mstar.
                hist,xedges,yedges = np.histogram2d(np.log10(mH_tot[mH_stellar > 0]),np.log10(mH_stellar/mH_tot[mH_stellar > 0]),bins=(Nbins1,Nbins2))
                whist,xedges,yedges = np.histogram2d(np.log10(mH_tot[mH_stellar > 0]),np.log10(mH_stellar/mH_tot[mH_stellar > 0]),bins=(Nbins1, Nbins2),weights=sSFR)

            
            else: 
                # hist will be an array with size Nbins1 X Nbins2; 
                # each element contains the number of haloes falling in the corresponding 2D bin of Mhalo and Mstar.
                hist,xedges,yedges = np.histogram2d(np.log10(mH_tot[mH_stellar > 0]),np.log10(mH_stellar/mH_tot[mH_stellar > 0]), bins=(xedges,yedges))
                whist,xedges,yedges = np.histogram2d(np.log10(mH_tot[mH_stellar > 0]),np.log10(mH_stellar/mH_tot[mH_stellar > 0]),bins=(xedges,yedges),weights=sSFR)
            
            #2D array of size Nbins1 X Nbins2, where each element contains by construction the MEAN sSFR of the haloes in the corresponding 2D bin
            r = whist/hist
    
    
            extent = [min(np.log10(mH_tot[mH_stellar > 0])) , max(np.log10(mH_tot[mH_stellar > 0])), min(np.log10(mH_stellar/mH_tot[mH_stellar > 0])) , max(np.log10(mH_stellar/mH_tot[mH_stellar > 0]))] 
            extent = [10,14,-3.5,-0.5]
            
            #im=plt.imshow(r.T, norm=mlp.colors.LogNorm(vmin=sSFRmin, vmax=sSFRmax), extent = extent, interpolation=None,origin='lower',aspect='auto')
    
        
            im = ax.flat[j].imshow(r.T, norm=mlp.colors.LogNorm(vmin=sSFRmin, vmax=sSFRmax), extent = extent, interpolation=None,origin='lower',aspect='auto',cmap="jet_r")
            
            #cbar = plt.colorbar()
            #cbar.set_label('sSFR [$yr^{-1}$]')
            
            ax.flat[j].set_ylim(-3.6,-0.6)
            ax.flat[j].set_xlim(9.8,14.2)
            
            #ax.flat[j].annotate('%s' %(fb_title),xycoords='data', xy=(-2,12),xytext=(0.5, 0.5), textcoords='axes fraction', color ='k')
            ax.flat[j].tick_params(axis="y",direction="inout")
            ax.flat[j].tick_params(axis="x",direction="inout")
            
            
            
            if j == 0 or j==1:
                ax.flat[j].set_xticklabels([])
            
            if j == 1 or j == 3:    
                ax.flat[j].set_yticklabels([])
            
            #fig, ax = plt.subplots()
            at = AnchoredText('%s' %(fb_title), prop=dict(size=12), frameon=True, loc='upper right')
            at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
            ax.flat[j].add_artist(at)
        
        ax.flat[2].set_ylabel('$log(M_{\star}/M_{H})$')
        ax.flat[0].set_ylabel('$log(M_{\star}/M_{H})$')
        ax.flat[3].set_xlabel('$log(M_{H}$ / [$M_{\odot}$])')
        ax.flat[2].set_xlabel('$log(M_{H}$ / [$M_{\odot}$])')
        
    
        
        cbar_ax = fig.colorbar(im, ax=ax.ravel().tolist(), shrink =0.5)
        cbar_ax.set_label('sSFR [$yr^{-1}$]')
        
        plt.gcf().text(0.76, 0.8, 'z = %s' %(z_round), fontsize=25)
        plt.savefig('%02d_%s_%s_SMR_SFR.png' %(count,z_round, fb_type), bbox_inches='tight',dpi=500)
    
        plt.show()


SMH_ratio(model, size, fb_types, fb_titles, snaps)        


#%% SMHM relation
                
def averaged_SMHM(model,snaps,fb_types, fb_names):
    
    linestyles = ['dotted','dashed','solid','dashdot',(0, (3, 10, 1, 10, 1, 10))]
    colors = ['gold','orange','darkorange','orangered','plum','mediumorchid','darkorchid','indigo', 'lime','limegreen','darkgreen','darkolivegreen']
   
    fig, ax2 = plt.subplots(1)
    n = 0
    for i in range(len(fb_types)):
        fb_type = fb_types[i]
        fb_name = fb_names[i]
        linestyle = linestyles[i] 
        #color = ccolors[i] 
        
        print(i)
        for j in range(len(snaps)):
        #################################### LOADING CATALOGUES & PRINTING IMPORTANT INFO    
            snap = snaps[j]
            color = colors[n]
            infile = '/Users/luciescharre/Desktop/Uni/MPhys/data_scripts/cats/%s/%s_%03d.hdf5' %(fb_type,model,snap)
            sim = caesar.load(infile)
            
            h = sim.simulation.hubble_constant    
            z = sim.simulation.redshift

            z_round =int(round(z))    
            print(z)
            print(' ')
            print(' ')

            #print(mH_tot)
            #########################################################
            #mH_tot, mH_stellar, SFR, SFR_100, bhmdot, bh_fedd = read_cat(sim)
            mH_tot, mH_stellar, SFR, SFR_100, bhmdot, bh_fedd, mH_BH = read_cat(sim, centralOnly=True)
            
        
            ##################################################    binning
            Nbin = 10
            
            #M_data, bin_edges, binnumber = stats.binned_statistic(np.log10(mH_tot), mH_stellar,statistic='mean', bins=Nbin+1)
            #M_std = stats.binned_statistic(np.log10(mH_tot), mH_stellar,statistic='std', bins=Nbin+1)[0]  
           
            #R_means = stats.binned_statistic(np.log10(mH_tot), mH_ratio,statistic='mean', bins=Nbin+1)[0]
            #R_std = stats.binned_statistic(np.log10(mH_tot), mH_ratio,statistic='std', bins=Nbin+1)[0]  
            
            M_data, bin_edges, binnumber = stats.binned_statistic(np.log10(mH_tot), mH_stellar,statistic='median', bins=Nbin+1)
            M_upper= stats.binned_statistic(np.log10(mH_tot), mH_stellar,statistic=percentile84, bins=Nbin+1)[0]  
            M_lower = stats.binned_statistic(np.log10(mH_tot), mH_stellar,statistic=percentile16, bins=Nbin+1)[0]  

            
            # find the bincentres 
            bincen = np.zeros(len(bin_edges)-1)
            for i in range(len(bin_edges)-1):
                bincen[i] = 0.5*(bin_edges[i]+bin_edges[i+1])

            
            #ax1.errorbar(10**bincen,M_means,yerr=0, label = 'z=%s %s' %(z_round,fb_type),linewidth = 1.0, linestyle=linestyle, color= color)
            #ax2.errorbar(10**bincen,R_means, yerr=0,label = 'z=%s %s' %(z_round,fb_type), linewidth = 1.0, linestyle=linestyle, color= color)
            
            
            ax2.errorbar(10**bincen,M_data, label = '%s z=%s' %(fb_name,z_round), linewidth = 1.0, linestyle=linestyle,color= color)
            ax2.fill_between(10**bincen, M_data-M_lower, M_data+M_upper,facecolor= color, alpha = 0.05)
            n+=1
            
    ax2.hlines(5.2e9,10**10.5,1e14,color ='k', label = '$M_{0}$')
    ax2.hlines(10**9.5,10**10.5,1e14,color ='grey', label = '$M_{\star,BH}$')
    ax2.loglog()
    ax2.set_xlabel('$M_{H}$ [$M_{\odot}]$')
    ax2.set_ylabel('$M_{\star}$ [$M_{\odot}]$')
    ax2.legend(bbox_to_anchor=(1.05, 1))
    plt.xlim(10**10.5,1e14)
    plt.ylim(5e8,1e12)
    plt.savefig('SMHM.pdf', bbox_inches='tight')
    plt.show()
    

            
snaps = [36,51,78,151]
fb_types = ['nofb','noagn','7jk']           
fb_names = ['no fb','stellar','x-ray']          
averaged_SMHM(model,snaps,fb_types, fb_names)          
        

    
    
    