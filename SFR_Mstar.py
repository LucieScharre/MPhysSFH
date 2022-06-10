#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 12:18:08 2022

@author: luciescharre
"""

import caesar
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mlp

#%% function to read in quantities in bulk
def read_cat(sim, centralOnly = False):
    if centralOnly == True:
        mH_tot = np.array([i.halo.masses['total'] for i in sim.galaxies if i.central==1])
        #i.masses only instead to do only the mass of the central galaxy
        # do i.halo.masses to do total mass of halo that has a central galaxy
        mH_stellar = np.array([i.masses['stellar'] for i in sim.galaxies if i.central==1])
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


#%% scatter plots of sSFR versus stellar mass coloured by various quantities (z_data)

def sSFR_mStar_subplots(model, snaps,size,zdata_s):
    for zdata in zdata_s:
        # loop calling the master_SMF function
        if zdata == 'gas_frac':
            fb_fols = ['nox','7jk']  
            fb_types = ['jets','x-ray',]  
            fig, ax = plt.subplots(nrows=2, ncols=2)
            plt.subplots_adjust(wspace=0, hspace=0)
            fig.set_size_inches(11,9)
        
        else:
            fb_fols  = ['noagn','nojet','nox','7jk']  
            fb_types = ['stellar','AGN winds','jets','x-ray'] 

            fig, ax = plt.subplots(nrows=4, ncols=2)
            plt.subplots_adjust(wspace=0, hspace=0)
            fig.set_size_inches(9,13)
        
        n=0
        
        for i in range(len(fb_types)):
            for j in range(len(snaps)):
                
                fb_type = fb_types[i]
                fb_fol = fb_fols[i]
    
                snap = snaps[j]
                infile = '/Users/luciescharre/Desktop/Uni/MPhys/data_scripts/cats/%s/%s_%03d.hdf5' %(fb_fol,model,snap)
                sim = caesar.load(infile)
                
                h = sim.simulation.hubble_constant    
                z = sim.simulation.redshift
                if z > 1.99:
                    z_round =round(z) 
                elif 2 > z > 0:
                    z_round =round(z)
                print(z)
                print(' ')
                print(' ')
                
                
                mH_tot, mH_stellar, SFR, SFR_100, bhmdot, bh_fedd,mH_BH = read_cat(sim, centralOnly=True)
                
                mH_gas = np.array([i.masses['gas'] for i in sim.galaxies if i.central==1])

                
                BHm_frac = mH_BH/mH_stellar
                gas_frac = mH_gas/mH_stellar
                sSFR = SFR/mH_stellar
                
                #print(max(np.log10(mH_BH)))

                    
                if zdata == 'bh_fedd':
                    vmin = 1e-4
                    vmax = 1e1
                    ax.flat[n].set_ylim(5e-13,1e-8)
                    ax.flat[n].set_xlim(2e9,1e12)
                    
                    plot = ax.flat[n].scatter(mH_stellar,sSFR, c= bh_fedd, norm=mlp.colors.LogNorm(vmin,vmax), s=1,cmap="jet_r")

                    clabel = '$f_{edd}$'
                    
                    
                elif zdata == 'BHm_frac':      
                    mH_BH_frac_below = BHm_frac[mH_BH<10**7.5]
                    mH_BH_frac_above = BHm_frac[mH_BH>10**7.5]
                    
                    
                    vmin = 1e-3
                    vmax = 1e-2
                    ax.flat[n].set_ylim(5e-13,1e-8)
                    ax.flat[n].set_xlim(2e9,1e12)
                    
                    plot = ax.flat[n].scatter(mH_stellar, sSFR, c= BHm_frac, norm=mlp.colors.LogNorm(vmin,vmax), s=1,cmap="jet_r")

                    clabel = '$M_{BH}/M_{\star}$'
                    
                    
                elif zdata == 'BHm':     
                    mH_stellar_below = mH_stellar[mH_BH<10**7.5]
                    mH_stellar_above = mH_stellar[mH_BH>10**7.5]
                    
                    sSFR_below = sSFR[mH_BH<10**7.5]
                    sSFR_above = sSFR[mH_BH>10**7.5]
                    
                    mH_BH_below = mH_BH[mH_BH<10**7.5]
                    mH_BH_above = mH_BH[mH_BH>10**7.5]
                    
                    
                    vmin = 10**7.5
                    vmax = 10**9.5
                    ax.flat[n].set_ylim(5e-13,1e-8)
                    ax.flat[n].set_xlim(2e9,1e12)
                    
                    
                    plot = ax.flat[n].scatter(mH_stellar, sSFR, c= mH_BH, norm=mlp.colors.LogNorm(vmin,vmax), s=2,cmap="jet_r")
                    #plot = ax.flat[n].scatter(mH_stellar_below, sSFR_below, c= mH_BH_below, norm=mlp.colors.LogNorm(vmin,vmax), s=1,cmap="jet_r",marker = '.')
                    #plot = ax.flat[n].scatter(mH_stellar_above, sSFR_above, c= mH_BH_above, norm=mlp.colors.LogNorm(vmin,vmax), s=2,cmap="jet_r", edgecolor='k')
                    #plot.cmap.set_under('darkred')
                    clabel = '$M_{BH} [M_{\odot}]$'
                    
                    
                    
                    
                elif zdata == 'gas_frac':      
                    vmin = 0.001
                    vmax = 20
                    ax.flat[n].set_ylim(5e-13,1e-8)
                    ax.flat[n].set_xlim(2e8,1e12)
                    
                    plot = ax.flat[n].scatter(mH_stellar, sSFR, c= gas_frac, norm=mlp.colors.LogNorm(vmin,vmax), s=1,cmap="jet_r")
                    
                    clabel = '$f_{gas}$'
                
                    
                
                
                ax.flat[n].tick_params(right=True, top=True)   
                ax.flat[n].tick_params(axis="y",direction="inout")
                ax.flat[n].tick_params(axis="x",direction="inout")
                ax.flat[n].loglog()
                n +=1
                    
                    
        if zdata == 'gas_frac':
            ax.flat[0].set_xticklabels([], direction='inout')
            ax.flat[1].set_xticklabels([], direction='inout') 
            ax.flat[2].set_xlabel('$M_{\star}$ [$M_{\odot}]$')
            ax.flat[3].set_xlabel('$M_{\star}$ [$M_{\odot}]$')
            
            ax_up0 = ax.flat[0].twiny()
            ax_up0.set_xscale('log')
            ax_up0.set_xlabel('z = 2')
            ax_up0.set_xlim(ax.flat[0].get_xlim())
            ax_up0.set_xticklabels([])
            
            ax_up1 = ax.flat[1].twiny()
            ax_up1.set_xscale('log')
            ax_up1.set_xlabel('z = 0')
            ax_up1.set_xlim(ax.flat[1].get_xlim())
            ax_up1.set_xticklabels([])
            
            ax.flat[0].set_ylabel('jets \nsSFR [$M_{\odot}$/yr]')
            ax.flat[1].set_yticklabels([], direction='inout')
            ax.flat[2].set_ylabel('x-ray \nsSFR [$M_{\odot}$/yr]')
            ax.flat[3].set_yticklabels([], direction='inout')

            
                    
        else:
            ax.flat[0].set_xticklabels([], direction='inout')
            ax.flat[1].set_xticklabels([], direction='inout')   
            ax.flat[2].set_xticklabels([], direction='inout')
            ax.flat[3].set_xticklabels([], direction='inout')      
            ax.flat[4].set_xticklabels([], direction='inout')
            ax.flat[5].set_xticklabels([], direction='inout') 
            ax.flat[6].set_xlabel('$M_{\star}$ [$M_{\odot}]$')
            ax.flat[7].set_xlabel('$M_{\star}$ [$M_{\odot}]$')
            
            ax_up0 = ax.flat[0].twiny()
            ax_up0.set_xscale('log')
            ax_up0.set_xlabel('z = 2')
            ax_up0.set_xlim(ax.flat[0].get_xlim())
            ax_up0.set_xticklabels([])
            
            ax_up1 = ax.flat[1].twiny()
            ax_up1.set_xscale('log')
            ax_up1.set_xlabel('z = 0')
            ax_up1.set_xlim(ax.flat[1].get_xlim())
            ax_up1.set_xticklabels([])
            
            

            ax.flat[0].set_ylabel('stellar \nsSFR [$M_{\odot}$/yr]')
            ax.flat[1].set_yticklabels([], direction='inout')
            ax.flat[2].set_ylabel('AGN winds \nsSFR [$M_{\odot}$/yr]')
            ax.flat[3].set_yticklabels([], direction='inout')
            ax.flat[4].set_ylabel('jets \nsSFR [$M_{\odot}$/yr]')
            ax.flat[5].set_yticklabels([], direction='inout')
            ax.flat[6].set_ylabel('x-ray\nsSFR [$M_{\odot}$/yr]')
            ax.flat[7].set_yticklabels([], direction='inout')
            
    
  
        cbar_ax = fig.colorbar(plot, ax=ax.ravel().tolist(), shrink =0.5)
        cbar_ax.set_label(clabel)
        plt.show()

        fig.savefig('%anew_s_sSFR.pdf' %(zdata),bbox_inches='tight')
        #fig.savefig('edd.pdf')

            
snaps = [78,151]   
zdata_s = ['BHm_frac','bh_fedd','BHm','gas_frac']

model = 'm50n512'
size = 50

#sSFR_mStar_subplots(model, snaps, size,zdata_s) 

#%%  scatter plots with separated galaxies by feedback type according to thresholds 

def sSFR_mStar_subplots_2x1(model, snaps,size, fb_fols, fb_types,label_positions):
        
        for j in range(len(fb_types)):
            fig, ax = plt.subplots(nrows=2, ncols=1)
            plt.subplots_adjust(wspace=0, hspace=0)
            fig.set_size_inches(9,13)
            
            fb_type = fb_types[j]
            fb_fol = fb_fols[j]
            label_pos = label_positions[j]
            
            for i in range(len(snaps)):
                snap = snaps[i]
                infile = '/Users/luciescharre/Desktop/Uni/MPhys/data_scripts/cats/%s/%s_%03d.hdf5' %(fb_fol,model,snap)
                sim = caesar.load(infile)
                
                h = sim.simulation.hubble_constant    
                z = sim.simulation.redshift
                z_round =round(z)
                print(z)
                print(' ')
                print(' ')

                #i.masses only instead to do only the mass of the central galaxy
                # do i.halo.masses to do total mass of halo that has a central galaxy
                mH_tot = np.array([i.halo.masses['total'] for i in sim.galaxies if i.central==1])
                mH_stellar = np.array([i.masses['stellar'] for i in sim.galaxies if i.central==1])
                mH_BH = np.array([i.masses['bh'] for i in sim.galaxies if i.central==1])
                
                SFR = np.array([i.sfr for i in sim.galaxies if i.central==1])
                bh_fedd= np.array([i.bh_fedd for i in sim.galaxies if i.central==1])  
                mH_gas = np.array([i.masses['gas'] for i in sim.galaxies if i.central==1])

                
                mH_tot = mH_tot[mH_BH>0]
                mH_stellar = mH_stellar[mH_BH>0]
                
                SFR =SFR[mH_BH>0]
                bh_fedd= bh_fedd[mH_BH>0]
                mH_gas = mH_gas[mH_BH>0]
                mH_BH = mH_BH[mH_BH>0]
                
                
                BHm_frac = mH_BH/mH_stellar
                gas_frac = mH_gas/mH_stellar
                sSFR = SFR/mH_stellar
                
                separate = True
                # separate out different parts
                if separate == True:
                    # black hole mass

                    mH_stellar_small = mH_stellar[np.ma.masked_inside(mH_BH, 0, 10**7.5).mask]
                    mH_stellar_8 = mH_stellar[np.ma.masked_inside(mH_BH, 10**7.5, 10**8.5).mask]
                    mH_stellar_9 = mH_stellar[np.ma.masked_inside(mH_BH, 10**8.5, 10**9.5).mask]
                    mH_stellar_10 = mH_stellar[np.ma.masked_inside(mH_BH, 10**9.5, 10**11).mask]
                    
                    sSFR_small = sSFR[np.ma.masked_inside(mH_BH, 0, 10**7.5).mask]
                    sSFR_8 = sSFR[np.ma.masked_inside(mH_BH, 10**7.5, 10**8.5).mask]
                    sSFR_9 = sSFR[np.ma.masked_inside(mH_BH, 10**8.5, 10**9.5).mask]
                    sSFR_10 = sSFR[np.ma.masked_inside(mH_BH, 10**9.5, 10**11).mask]
                     
                    mH_BH_small = mH_BH[np.ma.masked_inside(mH_BH, 0, 10**7.5).mask]
                    mH_BH_8 = mH_BH[np.ma.masked_inside(mH_BH, 10**7.5, 10**8.5).mask]
                    mH_BH_9 = mH_BH[np.ma.masked_inside(mH_BH, 10**8.5, 10**9.5).mask]
                    mH_BH_10 = mH_BH[np.ma.masked_inside(mH_BH, 10**9.5, 10**11).mask]
                    
                    gas_frac_small = gas_frac[np.ma.masked_inside(mH_BH, 0, 10**7.5).mask]
                    gas_frac_8 = gas_frac[np.ma.masked_inside(mH_BH, 10**7.5, 10**8.5).mask]
                    gas_frac_9 = gas_frac[np.ma.masked_inside(mH_BH, 10**8.5, 10**9.5).mask]
                    gas_frac_10 = gas_frac[np.ma.masked_inside(mH_BH, 10**9.5, 10**11).mask]
                    
                    bh_fedd_small = bh_fedd[np.ma.masked_inside(mH_BH, 0, 10**7.5).mask]
                    bh_fedd_8 = bh_fedd[np.ma.masked_inside(mH_BH, 10**7.5, 10**8.5).mask]
                    bh_fedd_9 = bh_fedd[np.ma.masked_inside(mH_BH, 10**8.5, 10**9.5).mask]
                    bh_fedd_10 = bh_fedd[np.ma.masked_inside(mH_BH, 10**9.5, 10**11).mask]
                    
                    SFR_small = bh_fedd[np.ma.masked_inside(mH_BH, 0, 10**7.5).mask]
                    SFR_8 = bh_fedd[np.ma.masked_inside(mH_BH, 10**7.5, 10**8.5).mask]
                    SFR_9 = bh_fedd[np.ma.masked_inside(mH_BH, 10**8.5, 10**9.5).mask]
                    SFR_10 = bh_fedd[np.ma.masked_inside(mH_BH, 10**9.5, 10**11).mask]
                    
                    
        
     # eddington fraction              
    # ################################################################### WIND
                    # above 0.2 (no jets) plotted as circles, wind galaxies 
                    mH_stellar_small_fedd_above = mH_stellar_small[bh_fedd_small > 0.2]
                    mH_stellar_8_fedd_above = mH_stellar_8[bh_fedd_8 > 0.2]
                    mH_stellar_9_fedd_above = mH_stellar_9[bh_fedd_9 > 0.2]
                    mH_stellar_10_fedd_above = mH_stellar_10[bh_fedd_10 > 0.2]
                    
                    sSFR_small_fedd_above= sSFR_small[bh_fedd_small > 0.2]
                    sSFR_8_fedd_above = sSFR_8[bh_fedd_8 > 0.2]
                    sSFR_9_fedd_above = sSFR_9[bh_fedd_9 > 0.2]
                    sSFR_10_fedd_above = sSFR_10[bh_fedd_10 > 0.2]
                     
                    mH_BH_small_fedd_above= mH_BH_small[bh_fedd_small > 0.2]
                    mH_BH_8_fedd_above = mH_BH_8[bh_fedd_8 > 0.2]
                    mH_BH_9_fedd_above = mH_BH_9[bh_fedd_9 > 0.2]
                    mH_BH_10_fedd_above = mH_BH_10[bh_fedd_10 > 0.2]
                    
                    gas_frac_small_fedd_above = gas_frac_small[bh_fedd_small > 0.2]
                    gas_frac_8_fedd_above = gas_frac_8[bh_fedd_8 > 0.2]
                    gas_frac_9_fedd_above = gas_frac_9[bh_fedd_9 > 0.2]
                    gas_frac_10_fedd_above = gas_frac_10[bh_fedd_10 > 0.2]
                    
                    bh_fedd_small_fedd_above = bh_fedd_small[bh_fedd_small > 0.2]
                    bh_fedd_8_fedd_above = bh_fedd_8[bh_fedd_8 > 0.2]
                    bh_fedd_9_fedd_above = bh_fedd_9[bh_fedd_9 > 0.2]
                    bh_fedd_10_fedd_above = bh_fedd_10[bh_fedd_10 > 0.2]
                    
     # eddington fraction  
    # ################################################################### 
                    # below 0.2 (jets except for small black holes, both x-ray and no x-ray)
                   
                    mH_stellar_8_fedd_below = mH_stellar_8[bh_fedd_8 < 0.2]
                    mH_stellar_9_fedd_below = mH_stellar_9[bh_fedd_9 < 0.2]
                    mH_stellar_10_fedd_below = mH_stellar_10[bh_fedd_10 < 0.2]
                    
        
                    sSFR_8_fedd_below = sSFR_8[bh_fedd_8 < 0.2]
                    sSFR_9_fedd_below = sSFR_9[bh_fedd_9 < 0.2]
                    sSFR_10_fedd_below = sSFR_10[bh_fedd_10 < 0.2]
                     
                    
                    mH_BH_8_fedd_below = mH_BH_8[bh_fedd_8 < 0.2]
                    mH_BH_9_fedd_below = mH_BH_9[bh_fedd_9 < 0.2]
                    mH_BH_10_fedd_below = mH_BH_10[bh_fedd_10 < 0.2]
                    
                    
                    gas_frac_8_fedd_below = gas_frac_8[bh_fedd_8 < 0.2]
                    gas_frac_9_fedd_below = gas_frac_9[bh_fedd_9 < 0.2]
                    gas_frac_10_fedd_below = gas_frac_10[bh_fedd_10 < 0.2]
                    
                   
                    bh_fedd_8_fedd_below = bh_fedd_8[bh_fedd_8 < 0.2]
                    bh_fedd_9_fedd_below = bh_fedd_9[bh_fedd_9 < 0.2]
                    bh_fedd_10_fedd_below = bh_fedd_10[bh_fedd_10 < 0.2]
                    
                    SFR_8_fedd_below = SFR_8[bh_fedd_8 < 0.2]
                    SFR_9_fedd_below = SFR_9[bh_fedd_9 < 0.2]
                    SFR_10_fedd_below = SFR_10[bh_fedd_10 < 0.2]
                    
    
    
     # gas fraction                 
    #
                    ################################################################### X-ray
    
                    # gasfraction below 0.2 (winds + jets + x-ray)
                    
                    mH_stellar_8_fedd_below_fgas_below = mH_stellar_8_fedd_below[gas_frac_8_fedd_below < 0.2]
                    mH_stellar_9_fedd_below_fgas_below  = mH_stellar_9_fedd_below[gas_frac_9_fedd_below < 0.2]
                    mH_stellar_10_fedd_below_fgas_below  = mH_stellar_10_fedd_below[gas_frac_10_fedd_below < 0.2]
                    
                   
                    sSFR_8_fedd_below_fgas_below  = sSFR_8_fedd_below[gas_frac_8_fedd_below < 0.2]
                    sSFR_9_fedd_below_fgas_below  = sSFR_9_fedd_below[gas_frac_9_fedd_below < 0.2]
                    sSFR_10_fedd_below_fgas_below  = sSFR_10_fedd_below[gas_frac_10_fedd_below < 0.2]
                     
                    
                    mH_BH_8_fedd_below_fgas_below  = mH_BH_8_fedd_below[gas_frac_8_fedd_below < 0.2]
                    mH_BH_9_fedd_below_fgas_below  = mH_BH_9_fedd_below[gas_frac_9_fedd_below < 0.2]
                    mH_BH_10_fedd_below_fgas_below  = mH_BH_10_fedd_below[gas_frac_10_fedd_below < 0.2]
                    
                    
                    gas_frac_8_fedd_below_fgas_below  = gas_frac_8_fedd_below[gas_frac_8_fedd_below < 0.2]
                    gas_frac_9_fedd_below_fgas_below  = gas_frac_9_fedd_below[gas_frac_9_fedd_below < 0.2]
                    gas_frac_10_fedd_below_fgas_below  = gas_frac_10_fedd_below[gas_frac_10_fedd_below < 0.2]
                    
                    
                    bh_fedd_8_fedd_below_fgas_below  = bh_fedd_8_fedd_below[gas_frac_8_fedd_below < 0.2]
                    bh_fedd_9_fedd_below_fgas_below  = bh_fedd_9_fedd_below[gas_frac_9_fedd_below < 0.2]
                    bh_fedd_10_fedd_below_fgas_below  = bh_fedd_10_fedd_below[gas_frac_10_fedd_below < 0.2]
                    
                    SFR_8_fedd_below_fgas_below  = SFR_8_fedd_below[gas_frac_8_fedd_below < 0.2]
                    SFR_9_fedd_below_fgas_below  = SFR_9_fedd_below[gas_frac_9_fedd_below < 0.2]
                    SFR_10_fedd_below_fgas_below  = SFR_10_fedd_below[gas_frac_10_fedd_below < 0.2]
    
                    
      # gas fraction                 
    #               
                    ################################################################### JETS
    
                    # gasfraction above 0.2 (winds + jets)
                    
                    mH_stellar_8_fedd_below_fgas_above = mH_stellar_8_fedd_below[gas_frac_8_fedd_below > 0.2]
                    mH_stellar_9_fedd_below_fgas_above  = mH_stellar_9_fedd_below[gas_frac_9_fedd_below > 0.2]
                    mH_stellar_10_fedd_below_fgas_above  = mH_stellar_10_fedd_below[gas_frac_10_fedd_below > 0.2]
                    
                   
                    sSFR_8_fedd_below_fgas_above  = sSFR_8_fedd_below[gas_frac_8_fedd_below > 0.2]
                    sSFR_9_fedd_below_fgas_above  = sSFR_9_fedd_below[gas_frac_9_fedd_below > 0.2]
                    sSFR_10_fedd_below_fgas_above = sSFR_10_fedd_below[gas_frac_10_fedd_below > 0.2]
                     
                    
                    mH_BH_8_fedd_below_fgas_above  = mH_BH_8_fedd_below[gas_frac_8_fedd_below > 0.2]
                    mH_BH_9_fedd_below_fgas_above  = mH_BH_9_fedd_below[gas_frac_9_fedd_below > 0.2]
                    mH_BH_10_fedd_below_fgas_above  = mH_BH_10_fedd_below[gas_frac_10_fedd_below > 0.2]
                    
                    
                    gas_frac_8_fedd_below_fgas_above  = gas_frac_8_fedd_below[gas_frac_8_fedd_below > 0.2]
                    gas_frac_9_fedd_below_fgas_above  = gas_frac_9_fedd_below[gas_frac_9_fedd_below > 0.2]
                    gas_frac_10_fedd_below_fgas_above  = gas_frac_10_fedd_below[gas_frac_10_fedd_below > 0.2]
                    
                    
                    bh_fedd_8_fedd_below_fgas_above = bh_fedd_8_fedd_below[gas_frac_8_fedd_below > 0.2]
                    bh_fedd_9_fedd_below_fgas_above  = bh_fedd_9_fedd_below[gas_frac_9_fedd_below > 0.2]
                    bh_fedd_10_fedd_below_fgas_above  = bh_fedd_10_fedd_below[gas_frac_10_fedd_below > 0.2]
                    
                    
                    
                    
    #  BH fraction                
                    
                    fBH_small = mH_BH_small/mH_stellar_small
                    fBH_8_fedd_above = mH_BH_8_fedd_above/mH_stellar_8_fedd_above
                    fBH_9_fedd_above = mH_BH_9_fedd_above/mH_stellar_9_fedd_above
                    fBH_10_fedd_above = mH_BH_10_fedd_above/mH_stellar_10_fedd_above
                    
                    
                    fBH_8_fedd_below_fgas_below = mH_BH_8_fedd_below_fgas_below/mH_stellar_8_fedd_below_fgas_below
                    fBH_9_fedd_below_fgas_below = mH_BH_9_fedd_below_fgas_below/mH_stellar_9_fedd_below_fgas_below
                    fBH_10_fedd_below_fgas_below = mH_BH_10_fedd_below_fgas_below/mH_stellar_10_fedd_below_fgas_below
                    
                    fBH_8_fedd_below_fgas_above = mH_BH_8_fedd_below_fgas_above/mH_stellar_8_fedd_below_fgas_above
                    fBH_9_fedd_below_fgas_above = mH_BH_9_fedd_below_fgas_above/mH_stellar_9_fedd_below_fgas_above
                    fBH_10_fedd_below_fgas_above = mH_BH_10_fedd_below_fgas_above/mH_stellar_10_fedd_below_fgas_above


#             PLOT                

                ax.flat[i].set_ylim(5e-13,1e-8)
                ax.flat[i].set_xlim(2e9,1e12)
                
                
                # gas fraction as colour bar
                vmin = 0.1
                vmax = 20
                cmap = 'jet_r'
                clabel = '$f_{gas}$'
                
                # eddington fraction as colour bar
                vmin = 1e-4
                vmax = 1
                cmap = 'jet_r'
                clabel = '$f_{edd}$'
                
                # bh mass fraction as colour bar
                """
                vmin = 1e-3
                vmax = 1e-2
                cmap = 'jet'
                clabel = '$M_{BH}/M_{\star}$'
                """
                
                
                ###
                s_small = 5
                s_8 = 20
                s_9 = 40
                s_10 = 80
                
                # AGN winds of 4 sizes, plotted as circles
                plot = ax.flat[i].scatter(mH_stellar_small,sSFR_small, c= bh_fedd_small, norm=mlp.colors.LogNorm(vmin,vmax),s=s_small,cmap=cmap)
                plot = ax.flat[i].scatter(mH_stellar_8_fedd_above,sSFR_8_fedd_above, c= bh_fedd_8_fedd_above, norm=mlp.colors.LogNorm(vmin,vmax), s=s_8,cmap=cmap)
                plot = ax.flat[i].scatter(mH_stellar_9_fedd_above,sSFR_9_fedd_above, c= bh_fedd_9_fedd_above, norm=mlp.colors.LogNorm(vmin,vmax), s=s_9,cmap=cmap)
                plot = ax.flat[i].scatter(mH_stellar_10_fedd_above,sSFR_10_fedd_above, c= bh_fedd_10_fedd_above, norm=mlp.colors.LogNorm(vmin,vmax), s=s_10,cmap=cmap)
                
                # JETS of three sizes
                plot = ax.flat[i].scatter(mH_stellar_8_fedd_below_fgas_above,sSFR_8_fedd_below_fgas_above, c= bh_fedd_8_fedd_below_fgas_above, norm=mlp.colors.LogNorm(vmin,vmax), s=s_8,cmap=cmap, marker = 's')
                plot = ax.flat[i].scatter(mH_stellar_9_fedd_below_fgas_above,sSFR_9_fedd_below_fgas_above, c= bh_fedd_9_fedd_below_fgas_above, norm=mlp.colors.LogNorm(vmin,vmax), s=s_9,cmap=cmap, marker = 's')
                plot = ax.flat[i].scatter(mH_stellar_10_fedd_below_fgas_above,sSFR_10_fedd_below_fgas_above, c= bh_fedd_10_fedd_below_fgas_above, norm=mlp.colors.LogNorm(vmin,vmax), s=s_10,cmap=cmap, marker = 's')

                # X-RAY of three sizes
                plot = ax.flat[i].scatter(mH_stellar_8_fedd_below_fgas_below,sSFR_8_fedd_below_fgas_below, c= bh_fedd_8_fedd_below_fgas_below, norm=mlp.colors.LogNorm(vmin,vmax), s=s_8,cmap=cmap, marker = 'v')
                plot = ax.flat[i].scatter(mH_stellar_9_fedd_below_fgas_below,sSFR_9_fedd_below_fgas_below, c= bh_fedd_9_fedd_below_fgas_below, norm=mlp.colors.LogNorm(vmin,vmax), s=s_9,cmap=cmap, marker = 'v')
                plot = ax.flat[i].scatter(mH_stellar_10_fedd_below_fgas_below,sSFR_10_fedd_below_fgas_below, c= bh_fedd_10_fedd_below_fgas_below, norm=mlp.colors.LogNorm(vmin,vmax), s=s_10,cmap=cmap, marker = 'v')

                print(len(mH_stellar_10_fedd_below_fgas_below))
                
                
                ax.flat[i].tick_params(right=True, top=True)   
                ax.flat[i].tick_params(axis="y",direction="inout")
                ax.flat[i].tick_params(axis="x",direction="inout")
                ax.flat[i].loglog()
                
                ax.flat[i].set_ylabel('z = %s \nsSFR [$M_{\odot}$/yr]' %int(z_round))
               
                
                    
                    
   
            ax.flat[0].set_xticklabels([], direction='inout')
            #ax.flat[1].set_xticklabels([], direction='inout') 
            ax.flat[1].set_xlabel('$M_{\star}$ [$M_{\odot}]$')
            #ax.flat[3].set_xlabel('$M_{\star}$ [$M_{\odot}]$')
            
            ax_up0 = ax.flat[0].twiny()
            ax_up0.set_xscale('log')
            #ax_up0.set_xlabel('z = 2')
            ax_up0.set_xlim(ax.flat[0].get_xlim())
            ax_up0.set_xticklabels([])
            
            ax_up1 = ax.flat[1].twiny()
            ax_up1.set_xscale('log')
            #ax_up1.set_xlabel('z = 0')
            ax_up1.set_xlim(ax.flat[1].get_xlim())
            ax_up1.set_xticklabels([])
            
            

            

            cbar_ax = fig.colorbar(plot, ax=ax.ravel().tolist(), shrink =0.5)
            cbar_ax.set_label(clabel)
            #fig.suptitle('%s' %fb_type)
            
            
            plt.gcf().text(label_pos, 0.8, '%s' %(fb_type), fontsize=15)
            plt.show()
    
            fig.savefig('%s_sSFR_fedd.png' %fb_type,bbox_inches='tight')
            fig.savefig('%s_sSFR_fedd.pdf' %fb_type,bbox_inches='tight')
        #fig.savefig('edd.pdf')
            


fb_fols  = ['noagn','nojet','nox','7jk']  
fb_types = ['stellar','AGN winds','jets','x-ray'] 

label_positions = [0.784,0.76,0.796,0.786]
snaps = [78,151]  

#sSFR_mStar_subplots_2x1(model, snaps,size, fb_fols, fb_types,label_positions)





def sSFR_mStar_subplots_8x2(model, snaps,size, fb_fols, fb_types):
        fig, ax = plt.subplots(nrows=4, ncols=2)
        plt.subplots_adjust(wspace=0, hspace=0)
        fig.set_size_inches(10,13)
        n =0
        for i in range(len(fb_types)):
            fb_type = fb_types[i]
            fb_fol = fb_fols[i]
            for j in range(len(snaps)):
    
                snap = snaps[j]
                
                infile = '/Users/luciescharre/Desktop/Uni/MPhys/data_scripts/cats/%s/%s_%03d.hdf5' %(fb_fol,model,snap)
                sim = caesar.load(infile)
                
                h = sim.simulation.hubble_constant    
                z = sim.simulation.redshift
                z_round =round(z)
                print(z)
                print(' ')
                print(' ')

                #i.masses only instead to do only the mass of the central galaxy
                # do i.halo.masses to do total mass of halo that has a central galaxy
                mH_tot = np.array([i.halo.masses['total'] for i in sim.galaxies if i.central==1])
                mH_stellar = np.array([i.masses['stellar'] for i in sim.galaxies if i.central==1])
                mH_BH = np.array([i.masses['bh'] for i in sim.galaxies if i.central==1])
                
                SFR = np.array([i.halo.sfr for i in sim.galaxies if i.central==1])
                bh_fedd= np.array([i.halo.bh_fedd for i in sim.galaxies if i.central==1])  
                mH_gas = np.array([i.masses['gas'] for i in sim.galaxies if i.central==1])

                
                mH_tot = mH_tot[mH_BH>0]
                mH_stellar = mH_stellar[mH_BH>0]
                
                SFR =SFR[mH_BH>0]
                bh_fedd= bh_fedd[mH_BH>0]
                mH_gas = mH_gas[mH_BH>0]
                mH_BH = mH_BH[mH_BH>0]
                
                """
                mH_tot = mH_tot[bh_fedd>1e-3]
                mH_stellar = mH_stellar[bh_fedd>1e-3]
                SFR =SFR[bh_fedd>1e-3]
                mH_gas = mH_gas[bh_fedd>1e-3]
                mH_BH = mH_BH[bh_fedd>1e-3]
                bh_fedd= bh_fedd[bh_fedd>1e-3]
                """
                
                BHm_frac = mH_BH/mH_stellar
                gas_frac = mH_gas/mH_stellar
                sSFR = SFR/mH_stellar
                
                
                #print(max(np.log10(mH_BH)))
                #print(np.log10(mH_BH[sSFR <2e-11]))
                separate = True
                # separate out different parts
                if separate == True:
                
                    # black hole mass
                    mH_stellar_small = mH_stellar[np.ma.masked_inside(mH_BH, 0, 10**7.5).mask]
                    mH_stellar_8 = mH_stellar[np.ma.masked_inside(mH_BH, 10**7.5, 10**8.5).mask]
                    mH_stellar_9 = mH_stellar[np.ma.masked_inside(mH_BH, 10**8.5, 10**9.5).mask]
                    mH_stellar_10 = mH_stellar[np.ma.masked_inside(mH_BH, 10**9.5, 10**11).mask]
                    
                    sSFR_small = sSFR[np.ma.masked_inside(mH_BH, 0, 10**7.5).mask]
                    sSFR_8 = sSFR[np.ma.masked_inside(mH_BH, 10**7.5, 10**8.5).mask]
                    sSFR_9 = sSFR[np.ma.masked_inside(mH_BH, 10**8.5, 10**9.5).mask]
                    sSFR_10 = sSFR[np.ma.masked_inside(mH_BH, 10**9.5, 10**11).mask]
                     
                    mH_BH_small = mH_BH[np.ma.masked_inside(mH_BH, 0, 10**7.5).mask]
                    mH_BH_8 = mH_BH[np.ma.masked_inside(mH_BH, 10**7.5, 10**8.5).mask]
                    mH_BH_9 = mH_BH[np.ma.masked_inside(mH_BH, 10**8.5, 10**9.5).mask]
                    mH_BH_10 = mH_BH[np.ma.masked_inside(mH_BH, 10**9.5, 10**11).mask]
                    
                    gas_frac_small = gas_frac[np.ma.masked_inside(mH_BH, 0, 10**7.5).mask]
                    gas_frac_8 = gas_frac[np.ma.masked_inside(mH_BH, 10**7.5, 10**8.5).mask]
                    gas_frac_9 = gas_frac[np.ma.masked_inside(mH_BH, 10**8.5, 10**9.5).mask]
                    gas_frac_10 = gas_frac[np.ma.masked_inside(mH_BH, 10**9.5, 10**11).mask]
                    
                    bh_fedd_small = bh_fedd[np.ma.masked_inside(mH_BH, 0, 10**7.5).mask]
                    bh_fedd_8 = bh_fedd[np.ma.masked_inside(mH_BH, 10**7.5, 10**8.5).mask]
                    bh_fedd_9 = bh_fedd[np.ma.masked_inside(mH_BH, 10**8.5, 10**9.5).mask]
                    bh_fedd_10 = bh_fedd[np.ma.masked_inside(mH_BH, 10**9.5, 10**11).mask]
                    
                    SFR_small = bh_fedd[np.ma.masked_inside(mH_BH, 0, 10**7.5).mask]
                    SFR_8 = bh_fedd[np.ma.masked_inside(mH_BH, 10**7.5, 10**8.5).mask]
                    SFR_9 = bh_fedd[np.ma.masked_inside(mH_BH, 10**8.5, 10**9.5).mask]
                    SFR_10 = bh_fedd[np.ma.masked_inside(mH_BH, 10**9.5, 10**11).mask]
                    
        
     # eddington fraction              
    # ################################################################### WIND
                    # above 0.2 (no jets) plotted as circles, wind galaxies 
                    mH_stellar_small_fedd_above = mH_stellar_small[bh_fedd_small > 0.2]
                    mH_stellar_8_fedd_above = mH_stellar_8[bh_fedd_8 > 0.2]
                    mH_stellar_9_fedd_above = mH_stellar_9[bh_fedd_9 > 0.2]
                    mH_stellar_10_fedd_above = mH_stellar_10[bh_fedd_10 > 0.2]
                    
                    sSFR_small_fedd_above= sSFR_small[bh_fedd_small > 0.2]
                    sSFR_8_fedd_above = sSFR_8[bh_fedd_8 > 0.2]
                    sSFR_9_fedd_above = sSFR_9[bh_fedd_9 > 0.2]
                    sSFR_10_fedd_above = sSFR_10[bh_fedd_10 > 0.2]
                     
                    mH_BH_small_fedd_above= mH_BH_small[bh_fedd_small > 0.2]
                    mH_BH_8_fedd_above = mH_BH_8[bh_fedd_8 > 0.2]
                    mH_BH_9_fedd_above = mH_BH_9[bh_fedd_9 > 0.2]
                    mH_BH_10_fedd_above = mH_BH_10[bh_fedd_10 > 0.2]
                    
                    gas_frac_small_fedd_above = gas_frac_small[bh_fedd_small > 0.2]
                    gas_frac_8_fedd_above = gas_frac_8[bh_fedd_8 > 0.2]
                    gas_frac_9_fedd_above = gas_frac_9[bh_fedd_9 > 0.2]
                    gas_frac_10_fedd_above = gas_frac_10[bh_fedd_10 > 0.2]
                    
                    bh_fedd_small_fedd_above = bh_fedd_small[bh_fedd_small > 0.2]
                    bh_fedd_8_fedd_above = bh_fedd_8[bh_fedd_8 > 0.2]
                    bh_fedd_9_fedd_above = bh_fedd_9[bh_fedd_9 > 0.2]
                    bh_fedd_10_fedd_above = bh_fedd_10[bh_fedd_10 > 0.2]
                    
     # eddington fraction  
    # ################################################################### 
                    # below 0.2 (jets except for small black holes, both x-ray and no x-ray)
                   
                    mH_stellar_8_fedd_below = mH_stellar_8[bh_fedd_8 < 0.2]
                    mH_stellar_9_fedd_below = mH_stellar_9[bh_fedd_9 < 0.2]
                    mH_stellar_10_fedd_below = mH_stellar_10[bh_fedd_10 < 0.2]
                    
        
                    sSFR_8_fedd_below = sSFR_8[bh_fedd_8 < 0.2]
                    sSFR_9_fedd_below = sSFR_9[bh_fedd_9 < 0.2]
                    sSFR_10_fedd_below = sSFR_10[bh_fedd_10 < 0.2]
                     
                    
                    mH_BH_8_fedd_below = mH_BH_8[bh_fedd_8 < 0.2]
                    mH_BH_9_fedd_below = mH_BH_9[bh_fedd_9 < 0.2]
                    mH_BH_10_fedd_below = mH_BH_10[bh_fedd_10 < 0.2]
                    
                    
                    gas_frac_8_fedd_below = gas_frac_8[bh_fedd_8 < 0.2]
                    gas_frac_9_fedd_below = gas_frac_9[bh_fedd_9 < 0.2]
                    gas_frac_10_fedd_below = gas_frac_10[bh_fedd_10 < 0.2]
                    
                   
                    bh_fedd_8_fedd_below = bh_fedd_8[bh_fedd_8 < 0.2]
                    bh_fedd_9_fedd_below = bh_fedd_9[bh_fedd_9 < 0.2]
                    bh_fedd_10_fedd_below = bh_fedd_10[bh_fedd_10 < 0.2]
                    
                    SFR_8_fedd_below = SFR_8[bh_fedd_8 < 0.2]
                    SFR_9_fedd_below = SFR_9[bh_fedd_9 < 0.2]
                    SFR_10_fedd_below = SFR_10[bh_fedd_10 < 0.2]
                    
    
    
     # gas fraction                 
    
                    ################################################################### X-ray
    
                    # gasfraction below 0.2 (winds + jets + x-ray)
                    
                    mH_stellar_8_fedd_below_fgas_below = mH_stellar_8_fedd_below[gas_frac_8_fedd_below < 0.2]
                    mH_stellar_9_fedd_below_fgas_below  = mH_stellar_9_fedd_below[gas_frac_9_fedd_below < 0.2]
                    mH_stellar_10_fedd_below_fgas_below  = mH_stellar_10_fedd_below[gas_frac_10_fedd_below < 0.2]
                    
                   
                    sSFR_8_fedd_below_fgas_below  = sSFR_8_fedd_below[gas_frac_8_fedd_below < 0.2]
                    sSFR_9_fedd_below_fgas_below  = sSFR_9_fedd_below[gas_frac_9_fedd_below < 0.2]
                    sSFR_10_fedd_below_fgas_below  = sSFR_10_fedd_below[gas_frac_10_fedd_below < 0.2]
                     
                    
                    mH_BH_8_fedd_below_fgas_below  = mH_BH_8_fedd_below[gas_frac_8_fedd_below < 0.2]
                    mH_BH_9_fedd_below_fgas_below  = mH_BH_9_fedd_below[gas_frac_9_fedd_below < 0.2]
                    mH_BH_10_fedd_below_fgas_below  = mH_BH_10_fedd_below[gas_frac_10_fedd_below < 0.2]
                    
                    
                    gas_frac_8_fedd_below_fgas_below  = gas_frac_8_fedd_below[gas_frac_8_fedd_below < 0.2]
                    gas_frac_9_fedd_below_fgas_below  = gas_frac_9_fedd_below[gas_frac_9_fedd_below < 0.2]
                    gas_frac_10_fedd_below_fgas_below  = gas_frac_10_fedd_below[gas_frac_10_fedd_below < 0.2]
                    
                    
                    bh_fedd_8_fedd_below_fgas_below  = bh_fedd_8_fedd_below[gas_frac_8_fedd_below < 0.2]
                    bh_fedd_9_fedd_below_fgas_below  = bh_fedd_9_fedd_below[gas_frac_9_fedd_below < 0.2]
                    bh_fedd_10_fedd_below_fgas_below  = bh_fedd_10_fedd_below[gas_frac_10_fedd_below < 0.2]
                    
                    SFR_8_fedd_below_fgas_below  = SFR_8_fedd_below[gas_frac_8_fedd_below < 0.2]
                    SFR_9_fedd_below_fgas_below  = SFR_9_fedd_below[gas_frac_9_fedd_below < 0.2]
                    SFR_10_fedd_below_fgas_below  = SFR_10_fedd_below[gas_frac_10_fedd_below < 0.2]
    
                    
      # gas fraction                 
                
                    ################################################################### JETS
    
                    # gasfraction above 0.2 (winds + jets)
                    
                    mH_stellar_8_fedd_below_fgas_above = mH_stellar_8_fedd_below[gas_frac_8_fedd_below > 0.2]
                    mH_stellar_9_fedd_below_fgas_above  = mH_stellar_9_fedd_below[gas_frac_9_fedd_below > 0.2]
                    mH_stellar_10_fedd_below_fgas_above  = mH_stellar_10_fedd_below[gas_frac_10_fedd_below > 0.2]
                    
                   
                    sSFR_8_fedd_below_fgas_above  = sSFR_8_fedd_below[gas_frac_8_fedd_below > 0.2]
                    sSFR_9_fedd_below_fgas_above  = sSFR_9_fedd_below[gas_frac_9_fedd_below > 0.2]
                    sSFR_10_fedd_below_fgas_above = sSFR_10_fedd_below[gas_frac_10_fedd_below > 0.2]
                     
                    
                    mH_BH_8_fedd_below_fgas_above  = mH_BH_8_fedd_below[gas_frac_8_fedd_below > 0.2]
                    mH_BH_9_fedd_below_fgas_above  = mH_BH_9_fedd_below[gas_frac_9_fedd_below > 0.2]
                    mH_BH_10_fedd_below_fgas_above  = mH_BH_10_fedd_below[gas_frac_10_fedd_below > 0.2]
                    
                    
                    gas_frac_8_fedd_below_fgas_above  = gas_frac_8_fedd_below[gas_frac_8_fedd_below > 0.2]
                    gas_frac_9_fedd_below_fgas_above  = gas_frac_9_fedd_below[gas_frac_9_fedd_below > 0.2]
                    gas_frac_10_fedd_below_fgas_above  = gas_frac_10_fedd_below[gas_frac_10_fedd_below > 0.2]
                    
                    
                    bh_fedd_8_fedd_below_fgas_above = bh_fedd_8_fedd_below[gas_frac_8_fedd_below > 0.2]
                    bh_fedd_9_fedd_below_fgas_above  = bh_fedd_9_fedd_below[gas_frac_9_fedd_below > 0.2]
                    bh_fedd_10_fedd_below_fgas_above  = bh_fedd_10_fedd_below[gas_frac_10_fedd_below > 0.2]
                    
                    
                    
                    
    #  BH fraction                
                    
                    fBH_small = mH_BH_small/mH_stellar_small
                    fBH_8_fedd_above = mH_BH_8_fedd_above/mH_stellar_8_fedd_above
                    fBH_9_fedd_above = mH_BH_9_fedd_above/mH_stellar_9_fedd_above
                    fBH_10_fedd_above = mH_BH_10_fedd_above/mH_stellar_10_fedd_above
                    
                    
                    fBH_8_fedd_below_fgas_below = mH_BH_8_fedd_below_fgas_below/mH_stellar_8_fedd_below_fgas_below
                    fBH_9_fedd_below_fgas_below = mH_BH_9_fedd_below_fgas_below/mH_stellar_9_fedd_below_fgas_below
                    fBH_10_fedd_below_fgas_below = mH_BH_10_fedd_below_fgas_below/mH_stellar_10_fedd_below_fgas_below
                    
                    fBH_8_fedd_below_fgas_above = mH_BH_8_fedd_below_fgas_above/mH_stellar_8_fedd_below_fgas_above
                    fBH_9_fedd_below_fgas_above = mH_BH_9_fedd_below_fgas_above/mH_stellar_9_fedd_below_fgas_above
                    fBH_10_fedd_below_fgas_above = mH_BH_10_fedd_below_fgas_above/mH_stellar_10_fedd_below_fgas_above


#             PLOT                

                ax.flat[n].set_ylim(5e-13,1e-8)
                ax.flat[n].set_xlim(2e9,1e12)
                
                
                # gas fraction as colour bar
                vmin = 0.1
                vmax = 20
                cmap = 'jet_r'
                clabel = '$f_{gas}$'
                
                # eddington fraction as colour bar
                vmin = 1e-4
                vmax = 1
                cmap = 'jet_r'
                clabel = '$f_{edd}$'
                
                # bh mass fraction as colour bar
                """
                vmin = 1e-3
                vmax = 1e-2
                cmap = 'jet'
                clabel = '$M_{BH}/M_{\star}$'
                """
                
                
                ###
                s_small = 5
                s_8 = 10
                s_9 = 15
                s_10 = 30
                
                s_small = 1
                s_8 = 5
                s_9 = 10
                s_10 = 20
                
                # AGN winds of 4 sizes, plotted as circles
                plot = ax.flat[n].scatter(mH_stellar_small,sSFR_small, c= bh_fedd_small, norm=mlp.colors.LogNorm(vmin,vmax),s=s_small,cmap=cmap)
                plot = ax.flat[n].scatter(mH_stellar_8_fedd_above,sSFR_8_fedd_above, c= bh_fedd_8_fedd_above, norm=mlp.colors.LogNorm(vmin,vmax), s=s_8,cmap=cmap)
                plot = ax.flat[n].scatter(mH_stellar_9_fedd_above,sSFR_9_fedd_above, c= bh_fedd_9_fedd_above, norm=mlp.colors.LogNorm(vmin,vmax), s=s_9,cmap=cmap)
                plot = ax.flat[n].scatter(mH_stellar_10_fedd_above,sSFR_10_fedd_above, c= bh_fedd_10_fedd_above, norm=mlp.colors.LogNorm(vmin,vmax), s=s_10,cmap=cmap)
                
                # JETS of three sizes
                plot = ax.flat[n].scatter(mH_stellar_8_fedd_below_fgas_above,sSFR_8_fedd_below_fgas_above, c= bh_fedd_8_fedd_below_fgas_above, norm=mlp.colors.LogNorm(vmin,vmax), s=s_8,cmap=cmap, marker = 's')
                plot = ax.flat[n].scatter(mH_stellar_9_fedd_below_fgas_above,sSFR_9_fedd_below_fgas_above, c= bh_fedd_9_fedd_below_fgas_above, norm=mlp.colors.LogNorm(vmin,vmax), s=s_9,cmap=cmap, marker = 's')
                plot = ax.flat[n].scatter(mH_stellar_10_fedd_below_fgas_above,sSFR_10_fedd_below_fgas_above, c= bh_fedd_10_fedd_below_fgas_above, norm=mlp.colors.LogNorm(vmin,vmax), s=s_10,cmap=cmap, marker = 's')

                # X-RAY of three sizes
                plot = ax.flat[n].scatter(mH_stellar_8_fedd_below_fgas_below,sSFR_8_fedd_below_fgas_below, c= bh_fedd_8_fedd_below_fgas_below, norm=mlp.colors.LogNorm(vmin,vmax), s=s_8,cmap=cmap, marker = 'v')
                plot = ax.flat[n].scatter(mH_stellar_9_fedd_below_fgas_below,sSFR_9_fedd_below_fgas_below, c= bh_fedd_9_fedd_below_fgas_below, norm=mlp.colors.LogNorm(vmin,vmax), s=s_9,cmap=cmap, marker = 'v')
                plot = ax.flat[n].scatter(mH_stellar_10_fedd_below_fgas_below,sSFR_10_fedd_below_fgas_below, c= bh_fedd_10_fedd_below_fgas_below, norm=mlp.colors.LogNorm(vmin,vmax), s=s_10,cmap=cmap, marker = 'v')

                #print(len(mH_stellar_10_fedd_below_fgas_below))
                
                
                ax.flat[n].tick_params(right=True, top=True)   
                ax.flat[n].tick_params(axis="y",direction="inout")
                ax.flat[n].tick_params(axis="x",direction="inout")
                ax.flat[n].loglog()
                n +=1
                print(n)    
   
        ax.flat[0].set_xticklabels([], direction='inout')
        ax.flat[1].set_xticklabels([], direction='inout')   
        ax.flat[2].set_xticklabels([], direction='inout')
        ax.flat[3].set_xticklabels([], direction='inout')      
        ax.flat[4].set_xticklabels([], direction='inout')
        ax.flat[5].set_xticklabels([], direction='inout') 
        ax.flat[6].set_xlabel('$M_{\star}$ [$M_{\odot}]$')
        ax.flat[7].set_xlabel('$M_{\star}$ [$M_{\odot}]$')
        
        ax_up0 = ax.flat[0].twiny()
        ax_up0.set_xscale('log')
        ax_up0.set_xlabel('z = 2')
        ax_up0.set_xlim(ax.flat[0].get_xlim())
        ax_up0.set_xticklabels([])
        
        ax_up1 = ax.flat[1].twiny()
        ax_up1.set_xscale('log')
        ax_up1.set_xlabel('z = 0')
        ax_up1.set_xlim(ax.flat[1].get_xlim())
        ax_up1.set_xticklabels([])
        
        

        ax.flat[0].set_ylabel('stellar \nsSFR [$M_{\odot}$/yr]')
        ax.flat[1].set_yticklabels([], direction='inout')
        ax.flat[2].set_ylabel('AGN winds \nsSFR [$M_{\odot}$/yr]')
        ax.flat[3].set_yticklabels([], direction='inout')
        ax.flat[4].set_ylabel('jets \nsSFR [$M_{\odot}$/yr]')
        ax.flat[5].set_yticklabels([], direction='inout')
        ax.flat[6].set_ylabel('x-ray\nsSFR [$M_{\odot}$/yr]')
        ax.flat[7].set_yticklabels([], direction='inout')
        


        cbar_ax = fig.colorbar(plot, ax=ax.ravel().tolist(), shrink =0.5)
        cbar_ax.set_label(clabel)
        #fig.suptitle('%s' %fb_type)
        
        
        #plt.gcf().text(label_pos, 0.8, '%s' %(fb_type), fontsize=15)
        plt.show()

        #fig.savefig('all_sSFR_fedd.png' ,bbox_inches='tight')
        #fig.savefig('all_sSFR_fedd.pdf' ,bbox_inches='tight')
            
            
#sSFR_mStar_subplots_8x2(model, snaps,size, fb_fols, fb_types)        
