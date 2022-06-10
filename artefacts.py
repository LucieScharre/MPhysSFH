import sys
#import pygad as pg
import yt
import math
from yt.units.yt_array import YTQuantity
import caesar
import h5py
import os
import numpy as np
from pygadgetreader import *
import matplotlib.pyplot as plt
import matplotlib as mlp
#import OBSSMF as obs
from random import sample
from scipy import stats
from matplotlib.offsetbox import AnchoredText
from matplotlib.offsetbox import (TextArea, DrawingArea, OffsetImage,
                                  AnnotationBbox)

model = "m50n512"  
size=50   

fb_types = ['noagn','nojet']   

fb_type = 'noagn'
fb_type = 'nojet'
fb_type = '7jk'
snap = 151

fig, ax = plt.subplots(nrows=2, ncols=2)
plt.subplots_adjust(wspace=0, hspace=0)
fig.set_size_inches(11,9)

n =0
for fb_type in fb_types:
    infile = '/Users/luciescharre/Desktop/Uni/MPhys/data_scripts/cats/%s/%s_%03d.hdf5' %(fb_type,model,snap)
    sim = caesar.load(infile)
    
    mH_tot = np.array([i.halo.masses['total'] for i in sim.galaxies if i.central==1])
    #i.masses only instead to do only the mass of the central galaxy
    # do i.halo.masses to do total mass of halo that has a central galaxy
    mH_stellar = np.array([i.masses['stellar'] for i in sim.galaxies if i.central==1])
    SFR = np.array([i.halo.sfr for i in sim.galaxies if i.central==1])
    
    halo_pos = np.array([i.halo.pos for i in sim.galaxies if i.central == 1])
    r200_halo = np.array([i.halo.virial_quantities['r200c']for i in sim.galaxies if i.central==1])
    r500_halo = np.array([i.halo.virial_quantities['r500c']for i in sim.galaxies if i.central==1])
    r2500_halo = np.array([i.halo.virial_quantities['r2500c']for i in sim.galaxies if i.central==1])
    
    SFR = SFR[mH_stellar > 1e10]
    mH_tot = mH_tot[mH_stellar > 1e10]
    halo_pos = halo_pos[mH_stellar > 1e10]
    r200_halo = r200_halo[mH_stellar > 1e10]
    mH_stellar = mH_stellar[mH_stellar > 1e10]
    
    SFR = SFR[r200_halo > 0]
    mH_tot = mH_tot[r200_halo > 0]
    halo_pos = halo_pos[r200_halo > 0]
    mH_stellar = mH_stellar[r200_halo > 0]
    r200_halo = r200_halo[r200_halo > 0]
    
    
    sSFR = SFR/mH_stellar
    mH_ratio = np.log10(mH_stellar/mH_tot)
    
    art_thresh = -1.0
    # artefact galaxies above -1.1
    mH_ratio_art = mH_ratio[mH_ratio > art_thresh]
    halo_pos_art = halo_pos[mH_ratio > art_thresh]
    r200_halo_art = r200_halo[mH_ratio > art_thresh]
    
    halo_pos_others = halo_pos[mH_ratio < art_thresh]
    
    #r500_halo_art = r500_halo[mH_ratio > art_thresh]
    #r2500_halo_art = r2500_halo[mH_ratio > art_thresh]
    print(len(mH_ratio_art))
    #print(galaxy_pos_art)
    
    
    # find distance to closest galaxy for all of them

    #print(av_dist)
    min_dist_all = []
    
    for pos1 in halo_pos:
        distances = []
        for pos2 in halo_pos:
            vector =  pos1 - pos2
            distance = np.linalg.norm(vector)
            distances.append(distance)
        
        distances = np.array(distances)
        min_dist_all.append(min(distances[distances>0]))
        
    

    min_dist_art = np.array(min_dist_all)[mH_ratio > art_thresh]
    r200_halo_art
    
    min_dist_others = np.array(min_dist_all)[mH_ratio < art_thresh]
    r200_halo_others = r200_halo[mH_ratio < art_thresh]
    
    min_dist_close_to_art= []
    r200_halo_close_to_art= []
    


    
    min_dist_close_to_art= []
    r200_halo_close_to_art= []
    
    for i in range(len(min_dist_art)):
        found = 'no'
        for j in range(len(halo_pos_others)):
            vector =  halo_pos_art[i] - halo_pos_others[j]
            distance = np.linalg.norm(vector)
        
            if distance == min_dist_art[i]:
                min_dist_close_to_art.append(min_dist_others[i])
                r200_halo_close_to_art.append(r200_halo_others[i])
                found = 'yes'
            
        if found == 'no':
            min_dist_close_to_art.append(np.nan)
            r200_halo_close_to_art.append(np.nan)
            
    
    print(len(min_dist_close_to_art),len(min_dist_art))
    
    
    
    s =20

    
    ax.flat[n].scatter(min_dist_all,r200_halo, s=s, label='other halos',color = 'teal')
    ax.flat[n].scatter(min_dist_art,r200_halo_art,color ='r' , s=s, label='$M_{\star}/M_{H} > 10^{-1}$')
    ax.flat[n].scatter(min_dist_close_to_art,r200_halo_close_to_art,color ='orange' , s=s, label='closest to halos with $M_{\star}/M_{H} > 10^{-1}$')
    ax.flat[n].set_xlabel('Distance to closest halo [kpc/h]')
    
    
    ax.flat[n].set_xlim(50,11000)
    ax.flat[n].set_ylim(50,1000)
    if n==0:
       
       ax.flat[n].set_ylabel('R200 [kpc/h]')
    else:
          ax.flat[n].legend(prop={'size': 9})
    #plt.savefig('mindist_R200_%s.png' %(fb_type))
    #plt.show()
    
    
    #plt.scatter(min_dist_all,min_dist_all/r200_halo, s=s)
    ax.flat[n+2].scatter(min_dist_art,min_dist_art/r200_halo_close_to_art,color ='r' , s=s, label='$M_{\star}/M_{H} > 10^{-1}$')
    #plt.scatter(min_dist_close_to_art,np.array(min_dist_close_to_art)/np.array(r200_halo_close_to_art),color ='orange' , s=s)
    ax.flat[n+2].set_xlabel('Distance to closest halo [kpc/h]')
    if n == 0:
        ax.flat[n+2].set_ylabel('Distance/R200 of closest halo')
    ax.flat[n+2].set_xlim(50,11000)
    ax.flat[n+2].set_ylim(0,19)
    ax_up1 = ax.flat[n+2].twiny()
    #ax_up1.set_xscale('log')
    #ax_up1.set_xlabel('z = 0')
    ax_up1.set_xlim(ax.flat[n+2].get_xlim())
    ax_up1.set_xticklabels([])
    #plt.ylim(0,50)
    #ax.flat[n+2].legend()
    n+=1
    
    ax_up0 = ax.flat[0].twiny()
    ax_up0.set_xlabel('Stellar')
    ax_up0.set_xlim(ax.flat[0].get_xlim())
    ax_up0.set_xticklabels([])
            
    ax_up1 = ax.flat[1].twiny()
    ax_up1.set_xlabel('AGN winds')
    ax_up1.set_xlim(ax.flat[1].get_xlim())
    ax_up1.set_xticklabels([])

  #%% 
ax.flat[0].set_xticklabels([], direction='inout')
ax.flat[1].set_xticklabels([], direction='inout')     
ax.flat[1].set_yticklabels([], direction='inout')
ax.flat[3].set_yticklabels([], direction='inout')        
 
plt.savefig('mindist_R200_closest_halo_ratio_2x2.pdf',bbox_inches='tight')
plt.show()



