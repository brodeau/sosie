#!/usr/bin/env python

#    B a r a K u d a
#
#    L. Brodeau, 2017

import sys
import os
from string import replace
import numpy as nmp

from netCDF4 import Dataset

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#import matplotlib.font_manager as font_manager
import matplotlib.patches as patches

import warnings
warnings.filterwarnings("ignore")

import matplotlib.dates as mdates

#import barakuda_orca   as bo
#import barakuda_ncio   as bnc
#import barakuda_colmap as bcm
#import barakuda_plot   as bp

import barakuda_tool as bt


color_dark_blue = '#091459'
color_red = '#AD0000'
b_gre = '#3A783E'
b_prp = '#8C008C'
b_org = '#ED7C4C'

color_text_colorbar = 'k'
color_stff_colorbar = 'k'
#color_continents    = '#9C5536'
#color_continents    = '#EDD0AB'
color_continents    = '0.75'



rDPI=100.



fig_ext='png'
#fig_ext='svg'

jt1=0 ; jt2=0

narg = len(sys.argv)

print narg

if narg != 4 and narg != 6 :
    print 'Usage: '+sys.argv[0]+' <file> <variable model> <variable ephem> (<jt1> <jt2>)'; sys.exit(0)

cf_in = sys.argv[1] ; cv_mdl=sys.argv[2] ; cv_eph=sys.argv[3]
if narg == 6:
    jt1=int(sys.argv[4]) ; jt2=int(sys.argv[5])

bt.chck4f(cf_in)

id_in    = Dataset(cf_in)
vt_epoch = id_in.variables['time'][:]
vmodel   = id_in.variables[cv_mdl][:]
vephem   = id_in.variables[cv_eph][:]
vlat     = id_in.variables['latitude'][:]
id_in.close()
print "  => READ!"





font_corr = 1.2
params = { 'font.family':'Ubuntu',
           'font.weight':    'normal',
           'font.size':       int(12*font_corr),
           'legend.fontsize': int(12*font_corr),
           'xtick.labelsize': int(11*font_corr),
           'ytick.labelsize': int(12*font_corr),
           'axes.labelsize':  int(12*font_corr),
           'axes.titlesize' : 20,
           'lines.linewidth' : 8,
           'lines.markersize' : 10 }
#           'figure.facecolor': 'w' }
mpl.rcParams.update(params)




nbr = len(vt_epoch)




# Initial raw plot:
# Create Matplotlib time array:

vtime = nmp.zeros(nbr)
for jt in range(nbr): vtime[jt] = mdates.epoch2num(vt_epoch[jt])




ii=nbr/300
ib=max(ii-ii%10,1)
print ' ii , ib =', ii, ib

xticks_d=30.*ib

fig = plt.figure(num = 1, figsize=(12,7), facecolor='w', edgecolor='k')
ax1 = plt.axes([0.07, 0.24, 0.9, 0.75])

ax1.set_xticks(vtime[::xticks_d])
ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
plt.xticks(rotation='60')


plt.plot(vtime, vephem, '-', color=color_dark_blue, linewidth=2, label='Satellite ("'+cv_eph+'")', zorder=10)
plt.plot(vtime, vmodel, '-', color=b_org, linewidth=2,  label='Model ("'+cv_mdl+'")', zorder=15)
ax1.set_ylim(-0.68,0.68) ; ax1.set_xlim(vtime[0],vtime[nbr-1])
plt.xlabel('Time [seconds since 1970]')
plt.ylabel('SSH [m]')
#cstep = '%5.5i'%(jpnij)
ax1.grid(color='k', linestyle='-', linewidth=0.3)

plt.legend(bbox_to_anchor=(0.55, 0.98), ncol=1, shadow=True, fancybox=True)


#ax2 = ax1.twinx()
#ax2.set_ylim(25.,68.)
#plt.plot(vtime, vlat, '--', color='0.4', linewidth=1.5, label='latitude', zorder=0.1)
#ax2.set_ylabel(r'Latitude [$^\circ$North]', color='0.4')
#[t.set_color('0.4') for t in ax2.yaxis.get_ticklabels()]

plt.savefig('fig_raw_data.'+fig_ext, dpi=120, transparent=True)
plt.close(1)






# Now 1 figure per sequence !

vmask = vmodel.mask

(idx_ok,) = nmp.where(vmask==False) # indexes with valid values!


print idx_ok

nbr_v = len(idx_ok)

print ' *** '+str(nbr_v)+' valid points out of '+str(nbr)+' !'

print vmask

# Will extract the N valid data sequences:
nb_seq=0
idx_seq_start = [] ; # index of first valid point of the sequence
idx_seq_stop  = [] ; # index of last  valid point of the sequence

jr=0
while jr < nbr:
    # Ignoring masked values and zeros...        
    if (not vmask[jr]) and (vmodel[jr]!=0.0):
        nb_seq = nb_seq + 1
        print '\n --- found seq #'+str(nb_seq)+' !'
        np_s = 1
        idx_seq_start.append(jr)
        print ' => starting at jt='+str(jr)
        while (not vmask[jr+1]) and (vmodel[jr+1]!=0.0) :
            jr = jr+1
            np_s = np_s+1
        idx_seq_stop.append(jr)
        print ' => and stoping at jt='+str(jr)
        
    jr = jr+1


if len(idx_seq_start) != nb_seq: print ' ERROR #1!'; sys.exit(1)

print '\n idx_seq_start =', idx_seq_start
print '\n idx_seq_stop =', idx_seq_stop

for js in range(nb_seq):
    print '\n\n ###################################'
    print '  *** Seq #'+str(js+1)+':'
    it1 = idx_seq_start[js]
    it2 = idx_seq_stop[js]
    #print vmodel[it1:it2+1]
    nbp = it2-it1+1

    print ' *** Considering '+str(nbp)+' points!\n'

    # Create Matplotlib time array:
    vtime = nmp.zeros(nbp)
    for jt in range(nbp): vtime[jt] = mdates.epoch2num(vt_epoch[it1+jt])




    ii=nbp/300
    ib=max(ii-ii%10,1)
    print ' ii , ib =', ii, ib

    xticks_d=30.*ib

    fig = plt.figure(num = 1, figsize=(12,7), facecolor='w', edgecolor='k')
    ax1 = plt.axes([0.07, 0.24, 0.9, 0.75])

    ax1.set_xticks(vtime[::xticks_d])
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
    plt.xticks(rotation='60')


    plt.plot(vtime, vephem[it1:it2+1], '-', color=color_dark_blue, linewidth=2, label='Satellite ("'+cv_eph+'")', zorder=10)
    plt.plot(vtime, vmodel[it1:it2+1], '-', color=b_org, linewidth=2,  label='Model ("'+cv_mdl+'")', zorder=15)
    ax1.set_ylim(-0.68,0.68) ; ax1.set_xlim(vtime[0],vtime[nbp-1])
    plt.xlabel('Time [seconds since 1970]')
    plt.ylabel('SSH [m]')
    #cstep = '%5.5i'%(jpnij)
    ax1.grid(color='k', linestyle='-', linewidth=0.3)
    
    plt.legend(bbox_to_anchor=(0.55, 0.98), ncol=1, shadow=True, fancybox=True)


    #ax2 = ax1.twinx()
    #ax2.set_ylim(25.,68.)
    #plt.plot(vtime, vlat, '--', color='0.4', linewidth=1.5, label='latitude', zorder=0.1)
    #ax2.set_ylabel(r'Latitude [$^\circ$North]', color='0.4')
    #[t.set_color('0.4') for t in ax2.yaxis.get_ticklabels()]

    plt.savefig('fig_seq_'+str(js+1)+'.'+fig_ext, dpi=120, transparent=True)
    plt.close(1)




