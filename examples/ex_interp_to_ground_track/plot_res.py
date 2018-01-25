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
id_in.close()
print "  => READ!"


nbr = len(vt_epoch)

if jt2 == 0: jt2 = nbr-1

if jt1 >= nbr or jt2 >= nbr:
    print 'ERROR: the file contains only '+str(nbr)+' time records!' ; sys.exit(0)
if jt1 >= jt2:
    print 'ERROR: jt2 must be > jt1!' ; sys.exit(0)




nbp = jt2-jt1+1

print ' *** Considering '+str(nbp)+' points!\n'

# Create Matplotlib time array:
vtime = nmp.zeros(nbp)
for jt in range(nbp): vtime[jt] = mdates.epoch2num(vt_epoch[jt1+jt])



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

ii=nbp/300
ib=max(ii-ii%10,1)
print ' ii , ib =', ii, ib

xticks_d=30.*ib

fig = plt.figure(num = 1, figsize=(12,7), facecolor='w', edgecolor='k')
ax1 = plt.axes([0.07, 0.24, 0.9, 0.75])

ax1.set_xticks(vtime[::xticks_d])
ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
plt.xticks(rotation='60')


plt.plot(vtime, vephem[jt1:jt2+1], '-', color=color_dark_blue, linewidth=2, label='Satellite ("'+cv_eph+'")', zorder=10)
plt.plot(vtime, vmodel[jt1:jt2+1], '-', color=b_org, linewidth=2,  label='Model ("'+cv_mdl+'")', zorder=15)
ax1.set_ylim(-0.68,0.68) ; ax1.set_xlim(vtime[0],vtime[nbp-1])
plt.xlabel('Time [seconds since 1970]')
plt.ylabel('SSH [m]')
#cstep = '%5.5i'%(jpnij)
ax1.grid(color='k', linestyle='-', linewidth=0.3)

plt.legend(bbox_to_anchor=(0.55, 0.98), ncol=1, shadow=True, fancybox=True)


#ax2 = ax1.twinx()
#ax2.set_ylim(25.,68.) ; ax2.set_xlim(0.,13.2)
#plt.plot(vtime, vlat, '--', color='0.4', linewidth=1.5, label='latitude', zorder=0.1)
#ax2.set_ylabel(r'Latitude [$^\circ$North]', color='0.4')
#[t.set_color('0.4') for t in ax2.yaxis.get_ticklabels()]

plt.savefig('fig.'+fig_ext, dpi=120, transparent=True)
plt.close(1)
