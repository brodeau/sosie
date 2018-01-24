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



#fig_ext='png'
fig_ext='svg'

#narg = len(sys.argv)
#if narg < 4: print 'Usage: '+sys.argv[0]+' <file> <variable> <# snapshot>'; sys.exit(0)
#cf_in = sys.argv[1] ; cv_in=sys.argv[2] ; jt=int(sys.argv[3])-1

cf_natl  = 'result.nc'
cf_saral = 'data_ephem.nc'
cv_in = 'sossheig'
cf_lat = 'lat_ephem.nc'
    
bt.chck4f(cf_natl)
bt.chck4f(cf_saral)


id_saral = Dataset(cf_saral)
vtime  =  id_saral.variables['time'][:]
vsaral =  id_saral.variables[cv_in][:]
id_saral.close()

id_lat = Dataset(cf_lat)
vlat  =  id_lat.variables['latitude'][:]
id_lat.close()



id_natl = Dataset(cf_natl)
vnatl = id_natl.variables[cv_in][:]
id_natl.close()
print "  => READ!"

vtime = vtime/60.

font_corr = 1.8


#params = { 'font.family':'Helvetica Neue',
#           'font.weight':    'light',
params = { 'font.family':'Ubuntu',
           'font.weight':    'normal',
           'font.size':       int(12*font_corr),
           'legend.fontsize': int(12*font_corr),
           'xtick.labelsize': int(12*font_corr),
           'ytick.labelsize': int(12*font_corr),
           'axes.labelsize':  int(12*font_corr),
           'axes.titlesize' : 20,
           'lines.linewidth' : 8,
           'lines.markersize' : 10 }
#           'figure.facecolor': 'w' }
mpl.rcParams.update(params)


fig = plt.figure(num = 1, figsize=(12,6), facecolor='w', edgecolor='k')
ax1 = plt.axes([0.105, 0.13, 0.81, 0.85])
plt.plot(vtime, vsaral, '-', color=color_dark_blue, linewidth=2, label='SARAL', zorder=10)
plt.plot(vtime, vnatl, '-', color=b_org, linewidth=2,  label='NATL60', zorder=15)
ax1.set_ylim(-0.68,0.68) ; ax1.set_xlim(0.,13.2)
plt.xlabel('Time [minutes]')
plt.ylabel('SSH [m]')
#cstep = '%5.5i'%(jpnij)
ax1.grid(color='k', linestyle='-', linewidth=0.3)

plt.legend(bbox_to_anchor=(0.55, 0.98), ncol=1, shadow=True, fancybox=True)


ax2 = ax1.twinx()
ax2.set_ylim(25.,68.) ; ax2.set_xlim(0.,13.2)
plt.plot(vtime, vlat, '--', color='0.4', linewidth=1.5, label='latitude', zorder=0.1)
ax2.set_ylabel(r'Latitude [$^\circ$North]', color='0.4')
[t.set_color('0.4') for t in ax2.yaxis.get_ticklabels()]

plt.savefig('fig.'+fig_ext, dpi=120, transparent=True)
plt.close(1)
