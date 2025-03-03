#!/home/ahutko/anaconda3/bin/python

import glob
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from cartopy import crs
from cartopy import feature
import geopy
import obspy
from obspy import read
from obspy.geodetics.base import gps2dist_azimuth
import datetime
from datetime import timedelta
import matplotlib.dates as md

vel = 5
km_miles = 1.0

evlat = 48.607
evlon = -122.804
otime = 120

#----- make dictionary of station Lat/Lon/startdate given a sta code
f = open('channels_squacids_west_coast')
lines = f.readlines()
f.close()
sta_dict = {}
for line in lines:
    if ( "#" not in line ):
        try:
            sta = line.split('.')[1]
            lat = float(line.split()[4])
            lon = float(line.split()[5])
            sta_dict[sta] = [lat,lon]
        except:
            pass

#----- Read in SAC files, for 6 ch sites, select only the HHZ
#      This list will need to be hand culled, but UW.A* - UW.C* is just a start.

saclist = glob.glob("UW.[A-C]*.SAC")
sncls, stations, sacfiles = [],[],[]
for sacfile in saclist:
     sta = sacfile.split('.')[1]
     cha = sacfile.split('.')[3]
     sncl = sacfile.split('.M.2025')[0]
     print(sta,cha, sncl)
     if sta in stations and cha == 'HHZ':
          sncls[:] = [x for x in sncls if sta not in x]
          sncls.append(sncl)
          stations.append(sta)
          sacfiles[:] = [x for x in sacfiles if sta not in x]
          sacfiles.append(sacfile)
     else:
          sncls.append(sncl)
          sacfiles.append(sacfile)

for sacfile in sacfiles:
    print(sacfile)

#----- Read in SAC files
stalist = []
data, srdict, distdict = {},{},{}
for sacfile in sacfiles:
    sta = sacfile.split('.')[1]
    stalist.append(sta)
    slat = sta_dict[sta][0]
    slon = sta_dict[sta][1]
    gps_1_to_2 = gps2dist_azimuth( evlat,evlon,slat,slon)
    distance_in_km = int(10.*(gps_1_to_2[0] / 1000.)) / 10.
    distdict[sta] = distance_in_km / km_miles
    st = read(sacfile)
    st.resample(100.0)
    tr = st[0]
    sr = 1.*tr.stats.sampling_rate
    srdict[sta] = sr
    it1 = int((otime-10)*sr)  # OT minus 10sec
    it2 = int((otime+300)*sr) # OT plus 5min, maybe make this 2min?
    tr.filter("bandpass", freqmin = 0.5, freqmax=20, corners=2,zerophase = True)  # Play around with this
    tr.detrend(type='demean')
    ampmax = max(abs(tr.data[it1:it2]))
    data[sta] = tr.data[it1:it2]/ampmax

#------ plot wiggles by distance
fig = plt.figure(figsize=(5,5),dpi=250)
gs1 = gridspec.GridSpec(1,1)
gs1.update(left=0.15,right=0.85,bottom=0.15,top=0.9,wspace=0.01)
ax = plt.subplot(gs1[:,:])

ymaxline = 100/km_miles
x2 = ymaxline/(vel*60.)
ax.plot([0,x2],[0,ymaxline],color='k',linewidth=0.15,alpha=0.25)
ax.plot([20,x2+20],[0,ymaxline],color='gray',linewidth=0.15,alpha=0.25)

index = 0
for sta in stalist:
    print(sta)
    index = distdict[sta]
    sr = srdict[sta]
    start_date = tr.stats.starttime
    datalen = len(data[sta])
#    x = [start_date + datetime.timedelta(seconds = idx/sr) for idx in range(datalen)]
    #x =  np.arange(20.5,48,1./(60*sr)) - otime
    x =  np.arange(otime-2,otime+28,1./(60*sr)) - otime #2024
    y = index + (1.899*data[sta])
    ax.plot(x,y,linewidth=0.1,color='r')
#    ax.text(28.6,index-(0.5/km_miles),sta,fontsize=4)


staxtext = 20
#ax.set_title('Seismic recordings from Sep 16, 2023 \n across the Puget Sound (filtered 40-50Hz)')
ax.set_title('Seismic recordings from Sep 28, 2024 \n across the Puget Sound (filtered 40-50Hz)')
ax.set_ylabel('Distance from ~2 km W of Browns Point, Tacoma (km)')
ax.set_xlabel('Minutes after 8:59pm')
ax.text(3.6,65/km_miles,'Booms begin',fontsize=6)
ax.text(23.5,65/km_miles,'Booms end',fontsize=6)
ax.text(23.5,63/km_miles,'20.0 min later',fontsize=6)
ax.text(1.6,25/km_miles,'(speed of sound in air = 340m/s)',fontsize=6)
#ax.text(17.8,65,'Booms end',fontsize=6,color='grey')
#ax.text(17.8,62,'20.0 min later',fontsize=6,color='grey')

#ax.text(staxtext,23,'Hansville',fontsize=6)
#ax.text(staxtext,45.5,'Port Townsend',fontsize=6)
#ax.text(staxtext,14,'Alki',fontsize=6)
#ax.text(staxtext,0.5,'North Bainbridge',fontsize=6)


ax.text(staxtext,70.5/km_miles,'Tolt Reservoir',fontsize=6)
ax.text(staxtext,58.5/km_miles,'Carnation',fontsize=6)
ax.text(staxtext,47.6/km_miles,'Near Bangor',fontsize=6)
ax.text(staxtext,41.5/km_miles,'N Bainbridge',fontsize=6)
ax.text(staxtext,27.5/km_miles,'Bremerton + W Seattle',fontsize=6)
ax.text(staxtext,21.5/km_miles,'Burley',fontsize=6)
ax.text(staxtext,10./km_miles,'Federal Way',fontsize=6)
ax.text(staxtext,8/km_miles,'University Place',fontsize=6)
ax.text(staxtext,3/km_miles,'Browns Point',fontsize=6)



ax.set_ylim([0,75/km_miles])  #,ymaxline])
ax.set_xlim([-2.5,30.5])

#logo
newax = fig.add_axes([0.75, 0.01, 0.2, 0.2], anchor='SE', zorder=-1)
im = plt.imread('PNSN_Small_RGB.png')
newax.imshow(im)
newax.axis('off')
plt.savefig("blah.png")



