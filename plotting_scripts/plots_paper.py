import fnmatch
import os
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib as mpl
import matplotlib.patches as mpatches
from matplotlib.offsetbox import AnchoredText
import numpy as np
import xarray as xr
from scipy.interpolate import Rbf
from scipy.ndimage import gaussian_filter
from matplotlib.patches import Rectangle
from matplotlib.legend_handler import HandlerBase
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

out_path = '/<your>/<output>/<path>/'

size = 16
font = {'family' : 'normal', 'weight' : 'normal', 'size' : size}
plt.rc('font', **font)


class HandlerColormap(HandlerBase):
    def __init__(self, cmap, num_stripes=8, **kw):
        HandlerBase.__init__(self, **kw)
        self.cmap = cmap
        self.num_stripes = num_stripes
    def create_artists(self, legend, orig_handle, 
                       xdescent, ydescent, width, height, fontsize, trans):
        stripes = []
        for i in range(self.num_stripes):
            s = Rectangle([xdescent + i * width / self.num_stripes, ydescent], 
                          width / self.num_stripes, 
                          height, 
                          fc=self.cmap((2 * i + 1) / (2 * self.num_stripes)), 
                          transform=trans)
            stripes.append(s)
        return stripes

cmaps = [plt.cm.Reds,plt.cm.Greens,plt.cm.Blues,plt.cm.Oranges,plt.cm.Greys]
cmap_handles = [Rectangle((0, 0), 1, 1) for _ in cmaps]
handler_map = dict(zip(cmap_handles, 
                       [HandlerColormap(cm, num_stripes=20) for cm in cmaps]))


#%%
###############################################################################
###CYCLONE TRACKS##############################################################
###############################################################################
in_path = '<your>/<input>/<path>/'

files = sorted(fnmatch.filter(os.listdir(in_path), '*.csv'))


track = []
for i in files:
    track.append(pd.read_csv(f'{in_path}{i}',sep=';'))

lstyle = ['-',
          '-',
          '-','--','-.',':',
          '-','--','-.',
          '-','--',
          '-','--']
col    = ['purple',
          '#075c07',
          '#e41a1c','#e41a1c','#e41a1c','#e41a1c',\
          '#4daf4a','#4daf4a','#4daf4a',\
          '#ffa100','#ffa100',\
          'b','b'] 
    
fig = plt.figure(figsize=(13,7.5))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-17.5,20,38.5,60])
ax.coastlines(resolution='50m',edgecolor='grey')
ax.gridlines()
ax.add_feature(cfeature.BORDERS, edgecolor='grey')
for i in range(len(files)):
    storm = files[i].split('_')
    plt.plot(track[i].lon,track[i].lat,label=storm[1][:-4],
             c=col[i],linestyle=lstyle[i],linewidth=2.5)
    plt.scatter(track[i].lon,track[i].lat,
                s=(abs(track[i].mslp-1013)/10)**3,marker='o',c=col[i])
l = plt.legend(ncol=2,loc='lower left')
ax.add_artist(l)

trick = plt.scatter(np.linspace(5,10,6),np.linspace(45,50,6),
                    s=(abs(np.linspace(950,1000,6)-1013)/10)**3,alpha=0)
handles, labels = trick.legend_elements(prop="sizes", alpha=0.5)
l2 = ax.legend(handles, ['1000','990','980','970','960','950'], 
               title = r'$p$ in hPa', loc="lower right")

plt.savefig(f'{out_path}fig01.pdf',bbox_inches="tight")
plt.close()

#%%
###############################################################################
###BURGLIND - OBS 6Z###########################################################
###############################################################################
in_path = '<your>/<input>/<path>/'
cf = xr.open_dataset(f"{in_path}Burglind_Jan18.nc").to_dataframe().reset_index()
cf = cf.loc[cf.daytime == 5]

xi,yi = np.meshgrid(np.arange(-5,20,.75),np.arange(42.5,57.5,.75)) #.33409395))
rbf = Rbf(cf.dropna(subset=['vgl']).lonp, cf.dropna(subset=['vgl']).latp, 
          cf.dropna(subset=['vgl']).vgl, function='inverse',epsilon=.1)
zi = rbf(xi,yi)

sig = 1 #gaussian_filter

lab = xr.open_dataset(f'{in_path}/labels/20180103_05.nc',
                      drop_variables=['wind','dwddt','wd','th','dthdt',
                                      'mslp','alt','dpdt','prec'])
lab.lab.values[lab.lab.values==5] = 4
lab.lab.values[lab.vgl < .8] = np.nan

lab = lab.to_dataframe()
lab = lab.loc[lab.vgl >= .8]

track = pd.read_csv('{in_path}/20180103_Burglind_Jan18.csv',sep=';')
track = track.loc[track.date == 201801030600]

my_cmap = mpl.colors.ListedColormap(['gray','red','green','blue','gold'])
t = ['NF','WJ','CFC','CJ','CS']

fig = plt.subplots(nrows=5,ncols=1,figsize=(10,20 ))
plt.subplots_adjust(wspace=.05,hspace=.05)
ax = plt.subplot(411,projection=ccrs.PlateCarree())
ax.set_extent([-5,20,42.5,57.5])
ax.coastlines(resolution='50m',zorder=4)
ax.gridlines(zorder=4)
ax.add_feature(cfeature.BORDERS, edgecolor='grey')
ax.add_feature(cfeature.OCEAN, color = 'white', edgecolor='black', zorder=3)
c = plt.scatter(cf.lonp,cf.latp,c=cf.dpdt,cmap='bwr')
plt.contour(xi,yi,gaussian_filter(zi,sigma=sig),levels=[.8,1,1.2],colors='k',linestyles=[':','--','-'])
plt.scatter(track.lon,track.lat,s=200,c='r', marker='x',
                linewidth=5,zorder=4)
c.set_clim([-3,3])        
cbar = plt.colorbar(c, extend='both')
cbar.set_label(r'$\Delta p$ in hPa$\,$h$^{-1}$')
at = AnchoredText("(a)", prop=dict(size=size), frameon=True, loc='upper left')
ax.add_artist(at)

ax = plt.subplot(412,projection=ccrs.PlateCarree())
ax.set_extent([-5,20,42.5,57.5])
ax.coastlines(resolution='50m',zorder=4)
ax.gridlines(zorder=4)
ax.add_feature(cfeature.BORDERS, edgecolor='grey')
ax.add_feature(cfeature.OCEAN, color = 'white', edgecolor='black', zorder=3)
c = plt.scatter(cf.lonp,cf.latp,c=cf.prec,cmap='Blues')
plt.contour(xi,yi,gaussian_filter(zi,sigma=sig),levels=[.8,1,1.2],colors='k',linestyles=[':','--','-'])
plt.scatter(track.lon,track.lat,s=200,c='r', marker='x',
                linewidth=5,zorder=4)
c.set_clim([0,8])        
cbar = plt.colorbar(c, extend='both')
cbar.set_label(r'$RR$ in mm$\,$h$^{-1}$')
at = AnchoredText("(b)", prop=dict(size=size), frameon=True, loc='upper left')
ax.add_artist(at)

ax = plt.subplot(413,projection=ccrs.PlateCarree())
ax.set_extent([-5,20,42.5,57.5])
ax.coastlines(resolution='50m',zorder=4)
ax.gridlines(zorder=4)
ax.add_feature(cfeature.BORDERS, edgecolor='grey')
ax.add_feature(cfeature.OCEAN, color = 'white', edgecolor='black', zorder=3)
c = plt.scatter(cf.lonp,cf.latp,c=cf.dthdt*10**2,cmap='bwr')
plt.contour(xi,yi,gaussian_filter(zi,sigma=sig),levels=[.8,1,1.2],colors='k',linestyles=[':','--','-'])
plt.scatter(track.lon,track.lat,s=200,c='r', marker='x',
                linewidth=5,zorder=4)
c.set_clim([-0.5,.5])        
cbar = plt.colorbar(c, extend='both')
cbar.set_label(r'$\Delta \tilde{\theta}$ in $10^2\times$h$^{-1}$')
at = AnchoredText("(c)", prop=dict(size=size), frameon=True, loc='upper left')
ax.add_artist(at)

ax = plt.subplot(414,projection=ccrs.PlateCarree())
ax.set_extent([-5,20,42.5,57.5])
ax.coastlines(resolution='50m',zorder=4)
ax.gridlines(zorder=4)
ax.add_feature(cfeature.BORDERS, edgecolor='grey')
ax.add_feature(cfeature.OCEAN, color = 'white', edgecolor='black', zorder=3)
plt.contour(xi,yi,gaussian_filter(zi,sigma=sig),levels=[.8,1,1.2],colors='k',linestyles=[':','--','-'])
plt.scatter(track.lon,track.lat,s=200,c='r', marker='x',
                linewidth=5,zorder=4)
c = plt.scatter(lab.lon,lab.lat,c=lab.lab,cmap=my_cmap)     
c.set_clim([-.5,4.5])
cbar = plt.colorbar(c,ticks=np.arange(0,4.5,1))
cbar.set_ticklabels(t)
at = AnchoredText("(d)", prop=dict(size=size), frameon=True, loc='upper left')
ax.add_artist(at)

plt.savefig(f'{out_path}fig03.pdf',bbox_inches="tight")
plt.close()

#%%
###############################################################################
###BURGLIND - KRIGING##########################################################
###############################################################################
in_path = '<your>/<input>/<path>/'

files = ['univ_normal_Burglind_2_full_preds.csv',
         'univ_normal_Burglind_5_full_preds.csv',
         'univ_normal_Burglind_8_full_preds.csv',
         'univ_normal_Burglind_11_full_preds.csv',
         'univ_normal_Burglind_14_full_preds.csv',
         'univ_normal_Burglind_17_full_preds.csv']

cmaps = [plt.cm.Reds,plt.cm.Greens,plt.cm.Blues,plt.cm.Oranges,plt.cm.Greys]

sb = [ '(a)','(b)','(c)','(d)','(e)','(f)']

track = pd.read_csv('{in_path}20180103_Burglind.csv',sep=';')
track = track.loc[(track.date == 201801030300) | (track.date == 201801030600) | 
                  (track.date == 201801030900) | (track.date == 201801031200) | 
                  (track.date == 201801031500) | (track.date == 201801031800) ]

xi,yi = np.meshgrid(np.arange(-5,20,.75),np.arange(42.5,57.5,.75))
sig = 1 #gaussian_filter

fig = plt.subplots(nrows=3,ncols=2,figsize=(20,20))
plt.subplots_adjust(wspace=.05,hspace=.05)
j=0
for i,k in zip(files,sb):
    j+=1
    df = pd.read_csv(f'{in_path}{i}', sep=',')
    #set neg probabilities to 0 and >1 to 1
    df.loc[:,'P0':'P5'] = df.loc[:,'P0':'P5'].where(df.loc[:,'P0':'P5']>=0,0)
    df.loc[:,'P0':'P5'] = df.loc[:,'P0':'P5'].where(df.loc[:,'P0':'P5']<=1,1)
    df = df.rename(columns={'label':'labs'})
    x = i.split('_')
    if i == files[0]:
        storm = x[2]
        date = '20180103'
    time = str(int(x[3])+1).zfill(2)
        
    ds = xr.open_dataset(f'{mntpnt}wo4u6/{date}_{str(int(x[3])).zfill(2)}.nc', 
                         drop_variables=('dpdt','dthdt','dwddt','th','mslp','wd','prec','lab')).to_dataframe()
    ds = ds.dropna(subset=['vgl'])
    rbf = Rbf(ds.lon, ds.lat, ds.vgl, function='inverse',epsilon=.1)
    zi = rbf(xi,yi)

    c = [df.labs.values,df.P0.values,df.P1.values,
         df.P2.values,df.P3.values,df.P5.values]
    
    ax = plt.subplot(3,2,j,projection=ccrs.PlateCarree())
    ax.set_extent([-5,20,42.5,57.5])
    ax.coastlines(resolution='50m',zorder=4)
    ax.gridlines(zorder=4)
    ax.add_feature(cfeature.BORDERS, edgecolor='grey')
    ax.add_feature(cfeature.OCEAN, color = 'white', edgecolor='black', zorder=3)
    if df.loc[df.labs==0].lon.any():
        cb0=plt.scatter(df.loc[df.labs==0].lon.values,df.loc[df.labs==0].lat.values,
                       c=df.loc[df.labs==0].P0.values,
                       cmap='Greys',alpha=df.loc[df.labs==0].P0.values,
                       marker='s', edgecolor = None)
        cb0.set_clim(0,1)
    if df.loc[df.labs==1].lon.any():
        cb=plt.scatter(df.loc[df.labs==1].lon.values,df.loc[df.labs==1].lat.values,
                       c=df.loc[df.labs==1].P1.values,
                       cmap='Reds',alpha=df.loc[df.labs==1].P1.values,
                       marker='s', edgecolor = None)
        cb.set_clim(0,1)
    if df.loc[df.labs==2].lon.any():
        cb=plt.scatter(df.loc[df.labs==2].lon.values,df.loc[df.labs==2].lat.values,
                       c=df.loc[df.labs==2].P2.values,
                       cmap='Greens',alpha=df.loc[df.labs==2].P2.values,
                       marker='s', edgecolor = None)
        cb.set_clim(0,1)
    if df.loc[df.labs==3].lon.any():
        cb=plt.scatter(df.loc[df.labs==3].lon.values,df.loc[df.labs==3].lat.values,
                       c=df.loc[df.labs==3].P3.values,
                       cmap='Blues',alpha=df.loc[df.labs==3].P3.values,
                       marker='s', edgecolor = None)
        cb.set_clim(0,1)
    if df.loc[df.labs==5].lon.any():
        cb=plt.scatter(df.loc[df.labs==5].lon.values,df.loc[df.labs==5].lat.values,
                       c=df.loc[df.labs==5].P5.values,
                       cmap='Oranges',alpha=df.loc[df.labs==5].P5.values,
                       marker='s',edgecolor = None)
        cb.set_clim(0,1)
    plt.contour(xi,yi,gaussian_filter(zi,sigma=sig),levels=[1,1.2],colors='k',
                linestyles=['--','-'])
    plt.scatter(track.iloc[j-1].lon,track.iloc[j-1].lat,s=200,c='r',marker='x',
                linewidth=5, zorder=4)
    if i == files[0]:
        plt.legend(handles=cmap_handles,labels=['WJ','CFC','CJ','CS','NF'],
                    handler_map=handler_map, loc='lower right',
                    handlelength = 3,title = '',
                    bbox_to_anchor=(0, 0.08, 1, 1),
                    bbox_transform=ax.transAxes, edgecolor='w') #r'$f \in [0,1]$')
        p=mpatches.Rectangle((14.7, 42.51), 5, 1.8,alpha=.9,color='w',zorder=4)
        ax.add_patch(p)
        
        cbaxes = inset_axes(ax, width="9%", height="3.8%", loc=4,
                           bbox_to_anchor=(.185, -.435+0.0875, 1, .8),
                           bbox_transform=ax.transAxes,borderpad=10)
        cb0_ = plt.colorbar(cb0,cax=cbaxes,ticks=[0,1], label=r'$f$ in %',
                            orientation='horizontal',
                           shrink=.1,aspect=5,drawedges=False)
        cb0_.set_ticklabels([0,100])
        cb0_.outline.set_color(None)
    
    at = AnchoredText(
    f"{k} {date}_{time}", prop=dict(size=size), frameon=True, loc='upper left')
    ax.add_artist(at)

plt.savefig(f'{out_path}fig04.pdf',bbox_inches="tight")
plt.close()


#%%
###############################################################################
###############################################################################
###############################################################################
in_path = '<your>/<input>/<path>/'

files = ['univ_normal_Bennet_11_full_preds.csv',
         'univ_normal_Herwart_8_full_preds.csv',
         'univ_normal_Friederike_11_full_preds.csv',
         'univ_normal_Xavier_11_full_preds.csv']

sig = 1 #gaussian_filter

cmaps = [plt.cm.Reds,plt.cm.Greens,plt.cm.Blues,plt.cm.Oranges,plt.cm.Greys]

track = pd.read_csv('{in_path}20190304_Bennet_Mar19.csv',sep=';')
track = track.loc[(track.date == 201903041200)]
dummy = pd.read_csv('{in_path}20171029_Herwart.csv',sep=';')
track = track.append(dummy.loc[dummy.date == 201710290900])
dummy = pd.read_csv('{in_path}20180118_Friederike_Jan18.csv',sep=';')
track = track.append(dummy.loc[dummy.date == 201801181200]).reset_index()
dummy = pd.read_csv('{in_path}20171005_Xavier_Oct17.csv',sep=';')
track = track.append(dummy.loc[dummy.date == 201710051200]).reset_index()

xi,yi = np.meshgrid(np.arange(-5,20,.75),np.arange(42.5,57.5,.75))
    
sb = ['(a)','(b)','(c)','(d)']

fig = plt.subplots(nrows=2,ncols=2,figsize=(20,13))
plt.subplots_adjust(wspace=.05,hspace=.05)
j=0
for i,k in zip(files,sb):
    j+=1
    df = pd.read_csv(f'{in_path}{i}', sep=',')
    #set neg probabilities to 0 and >1 to 1
    df.loc[:,'P0':'P5'] = df.loc[:,'P0':'P5'].where(df.loc[:,'P0':'P5']>=0,0)
    df.loc[:,'P0':'P5'] = df.loc[:,'P0':'P5'].where(df.loc[:,'P0':'P5']<=1,1)
    df = df.rename(columns={'label':'labs'})
    x = i.split('_')
    date = str(track.iloc[j-1].date)[:8]
    time = str(int(x[3])+1).zfill(2)
        
    ds = xr.open_dataset(f'{mntpnt}wo4u6/{date}_{str(int(x[3])).zfill(2)}.nc', 
                         drop_variables=('dpdt','dthdt','dwddt','th','mslp','wd','prec','lab')).to_dataframe()
    ds = ds.dropna(subset=['vgl'])
    rbf = Rbf(ds.lon, ds.lat, ds.vgl, function='inverse',epsilon=.1)
    zi = rbf(xi,yi)

    c = [df.labs.values,df.P0.values,df.P1.values,
         df.P2.values,df.P3.values,df.P5.values]
    
    ax = plt.subplot(2,2,j,projection=ccrs.PlateCarree())
    ax.set_extent([-5,20,42.5,57.5])
    ax.coastlines(resolution='50m',zorder=4)
    ax.gridlines(zorder=4)
    ax.add_feature(cfeature.BORDERS, edgecolor='grey')
    ax.add_feature(cfeature.OCEAN, color = 'white', edgecolor='black', zorder=3)
    if df.loc[df.labs==0].lon.any():
        cb0=plt.scatter(df.loc[df.labs==0].lon.values,df.loc[df.labs==0].lat.values,
                       c=df.loc[df.labs==0].P0.values,
                       cmap='Greys',alpha=df.loc[df.labs==0].P0.values,
                       marker='s', edgecolor = None)
        cb0.set_clim(0,1)
    if df.loc[df.labs==1].lon.any():
        cb=plt.scatter(df.loc[df.labs==1].lon.values,df.loc[df.labs==1].lat.values,
                       c=df.loc[df.labs==1].P1.values,
                       cmap='Reds',alpha=df.loc[df.labs==1].P1.values,
                       marker='s', edgecolor = None)
        cb.set_clim(0,1)
    if df.loc[df.labs==2].lon.any():
        cb=plt.scatter(df.loc[df.labs==2].lon.values,df.loc[df.labs==2].lat.values,
                       c=df.loc[df.labs==2].P2.values,
                       cmap='Greens',alpha=df.loc[df.labs==2].P2.values,
                       marker='s', edgecolor = None)
        cb.set_clim(0,1)
    if df.loc[df.labs==3].lon.any():
        cb=plt.scatter(df.loc[df.labs==3].lon.values,df.loc[df.labs==3].lat.values,
                       c=df.loc[df.labs==3].P3.values,
                       cmap='Blues',alpha=df.loc[df.labs==3].P3.values,
                       marker='s', edgecolor = None)
        cb.set_clim(0,1)
    if df.loc[df.labs==5].lon.any():
        cb=plt.scatter(df.loc[df.labs==5].lon.values,df.loc[df.labs==5].lat.values,
                       c=df.loc[df.labs==5].P5.values,
                       cmap='Oranges',alpha=df.loc[df.labs==5].P5.values,
                       marker='s',edgecolor = None)
        cb.set_clim(0,1)
    plt.contour(xi,yi,gaussian_filter(zi,sigma=sig),levels=[1,1.2],colors='k',
                linestyles=['--','-'])
    plt.scatter(track.iloc[j-1].lon,track.iloc[j-1].lat,s=200,c='r',marker='x',
                linewidth=5,zorder=4)
    if i == files[0]:
        plt.legend(handles=cmap_handles,labels=['WJ','CFC','CJ','CS','NF'],
                    handler_map=handler_map, loc='lower right',
                    handlelength = 3,title = '',
                    bbox_to_anchor=(0, 0.08, 1, 1),
                    bbox_transform=ax.transAxes, edgecolor='w')
        p=mpatches.Rectangle((14.7, 42.51), 5, 1.8,alpha=.9,color='w',zorder=4)
        ax.add_patch(p)
        
        cbaxes = inset_axes(ax, width="9%", height="3.8%", loc=4,
                           bbox_to_anchor=(.185, -.435+0.0875, 1, .8),
                           bbox_transform=ax.transAxes,borderpad=10)
        cb0_ = plt.colorbar(cb0,cax=cbaxes,ticks=[0,1], label=r'$f$ in %',
                            orientation='horizontal',
                           shrink=.1,aspect=5,drawedges=False)
        cb0_.set_ticklabels([0,100])
        cb0_.outline.set_color(None)
    
    at = AnchoredText(f"{k} {date}_{time}", prop=dict(size=size), frameon=True, 
                      loc='upper left')
    ax.add_artist(at)

plt.savefig(f'{out_path}fig11.pdf',bbox_inches="tight")
plt.close()

#%%
###############################################################################
#CREA##########################################################################
###############################################################################
in_path = '<your>/<input>/<path>/'
cf = xr.open_dataset(f"{in_path}20180103.nc") 
sig = 2 

df = pd.read_csv(f'{in_path}cosmo-rea_preds_Burglind.csv', sep=',')
df = df.loc[df.time == 5]
df.loc[:,'p0':'p5'] = df.loc[:,'p0':'p5'].where(df.loc[:,'p0':'p5']>=0,0)
df.loc[:,'p0':'p5'] = df.loc[:,'p0':'p5'].where(df.loc[:,'p0':'p5']<=1,1)
storm = 'Burglind'; date = '20180103'

df = df.rename(columns={'p0':'P0','p1':'P1','p2':'P2','p3':'P3','p5':'P5'})
df['labs'] = df.loc[:,'P0':'P5'].idxmax(axis=1)

df.labs = df.labs.where(df.labs != 'P0',0)
df.labs = df.labs.where(df.labs != 'P1',1)
df.labs = df.labs.where(df.labs != 'P2',2)
df.labs = df.labs.where(df.labs != 'P3',3)
df.labs = df.labs.where(df.labs != 'P5',4)

track = pd.read_csv('{in_path}20180103_Burglind_Jan18.csv',sep=';')
track = track.loc[track.date == 201801030600]

cmaps = [plt.cm.Reds,plt.cm.Greens,plt.cm.Blues,plt.cm.Oranges,plt.cm.Greys]
    
fig = plt.figure(figsize=(10,6))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-5,20,42.5,57.5])
ax.coastlines(resolution='50m',zorder=4)
ax.gridlines(zorder=4)
ax.add_feature(cfeature.BORDERS, edgecolor='grey')
ax.add_feature(cfeature.OCEAN, color = 'white', edgecolor='black', zorder=3)
if df.loc[df.labs==0].lon.any():
    cb0=plt.scatter(df.loc[df.labs==0].lon.values,df.loc[df.labs==0].lat.values,
                    c=df.loc[df.labs==0].P0.values,
                    cmap='Greys',alpha=df.loc[df.labs==0].P0.values,
                    marker='s', edgecolor = None)
    cb0.set_clim(0,1)
if df.loc[df.labs==1].lon.any():
    cb1=plt.scatter(df.loc[df.labs==1].lon.values,df.loc[df.labs==1].lat.values,
                    c=df.loc[df.labs==1].P1.values,
                    cmap='Reds',alpha=df.loc[df.labs==1].P1.values,
                    marker='s', edgecolor = None)
    cb1.set_clim(0,1)
if df.loc[df.labs==2].lon.any():
    cb2=plt.scatter(df.loc[df.labs==2].lon.values,df.loc[df.labs==2].lat.values,
                    c=df.loc[df.labs==2].P2.values,
                    cmap='Greens',alpha=df.loc[df.labs==2].P2.values,
                    marker='s', edgecolor = None)
    cb2.set_clim(0,1)
if df.loc[df.labs==3].lon.any():
    cb3=plt.scatter(df.loc[df.labs==3].lon.values,df.loc[df.labs==3].lat.values,
                    c=df.loc[df.labs==3].P3.values,
                    cmap='Blues',alpha=df.loc[df.labs==3].P3.values,
                    marker='s', edgecolor = None)
    cb3.set_clim(0,1)
if df.loc[df.labs==4].lon.any():
    cb4=plt.scatter(df.loc[df.labs==4].lon.values,df.loc[df.labs==4].lat.values,
                    c=df.loc[df.labs==4].P5.values,
                    cmap='Oranges',alpha=df.loc[df.labs==4].P5.values,
                    marker='s',edgecolor = None)
    cb4.set_clim(0,1)
plt.scatter(track.lon,track.lat,s=200,c='r', marker='x',linewidth=5,zorder=4)
plt.legend(handles=cmap_handles,labels=['WJ','CFC','CJ','CS','NF'],
            handler_map=handler_map, loc='lower right',
            handlelength = 3,title = '',
            bbox_to_anchor=(0, 0.08, 1, 1),
            bbox_transform=ax.transAxes, edgecolor='w')

p=mpatches.Rectangle((14.7, 42.5), 5, 1.8,alpha=.9,color='w',zorder=4)
ax.add_patch(p)

cbaxes = inset_axes(ax, width="9%", height="4%", loc=4,
                   bbox_to_anchor=(.19, -.435+0.08, 1, .8),
                   bbox_transform=ax.transAxes,borderpad=10)
cb0_ = plt.colorbar(cb0,cax=cbaxes,ticks=[0,1], label=r'$f$ in %',orientation='horizontal',
                   shrink=.1,aspect=5,drawedges=False)
cb0_.set_ticklabels([0,100])
cb0_.outline.set_color(None)

at = AnchoredText(f"{date}_6", prop=dict(size=size), frameon=True, loc='upper left')
ax.add_artist(at)

plt.savefig(f'{out_path}fig05.pdf',bbox_inches="tight")
plt.close()
