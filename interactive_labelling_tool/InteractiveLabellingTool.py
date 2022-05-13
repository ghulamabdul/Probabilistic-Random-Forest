#!/usr/bin/env python
########################################################################################################################
# RAMEFI - INTERACTIVE LABELLING TOOL
# VISUALIZE, EXPLORE AND LABEL STATION OR GRIDDED DATA
# -----------------------------------------------------------------------------------------------------------------------
# Author:      Lea Eisenstein, lea.eisenstein@kit.edu
# Institute:   IMK-TRO, KIT, Germany
# Date:        May 2022
########################################################################################################################
'''
Use the ``bokeh serve`` command for the directory application by executing:
    bokeh serve InteractiveLabellingTool --show --port 4009 
at your command prompt. Then navigate to the URL
    http://localhost:4009/InteractiveLabellingTool
in your browser.
'''
#%%
### LOAD PACKAGES ######################################################################################################
import xarray as xr
import numpy as np
import pandas as pd
import copy
import os
from bokeh.io import curdoc
from bokeh.layouts import row, column, gridplot, Spacer
from bokeh.models import ColumnDataSource, CustomJS, ColorBar, LinearColorMapper, BasicTicker, Panel, Tabs, Select
from bokeh.models.widgets import Button, Slider, Div, DataTable, TableColumn
from bokeh.plotting import figure
from bokeh.palettes import Reds9, RdYlBu11, Blues9, Viridis11

#%%
### ADAPT ACCORDINGLY ###################################################################################
# Input path
in_path = '/<your>/<input>/<path>/' 
# Output path
out_path = '/<your>/<output>/<path>/'
in_path = '/home/lea/db/WP1/Paper/'
# Which storms do you want to read in?
# NOTE: The more storms or the higher the resolution, the slower the tool gets!
#Example:
stormlist = ['example_obs']
stormdates = ['20220513']

"""
Tools to interact with the plots. 
For more information, see https://docs.bokeh.org/en/latest/docs/user_guide/tools.html

    lasso_select: Define an arbitrary region for selection by left-dragging a 
                mouse of dragging a finger across the plot area.
    poly_select:  Define an arbitrary polygonal region for selection by left-
                clicking a mouse, or tapping a finger at different locations.
    tap:          Select single points by clicking the left mouse button, or 
                tapping with a finger
    box_zoom:     Define a rectangular region to zoom the plot bounds to by 
                left-dragging a mouse, or dragging a finger across the plot area.
    hover:        passive inspector tool
    reset:        Restores the plot ranges to their original values.
"""
tools = ['lasso_select, poly_select,tap,box_zoom, hover, reset']
tooltips = [("dp","@dp{0.00}"),("RR","@rr{0.00}")]

#%% Load data
select = Select(title='Storm',value=stormlist[0],options=stormlist)
stormdate = stormdates[0]

"""
Used variables and their names in our original netCDF file and here in the tool.
For more details see README.md.
name nc-file   name tool   description
vgl            vgl         normalised wind speed vgl = v/v98
prec           rr          precipitation amount
mslp           p           mean-sea level pressure
dpdt           dp          mean-sea level pressure tendency
th             th          normalised potential temperature
dthdt          dth         normalised pot. temperature tendency
wd             d           wind direction
dwddt          dd          wind direction tendency 
label          lab         set labels (0 in beginning)
"""

df = pd.DataFrame()
for i in range(len(stormlist)):
    ds = xr.open_mfdataset(f'{in_path}{stormlist[i]}.nc', combine='by_coords')
    
    df_ = ds.transpose().to_dataframe().reset_index(level=['index', 'time'])
    #drop rows with too many NaN values (thresh = number of non-NaN values)
    df_ = df_.drop(columns='index').dropna(subset=['vgl','dpdt','dthdt','prec','dwddt'],thresh=3).reset_index(drop=True)
    ds.close()
    
    
    #add column called "label" for the feature:
    # 0: None, 1: WJ, 2: CFC, 3: CJ, 4: SJ, 5: CS
    df_.insert(len(df_.columns), 'label', 0)
    df_.insert(len(df_.columns), 'storm', stormlist[i])
    
    df = pd.concat([df,df_],ignore_index=True,sort=True)
    del df_

''' 
Only stations where the normalisied wind (vgl = v/v98)
is over 0.8 are labelled. All other stations are labelled -1
'''
df.loc[df.vgl <.8,'label'] = -1

daytimes = np.unique(df.loc[df['storm']==select.value].time)
initdata = df.loc[(df['storm']==select.value) & (df['time'] == daytimes[0])]

source = ColumnDataSource(data=dict(dt=initdata['time'].values, 
                    lon=initdata['lon'].values, lat=initdata['lat'].values, 
                    vgl=initdata['vgl'].values, rr=initdata['prec'].values, 
                    p=initdata['mslp'].values, dp=initdata['dpdt'].values, 
                    th=initdata['th'].values, dth=initdata['dthdt'].values, 
                    d=initdata['wd'].values, dd=initdata['dwddt'].values,
                    lab=initdata['label'].values))

num_stat = len(source.data['lon']) #number of current stations

lon_max = 20;     lon_min = -10
lat_max = 60;     lat_min = 40
dp_max  = 5;      dp_min  = -5
dth_max = .006;   dth_min = -.006
th_min  = .95;    th_max  = 1.05
rr_max  = 5;      rr_min  = 0
dd_max  = 150;    dd_min  = -75 
p_max   = 1025;   p_min   = 970 
d_max   = 359;    d_min   = 0

#%%
#CHANGE TIMESTEP
def dt_prior():
    if t_slider.value > t_slider.start: t_slider.value -= 1
def dt_next():
    if t_slider.value < t_slider.end: t_slider.value += 1
    
#UPDATE DATA - TIMESTEP
def update_data(attrname, old, new):
    '''Update the data source of the plots'''
    # Generate the new subset of data
    updateddata = df.loc[(df['storm']==select.value) & (df['time'] == t_slider.value)]
    source.data = dict(dt=updateddata['time'].values, 
                       lon=updateddata['lon'].values, lat=updateddata['lat'].values, 
                       vgl=updateddata['vgl'].values, rr=updateddata['prec'].values, 
                       p=updateddata['mslp'].values, dp=updateddata['dpdt'].values, 
                       dth=updateddata['dthdt'].values, th=updateddata['th'].values,
                       dd=updateddata['dwddt'].values, d=updateddata['wd'].values, 
                       lab=updateddata['label'])
    update_count()
    
    #Update the histograms
    bins=40
    hh_dp2, _  = np.histogram(source.data['dp'], bins=bins,range=[dp_min,dp_max])
    hh_dth2, _ = np.histogram(source.data['dth'], bins=bins,range=[dth_min,dth_max])
    hh_th2, _  = np.histogram(source.data['th'], bins=bins,range=[th_min,th_max])
    hh_rr2, _  = np.histogram(source.data['rr'], bins=bins,range=[rr_min,rr_max])
    hh_p2, _   = np.histogram(source.data['p'], bins=bins,range=[p_min,p_max])
    hh_d2, _   = np.histogram(source.data['d'], bins=int(bins/2),range=[d_min,d_max])
    hh_dd2, _  = np.histogram(source.data['dd'], bins=bins,range=[dd_min,dd_max])
    h_dp2.data_source.data["top"]  = hh_dp2
    h_dth2.data_source.data["top"] = hh_dth2
    h_th2.data_source.data["top"]  = hh_th2
    h_rr2.data_source.data["top"]  = hh_rr2
    h_p2.data_source.data["top"]   = hh_p2
    h_d2.data_source.data["top"]   = hh_d2
    h_dd2.data_source.data["top"]  = hh_dd2
    
    if len(source.selected.indices) > 0:
        hh_dp1, _  = np.histogram(source.data['dp'][source.selected.indices], 
                                    bins=bins,range=[dp_min,dp_max])
        hh_dth1, _ = np.histogram(source.data['dth'][source.selected.indices], 
                                    bins=bins,range=[dth_min,dth_max])
        hh_th1, _  = np.histogram(source.data['th'][source.selected.indices], 
                                    bins=bins,range=[th_min,th_max])
        hh_rr1, _  = np.histogram(source.data['rr'][source.selected.indices], 
                                    bins=bins,range=[rr_min,rr_max])
        hh_p1, _   = np.histogram(source.data['p'][source.selected.indices], 
                                    bins=bins,range=[p_min,p_max])
        hh_d1, _   = np.histogram(source.data['d'][source.selected.indices], 
                                    bins=int(bins/2),range=[d_min,d_max])
        hh_dd1, _  = np.histogram(source.data['dd'][source.selected.indices], 
                                    bins=bins,range=[dd_min,dd_max])
        h_dp1.data_source.data["top"]  = hh_dp1
        h_dth1.data_source.data["top"] = hh_dth1
        h_th1.data_source.data["top"]  = hh_th1
        h_rr1.data_source.data["top"]  = hh_rr1
        h_p1.data_source.data["top"]   = hh_p1
        h_d1.data_source.data["top"]   = hh_d1
        h_dd1.data_source.data["top"]  = hh_dd1
        
        num_stat = len(source.selected.indices)
        count.text= f'<p style="font-size:20px; color:#000000">Stations selected: {num_stat} </font>'
    else:
        h_dp1.data_source.data["top"]  = hh_dp2
        h_dth1.data_source.data["top"] = hh_dth2
        h_th1.data_source.data["top"]  = hh_th2
        h_rr1.data_source.data["top"]  = hh_rr2
        h_p1.data_source.data["top"]   = hh_p2
        h_d1.data_source.data["top"]   = hh_d2
        h_dd1.data_source.data["top"]  = hh_dd2
        
        num_stat = len(source.data['lon'])
        count.text= f'<p style="font-size:20px; color:#000000">Stations selected: {num_stat} </font>'

#UPDATE DATA - STORM
def sel_update(attrname,old,new):
    print(f'Change from {old} to {new}')
    
    stormdate = stormdates[stormlist.index(new)]    
    daytimes = np.unique(df.loc[df['storm'] == new].time)
    t_slider.start = daytimes[0]
    t_slider.value = daytimes[0]
    t_slider.end = daytimes[-1]
    initdata = df.loc[(df['storm'] == new) & (df['time'] == daytimes[0])]
    source.data = dict(dt=initdata['time'].values, 
                       lon=initdata['lon'].values, lat=initdata['lat'].values, 
                       vgl=initdata['vgl'].values, rr=initdata['prec'].values,
                       p=initdata['mslp'].values, dp=initdata['dpdt'].values, 
                       th=initdata['th'].values, dth=initdata['dthdt'].values, 
                       d=initdata['wd'].values, dd=initdata['dwddt'].values,
                       lab=initdata['label'])
    
    storminfo.text = f'<p style="font-size:20px; color:#000000">{new[:-6]} ({stormdate[0][6:8]}.{stormdate[0][4:6]}.{stormdate[0][:4]})</font>'
        
    num_stat = len(source.data['lon'])
    count.text= f'<p style="font-size:20px; color:#000000">Stations selected: {num_stat} </font>'

#UPDATE DATA - SELECTION
def update(attr,old,new):
    new = new
    #Update histograms
    if len(new) > 0:
        hh_dp1, _  = np.histogram(source.data['dp'][new],bins = dp_edges)
        hh_dth1, _ = np.histogram(source.data['dth'][new],bins = dth_edges)
        hh_th1, _  = np.histogram(source.data['th'][new],bins = th_edges)
        hh_rr1, _  = np.histogram(source.data['rr'][new],bins = rr_edges)
        hh_p1, _   = np.histogram(source.data['p'][new],bins = p_edges)
        hh_d1, _   = np.histogram(source.data['d'][new],bins = d_edges)
        hh_dd1, _  = np.histogram(source.data['dd'][new],bins = dd_edges)
    else:
        hh_dp1, _  = np.histogram(source.data['dp'],bins = dp_edges)
        hh_dth1, _ = np.histogram(source.data['dth'],bins = dth_edges)
        hh_th1, _  = np.histogram(source.data['th'],bins = th_edges)
        hh_rr1, _  = np.histogram(source.data['rr'],bins = rr_edges)
        hh_p1, _   = np.histogram(source.data['p'],bins = p_edges)
        hh_d1, _   = np.histogram(source.data['d'],bins = d_edges)
        hh_dd1, _  = np.histogram(source.data['dd'],bins = dd_edges)
    h_dp1.data_source.data["top"]  = hh_dp1
    h_dth1.data_source.data["top"] = hh_dth1
    h_th1.data_source.data["top"]  = hh_th1
    h_rr1.data_source.data["top"]  = hh_rr1
    h_p1.data_source.data["top"]   = hh_p1
    h_d1.data_source.data["top"]   = hh_d1
    h_dd1.data_source.data["top"]  = hh_dd1
    
    count.text= f'<p style="font-size:20px; color:#000000">Stations selected: {len(new)} </font>'

#SETTING LABELS
def set_label(val):
    datas = copy.copy(source.data['lab'])
    boolean = np.full(len(datas),False)
    boolean[source.selected.indices] = True
    datas[boolean] = val 
    datas[(source.data['vgl']<.8)] = -1
    source.data['lab'] = datas
    update_count()
    
def update_count():
    c_WJ = sum(source.data['lab']==1)
    count_WJ.text=f'<p style="font-size:18px; color:#000000">labeled WJ: {c_WJ} </font>'
    c_CFC = sum(source.data['lab']==2)
    count_CFC.text=f'<p style="font-size:18px; color:#000000">labeled CFC: {c_CFC} </font>'
    c_CJ = sum(source.data['lab']==3)
    count_CJ.text=f'<p style="font-size:18px; color:#000000">labeled CJ: {c_CJ} </font>'
    c_SJ = sum(source.data['lab']==4)
    count_SJ.text=f'<p style="font-size:18px; color:#000000">labeled SJ: {c_SJ} </font>'
    c_CS = sum(source.data['lab']==5)
    count_CS.text=f'<p style="font-size:18px; color:#000000">labeled CS: {c_CS} </font>'
    
def set_label_WJ():
    WJ_button.button_type = 'danger'
    set_label(1)
    WJ_button.button_type = 'default'
    print(f'SET {len(source.selected.indices)} WJ LABELS')
def set_label_CFC():
    CFC_button.button_type = 'danger'
    set_label(2)
    CFC_button.button_type = 'default'
    print(f'SET {len(source.selected.indices)} CFC LABELS')
def set_label_CJ():
    CJ_button.button_type = 'danger'
    set_label(3)
    CJ_button.button_type = 'default'
    print(f'SET {len(source.selected.indices)} CJ LABELS')
def set_label_SJ():
    SJ_button.button_type = 'danger'
    set_label(4)
    SJ_button.button_type = 'default'
    print(f'SET {len(source.selected.indices)} SJ LABELS')
def set_label_CS():
    CS_button.button_type = 'danger'
    set_label(5)
    CS_button.button_type = 'default'
    print(f'SET {len(source.selected.indices)} CS LABELS')
def clear_labels():
    clabels_button.button_type = 'danger'
    set_label(0)
    clabels_button.button_type = 'default'
    print(f'CLEARED {len(source.selected.indices)} LABELS')
def clear_all_labels():
    calabels_button.button_type = 'danger'
    source.data['lab'][:] = 0
    update_count()
    calabels_button.button_type = 'default'
    print('CLEARED ALL LABELS')

#Save set labels and all parameters for the selected time step
def save_timestep():
    save_button.button_type = 'danger'
    
    c_WJ = sum(source.data['lab']==1)
    count_WJ.text=f'<p style="font-size:18px; color:#000000">labeled WJ: {c_WJ} </font>'
    c_CFC = sum(source.data['lab']==2)
    count_CFC.text=f'<p style="font-size:18px; color:#000000">labeled CFC: {c_CFC} </font>'
    c_CJ = sum(source.data['lab']==3)
    count_CJ.text=f'<p style="font-size:18px; color:#000000">labeled CJ: {c_CJ} </font>'
    c_SJ = sum(source.data['lab']==4)
    count_SJ.text=f'<p style="font-size:18px; color:#000000">labeled SJ: {c_SJ} </font>'
    c_CS = sum(source.data['lab']==5)
    count_CS.text=f'<p style="font-size:18px; color:#000000">labeled CS: {c_CS} </font>'
    
    df_end = pd.DataFrame(data={'time':source.data['dt'], 
                                'lon':source.data['lon'], 'lat':source.data['lat'], 
                                'vgl':source.data['vgl'], 'prec':source.data['rr'],
                                'mslp':source.data['p'], 'dpdt':source.data['dp'], 
                                'th':source.data['th'], 'dthdt':source.data['dth'], 
                                'wd':source.data['d'], 'dwddt':source.data['dd'],
                                'label':source.data['lab']})
    #Make sure just stations with v/v98 > 0.8 are labelled
    df_end.loc[df_end.vgl < .8,'lab'] =-1
    
    #Save labels also in loaded data
    global df
    df.loc[(df['storm']==select.value) & (df['time'] == t_slider.value),'label'] = df_end.label
    
    stormdate = stormdates[stormlist.index(select.value)]
    if not os.path.exists(out_path): os.makedirs(out_path)
    df_end = df_end.reset_index().set_index(['index','time']).to_xarray().to_netcdf(f"{out_path}{stormdate[0]}_{str(int(t_slider.value)).zfill(2)}.nc")
    save_button.button_type = 'default'
    print(f'SAVED {stormdate[0]}_{int(t_slider.value)}.nc')

#Histograms
bins=40
dp_hist, dp_edges   = np.histogram(source.data['dp'], bins=bins,range=[dp_min,dp_max])
dth_hist, dth_edges = np.histogram(source.data['dth'], bins=bins,range=[dth_min,dth_max])
th_hist, th_edges   = np.histogram(source.data['th'], bins=bins,range=[th_min,th_max])
rr_hist, rr_edges   = np.histogram(source.data['rr'], bins=bins,range=[rr_min,rr_max])
p_hist, p_edges     = np.histogram(source.data['p'], bins=bins,range=[p_min,p_max])
d_hist, d_edges     = np.histogram(source.data['d'], bins=int(bins/2),range=[d_min,d_max])
dd_hist, dd_edges   = np.histogram(source.data['dd'], bins=bins,range=[dd_min,dd_max])
hzeros = np.zeros(bins)

#%% PLOTTING
Reds = list(Reds9); RdYlBu = list(RdYlBu11); Blues = list(Blues9)
Viridis = list(Viridis11)
Reds.reverse(); Blues.reverse()
lab_cols = ['#888888','#ff0000','#008800','#0000ff','#8800ff','#ffff00']
plot_width = np.round(500*(lon_max-lon_min)/(lat_max-lat_min)).astype(int)


### PLOT MAPS
colormapper = LinearColorMapper(palette=Reds[1:-1], low=.8, high=1.7, 
                                 low_color= '#c9c9c9', high_color=Reds[-1])
colormapper_ns = LinearColorMapper(palette=Blues[1:-1], low=.8, high=1.5, 
                                 low_color= Blues[0], high_color=Blues[-1]) # not selected
map_vgl = figure(plot_height=500, plot_width=plot_width, title="Normalised wind speed v/v98", 
                   x_range=[lon_min, lon_max], y_range=[lat_min, lat_max], 
                   toolbar_location=None, tools=tools, tooltips=tooltips)
p1_vgl = map_vgl.scatter('lon', 'lat', source=source, 
                           color={'field': 'vgl', 'transform': colormapper}, 
                           line_color=None, alpha= 1.0 , 
                           nonselection_fill_color="fill_color", size = 8)
p1_vgl.nonselection_glyph.fill_alpha= 0.8
p1_vgl.nonselection_glyph.fill_color={'field': 'vgl', 'transform': colormapper_ns}
color_bar = ColorBar(color_mapper=colormapper, ticker=BasicTicker(), 
                     label_standoff=12, border_line_color=None, location=(0,0))
map_vgl.add_layout(color_bar,'right')
tab1 = Panel(child=map_vgl, title='v/v98')

colormapper = LinearColorMapper(palette=RdYlBu, low=-2, high=2)
map_dp = figure(plot_height=500, plot_width=plot_width, 
                title="Mean-sea level pressure tendency dp in hPa/h", 
                x_range=map_vgl.x_range, y_range=map_vgl.y_range, 
                toolbar_location=None, tools=tools, tooltips=tooltips)
p1_dp = map_dp.scatter('lon', 'lat', source=source, 
                       color={'field': 'dp', 'transform': colormapper}, 
                       line_color=None, alpha= 1.0, size = 8)
p1_dp.nonselection_glyph.fill_alpha= .5
p1_dp.nonselection_glyph.fill_color={'field': 'dp', 'transform': colormapper}
color_bar = ColorBar(color_mapper=colormapper, ticker=BasicTicker(), 
                     label_standoff=12, border_line_color=None, location=(0,0))
map_dp.add_layout(color_bar,'right')
tab2 = Panel(child=map_dp, title='dp')

colormapper = LinearColorMapper(palette=Viridis, low=960, high=1015, 
                                low_color=Viridis[0],high_color=Viridis[-1])
map_p = figure(plot_height=500, plot_width=plot_width, 
               title="Mean-sea level pressure p in hPa", 
               x_range=map_vgl.x_range, y_range=map_vgl.y_range, 
               toolbar_location=None, tools=tools, tooltips=tooltips)
p1_p = map_p.scatter('lon', 'lat', source=source, 
                     color={'field': 'p', 'transform': colormapper}, 
                     line_color=None, alpha= 1.0, size = 8)
p1_p.nonselection_glyph.fill_alpha= .5
p1_p.nonselection_glyph.fill_color={'field': 'p', 'transform': colormapper}
color_bar = ColorBar(color_mapper=colormapper, ticker=BasicTicker(), 
                     label_standoff=12, border_line_color=None, location=(0,0))
map_p.add_layout(color_bar,'right')
tab2_ = Panel(child=map_p, title='p')

colormapper = LinearColorMapper(palette=RdYlBu, low=-0.008, high=0.008)
map_dth = figure(plot_height=500, plot_width=plot_width, 
                 title=r"Normalised potential temperature tendency d(th/th50)", 
                 x_range=map_vgl.x_range, y_range=map_vgl.y_range, 
                 toolbar_location=None, tools=tools, tooltips=tooltips)
p1_dth = map_dth.scatter('lon', 'lat', source=source, 
                         color={'field': 'dth', 'transform': colormapper}, 
                         line_color=None, alpha= 1.0, size = 8)
p1_dth.nonselection_glyph.fill_alpha= .5
p1_dth.nonselection_glyph.fill_color={'field': 'dth', 'transform': colormapper}
color_bar = ColorBar(color_mapper=colormapper, ticker=BasicTicker(), 
                     label_standoff=12, border_line_color=None, location=(0,0))
map_dth.add_layout(color_bar,'right')
tab3 = Panel(child=map_dth, title='dth')

colormapper = LinearColorMapper(palette=RdYlBu, low=.975, high=1.025)
map_th = figure(plot_height=500, plot_width=plot_width, 
                title="Normalised potential temperature th/th50", 
                x_range=map_vgl.x_range, y_range=map_vgl.y_range, 
                toolbar_location=None, tools=tools, tooltips=tooltips)
p1_th = map_th.scatter('lon', 'lat', source=source, 
                       color={'field': 'th', 'transform': colormapper}, 
                       line_color=None, alpha= 1.0, size = 8)
p1_th.nonselection_glyph.fill_alpha= .5
p1_th.nonselection_glyph.fill_color={'field': 'th', 'transform': colormapper}
color_bar = ColorBar(color_mapper=colormapper, ticker=BasicTicker(), 
                     label_standoff=12, border_line_color=None, location=(0,0))
map_th.add_layout(color_bar,'right')
tab3_ = Panel(child=map_th, title='th/th50')

colormapper = LinearColorMapper(palette=Reds[1:-1], low=.2, high=3.6, 
                                high_color=Reds[-1], low_color='#ffffff')
colormapper_ns = LinearColorMapper(palette=Blues[1:-1], low=.2, high=3.6, 
                                   high_color=Blues[-1], low_color='#ffffff')
map_rr = figure(plot_height=500, plot_width=plot_width, 
                title="Precipitation amount RR in mm/h", 
                x_range=map_vgl.x_range, y_range=map_vgl.y_range, 
                toolbar_location=None, tools=tools, tooltips=tooltips)
p1_rr = map_rr.scatter('lon', 'lat', source=source, 
                        color={'field': 'rr', 'transform': colormapper}, 
                        line_color=None, alpha= 1.0, size = 8)
p1_rr.nonselection_glyph.fill_alpha= 0.8
p1_rr.nonselection_glyph.fill_color={'field': 'rr', 'transform': colormapper_ns}
color_bar = ColorBar(color_mapper=colormapper, ticker=BasicTicker(), 
                     label_standoff=12, border_line_color=None, location=(0,0))
map_rr.add_layout(color_bar,'right')
tab4 = Panel(child=map_rr, title='RR')

colormapper = LinearColorMapper(palette=RdYlBu, low=-50, high=50)
map_dd = figure(plot_height=500, plot_width=plot_width, 
                title="Wind direction tendency dd in deg/h", 
                x_range=map_vgl.x_range, y_range=map_vgl.y_range, 
                toolbar_location=None, tools=tools, tooltips=tooltips)
p1_dd = map_dd.scatter('lon', 'lat', source=source, 
                       color={'field': 'dd', 'transform': colormapper}, 
                       line_color=None, alpha= 1.0, size = 8)
p1_dd.nonselection_glyph.fill_alpha= .5
p1_dd.nonselection_glyph.fill_color={'field': 'dd', 'transform': colormapper}
color_bar = ColorBar(color_mapper=colormapper, ticker=BasicTicker(), 
                     label_standoff=12, border_line_color=None, location=(0,0))
map_dd.add_layout(color_bar,'right')
tab5 = Panel(child=map_dd, title='dd')

colormapper = LinearColorMapper(palette=Viridis, low=0, high=359)
map_d = figure(plot_height=500, plot_width=plot_width, title="Wind direction d in deg", 
                     x_range=map_vgl.x_range, y_range=map_vgl.y_range, 
                     toolbar_location=None, tools=tools, tooltips=tooltips)
p1_d = map_d.scatter('lon', 'lat', source=source, 
                         color={'field': 'd', 'transform': colormapper}, 
                         line_color=None, alpha= 1.0, size = 8)
p1_d.nonselection_glyph.fill_alpha= .5
p1_d.nonselection_glyph.fill_color={'field': 'd', 'transform': colormapper}
color_bar = ColorBar(color_mapper=colormapper, ticker=BasicTicker(), 
                     label_standoff=12, border_line_color=None, location=(0,0))
map_d.add_layout(color_bar,'right')
tab5_5 = Panel(child=map_d, title='d')

colormapper = LinearColorMapper(palette=lab_cols, low=-.5, high=5.5, 
                                low_color='#c9c9c9')
map_lab = figure(plot_height=500, plot_width=plot_width, title="Set labels", 
                   x_range=map_vgl.x_range, y_range=map_vgl.y_range, 
                   toolbar_location=None, tools=tools, tooltips=tooltips)
p1_lab = map_lab.scatter('lon', 'lat', source=source, 
                           color={'field': 'lab', 'transform': colormapper}, 
                           line_color=None, alpha= 1.0, size = 8)
p1_lab.nonselection_glyph.fill_alpha= .8
p1_lab.nonselection_glyph.fill_color={'field': 'lab', 'transform': colormapper}
color_bar = ColorBar(color_mapper=colormapper, ticker=BasicTicker(), 
                     label_standoff=12, border_line_color=None, location=(0,0),
                     major_label_overrides={0: 'NF', 1: 'WJ', 2: 'CFC',
                                     3: 'CJ', 4: 'SJ', 5: 'CS'})
map_lab.add_layout(color_bar,'right')
tab6 = Panel(child=map_lab, title='label')

### PLOT HISTOGRAMS
hg_dp = figure(plot_height=300, plot_width=600, toolbar_location=None, 
                   x_axis_label="dp in hPa/h", y_axis_label="", 
                   x_range=[dp_min, dp_max], tools=tools)
h_dp1 = hg_dp.quad(top=dp_hist, bottom=0, left=dp_edges[:-1], 
                   right=dp_edges[1:],fill_color="navy", line_color="white")
h_dp2 = hg_dp.quad(top=dp_hist, bottom=0, left=dp_edges[:-1], 
                   right=dp_edges[1:],fill_color="navy", line_color="white",alpha=.2)

hg_dth = figure(plot_height=300, plot_width=600, toolbar_location=None, 
                       x_axis_label="dth in 1/h", y_axis_label="", 
                       x_range=[dth_min, dth_max], tools=tools)
h_dth1 = hg_dth.quad(top=dth_hist, bottom=0, left=dth_edges[:-1], 
                     right=dth_edges[1:],fill_color="navy", line_color="white")
h_dth2 = hg_dth.quad(top=dth_hist, bottom=0, left=dth_edges[:-1], 
                     right=dth_edges[1:],fill_color="navy", line_color="white",alpha=.2)

hg_th = figure(plot_height=300, plot_width=600, toolbar_location=None, 
                    x_axis_label="th/th50 in K/K", y_axis_label="", 
                    x_range=[th_min, th_max], tools=tools)
h_th1 = hg_th.quad(top=th_hist, bottom=0, left=th_edges[:-1], 
                   right=th_edges[1:],fill_color="navy", line_color="white")
h_th2 = hg_th.quad(top=th_hist, bottom=0, left=th_edges[:-1], 
                   right=th_edges[1:],fill_color="navy", line_color="white",alpha=.2)

hg_rr = figure(plot_height=300, plot_width=600, toolbar_location=None, 
               x_axis_label="RR in mm/h", y_axis_label="", 
               x_range=[rr_min, rr_max], tools=tools)
h_rr1 = hg_rr.quad(top=rr_hist, bottom=0, left=rr_edges[:-1], 
                   right=rr_edges[1:],fill_color="navy", line_color="white")
h_rr2 = hg_rr.quad(top=rr_hist, bottom=0, left=rr_edges[:-1], 
                   right=rr_edges[1:],fill_color="navy", line_color="white",alpha=.2)

hg_p = figure(plot_height=300, plot_width=600, toolbar_location=None, 
               x_axis_label="p in hPa", y_axis_label="", 
               x_range=[p_min, p_max], tools=tools)
h_p1 = hg_p.quad(top=p_hist, bottom=0, left=p_edges[:-1], 
                   right=p_edges[1:],fill_color="navy", line_color="white")
h_p2 = hg_p.quad(top=p_hist, bottom=0, left=p_edges[:-1], 
                   right=p_edges[1:],fill_color="navy", line_color="white",alpha=.2)

hg_d = figure(plot_height=300, plot_width=600, toolbar_location=None, 
               x_axis_label="d in deg", y_axis_label="", 
               x_range=[d_min, d_max], tools=tools)
h_d1 = hg_d.quad(top=d_hist, bottom=0, left=d_edges[:-1], 
                   right=d_edges[1:],fill_color="navy", line_color="white")
h_d2 = hg_d.quad(top=d_hist, bottom=0, left=d_edges[:-1], 
                   right=d_edges[1:],fill_color="navy", line_color="white",alpha=.2)

hg_dd = figure(plot_height=300, plot_width=600, toolbar_location=None, 
               x_axis_label="dd in deg/h", y_axis_label="", 
               x_range=[dd_min, dd_max], tools=tools)
h_dd1 = hg_dd.quad(top=dd_hist, bottom=0, left=dd_edges[:-1], 
                   right=dd_edges[1:],fill_color="navy", line_color="white")
h_dd2 = hg_dd.quad(top=dd_hist, bottom=0, left=dd_edges[:-1], 
                   right=dd_edges[1:],fill_color="navy", line_color="white",alpha=.2)

#TABLE
columns = [TableColumn(field="dt", title="time"),
           TableColumn(field="lat", title="lat"),
           TableColumn(field="lon", title="lon"),
           TableColumn(field="vgl", title="v/v98"),
           TableColumn(field="dp", title="dp"),
           TableColumn(field="dth", title="d(th/th50)"),
           TableColumn(field="th", title="th/th50"),
           TableColumn(field="rr", title="RR"),
           TableColumn(field="dd", title="dd"),
           TableColumn(field="p", title="p"),
           TableColumn(field="d", title="d"),
           TableColumn(field="lab", title="label")]
data_table = DataTable(source=source,columns=columns, width=600, height=500)
    
p1_vgl.data_source.selected.on_change('indices',update)
p1_dp.data_source.selected.on_change('indices',update)
p1_p.data_source.selected.on_change('indices',update)
p1_dth.data_source.selected.on_change('indices',update)
p1_th.data_source.selected.on_change('indices',update)
p1_rr.data_source.selected.on_change('indices',update)
p1_dd.data_source.selected.on_change('indices',update)
p1_d.data_source.selected.on_change('indices',update)


# Set up the widgets and connect them to the function 'update_data()'
t_slider = Slider(start=daytimes[0], end=daytimes[-1], value=daytimes[0], step=1, title="Daytime", width=120)
t_slider.on_change('value', update_data)
daytime_prior = Button(label='<', width = 35)
daytime_next = Button(label='>', width = 35)
daytime_prior.on_click(dt_prior);   daytime_next.on_click(dt_next)
select.on_change('value', sel_update)


### DASHBOARD
count = Div(text=f'<p style="font-size:20px; color:#000000">Stations selected: {num_stat} </font>')
clear_button = Button(label='Clear selection', width=100)

clear_button.js_on_click(CustomJS(args=dict(source = source), code="""
                                  source.selected.indices = [];
                                  source.change.emit()
                                  """))
    

count_WJ = Div(text='<p style="font-size:18px; color:#000000">labeled WJ: 0 </font>')
count_CFC = Div(text='<p style="font-size:18px; color:#000000">labeled CFC: 0 </font>')
count_CJ = Div(text='<p style="font-size:18px; color:#000000">labeled CJ: 0 </font>')
count_SJ = Div(text='<p style="font-size:18px; color:#000000">labeled SJ: 0 </font>')
count_CS = Div(text='<p style="font-size:18px; color:#000000">labeled CS: 0 </font>')

WJ_button  = Button(label='Set WJ label', width=100)
CFC_button = Button(label='Set CFC label', width=100)
CJ_button  = Button(label='Set CJ label', width=100)
SJ_button  = Button(label='Set SJ label', width=100)
CS_button  = Button(label='Set CS label', width=100)
clabels_button  = Button(label='Clear', width=100)
calabels_button = Button(label='Clear all', width=100)
save_button = Button(label='Save timestep', width=100)

WJ_button.on_click(set_label_WJ)
CFC_button.on_click(set_label_CFC)
CJ_button.on_click(set_label_CJ)
SJ_button.on_click(set_label_SJ)
CS_button.on_click(set_label_CS)
clabels_button.on_click(clear_labels)
calabels_button.on_click(clear_all_labels)
save_button.on_click(save_timestep)


title = '<p style="font-size:25px; font-weight: bold; color:#000000">RAMEFI - data explorer</p>'
instruction = 'Select on maps and set labels according to features; clear selection to reset.'
storminfo = Div(text=f'<p style="font-size:20px; color:#000000">{select.value[:-6]} ({stormdate[6:8]}.{stormdate[4:6]}.{stormdate[:4]})</font>', width=500, height=20)
w2w_logo   = '<a href="https://www.wavestoweather.de/research_areas/phase2/c5/index.html"><img src="https://www.ipa.uni-mainz.de/files/2015/11/Logo_W2W_text.png" width=170></a>'

header = row(Div(text=w2w_logo, width=170), Spacer(width=20),
                column(Div(text=title, height=55, width=750), 
                       Div(text=instruction, width=500, height=20), storminfo),
                Spacer(width=20), column(count, clear_button, 
                                         row(daytime_prior, t_slider, daytime_next), select), 
                sizing_mode='fixed', css_classes=['scrollable'])

buttons = column(row(WJ_button, CFC_button, CJ_button, SJ_button, CS_button, 
                     clabels_button, calabels_button, save_button),
                 row(count_WJ, count_CFC, count_CJ, count_SJ, count_CS))

# Group the plots in a gridplot with a merged toolbar
body = gridplot([[Tabs(tabs=[tab1,tab2,tab2_,tab3,tab3_,tab4,tab5,tab5_5,tab6]), data_table],
                 [hg_dp, hg_rr], [hg_dth, hg_th],[hg_p, hg_d]],
                 toolbar_location='left')
body.children[1].css_classes = ['scrollable']
body.children[1].sizing_mode = 'fixed'
body.children[1].height = 2000
body.children[1].width = 1400

curdoc().add_root(column(header, Spacer(width=50), buttons, Spacer(width=50), body))
curdoc().title = "RAMEFI - Data Explorer"
