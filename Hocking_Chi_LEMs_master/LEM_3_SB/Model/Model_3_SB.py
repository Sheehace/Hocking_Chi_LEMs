# Hocking Hills Chi-Experiment Models
# Author: Chris Sheehan

#%% Imports and paths

# Print status
print('Importing libraries, capturing paths, importing input parameters, and creating output directories...')

# Libraries
import inspect
import os
import os.path
from os import path
import numpy as np
import pandas as pd
from landlab import RasterModelGrid, imshow_grid
from landlab.components import StreamPowerEroder, LinearDiffuser, FlowAccumulator, ChannelProfiler, PriorityFloodFlowRouter
from landlab.io.esri_ascii import write_esri_ascii
from matplotlib import pyplot as plt
from osgeo import gdal

# Capture paths
script_name = os.path.basename(__file__)
script_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
model_directory = script_directory.replace('\\model', '')
input_directory = model_directory + '\\input'

# Import inputs
su_Params = pd.read_csv(input_directory + '\Spin_up.csv')
mr_Params = pd.read_csv(input_directory + '\Main_run.csv')

# Create output directories
if path.exists(str(model_directory)+'/Output') == False:
    os.mkdir(str(model_directory)+'/Output')
if path.exists(str(model_directory)+'/Output/Spin_up') == False:
    os.mkdir(str(model_directory)+'/Output/Spin_up')    
if path.exists(str(model_directory)+'/Output/Main_run') == False:
    os.mkdir(str(model_directory)+'/Output/Main_run')       
if path.exists(str(model_directory)+'/Output/Spin_up/topographic__elevation__png') == False:
    os.mkdir(str(model_directory)+'/Output/Spin_up/topographic__elevation__png')
if path.exists(str(model_directory)+'/Output/Spin_up/topographic__elevation__pdf') == False:
    os.mkdir(str(model_directory)+'/Output/Spin_up/topographic__elevation__pdf')
if path.exists(str(model_directory)+'/Output/Spin_up/topographic__elevation__asc') == False:
    os.mkdir(str(model_directory)+'/Output/Spin_up/topographic__elevation__asc')
if path.exists(str(model_directory)+'/Output/Spin_up/topographic__elevation__tif') == False:
    os.mkdir(str(model_directory)+'/Output/Spin_up/topographic__elevation__tif')
if path.exists(str(model_directory)+'/Output/Spin_up/dzdt__png') == False:
    os.mkdir(str(model_directory)+'/Output/Spin_up/dzdt__png')
if path.exists(str(model_directory)+'/Output/Spin_up/dzdt__pdf') == False:
    os.mkdir(str(model_directory)+'/Output/Spin_up/dzdt__pdf')
if path.exists(str(model_directory)+'/Output/Spin_up/dzdt__asc') == False:
    os.mkdir(str(model_directory)+'/Output/Spin_up/dzdt__asc')
if path.exists(str(model_directory)+'/Output/Spin_up/dzdt__tif') == False:
    os.mkdir(str(model_directory)+'/Output/Spin_up/dzdt__tif')
if path.exists(str(model_directory)+'/Output/Spin_up/channel_map__png') == False:
    os.mkdir(str(model_directory)+'/Output/Spin_up/channel_map__png')
if path.exists(str(model_directory)+'/Output/Spin_up/channel_map__pdf') == False:
    os.mkdir(str(model_directory)+'/Output/Spin_up/channel_map__pdf')
if path.exists(str(model_directory)+'/Output/Main_run/topographic__elevation__png') == False:
    os.mkdir(str(model_directory)+'/Output/Main_run/topographic__elevation__png')
if path.exists(str(model_directory)+'/Output/Main_run/topographic__elevation__pdf') == False:
    os.mkdir(str(model_directory)+'/Output/Main_run/topographic__elevation__pdf')
if path.exists(str(model_directory)+'/Output/Main_run/topographic__elevation__asc') == False:
    os.mkdir(str(model_directory)+'/Output/Main_run/topographic__elevation__asc')
if path.exists(str(model_directory)+'/Output/Main_run/topographic__elevation__tif') == False:
    os.mkdir(str(model_directory)+'/Output/Main_run/topographic__elevation__tif')
if path.exists(str(model_directory)+'/Output/Main_run/dzdt__png') == False:
    os.mkdir(str(model_directory)+'/Output/Main_run/dzdt__png')
if path.exists(str(model_directory)+'/Output/Main_run/dzdt__pdf') == False:
    os.mkdir(str(model_directory)+'/Output/Main_run/dzdt__pdf')
if path.exists(str(model_directory)+'/Output/Main_run/dzdt__asc') == False:
    os.mkdir(str(model_directory)+'/Output/Main_run/dzdt__asc')
if path.exists(str(model_directory)+'/Output/Main_run/dzdt__tif') == False:
    os.mkdir(str(model_directory)+'/Output/Main_run/dzdt__tif')
if path.exists(str(model_directory)+'/Output/Main_run/channel_map__png') == False:
    os.mkdir(str(model_directory)+'/Output/Main_run/channel_map__png')
if path.exists(str(model_directory)+'/Output/Main_run/channel_map__pdf') == False:
    os.mkdir(str(model_directory)+'/Output/Main_run/channel_map__pdf')
    
# Print space
print(' ')
    
#%% Retrieve spin-up parameters

# Print status
print('Retrieving spin-up parameters...')

# Assign inputs
number_of_rows = su_Params.number_of_rows[0]         
number_of_columns = su_Params.number_of_columns[0]     
dxy = su_Params.dxy[0]
su_msp = su_Params.msp[0]
su_nsp = su_Params.nsp[0]
su_ksp_1 = su_Params.ksp_1[0]
su_ksp_2 = su_Params.ksp_2[0] 
su_contact_x = su_Params.contact_x[0] 
su_khs = su_Params.khs[0]
su_u = su_Params.u[0]
su_dt = su_Params.dt[0]     
su_tmax = su_Params.tmax[0]
su_export_interval = su_Params.export_interval[0]

# Print space
print(' ')

#%% Create grid

# Print status
print('Creating model grid...')

# Create grid
mg = RasterModelGrid((number_of_rows, number_of_columns), dxy)
mg.add_zeros('node', 'topographic__elevation')
mg.add_zeros('node', 'previous__topographic__elevation')
mg.add_zeros('node', 'dzdt')
np.random.seed(0)                                      
mg_noise = np.random.rand(mg.number_of_nodes)/1000.    
mg.at_node['topographic__elevation'] += mg_noise

# Print space
print(' ')

#%% Spin-up boundary conditions and parameter assignments

# Print status
print('Assigning spin-up boundary conditions and parameter values...')

# Boundary conditions
mg.set_status_at_node_on_edges(right=4, top=4, left=4, bottom=4)
mg.status_at_node[int(number_of_columns * (number_of_rows / 2))] = mg.BC_NODE_IS_FIXED_VALUE

# Create lithology
Ksp = np.ones(mg.number_of_nodes) * su_ksp_1 
Khs = np.ones(mg.number_of_nodes) * su_khs 

# Create uplift field
U = np.ones(mg.number_of_nodes) * su_u

# Create clock
t = np.arange(0, su_tmax, su_dt) 

# Set tickers
total_time = 0 
export_ticker = 0
lith_ticker = 1

# Print space
print(' ')

#%% Spin-up

# Print status
print('Performing model spin-up for ', su_tmax, ' years...')

# Initiate components
frr = FlowAccumulator(mg, flow_director='D8')
spr = StreamPowerEroder(mg, K_sp=Ksp, m_sp=su_msp, n_sp=su_nsp)
dfn = LinearDiffuser(mg, linear_diffusivity=su_khs)

# Initiate elevational tracking
mg.at_node['previous__topographic__elevation'][mg.core_nodes] = mg.at_node['topographic__elevation'][mg.core_nodes]

# Time Loop
for ti in t:
    
    # Set hard lithology
    if lith_ticker == 1:
        if ti >= 1E5:
            Ksp[np.where(mg.node_x > su_contact_x)] = su_ksp_2 
            spr = StreamPowerEroder(mg, K_sp=Ksp, m_sp=su_msp, n_sp=su_nsp)
            lith_ticker = 0
    
    # Uplift
    mg.at_node['topographic__elevation'][mg.core_nodes] += U[mg.core_nodes] * su_dt     
    
    # Erode
    dfn.run_one_step(su_dt)                                
    frr.run_one_step()                    
    spr.run_one_step(su_dt)

    # Update elevational tracking
    mg.at_node['dzdt'][mg.core_nodes] = (mg.at_node['topographic__elevation'][mg.core_nodes] - mg.at_node['previous__topographic__elevation'][mg.core_nodes]) / su_dt
    mg.at_node['previous__topographic__elevation'][mg.core_nodes] = mg.at_node['topographic__elevation'][mg.core_nodes]   
    
    # Advance clock
    total_time += su_dt     
    export_ticker += su_dt
    
    # Display status
    print('Year = ', total_time)
    
    # Export
    if export_ticker == su_export_interval:
        
        # Initiate ChannelProfiler during first export. Requires threshold drainage area during ititiation.
        if total_time == su_export_interval:
            prf = ChannelProfiler(mg, number_of_watersheds=1, main_channel_only=False, minimum_channel_threshold=100000)
        
        # Display status
        print('Exporting data... Please be patient!')
        
        # Plot topographic__elevation
        plt.ioff()
        fig = plt.figure()         
        imshow_grid(mg, 'topographic__elevation', grid_units=('m', 'm'), var_name="Elevation (m)", cmap='terrain', allow_colorbar=True)
        plt.title('$Year$=' + str(total_time))
        plt.tight_layout()
        fig.savefig(str(model_directory) + '/Output/Spin_up/topographic__elevation__png/' + str(total_time) + '.png',  format='png', dpi=300)
        fig.savefig(str(model_directory) + '/Output/Spin_up/topographic__elevation__pdf/' + str(total_time) + '.pdf',  format='pdf', dpi=300)
        write_esri_ascii(str(model_directory) + '/Output/Spin_up/topographic__elevation__asc/' + str(total_time) + '.asc', mg, names='topographic__elevation', clobber=True)
        ds = gdal.Open(str(model_directory) + '/Output/Spin_up/topographic__elevation__asc/' + str(total_time) + '.asc')
        ds = gdal.Warp(str(model_directory) + '/Output/Spin_up/topographic__elevation__tif/' + str(total_time) + '.tif', ds, format="GTiff", dstSRS="EPSG:32616")
        plt.close(fig)
        
        # Plot dzdt
        plt.ioff()
        fig = plt.figure()         
        imshow_grid(mg, 'dzdt', grid_units=('m', 'm'), var_name="Rate of elevational change ($\mathregular{m\ yr^{-1}}$)", cmap='seismic_r', allow_colorbar=True, symmetric_cbar = True)
        plt.title('$Year$=' + str(total_time))
        plt.tight_layout()
        fig.savefig(str(model_directory) + '/Output/Spin_up/dzdt__png/' + str(total_time) + '.png',  format='png', dpi=300)
        fig.savefig(str(model_directory) + '/Output/Spin_up/dzdt__pdf/' + str(total_time) + '.pdf',  format='pdf', dpi=300)
        write_esri_ascii(str(model_directory) + '/Output/Spin_up/dzdt__asc/' + str(total_time) + '.asc', mg, names='dzdt', clobber=True)
        ds = gdal.Open(str(model_directory) + '/Output/Spin_up/dzdt__asc/' + str(total_time) + '.asc')
        ds = gdal.Warp(str(model_directory) + '/Output/Spin_up/dzdt__tif/' + str(total_time) + '.tif', ds, format="GTiff", dstSRS="EPSG:32616")
        plt.close(fig)

        # Plot channel map
        plt.ioff()
        fig = plt.figure() 
        prf.run_one_step()        
        prf.plot_profiles_in_map_view()
        plt.title('$Year$=' + str(total_time))
        plt.tight_layout()
        fig.savefig(str(model_directory) + '/Output/Spin_up/channel_map__png/' + str(total_time) + '.png',  format='png', dpi=300)
        fig.savefig(str(model_directory) + '/Output/Spin_up/channel_map__pdf/' + str(total_time) + '.pdf',  format='pdf', dpi=300)
        plt.close(fig)
        
        # Reset ticker
        export_ticker = 0
        
# Print status / space
print('Spin-up complete!')
print(' ')
        
#%% Retrieve main run parameters

# Print status
print('Retrieving main_run parameters...')

# Assign inputs
breach_x = mr_Params.breach_x[0]
mr_msp = mr_Params.msp[0]
mr_nsp = mr_Params.nsp[0]
mr_ksp_1 = mr_Params.ksp_1[0]
mr_ksp_2 = mr_Params.ksp_2[0] 
mr_contact_x = mr_Params.contact_x[0] 
mr_khs = mr_Params.khs[0]
mr_u = mr_Params.u[0]
mr_dt = int(mr_Params.dt[0])     
mr_tmax = mr_Params.tmax[0]
mr_export_interval = mr_Params.export_interval[0]

# Print space
print(' ')
       
#%% Force boundary condition change with an ice sheet

# Print status
print('Forcing main run boundary condition with virtual ice sheet...')

# Close all boundaries
mg.set_status_at_node_on_edges(right=4, top=4, left=4, bottom=4)

# Open breach
mg.status_at_node[int(breach_x)] = mg.BC_NODE_IS_FIXED_VALUE

# Print space
print(' ')

#%% Main run parameter assignments

# Print status
print('Assigning main run parameter values...')

# Create clock
t = np.arange(0, mr_tmax, mr_dt) 

# Set tickers
total_time = 0 
export_ticker = 0

# Print space
print(' ')

#%% Main run

# Print status
print('Performing main run for ', mr_tmax, ' years...')

# Initiate components
pffr = PriorityFloodFlowRouter(mg, depression_handler = 'fill')

# Initiate elevational tracking
mg.at_node['previous__topographic__elevation'][mg.core_nodes] = mg.at_node['topographic__elevation'][mg.core_nodes]

# Time Loop
for ti in t:
    
    # Uplift
    mg.at_node['topographic__elevation'][mg.core_nodes] += U[mg.core_nodes] * mr_dt     
    
    # Erode
    dfn.run_one_step(mr_dt)                                
    pffr.run_one_step()                    
    spr.run_one_step(mr_dt)

    # Update elevational tracking
    mg.at_node['dzdt'][mg.core_nodes] = (mg.at_node['topographic__elevation'][mg.core_nodes] - mg.at_node['previous__topographic__elevation'][mg.core_nodes]) / mr_dt
    mg.at_node['previous__topographic__elevation'][mg.core_nodes] = mg.at_node['topographic__elevation'][mg.core_nodes]   
    
    # Advance clock
    total_time += mr_dt     
    export_ticker += mr_dt
    
    # Display status
    print('Year = ', total_time)
    
    # Export
    if export_ticker == mr_export_interval:
        
        # Initiate ChannelProfiler during first export. Requires threshold drainage area during ititiation.
        if total_time == mr_export_interval:
            prf = ChannelProfiler(mg, number_of_watersheds=1, main_channel_only=False, minimum_channel_threshold=100000)
        
        # Display status
        print('Exporting data... Please be patient!')
        
        # Plot topographic__elevation
        plt.ioff()
        fig = plt.figure()         
        imshow_grid(mg, 'topographic__elevation', grid_units=('m', 'm'), var_name="Elevation (m)", cmap='terrain', allow_colorbar=True)
        plt.title('$Year$=' + str(total_time))
        plt.tight_layout()
        fig.savefig(str(model_directory) + '/Output/Main_run/topographic__elevation__png/' + str(total_time) + '.png',  format='png', dpi=300)
        fig.savefig(str(model_directory) + '/Output/Main_run/topographic__elevation__pdf/' + str(total_time) + '.pdf',  format='pdf', dpi=300)
        write_esri_ascii(str(model_directory) + '/Output/Main_run/topographic__elevation__asc/' + str(total_time) + '.asc', mg, names='topographic__elevation', clobber=True)
        ds = gdal.Open(str(model_directory) + '/Output/Main_run/topographic__elevation__asc/' + str(total_time) + '.asc')
        ds = gdal.Warp(str(model_directory) + '/Output/Main_run/topographic__elevation__tif/' + str(total_time) + '.tif', ds, format="GTiff", dstSRS="EPSG:32616")
        plt.close(fig)
        
        # Plot dzdt
        plt.ioff()
        fig = plt.figure()         
        imshow_grid(mg, 'dzdt', grid_units=('m', 'm'), var_name="Rate of elevational change ($\mathregular{m\ yr^{-1}}$)", cmap='seismic_r', allow_colorbar=True, symmetric_cbar = True)
        plt.title('$Year$=' + str(total_time))
        plt.tight_layout()
        fig.savefig(str(model_directory) + '/Output/Main_run/dzdt__png/' + str(total_time) + '.png',  format='png', dpi=300)
        fig.savefig(str(model_directory) + '/Output/Main_run/dzdt__pdf/' + str(total_time) + '.pdf',  format='pdf', dpi=300)
        write_esri_ascii(str(model_directory) + '/Output/Main_run/dzdt__asc/' + str(total_time) + '.asc', mg, names='dzdt', clobber=True)
        ds = gdal.Open(str(model_directory) + '/Output/Main_run/dzdt__asc/' + str(total_time) + '.asc')
        ds = gdal.Warp(str(model_directory) + '/Output/Main_run/dzdt__tif/' + str(total_time) + '.tif', ds, format="GTiff", dstSRS="EPSG:32616")
        plt.close(fig)

        # Plot channel map
        plt.ioff()
        fig = plt.figure() 
        prf.run_one_step()        
        prf.plot_profiles_in_map_view()
        plt.title('$Year$=' + str(total_time))
        plt.tight_layout()
        fig.savefig(str(model_directory) + '/Output/Main_run/channel_map__png/' + str(total_time) + '.png',  format='png', dpi=300)
        fig.savefig(str(model_directory) + '/Output/Main_run/channel_map__pdf/' + str(total_time) + '.pdf',  format='pdf', dpi=300)
        plt.close(fig)
        
        # Reset ticker
        export_ticker = 0
        
# Print status / space
print('Main run complete!')
print(' ')
print('Landscape Evolution Model: "', script_name, '" complete. Happy trails!')