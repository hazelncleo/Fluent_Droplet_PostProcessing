import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import os

plt.style.use('ggplot')

df_1 = pd.read_csv(os.path.join('output','output_data_1000.csv'), index_col='timestep_number')
df_075 = pd.read_csv(os.path.join('output','output_data_750.csv'), index_col='timestep_number')
df_05 = pd.read_csv(os.path.join('output','output_data_500.csv'), index_col='timestep_number')

f,ax = plt.subplots(1,3, figsize = (10,4))

ax[0].plot(df_1['analysis_time'],df_1['max_shearrate'],'-k')
ax[0].plot(df_075['analysis_time'],df_075['max_shearrate'],'-r')
ax[0].plot(df_05['analysis_time'],df_05['max_shearrate'],'-g')
ax[0].set_yscale('log')

ax[1].plot(df_1['analysis_time'],df_1['total_volume_delivered'],'-k')
ax[1].plot(df_075['analysis_time'],df_075['total_volume_delivered'],'-r')
ax[1].plot(df_05['analysis_time'],df_05['total_volume_delivered'],'-g')

ax[2].plot(df_1['analysis_time'],df_1['volumetric_flowrate'],'-k')
ax[2].plot(df_075['analysis_time'],df_075['volumetric_flowrate'],'-r')
ax[2].plot(df_05['analysis_time'],df_05['volumetric_flowrate'],'-g')

plt.savefig(os.path.join('output','plot.png'), dpi=750)





class Data_Plotter:
    '''
    
    '''
    def __init__(self):
        
        pass


    def load_csv_data(self):
        pass


    def plot_flowrates(self):
        pass


    def plot_shearrates(self):
        pass


    def plot_droplet_sizing(self):
        pass


