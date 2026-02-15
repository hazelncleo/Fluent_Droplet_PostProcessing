import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.style.use('ggplot')


files = [
#    'output/250nm.csv',
    'output/500nm.csv',
    'output/750nm.csv',
    'output/1000nm.csv',
    'output/500nm_noisy.csv',
    'output/1000nm_noisy.csv'
]

data = {
#    '250nm' : pd.read_csv(files[0]),
    '500nm' : pd.read_csv(files[0]),
    '750nm' : pd.read_csv(files[1]),
    '1000nm' : pd.read_csv(files[2]),
    '500nm_noisy' : pd.read_csv(files[3]),
    '1000nm_noisy' : pd.read_csv(files[4])
}


f,(ax_a,ax_n) = plt.subplots(2,3, figsize=(18,12))
f.tight_layout(pad=5)

#ax_a[0].plot(data['250nm']['cycle_time'],data['250nm']['max_shearrate'],'-r')
ax_a[0].plot(data['500nm']['cycle_time'],data['500nm']['max_shearrate'],'-g')
ax_a[0].plot(data['750nm']['cycle_time'],data['750nm']['max_shearrate'],'-b')
ax_a[0].plot(data['1000nm']['cycle_time'],data['1000nm']['max_shearrate'],'-k')
ax_a[0].set_yscale('log')
ax_a[0].set(
    title = 'Max Shearrate for different amplitudes',
    xlim = [-0.1,60],
    ylim = [100,1e8],
    xlabel = 'N Cycles',
    ylabel = r'Max Shearrate $(\frac{1}{s})$'
)

#ax_a[1].plot(data['250nm']['cycle_time'],data['250nm']['total_volume_delivered'],'-r')
ax_a[1].plot(data['500nm']['cycle_time'],data['500nm']['total_volume_delivered'],'-g')
ax_a[1].plot(data['750nm']['cycle_time'],data['750nm']['total_volume_delivered'],'-b')
ax_a[1].plot(data['1000nm']['cycle_time'],data['1000nm']['total_volume_delivered'],'-k')
ax_a[1].set(
    title = 'Total volume delivered for different amplitudes',
    xlim = [-0.1,60],
    ylim = [-0.00001,0.00025],
    xlabel = 'N Cycles',
    ylabel = r'Total Volume Delivered $(\mu L)$'
)

#ax_a[2].plot(data['250nm']['cycle_time'],data['250nm']['volumetric_flowrate'],'-r')
ax_a[2].plot(data['500nm']['cycle_time'],data['500nm']['volumetric_flowrate'],'-g', label='500nm')
ax_a[2].plot(data['750nm']['cycle_time'],data['750nm']['volumetric_flowrate'],'-b', label='750nm')
ax_a[2].plot(data['1000nm']['cycle_time'],data['1000nm']['volumetric_flowrate'],'-k', label='1000nm')
ax_a[2].plot([-500,-500],[-250,-250],'-y', label='500nm noisy')
ax_a[2].plot([-500,-500],[-250,-250],'-r', label='1000nm noisy')
ax_a[2].set(
    title = 'Volumetric Flowrate for different amplitudes',
    xlim = [-0.1,60],
    ylim = [-0.1,7],
    xlabel = 'N Cycles',
    ylabel = r'Volumetric Flowrate $(\frac{\mu L}{s})$'
)
ax_a[2].legend(loc='center right', bbox_to_anchor=(0.4, 0.75), fancybox=True, shadow=True)

ax_n[0].plot(data['500nm']['cycle_time'],data['500nm']['max_shearrate'],'-g')
ax_n[0].plot(data['1000nm']['cycle_time'],data['1000nm']['max_shearrate'],'-k')
ax_n[0].plot(data['500nm_noisy']['cycle_time'],data['500nm_noisy']['max_shearrate'],'-y')
ax_n[0].plot(data['1000nm_noisy']['cycle_time'],data['1000nm_noisy']['max_shearrate'],'-r')
ax_n[0].set_yscale('log')
ax_n[0].set(
    title = 'Max shearrate with and without noise',
    xlim = [-0.1,60],
    ylim = [100,1e8],
    xlabel = 'N Cycles',
    ylabel = r'Max Shearrate $(\frac{1}{s})$'
)


ax_n[1].plot(data['1000nm']['cycle_time'],data['1000nm']['total_volume_delivered'],'-k')
ax_n[1].plot(data['1000nm_noisy']['cycle_time'],data['1000nm_noisy']['total_volume_delivered'],'-r')
ax_n[1].set(
    title = 'total volume delivered with and without noise',
    xlim = [-0.1,60],
    ylim = [-0.00001,0.00025],
    xlabel = 'N Cycles',
    ylabel = r'Total Volume Delivered $(\mu L)$'
)

ax_n[2].plot(data['1000nm']['cycle_time'],data['1000nm']['volumetric_flowrate'],'-k')
ax_n[2].plot(data['1000nm_noisy']['cycle_time'],data['1000nm_noisy']['volumetric_flowrate'],'-r')
ax_n[2].set(
    title = 'Volumetric flowrate with and without noise',
    xlim = [0,60],
    ylim = [-0.1,7],
    xlabel = 'N Cycles',
    ylabel = r'Volumetric Flowrate $(\frac{\mu L}{s})$'
)

plt.savefig('output/quickplot.png',dpi=750)

class Plotter:
    def __init__(self):
        pass

    def plot_data(self):
        pass
    