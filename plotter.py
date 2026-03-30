import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

plt.style.use('ggplot')


files = [
    '250nm',
    '500nm',
    '750nm',
    '1000nm',
    '500nm_noisy',
    '1000nm_noisy',
    '1000nm_fullduty',
    '1000nm_noisy_fullduty'
]

data = {}

for fpath in files:
    data[fpath] = pd.read_csv(os.path.join('output',fpath+'.csv'))


freq_1 = 1.63
freq_2 = 0.14
A = 1

t = np.linspace(0, 15/freq_1, 500)
x_t = A * (1 - np.cos(2 * np.pi * freq_1 * t))
noise_t = A / 10 * (1 - np.cos(2 * np.pi * freq_2 * t))

f,ax = plt.subplots(1,1, figsize = (9,5))
f.tight_layout(pad=5)

ax.plot(t,x_t,'-k', label='Rigid vibration')
ax.plot(t,x_t+noise_t,'-r', label='Rigid vibration & noise')

f.suptitle('Vibrations used')

ax.set(
    xlabel = r'time $(\mu s)$',
    ylabel = r'$z$ displacement $(\mu m)$',
    ylim = [-0.1,2.8]
)

ax.legend(loc='upper right', fancybox=True, shadow=True)

f.savefig('output/Vibrations.png',dpi=750)




f,ax = plt.subplots(1,3, figsize=(15,5))
f.tight_layout(pad=5)
f.suptitle('Vibration amplitudes effect on values of interest')

ax[0].plot(data['250nm']['cycle_time'],data['250nm']['max_shearrate'].rolling(window=15,min_periods=1,center=True).mean(),'-r')
ax[0].plot(data['500nm']['cycle_time'],data['500nm']['max_shearrate'].rolling(window=15,min_periods=1,center=True).mean(),'-g')
ax[0].plot(data['750nm']['cycle_time'],data['750nm']['max_shearrate'].rolling(window=15,min_periods=1,center=True).mean(),'-b')
ax[0].plot(data['1000nm_fullduty']['cycle_time'],data['1000nm_fullduty']['max_shearrate'].rolling(window=15,min_periods=1,center=True).mean(),'-k')
ax[0].set_yscale('log')
ax[0].set(
    title = 'Max Shearrate for different amplitudes',
    xlim = [-0.1,60],
    ylim = [1e5,1e8],
    xlabel = 'N Cycles',
    ylabel = r'Max Shearrate $(\frac{1}{s})$'
)

ax[1].plot(data['250nm']['cycle_time'],data['250nm']['total_volume_delivered'],'-r')
ax[1].plot(data['500nm']['cycle_time'],data['500nm']['total_volume_delivered'],'-g')
ax[1].plot(data['750nm']['cycle_time'],data['750nm']['total_volume_delivered'],'-b')
ax[1].plot(data['1000nm_fullduty']['cycle_time'],data['1000nm_fullduty']['total_volume_delivered'],'-k')
ax[1].set(
    title = 'Total volume delivered for different amplitudes',
    xlim = [-0.1,60],
    ylim = [-0.00001,0.0003],
    xlabel = 'N Cycles',
    ylabel = r'Total Volume Delivered $(\mu L)$'
)

ax[2].plot(data['250nm']['cycle_time'],data['250nm']['volumetric_flowrate'],'-r', label='250nm')
ax[2].plot(data['500nm']['cycle_time'],data['500nm']['volumetric_flowrate'],'-g', label='500nm')
ax[2].plot(data['750nm']['cycle_time'],data['750nm']['volumetric_flowrate'],'-b', label='750nm')
ax[2].plot(data['1000nm_fullduty']['cycle_time'],data['1000nm_fullduty']['volumetric_flowrate'],'-k', label='1000nm')
ax[2].set(
    title = 'Volumetric Flowrate for different amplitudes',
    xlim = [-0.1,60],
    ylim = [-0.1,6],
    xlabel = 'N Cycles',
    ylabel = r'Volumetric Flowrate $(\frac{\mu L}{s})$'
)
ax[2].legend(loc='upper left', fancybox=True, shadow=True)

f.savefig('output/Amplitudes.png',dpi=750)


f,ax = plt.subplots(1,3, figsize=(15,5))
f.tight_layout(pad=5)
f.suptitle('Effect of adding low amplitude 140kHz noise.')

ax[0].plot(data['500nm']['cycle_time'],data['500nm']['max_shearrate'].rolling(window=15,min_periods=1,center=True).mean(),'-r',label='500nm')
ax[0].plot(data['500nm_noisy']['cycle_time'],data['500nm_noisy']['max_shearrate'].rolling(window=15,min_periods=1,center=True).mean(),'-g',label='500nm with noise')
ax[0].plot(data['1000nm_fullduty']['cycle_time'],data['1000nm_fullduty']['max_shearrate'].rolling(window=15,min_periods=1,center=True).mean(),'-b',label='1000nm')
ax[0].plot(data['1000nm_noisy_fullduty']['cycle_time'],data['1000nm_noisy_fullduty']['max_shearrate'].rolling(window=15,min_periods=1,center=True).mean(),'-k',label='1000nm with noise')
ax[0].set_yscale('log')
ax[0].set(
    title = 'Max shearrate with and without noise',
    xlim = [-0.1,60],
    ylim = [1e5,1e8],
    xlabel = 'N Cycles',
    ylabel = r'Max Shearrate $(\frac{1}{s})$'
)

ax[1].plot(data['1000nm_fullduty']['cycle_time'],data['1000nm_fullduty']['total_volume_delivered'],'-b')
ax[1].plot(data['1000nm_noisy_fullduty']['cycle_time'],data['1000nm_noisy_fullduty']['total_volume_delivered'],'-k')
ax[1].set(
    title = 'total volume delivered with and without noise',
    xlim = [-0.1,60],
    ylim = [-0.00001,0.0003],
    xlabel = 'N Cycles',
    ylabel = r'Total Volume Delivered $(\mu L)$'
)

ax[2].plot(data['1000nm_fullduty']['cycle_time'],data['1000nm_fullduty']['volumetric_flowrate'],'-b')
ax[2].plot(data['1000nm_noisy_fullduty']['cycle_time'],data['1000nm_noisy_fullduty']['volumetric_flowrate'],'-k')
ax[2].set(
    title = 'Volumetric flowrate with and without noise',
    xlim = [-0.1,60],
    ylim = [-0.1,6],
    xlabel = 'N Cycles',
    ylabel = r'Volumetric Flowrate $(\frac{\mu L}{s})$'
)

ax[0].legend(loc='lower right', fancybox=True, shadow=True)

f.savefig('output/Noise.png',dpi=750)
















if False:












    f,(ax_a,ax_n) = plt.subplots(2,3, figsize=(18,12))
    f.tight_layout(pad=5)

    #ax_a[0].plot(data['250nm']['cycle_time'],data['250nm']['max_shearrate'],'-r')
    ax_a[0].plot(data['500nm']['cycle_time'],data['500nm']['max_shearrate'].rolling(window=15,min_periods=1,center=True).mean(),'-g')
    ax_a[0].plot(data['750nm']['cycle_time'],data['750nm']['max_shearrate'].rolling(window=15,min_periods=1,center=True).mean(),'-b')
    ax_a[0].plot(data['1000nm']['cycle_time'],data['1000nm']['max_shearrate'].rolling(window=15,min_periods=1,center=True).mean(),'-k')
    ax_a[0].plot(data['1000nm_fullduty']['cycle_time'],data['1000nm_fullduty']['max_shearrate'].rolling(window=15,min_periods=1,center=True).mean(),'-c')
    ax_a[0].set_yscale('log')
    ax_a[0].set(
        title = 'Max Shearrate for different amplitudes',
        xlim = [-0.1,200],
        ylim = [1e4,1e8],
        xlabel = 'N Cycles',
        ylabel = r'Max Shearrate $(\frac{1}{s})$'
    )

    #ax_a[1].plot(data['250nm']['cycle_time'],data['250nm']['total_volume_delivered'],'-r')
    ax_a[1].plot(data['500nm']['cycle_time'],data['500nm']['total_volume_delivered'],'-g')
    ax_a[1].plot(data['750nm']['cycle_time'],data['750nm']['total_volume_delivered'],'-b')
    ax_a[1].plot(data['1000nm']['cycle_time'],data['1000nm']['total_volume_delivered'],'-k')
    ax_a[1].plot(data['1000nm_fullduty']['cycle_time'],data['1000nm_fullduty']['total_volume_delivered'],'-c')
    ax_a[1].set(
        title = 'Total volume delivered for different amplitudes',
        xlim = [-0.1,200],
        ylim = [-0.00001,0.001],
        xlabel = 'N Cycles',
        ylabel = r'Total Volume Delivered $(\mu L)$'
    )

    #ax_a[2].plot(data['250nm']['cycle_time'],data['250nm']['volumetric_flowrate'],'-r')
    ax_a[2].plot(data['500nm']['cycle_time'],data['500nm']['volumetric_flowrate'],'-g', label='500nm')
    ax_a[2].plot(data['750nm']['cycle_time'],data['750nm']['volumetric_flowrate'],'-b', label='750nm')
    ax_a[2].plot(data['1000nm']['cycle_time'],data['1000nm']['volumetric_flowrate'],'-k', label='1000nm')
    ax_a[2].plot(data['1000nm_fullduty']['cycle_time'],data['1000nm_fullduty']['volumetric_flowrate'],'-c', label='1000nm full dutycycle')
    ax_a[2].plot([-500,-500],[-250,-250],'-y', label='500nm noisy')
    ax_a[2].plot([-500,-500],[-250,-250],'-r', label='1000nm noisy')
    ax_a[2].plot([-500,-500],[-250,-250],'-m', label='1000nm noisy fullduty')
    ax_a[2].set(
        title = 'Volumetric Flowrate for different amplitudes',
        xlim = [-0.1,200],
        ylim = [-0.1,8],
        xlabel = 'N Cycles',
        ylabel = r'Volumetric Flowrate $(\frac{\mu L}{s})$'
    )
    ax_a[2].legend(loc='lower right', fancybox=True, shadow=True)

    ax_n[0].plot(data['500nm']['cycle_time'],data['500nm']['max_shearrate'].rolling(window=15,min_periods=1,center=True).mean(),'-g')
    ax_n[0].plot(data['1000nm']['cycle_time'],data['1000nm']['max_shearrate'].rolling(window=15,min_periods=1,center=True).mean(),'-k')
    ax_n[0].plot(data['500nm_noisy']['cycle_time'],data['500nm_noisy']['max_shearrate'].rolling(window=15,min_periods=1,center=True).mean(),'-y')
    ax_n[0].plot(data['1000nm_noisy']['cycle_time'],data['1000nm_noisy']['max_shearrate'].rolling(window=15,min_periods=1,center=True).mean(),'-r')
    ax_n[0].plot(data['1000nm_fullduty']['cycle_time'],data['1000nm_fullduty']['max_shearrate'].rolling(window=15,min_periods=1,center=True).mean(),'-c')
    ax_n[0].plot(data['1000nm_noisy_fullduty']['cycle_time'],data['1000nm_noisy_fullduty']['max_shearrate'].rolling(window=15,min_periods=1,center=True).mean(),'-m')
    ax_n[0].set_yscale('log')
    ax_n[0].set(
        title = 'Max shearrate with and without noise',
        xlim = [-0.1,200],
        ylim = [1e4,1e8],
        xlabel = 'N Cycles',
        ylabel = r'Max Shearrate $(\frac{1}{s})$'
    )


    ax_n[1].plot(data['1000nm']['cycle_time'],data['1000nm']['total_volume_delivered'],'-k')
    ax_n[1].plot(data['1000nm_noisy']['cycle_time'],data['1000nm_noisy']['total_volume_delivered'],'-r')
    ax_n[1].plot(data['1000nm_fullduty']['cycle_time'],data['1000nm_fullduty']['total_volume_delivered'],'-c')
    ax_n[1].plot(data['1000nm_noisy_fullduty']['cycle_time'],data['1000nm_noisy_fullduty']['total_volume_delivered'],'-m')
    ax_n[1].set(
        title = 'total volume delivered with and without noise',
        xlim = [-0.1,200],
        ylim = [-0.00001,0.001],
        xlabel = 'N Cycles',
        ylabel = r'Total Volume Delivered $(\mu L)$'
    )

    ax_n[2].plot(data['1000nm']['cycle_time'],data['1000nm']['volumetric_flowrate'],'-k')
    ax_n[2].plot(data['1000nm_noisy']['cycle_time'],data['1000nm_noisy']['volumetric_flowrate'],'-r')
    ax_n[2].plot(data['1000nm_fullduty']['cycle_time'],data['1000nm_fullduty']['volumetric_flowrate'],'-c')
    ax_n[2].plot(data['1000nm_noisy_fullduty']['cycle_time'],data['1000nm_noisy_fullduty']['volumetric_flowrate'],'-m')
    ax_n[2].set(
        title = 'Volumetric flowrate with and without noise',
        xlim = [0,200],
        ylim = [-0.1,8],
        xlabel = 'N Cycles',
        ylabel = r'Volumetric Flowrate $(\frac{\mu L}{s})$'
    )

    plt.savefig('output/quickplot.png',dpi=750)

    class Plotter:
        def __init__(self):
            pass

        def plot_data(self):
            pass
        