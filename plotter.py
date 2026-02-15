import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.style.use('ggplot')


files = [
    'output/250nm.csv',
    'output/500nm.csv',
    'output/750nm.csv',
    'output/1000nm.csv',
    'output/500nm_noisy.csv',
    'output/1000nm_noisy.csv'
]

data = {
    '250nm' : pd.read_csv(file[0]),
    '500nm' : pd.read_csv(file[1]),
    '750nm' : pd.read_csv(file[2]),
    '1000nm' : pd.read_csv(file[3]),
    '500nm_noisy' : pd.read_csv(file[4]),
    '1000nm_noisy' : pd.read_csv(file[5]),
}


f,(ax_a,ax_n) = plt.subplots(2,3)

ax_a[0].plot(data['250nm']['cycle_time'],data['250nm']['max_shearrate'],'-r')
ax_a[0].plot(data['500nm']['cycle_time'],data['500nm']['max_shearrate'],'-g')
ax_a[0].plot(data['750nm']['cycle_time'],data['750nm']['max_shearrate'],'-b')
ax_a[0].plot(data['1000nm']['cycle_time'],data['1000nm']['max_shearrate'],'-k')

ax_a[1].plot(data['250nm']['cycle_time'],data['250nm']['total_volume_delivered'],'-r')
ax_a[1].plot(data['500nm']['cycle_time'],data['500nm']['total_volume_delivered'],'-g')
ax_a[1].plot(data['750nm']['cycle_time'],data['750nm']['total_volume_delivered'],'-b')
ax_a[1].plot(data['1000nm']['cycle_time'],data['1000nm']['total_volume_delivered'],'-k')

ax_a[2].plot(data['250nm']['cycle_time'],data['250nm']['volumetric_flowrate'],'-r')
ax_a[2].plot(data['500nm']['cycle_time'],data['500nm']['volumetric_flowrate'],'-g')
ax_a[2].plot(data['750nm']['cycle_time'],data['750nm']['volumetric_flowrate'],'-b')
ax_a[2].plot(data['1000nm']['cycle_time'],data['1000nm']['volumetric_flowrate'],'-k')

ax_n[0].plot(data['500nm']['cycle_time'],data['500nm']['max_shearrate'],'-g')
ax_n[0].plot(data['1000nm']['cycle_time'],data['1000nm']['max_shearrate'],'-k')
ax_n[0].plot(data['500nm_noisy']['cycle_time'],data['500nm_noisy']['max_shearrate'],'-y')
ax_n[0].plot(data['1000nm_noisy']['cycle_time'],data['1000nm_noisy']['max_shearrate'],'-c')

ax_n[1].plot(data['500nm']['cycle_time'],data['500nm']['total_volume_delivered'],'-g')
ax_n[1].plot(data['1000nm']['cycle_time'],data['1000nm']['total_volume_delivered'],'-k')
ax_n[1].plot(data['500nm_noisy']['cycle_time'],data['500nm_noisy']['total_volume_delivered'],'-y')
ax_n[1].plot(data['1000nm_noisy']['cycle_time'],data['1000nm_noisy']['total_volume_delivered'],'-c')



plt.savefig('quickplot.png',dpi=750)