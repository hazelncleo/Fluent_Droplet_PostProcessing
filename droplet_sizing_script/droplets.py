import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import glob
import os

files = glob.glob('droplet_sizing_script/droplets_*.csv')

for i,f in enumerate(files):
    df = pd.read_csv(f)
    df=df[df['volume'] < 1e-13]
    if i == 0:
        volume = df['volume']*1e18
        ux = df['ux']
        uy = df['uy']
        uz = df['uz']
        vx = df['vx']
        vy = df['vy']
        vz = df['vz']
    else:
        volume = np.append(volume, df['volume']*1e18)
        ux = np.append(ux, df['ux'])
        uy = np.append(ux, df['uy'])
        uz = np.append(ux, df['uz'])
        vx = np.append(ux, df['vx'])
        vy = np.append(ux, df['vy'])
        vz = np.append(ux, df['vz'])

velocity = np.sqrt(vx**2 + vy**2 + vz**2)
diameter = 2 * np.cbrt((3 / (4 * np.pi)) * volume)

f,ax = plt.subplots(1,1, figsize=(5,4.2))
f.tight_layout(pad=2.8)

ax.hist(diameter,weights=volume,bins=15,density=True,color='xkcd:sky blue',edgecolor='k')
sns.kdeplot(x=diameter, weights=volume, ax=ax,color='r')

ax.set(
    title = 'Droplet diameter distribution (Volume)',
    xlabel = r'Diameter $(\mu m)$',
    ylabel = 'Density'
)




f.savefig(f'droplet_sizing_script/test.png', dpi=750)