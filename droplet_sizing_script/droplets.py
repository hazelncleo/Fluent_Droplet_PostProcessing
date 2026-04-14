import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import glob
import os

if True:

    df1 = pd.read_csv('droplet_sizing_script/correct.csv')
    df2 = pd.read_csv('droplet_sizing_script/droplets_3378.csv')

    df1=df1[df1['volume'] < 1e-15]
    df2=df2[df2['volume'] < 1e-15]

    volume1 = df1['volume']*1e18
    vx1 = df1['vx']
    vy1 = df1['vy']
    vz1 = df1['vz']

    volume2 = df2['volume']*1e18
    print(max(volume1),max(volume2))
    vx2 = df2['vx']
    vy2 = df2['vy']
    vz2 = df2['vz']

    velocity1 = np.sqrt(vx1*vx1 + vy1*vy1 + vz1*vz1)
    velocity2 = np.sqrt(vx2*vx2 + vy2*vy2 + vz2*vz2)
    diameter1 = 2 * np.cbrt((3 / (4 * np.pi)) * volume1)
    diameter2 = 2 * np.cbrt((3 / (4 * np.pi)) * volume2)

    f,ax = plt.subplots(1,2, figsize=(8.8,4.2))
    f.tight_layout(pad=2.8)

    ax[0].hist(diameter1,weights=volume1,bins=15,density=True,color='xkcd:sky blue',edgecolor='k')
    sns.kdeplot(x=diameter1, weights=volume1, ax=ax[0],color='r')

    ax[0].set(
        title = 'Correct',
        xlabel = r'Diameter $(\mu m)$',
        ylabel = 'Density'
    )

    ax[1].hist(diameter2,weights=volume2,bins=15,density=True,color='xkcd:sky blue',edgecolor='k')
    sns.kdeplot(x=diameter2, weights=volume2, ax=ax[1],color='r')

    ax[1].set(
        title = 'Wrong',
        xlabel = r'Diameter $(\mu m)$',
        ylabel = 'Density'
    )


    f.savefig(f'droplet_sizing_script/test2.png', dpi=750)
else:

    files = glob.glob('droplet_sizing_script/droplet_data/droplets_*.csv')


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

    ax.hist(diameter,weights=volume,bins=30,density=True,color='xkcd:sky blue',edgecolor='k')
    sns.kdeplot(x=diameter, weights=volume, ax=ax,color='r')

    ax.set(
        title = 'Droplet diameter distribution (Volume)',
        xlabel = r'Diameter $(\mu m)$',
        ylabel = 'Density'
    )




    f.savefig(f'droplet_sizing_script/test.png', dpi=750)