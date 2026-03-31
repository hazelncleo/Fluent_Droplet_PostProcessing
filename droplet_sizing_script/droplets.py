import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv('droplet_sizing_script/droplets.csv')

droplet_id = df['droplet_id'][1:]
volume = df['volume'][1:]*1e18
ux = df['ux'][1:]
uy = df['uy'][1:]
uz = df['uz'][1:]
vx = df['vx'][1:]
vy = df['vy'][1:]
vz = df['vz'][1:]
velocity = np.sqrt(vx**2 + vy**2 + vz**2)
diameter = 2 * np.power((3 / (4 * np.pi)) * volume, 1/3)

f,ax = plt.subplots(1,4, figsize = (14,4))
ax[0].hist(diameter,bins=15)
ax[1].hist(diameter,weights=volume,bins=15)
ax[2].hist(velocity,bins=15)
ax[3].hist(velocity[uz < 500e-6])


f.savefig('test.png', dpi=750)