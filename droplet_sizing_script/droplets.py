import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv('droplet_sizing_script/droplets.csv')

droplet_id = df['droplet_id'][1:]
volume = df['volume'][1:]*1e18
ux = df['ux'][1:]
uy = df['uy'][1:]
uz = df['uz'][1:]
diameter = 2 * np.power((3 / (4 * np.pi)) * volume, 1/3)

f,ax = plt.subplots(2,3, figsize = (15,6))
ax[0,0].hist(volume*1e-9)
ax[0,1].hist(diameter)
ax[0,2].scatter(ux,uy,s=diameter,c='k')
ax[1,0].scatter(ux,uz,s=diameter,c='k')

f.savefig('test.png', dpi=750)