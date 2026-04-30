import pandas as pd
import glob
import os
import re
import numpy as np
import matplotlib.pyplot as plt

def tryint(s):
    try:
        return int(s)
    except:
        return s

def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]


class DropletSizingScriptController:
    def __init__(self, parameters, folder):
        self.parameters = parameters
        self.folder = folder


    def droplet_sizing_calculation(self, plot_results = False, animate_results = False): # TODO
        '''
        Docstring for droplet_sizing

        :param self: Description
        '''



        csv_files = glob.glob(os.path.join(self.folder, '*.csv'))
        csv_files.sort(key=alphanum_key)

        qs = np.zeros((len(csv_files), 3))
        max_diam = 6
        n_bins = 25

        for i,droplet_file in enumerate(csv_files):

            droplet_data = pd.read_csv(droplet_file, skiprows=8)

            droplet_data['diameter'] = 2 * np.cbrt((3 / (4 * np.pi)) * 1e18 * droplet_data.volume)

            no_large_diameters = droplet_data[droplet_data['diameter'] < max_diam]

            if no_large_diameters.shape[0] >= 1:
                quantiles = np.quantile(no_large_diameters['diameter'], [0.1,0.5,0.9], weights=no_large_diameters['volume'], method='inverted_cdf')
            else:
                quantiles = np.zeros(3)

            if quantiles[0] == quantiles[0]:
                qs[i,0] = quantiles[0]

            if quantiles[1] == quantiles[1]:
                qs[i,1] = quantiles[1]

            if quantiles[2] == quantiles[2]:
                qs[i,2] = quantiles[2]

            if i == 0:
                diam = no_large_diameters['diameter']
                volume = no_large_diameters['volume']
            else:
                diam = np.hstack((diam, no_large_diameters['diameter']))
                volume = np.hstack((volume, no_large_diameters['volume']))

        f,ax=plt.subplots(1,1)
        ax.set_xlim([400, 1500])
        ax.plot(qs[:,0],'-k', label=r'$D_v 10$')
        ax.plot(qs[:,1],'-r', label=r'$D_v 50$')
        ax.plot(qs[:,2],'-g', label=r'$D_v 90$')
        ax.legend(loc='upper right')
        f.savefig('test.png')

        f,ax=plt.subplots(1,1)

        logbins = np.logspace(np.log10(0.7), np.log10(max_diam), n_bins)
        ax.hist(diam, weights=volume, bins = logbins, density=True, color='#7CB9E2')
        ax.set_xscale('log')
        f.savefig('hist.png')












    def droplet_sizing_plot(self): # TODO
        pass


    def droplet_sizing_animation(self): # TODO
        pass



if __name__ == '__main__':
    sizer = DropletSizingScriptController(None, 'test_datasets/droplets')
    sizer.droplet_sizing_calculation()