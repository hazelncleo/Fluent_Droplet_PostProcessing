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
        
        
        for i,droplet_file in enumerate(csv_files):

            droplet_data = pd.read_csv(droplet_file, skiprows=8)

            droplet_data['diameter'] = 2 * np.cbrt((3 / (4 * np.pi)) * 1e18 * droplet_data.volume)

            no_large_diameters = droplet_data[droplet_data['diameter'] < 20]

            quantiles = no_large_diameters['diameter'].quantile(q = [0.1,0.5,0.9]).values

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

        plt.plot(qs[:,0],'-k')
        plt.plot(qs[:,1],'-r')
        plt.plot(qs[:,2],'-g')
        plt.savefig('test.png')

        f,ax=plt.subplots(1,1)
        ax.hist(diam, weights=volume)
        f.savefig('hist.png')

            



            






    def droplet_sizing_plot(self): # TODO
        pass


    def droplet_sizing_animation(self): # TODO
        pass



if __name__ == '__main__':
    sizer = DropletSizingScriptController(None, 'test_datasets/droplets')
    sizer.droplet_sizing_calculation()