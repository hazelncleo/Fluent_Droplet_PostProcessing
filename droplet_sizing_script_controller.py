import pandas as pd
import glob
import os
import re

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

        csv_files = glob.glob(os.path.join(self.folder, 'droplets_data', '*.csv'))
        csv_files.sort(key=alphanum_key)

        for droplet_file in csv_files:
            droplet_data = pd.read_csv(droplet_file, skiprows=8)
            print('bruh')




    def droplet_sizing_plot(self): # TODO
        pass


    def droplet_sizing_animation(self): # TODO
        pass



if __name__ == '__main__':
    sizer = DropletSizingScriptController(None, 'test_datasets/droplets')

    sizer.droplet_sizing_calculation()