import ansys.pyensight.core as ens
import numpy as np
import os
import matplotlib.pyplot as plt
import glob
from fft_iso import FFT_ISO
from HazelsAwesomeTheme import red_text,green_text,blue_text,yellow_text
import pandas as pd
from tkinter import Tk
from tkinter.filedialog import askdirectory


class DropletSizingScriptController:
    def __init__(self):
        pass

    def droplet_sizing_calculation(self, plot_results = False): # TODO
        '''
        Docstring for droplet_sizing
        
        :param self: Description
        '''
        self.droplets_calculated = True


    def droplet_sizing_plot(self): # TODO
        pass

