import ansys.pyensight.core as ens
import numpy as np
import os
import matplotlib.pyplot as plt
import glob
from HazelsAwesomeTheme import red_text,green_text,blue_text,yellow_text
import pandas as pd
from tkinter import Tk
from tkinter.filedialog import askdirectory

from droplet_sizing_script_controller import DropletSizingScriptController
from post_processor import PostProcessor
from ensight_controller import EnsightController
from fft_iso import FFT_ISO

class PostProcessor:
    def __init__(self, parameters):


        self.parameters = parameters

        self.droplet_sizer = DropletSizingScriptController()

        

    def post_process(self, options = None): # TODO
        '''
        Data columns:
            time,
            cycle_time,
            max_shearrate,
            total_volume_delivered,
            volumetric_flowrate,
            fpf,
            dv10,
            dv50,
            dv90
        '''

        self.set_iso_view()
        
        self.shearrate_calculation(plot_results = False, animate_results = False)
        
        self.flowrate_calculation(plot_results = False, animate_results = False)

        self.save_data_to_csv()

    
