from tkinter import Tk
from tkinter.filedialog import askdirectory
import pandas as pd
import os
import inquirer
from HazelsAwesomeTheme import HazelsAwesomeTheme, blue_text, yellow_text, red_text, green_text

from droplet_sizing_script_controller import DropletSizingScriptController
from ensight_controller import EnsightController

class PostProcessor:
    '''
    
    '''
    def __init__(self, parameters):
        '''
        Docstring for __init__
        
        :param self: Description
        :param parameters: Description
        '''
        self.parameters = parameters

        self.select_data_folder()

        self.droplet_sizer = DropletSizingScriptController(self.parameters, self.folder)

        self.ensight_controller = EnsightController(self.parameters, self.folder)


    def select_data_folder(self):
        '''
        ---------------------------------------------------
        Prompt the user to select a folder containing the fluent results data to be post-processed (.cas.h5 & .dat.h5)
        ---------------------------------------------------
        RETURNS
        ---------------------------------------------------
        file_path : str
            A string pointing to the folder the user selected. 
        ---------------------------------------------------
        NOTE: A FileExistsError will be raised if the dialog is cancelled
        ---------------------------------------------------
        '''

        # Create window and remove from view
        root = Tk()
        root.iconbitmap('cade.ico')
        root.overrideredirect(1)
        root.geometry('0x0+0+0')
        root.withdraw()
        root.lift()
        root.attributes("-topmost", True)
        
        self.folder = os.path.abspath(askdirectory(title = 'Select folder to read .cas.h5 & .dat.h5 files from: ', initialdir = os.path.abspath(os.getcwd())))
        root.destroy()

        if not self.folder: 
            raise FileExistsError(red_text('A folder was not selected.'))


    def prompt_for_options(self):

        base_choices = [
            'Droplet sizing',
            'General animation',
            'Shearrate',
            'Velocity animation',
            'Flowrate',
            'FFT'
        ]

        further_choices = [
            'Calculation',
            'Plot',
            'Animation'
        ]

        answers = inquirer.prompt(
            questions = [
                inquirer.Checkbox(
                    name     = 'options',
                    message  = 'Select the post-processing options to apply to the data selected',
                    choices  = base_choices,
                    carousel = True
                ),
                inquirer.List(
                    name     = 'droplet_options',
                    message  = 'Select the post-processing type for the ' + blue_text('droplet sizing') + ' script',
                    choices  = further_choices,
                    carousel = True,
                    ignore   = lambda x: 'Droplet sizing' not in x['options'] 
                ),
                inquirer.List(
                    name     = 'shearrate_options',
                    message  = 'Select the post-processing type for the ' + blue_text('shearrate') + ' calculations',
                    choices  = further_choices,
                    carousel = True,
                    ignore   = lambda x: 'Shearrate' not in x['options'] 
                ),
                inquirer.List(
                    name     = 'flowrate_options',
                    message  = 'Select the post-processing type for the ' + blue_text('flowrate') + ' calculations',
                    choices  = further_choices,
                    carousel = True,
                    ignore   = lambda x: 'Flowrate' not in x['options'] 
                )
            ],
            theme = HazelsAwesomeTheme()
        )

        if answers is None:
            raise ValueError(red_text('The options selected were not valid'))
        
        self.options = {choice : (choice in answers['options']) for choice in base_choices}

        if self.options['Droplet sizing']: self.options['Droplet sizing'] = answers['droplet_options']
        if self.options['Shearrate']:      self.options['Shearrate']      = answers['shearrate_options']
        if self.options['Flowrate']:       self.options['Flowrate']       = answers['flowrate_options']


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

        self.prompt_for_options()

        if not any(self.options.values()):
            print(yellow_text('No post-processing options selected, closing program.'))
            return
        
        os.makedirs(os.path.join(self.folder,'output'), exist_ok=True)

        # run droplet sizing fluent script
        if self.options['Droplet sizing']:
            self.droplet_sizer.droplet_sizing_calculation(
                plot_results    = ('Animate' in self.options['Droplet sizing']) or ('Plot' in self.options['Droplet sizing']),
                animate_results = ('Animate' in self.options['Droplet sizing'])
            )

        # Only boot ensight if required
        if any([value for key, value in self.options.items() if key not in 'Droplet sizing']):

            self.ensight_controller.start_ensight()

            self.ensight_controller.set_iso_view()

            # Run post-processing operations
            if self.options['General animation']:  self.ensight_controller.basic_animation()
            if self.options['Velocity animation']: self.ensight_controller.velocity_animation()
            if self.options['FFT']:                self.ensight_controller.fft_of_surface() # FFT should only be run for cases with no droplet formation

            if self.options['Shearrate']:
                self.ensight_controller.shearrate_calculation(
                    plot_results    = ('Animate' in self.options['Shearrate']) or ('Plot' in self.options['Shearrate']), 
                    animate_results = ('Animate' in self.options['Shearrate'])
                )

            if self.options['Flowrate']:
                self.ensight_controller.flowrate_calculation(
                    plot_results    = ('Animate' in self.options['Flowrate']) or ('Plot' in self.options['Flowrate']), 
                    animate_results = ('Animate' in self.options['Flowrate'])
                )

        self.save_data_to_csv()


    def save_data_to_csv(self):
        '''
        Docstring for save_data_to_csv
        
        :param self: Description
        '''

        
        fpath = os.path.join(self.folder, 'output', 'output_data.csv')
        
        if os.path.exists(fpath):
            print(yellow_text('Warning: The file "{}" was overwritten upon saving the .csv output file.'.format(os.path.join('output', 'output_data.csv'))))

        self.ensight_controller.results_data.to_csv(fpath)

    
