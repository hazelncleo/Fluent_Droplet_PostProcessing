from tkinter import Tk
from tkinter.filedialog import askdirectory
import pandas as pd
import os
import inquirer
from HazelsAwesomeTheme import HazelsAwesomeTheme, blue_text, yellow_text, red_text, green_text

from droplet_sizing_script_controller import DropletSizingScriptController
from ensight_controller import EnsightController

class ControlPostProcessor:
    '''

    '''
    def __init__(self, parameters, folder = None, options = None):
        '''
        Docstring for __init__

        :param self: Description
        :param parameters: Description
        '''
        self.parameters = parameters
        self.options    = options

        if folder:
            self.folder = folder
        else:
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

        self.folder = askdirectory(title = 'Select folder to read .cas.h5 & .dat.h5 files from: ', initialdir = os.path.abspath(os.getcwd()))
        root.destroy()

        if not self.folder:
            raise FileExistsError(red_text('A folder was not selected.'))

        self.folder = os.path.abspath(self.folder)


    def prompt_for_options(self):

        postprocessing_choices = [
            'Droplet sizing',
            'General animation',
            'Shearrate',
            'Velocity animation',
            'Flowrate',
            'FFT'
        ]

        type_of_postprocessing = [
            'Calculation',
            'Plot',
            'Animation'
        ]

        selected_answers = inquirer.prompt(
            questions = [
                inquirer.Checkbox(
                    name     = 'options',
                    message  = 'Select the post-processing options to apply to the data selected',
                    choices  = postprocessing_choices,
                    carousel = True
                ),
                inquirer.List(
                    name     = 'droplet_options',
                    message  = 'Select the post-processing type for the ' + blue_text('droplet sizing') + ' script',
                    choices  = type_of_postprocessing,
                    carousel = True,
                    ignore   = lambda x: 'Droplet sizing' not in x['options']
                ),
                inquirer.List(
                    name     = 'shearrate_options',
                    message  = 'Select the post-processing type for the ' + blue_text('shearrate') + ' calculations',
                    choices  = type_of_postprocessing,
                    carousel = True,
                    ignore   = lambda x: 'Shearrate' not in x['options']
                ),
                inquirer.List(
                    name     = 'flowrate_options',
                    message  = 'Select the post-processing type for the ' + blue_text('flowrate') + ' calculations',
                    choices  = type_of_postprocessing,
                    carousel = True,
                    ignore   = lambda x: 'Flowrate' not in x['options']
                )
            ],
            theme = HazelsAwesomeTheme()
        )

        if selected_answers is None:
            raise ValueError(red_text('The options selected were not valid'))

        self.options = {choice : (choice in selected_answers['options']) for choice in postprocessing_choices}

        if self.options['Droplet sizing']: self.options['Droplet sizing'] = selected_answers['droplet_options']
        if self.options['Shearrate']:      self.options['Shearrate']      = selected_answers['shearrate_options']
        if self.options['Flowrate']:       self.options['Flowrate']       = selected_answers['flowrate_options']


    def post_process(self):
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
        if self.options is None:
            self.prompt_for_options()
        else:
            print('Options specified by file.')

        if not any(self.options.values()):
            print(yellow_text('No post-processing options selected, closing program.'))
            return

        os.makedirs(os.path.join(self.folder,'output'), exist_ok = True)

        # run droplet sizing fluent script
        if self.options['Droplet sizing']:
            self.droplet_sizer.droplet_sizing_calculation(
                plot_results    = ('Animation' in self.options['Droplet sizing']) or ('Plot' in self.options['Droplet sizing']),
                animate_results = ('Animation' in self.options['Droplet sizing'])
            )

            self.save_droplet_data_to_csv()

        # Only boot ensight if required
        if any([value for key, value in self.options.items() if key not in 'Droplet sizing']):

            self.ensight_controller.start_ensight()

            self.ensight_controller.set_iso_view()

            # Run post-processing operations
            if self.options['General animation']:  self.ensight_controller.basic_animation()
            if self.options['Velocity animation']: self.ensight_controller.velocity_animation()
            if self.options['FFT']:                self.ensight_controller.fft_of_surface(plot_results=True) # FFT should only be run for cases with no droplet formation

            if self.options['Shearrate']:
                self.ensight_controller.shearrate_calculation(
                    plot_results    = ('Animation' in self.options['Shearrate']) or ('Plot' in self.options['Shearrate']),
                    animate_results = ('Animation' in self.options['Shearrate'])
                )

            if self.options['Flowrate']:
                self.ensight_controller.flowrate_calculation(
                    plot_results    = ('Animation' in self.options['Flowrate']) or ('Plot' in self.options['Flowrate']),
                    animate_results = ('Animation' in self.options['Flowrate'])
                )

            self.save_ensight_data_to_csv()


    def save_droplet_data_to_csv(self):
        '''
        Docstring for save_data_to_csv

        :param self: Description
        '''

        fpath = os.path.join(self.folder, 'output', 'cumulative_droplet_data.csv')

        if os.path.exists(fpath):
            print(yellow_text('Warning: The file "{}" was overwritten upon saving the .csv output file.'.format(os.path.join('output', 'cumulative_droplet_data.csv'))))

        self.droplet_sizer.results_data.to_csv(fpath)

        fpath = os.path.join(self.folder, 'output', 'individual_droplet_data.csv')

        if os.path.exists(fpath):
            print(yellow_text('Warning: The file "{}" was overwritten upon saving the .csv output file.'.format(os.path.join('output', 'individual_droplet_data.csv'))))

        self.droplet_sizer.individual_droplet_data.to_csv(fpath)


    def save_ensight_data_to_csv(self):
        '''
        Docstring for save_data_to_csv

        :param self: Description
        '''

        fpath = os.path.join(self.folder, 'output', 'ensight_data.csv')

        if os.path.exists(fpath):
            print(yellow_text('Warning: The file "{}" was overwritten upon saving the .csv output file.'.format(os.path.join('output', 'ensight_data.csv'))))

            temp_data = pd.read_csv(fpath, index_col = 'timestep_number')

            for column in self.ensight_controller.results_data:
                if column not in temp_data:
                    temp_data[column] = self.ensight_controller.results_data[column]

            temp_data.to_csv(fpath)


        else:
            self.ensight_controller.results_data.to_csv(fpath)


