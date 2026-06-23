import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fftfreq, fftshift, fftn
from scipy.interpolate import LinearNDInterpolator, make_interp_spline
from skimage.filters import window
from seaborn import color_palette, dark_palette
import pandas as pd
import os


def multiply_along_axis(A, B, axis):
    return np.swapaxes(np.swapaxes(A, axis, -1) * B, -1, axis)


class FFT_ISO:
    def __init__(
        self,
        parameters = {
            'vibration_frequency' : 1.63e6,
            'vibration_amplitude' : 1e-6,
            'n_cycles'            : 60,
            'n_elements'          : 35,
            'channel_width'       : 50,
            'grid_size'           : 500
        },
        n_timesteps = 40,
        times       = [0, 2.5e-5],
        use_ggplot_styling = True
    ):
        '''
        Docstring for __init__

        :param self: Description
        :param data: Description
        :param file: Description
        '''

        if use_ggplot_styling:
            plt.style.use('ggplot')

        self.cmap_psd    = color_palette('rocket', as_cmap = True)
        self.cmap_deform = dark_palette("#69d", as_cmap = True)
        self.parameters  = parameters
        self.times       = times

        self.all_raw_data_arrays = {}
        self.centred_data = {}
        self.interpolated_data = {}
        self.windowed_data = {}
        self.fft_data = {}
        self.PSD = {}
        self.wavelengths = {}
        self.temporal_data = {}

        # Calculate sampling data from provided parameters TODO Parameter handling
        self.sampling_data = {
            'elements_along_width'  : self.parameters['n_elements'],
            'elements_along_length' : int(np.ceil(self.parameters['n_elements'] * self.parameters['grid_size'] / self.parameters['channel_width'])),
            'n_samples_time'        : 1 * n_timesteps
        }

        self.sampling_data.update({
            'spatial_sample_spacing'  : self.parameters['grid_size'] / self.sampling_data['elements_along_length'],
            'temporal_sample_spacing' : (self.times[1] - self.times[0]) / (self.sampling_data['n_samples_time'] - 1),
            'n_samples_width'         : self.sampling_data['elements_along_width'] + 1,
            'n_samples_length'        : self.sampling_data['elements_along_length'] + 1
        })

        self.data_loaded = False

        self.create_meshes()

        self.create_window_functions()


    def calculate_theory_wavelength(self) -> float:
        ''' Calculate the theoretical surface wavelength in micron '''
        return 1e6 * np.power((2 * np.pi * 0.072) / (997 * (self.parameters['vibration_frequency'] / 2)**2), 1/3)


    def create_meshes(self):
        '''
        Create both spatial and frequency domain meshes
        '''

        # Get 1D spatial and frequency domain meshes
        self.meshes = {
            'channel_width' : np.linspace(
                start = -int(self.parameters['channel_width'] / 2),
                stop  =  int(self.parameters['channel_width'] / 2),
                num   =  self.sampling_data['n_samples_width']
            ),
            'channel_length' : np.linspace(
                start = -int(self.parameters['grid_size'] / 2),
                stop  =  int(self.parameters['grid_size'] / 2),
                num   =  self.sampling_data['n_samples_length']
            ),
            'time' : np.linspace(
                start = self.times[0],
                stop  = self.times[1],
                num   = self.sampling_data['n_samples_time']
            )
        }

        self.meshes['frequency']                = fftshift(fftfreq(self.sampling_data['n_samples_time'], self.sampling_data['temporal_sample_spacing']))
        self.meshes['channel_width_frequency']  = fftshift(fftfreq(self.sampling_data['n_samples_width'], self.sampling_data['spatial_sample_spacing']))
        self.meshes['channel_length_frequency'] = fftshift(fftfreq(self.sampling_data['n_samples_length'], self.sampling_data['spatial_sample_spacing']))

        self.meshes['x_mesh'], self.meshes['y_mesh'] = np.meshgrid(self.meshes['channel_width'],  self.meshes['channel_length'])

        self.meshes['x_frequency_mesh'], self.meshes['y_frequency_mesh'] = np.meshgrid(self.meshes['channel_width_frequency'],  self.meshes['channel_length_frequency'])


    def create_window_functions(self):
        '''
        Create 2D window functions
        '''

        self.spatial_window = window(
            window_type = ('tukey', 0.65),
            shape = (self.sampling_data['n_samples_length'], self.sampling_data['n_samples_width'])
        )

        self.temporal_window = window(('tukey', 0.4), self.sampling_data['n_samples_time'])


    def send_data(self, data, time, index):

        self.all_raw_data_arrays[index] = {
            'time'  : time,
            'index' : index,
            'data'  : data * 1e6
        }


    def solve(self):
        '''
        Solve the fft of the surface data
        '''

        self.centre_data()

        self.interpolate_data()

        self.calculate_ffts()

        self.calculate_PSDs()

        self.calculate_normed_wavelengths()

        self.data_loaded = True


    def centre_data(self):
        '''
        Docstring for centre_data

        :param self: Description
        '''

        for index in self.all_raw_data_arrays:

            # Centre z data
            self.centred_data[index] = {
                'time'  : self.all_raw_data_arrays[index]['time'],
                'index' : index,
                'mean'  : np.mean(self.all_raw_data_arrays[index]['data'][:,2]),
                'data'  : self.all_raw_data_arrays[index]['data']
            }

            self.centred_data[index]['data'][:,2] = self.centred_data[index]['data'][:,2] - self.centred_data[index]['mean']

        self.temporal_data['times'] = np.array([self.centred_data[index]['time'] for index in self.centred_data])
        self.temporal_data['sort_index'] = np.argsort(self.temporal_data['times'])
        self.temporal_data['times'] = self.temporal_data['times'][self.temporal_data['sort_index']]

        self.temporal_data['raw_data'] = np.array([self.centred_data[index]['mean'] for index in self.centred_data])[self.temporal_data['sort_index']]
        self.temporal_data['centred'] = self.temporal_data['raw_data'] - np.mean(self.temporal_data['raw_data'])


    def interpolate_data(self):
        '''
        Interpolate the raw data onto a cartesian grid.
        '''

        self.temporal_data['bulk_data'] = np.empty((*self.meshes['x_mesh'].shape, self.sampling_data['n_samples_time']))

        for index in self.centred_data:

            interpolator = LinearNDInterpolator(self.centred_data[index]['data'][:,:2], self.centred_data[index]['data'][:,2], fill_value = 0)

            self.interpolated_data[index] = {
                'time'  : self.centred_data[index]['time'],
                'index' : index,
                'data'  : interpolator(self.meshes['x_mesh'], self.meshes['y_mesh'])
            }

            self.temporal_data['bulk_data'][:,:,index - 1] = self.interpolated_data[index]['data']

            self.windowed_data[index] = {
                'time'  : self.centred_data[index]['time'],
                'index' : index,
                'data'  : self.interpolated_data[index]['data'] * self.spatial_window
            }

        self.temporal_data['interpolated'] = np.interp(self.meshes['time'], self.temporal_data['times'], self.temporal_data['centred'])

        temporal_interpolator = make_interp_spline(self.temporal_data['times'], self.temporal_data['bulk_data'], k = 1, axis = 2)
        self.temporal_data['bulk_data'] = temporal_interpolator(self.meshes['time'])

        self.temporal_data['bulk_data'] = self.temporal_data['bulk_data'] * self.temporal_window

        self.temporal_data['windowed'] = self.temporal_data['interpolated'] * self.temporal_window


    def calculate_ffts(self):
        '''

        '''

        for index in self.windowed_data:

            self.fft_data[index] = {
                'time'  : self.windowed_data[index]['time'],
                'index' : index,
                'data'  : fftshift(fftn(self.windowed_data[index]['data'], workers = -1))
            }
            self.fft_data[index]['data_magnitude'] = np.abs(self.fft_data[index]['data'])

        self.temporal_data['fft'] = fftshift(fftn(self.temporal_data['windowed']))

        self.temporal_data['fft_bulk'] = fftshift(fftn(self.temporal_data['bulk_data'], axes = 2))


    def calculate_PSDs(self):

        for index in self.fft_data:

            self.PSD[index] = {
                'time'  : self.fft_data[index]['time'],
                'index' : index,
                'data'  : np.log(np.abs(self.fft_data[index]['data'])**2)
            }

            self.PSD[index]['flat_data'] = self.PSD[index]['data'].flatten()

        self.mean_fft      = np.mean(np.stack([self.fft_data[index]['data_magnitude'] for index in self.fft_data], axis = 2), axis = 2)
        self.mean_PSD      = np.log(self.mean_fft**2)
        self.mean_flat_PSD = self.mean_PSD.flatten()

        self.temporal_data['PSD'] = np.log(np.abs(self.temporal_data['fft'])**2)
        self.temporal_data['mean_bulk'] = np.mean(np.abs(self.temporal_data['fft_bulk']), axis = (0,1))
        self.temporal_data['mean_PSD'] = np.log(np.abs(self.temporal_data['mean_bulk']**2))


    def calculate_normed_wavelengths(self):
        '''
        Calculate the normed wavelengths from the frequency data.
        '''

        self.norms = {
            'frequency' : np.sqrt(self.meshes['x_frequency_mesh']**2 + self.meshes['y_frequency_mesh']**2)
        }

        self.norms['flat_frequency'] = self.norms['frequency'].flatten()

        self.calculate_masks()

        self.wavelengths = np.reciprocal(self.norms['flat_frequency'][self.masks['remove_zero']])

        self.ranged_wavelengths = self.wavelengths[self.masks['range']]



    def calculate_masks(self):
        '''
        Calculate Masks to filter the data for plotting
        '''

        self.masks = {
            'remove_zero'        : self.norms['flat_frequency'] != 0,
            'remove_large_freqs' : self.norms['flat_frequency'] > 0.02
        }

        self.masks['range'] = self.masks['remove_large_freqs'][self.masks['remove_zero']]


    def angles_from_frequency_data(self):
        '''
        TODO: Look at angles that waves are occurring in
        '''
        pass

    '''
    ------------------
    PLOTTING & OUTPUT DATA
    ------------------
    '''

    def full_plot(self, title = None, file_name = None, display = False):
        '''
        Docstring for full_plot

        :param self: Description
        :param title: Description
        :param file_name: Description
        :param display: Description
        '''

        # Create plot objects
        f,ax = plt.subplots(
            nrows   = 2,
            ncols   = 4,
            layout  = 'constrained',
            figsize = (22,13)
        )

        f.get_layout_engine().set(
            w_pad  = 0.2,
            h_pad  = 0.1,
            hspace = 0.05,
            wspace = 0
        )

        if title:
            f.suptitle(title, size = 50)

        # Set Styling
        ax[0,0].set(
            xlim   = [-int(self.parameters['grid_size'] / 2), int(self.parameters['grid_size'] / 2)],
            ylim   = [-int(self.parameters['grid_size'] / 2), int(self.parameters['grid_size'] / 2)],
            aspect = 'equal'
        )

        ax[0,1].set(
            xlim   = [-int(self.parameters['grid_size'] / 2), int(self.parameters['grid_size'] / 2)],
            ylim   = [-int(self.parameters['grid_size'] / 2), int(self.parameters['grid_size'] / 2)],
            aspect = 'equal'
        )

        ax[0,2].set(
            xlim   = [-int(self.parameters['grid_size'] / 2), int(self.parameters['grid_size'] / 2)],
            ylim   = [-int(self.parameters['grid_size'] / 2), int(self.parameters['grid_size'] / 2)],
            aspect = 'equal'
        )

        ax[0,3].set(
            xlim   = [-int(self.parameters['grid_size'] / 2), int(self.parameters['grid_size'] / 2)],
            ylim   = [-int(self.parameters['grid_size'] / 2), int(self.parameters['grid_size'] / 2)],
            aspect = 'equal'
        )

        ax[1,0].set(
            xlim   = [np.min(self.meshes['vertical_y_frequency_mesh']), np.max(self.meshes['vertical_y_frequency_mesh'])],
            ylim   = [np.min(self.meshes['vertical_y_frequency_mesh']), np.max(self.meshes['vertical_y_frequency_mesh'])],
            aspect = 'equal'
        )

        ax[1,1].set(
            xlim   = [np.min(self.meshes['vertical_y_frequency_mesh']), np.max(self.meshes['vertical_y_frequency_mesh'])],
            ylim   = [np.min(self.meshes['vertical_y_frequency_mesh']), np.max(self.meshes['vertical_y_frequency_mesh'])],
            aspect = 'equal'
        )

        ax[1,2].set_aspect(1.75)
        ax[1,3].set_aspect(1.75)

        ax[0,0].set_title('Raw Data for Vertical Section', fontsize = 12)
        ax[0,1].set_title('Raw Data for Horizontal Section', fontsize = 12)
        ax[0,2].set_title('Windowed Vertical Data', fontsize = 12)
        ax[0,3].set_title('Windowed Horizontal Data', fontsize = 12)
        ax[1,0].set_title('Logged PSD, Vertical', fontsize = 12)
        ax[1,1].set_title('Logged PSD, Horizontal', fontsize = 12)
        ax[1,2].set_title('Normed Wavelength Powers, Vertical', fontsize = 12)
        ax[1,3].set_title('Normed Wavelength Powers, Horizontal', fontsize = 12)

        ax[0,0].set_xlabel(r'Position $(\mu m)$', fontsize = 9)
        ax[0,1].set_xlabel(r'Position $(\mu m)$', fontsize = 9)
        ax[0,2].set_xlabel(r'Position $(\mu m)$', fontsize = 9)
        ax[0,3].set_xlabel(r'Position $(\mu m)$', fontsize = 9)
        ax[1,0].set_xlabel(r'Spatial Frequency X $(\frac{1}{\mu m})$', fontsize = 9)
        ax[1,1].set_xlabel(r'Spatial Frequency X $(\frac{1}{\mu m})$', fontsize = 9)
        ax[1,2].set_xlabel(r'Wavelength $(\mu m)$', fontsize = 9)
        ax[1,3].set_xlabel(r'Wavelength $(\mu m)$', fontsize = 9)

        ax[0,0].set_ylabel(r'Position $(\mu m)$', fontsize = 9)
        ax[0,1].set_ylabel(r'Position $(\mu m)$', fontsize = 9)
        ax[0,2].set_ylabel(r'Position $(\mu m)$', fontsize = 9)
        ax[0,3].set_ylabel(r'Position $(\mu m)$', fontsize = 9)
        ax[1,0].set_ylabel(r'Spatial Frequency Y $(\frac{1}{\mu m})$', fontsize = 9)
        ax[1,1].set_ylabel(r'Spatial Frequency Y $(\frac{1}{\mu m})$', fontsize = 9)
        ax[1,2].set_ylabel(r'Logged PSD', fontsize = 9)
        ax[1,3].set_ylabel(r'Logged PSD', fontsize = 9)

        # Calculate colorbar ranges
        height_cmap_min = min(np.min(self.data['vertical']), np.min(self.data['horizontal']))
        height_cmap_max = max(np.max(self.data['vertical']), np.max(self.data['horizontal']))
        PSD_cmap_min    = min(np.min(self.PSD['vertical_windowed']), np.min(self.PSD['horizontal_windowed']))
        PSD_cmap_max    = max(np.max(self.PSD['vertical_windowed']), np.max(self.PSD['horizontal_windowed']))


        ax[0,0].pcolormesh(
            self.meshes['vertical_x_mesh'],
            self.meshes['vertical_y_mesh'],
            self.data['vertical'],
            cmap = self.cmap_deform,
            vmin = height_cmap_min,
            vmax = height_cmap_max
        )

        ax[0,1].pcolormesh(
            self.meshes['horizontal_x_mesh'],
            self.meshes['horizontal_y_mesh'],
            self.data['horizontal'],
            cmap = self.cmap_deform,
            vmin = height_cmap_min,
            vmax = height_cmap_max
        )

        ax[0,2].pcolormesh(
            self.meshes['vertical_x_mesh'],
            self.meshes['vertical_y_mesh'],
            self.data['vertical_windowed'],
            cmap = self.cmap_deform,
            vmin = height_cmap_min,
            vmax = height_cmap_max
        )

        display_height_cmap = ax[0,3].pcolormesh(
            self.meshes['horizontal_x_mesh'],
            self.meshes['horizontal_y_mesh'],
            self.data['horizontal_windowed'],
            cmap = self.cmap_deform,
            vmin = height_cmap_min,
            vmax = height_cmap_max
        )

        ax[1,0].pcolormesh(
            self.meshes['vertical_x_frequency_mesh'],
            self.meshes['vertical_y_frequency_mesh'],
            self.PSD['vertical_windowed'],
            cmap = self.cmap_psd,
            vmin = PSD_cmap_min,
            vmax = PSD_cmap_max
        )

        display_PSD_cmap = ax[1,1].pcolormesh(
            self.meshes['horizontal_x_frequency_mesh'],
            self.meshes['horizontal_y_frequency_mesh'],
            self.PSD['horizontal_windowed'],
            cmap = self.cmap_psd,
            vmin = PSD_cmap_min,
            vmax = PSD_cmap_max
        )

        ax[1,2].scatter(
            x     = self.norms['vertical_wavelength_flat'][self.masks['vertical_zero-range']],
            y     = self.PSD['vertical_windowed_flat'][self.masks['vertical_range']],
            color = 'k',
            s     = 0.75
        )

        ax[1,3].scatter(
            x     = self.norms['horizontal_wavelength_flat'][self.masks['horizontal_zero-range']],
            y     = self.PSD['horizontal_windowed_flat'][self.masks['horizontal_range']],
            color = 'k',
            s     = 0.75
        )

        # Finalizing Styling TODO: Add calculation for predicted wavelength from theory
        for mult in self.frequency_multiplier:
            ax[1,2].axvline(
                x      = 8.8 * mult,
                ymin   = 0,
                ymax   = 1,
                color  = 'r',
                lw     = 1,
                alpha  = 0.75,
                ls     = '--',
                dashes = (4, mult * 3),
                label  = r'{:.1f}$\mu m$'.format(8.8 * mult)
            )

            ax[1,3].axvline(
                x      = 8.8 * mult,
                ymin   = 0,
                ymax   = 1,
                color  = 'r',
                lw     = 1,
                alpha  = 0.75,
                ls     = '--',
                dashes = (4, mult * 3),
                label  = r'{:.1f}$\mu m$'.format(8.8 * mult)
            )

        f.colorbar(
            display_height_cmap,
            ax           = ax[0, 1:3],
            orientation  = 'horizontal',
            shrink       = 0.6,
            pad          = 0.06,
            label        = r'$z$ Displacement $(\mu m)$',
            ticklocation = 'top'
        )

        f.colorbar(
            display_PSD_cmap,
            ax           = ax[1, :2],
            orientation  = 'horizontal',
            shrink       = 0.6,
            pad          = 0.06,
            label        = 'Power',
            ticklocation = 'top'
        )

        ax[1,2].legend(
            loc            = 'lower right',
            title          = 'Theory Wavelengths',
            fancybox       = True,
            fontsize       = 'small',
            title_fontsize = 'small'
        )

        ax[1,3].legend(
            loc            = 'lower right',
            title          = 'Theory Wavelengths',
            fancybox       = True,
            fontsize       = 'small',
            title_fontsize = 'small'
        )

        if display:
            plt.show()

        if file_name:
            f.savefig(file_name, dpi = 750)


    def small_plot(self, title = None, file_name = None, display = False, index = 0): # TODO


        # Create plot objects
        fig,ax = plt.subplots(
            nrows   = 1,
            ncols   = 3,
            layout  = 'constrained',
            figsize = (11,4.5)
        )

        fig.get_layout_engine().set(
            w_pad  = 0.15,
            h_pad  = 0.075,
            hspace = 0.05,
            wspace = 0.05
        )

        if title:
            fig.suptitle(title, size = 16)

        # Set Styling
        ax[0].set(
            xlim   = [-int(self.parameters['grid_size'] / 2), int(self.parameters['grid_size'] / 2)],
            ylim   = [-int(self.parameters['grid_size'] / 2), int(self.parameters['grid_size'] / 2)],
            aspect = 'equal'
        )

        ax[1].set(
            xlim   = [np.min(self.meshes['y_frequency_mesh']), np.max(self.meshes['y_frequency_mesh'])],
            ylim   = [np.min(self.meshes['y_frequency_mesh']), np.max(self.meshes['y_frequency_mesh'])],
            aspect = 'equal'
        )

        ax[0].set_title('Raw displacement data of surface.', fontsize = 10)
        ax[1].set_title('Logged PSD', fontsize = 10)
        ax[2].set_title('Normed Wavelength Powers', fontsize = 10)

        ax[0].set_xlabel(r'$x$ Position $(\mu m)$', fontsize = 9)
        ax[1].set_xlabel(r'Spatial Frequency X $(\frac{1}{\mu m})$', fontsize = 9)
        ax[2].set_xlabel(r'Wavelength $(\mu m)$', fontsize = 9)

        ax[0].set_ylabel(r'$y$ Position $(\mu m)$', fontsize = 9)
        ax[1].set_ylabel(r'Spatial Frequency Y $(\frac{1}{\mu m})$', fontsize = 9)
        ax[2].set_ylabel(r'Power $(\mu m^4)$', fontsize = 9)

        # Calculate colorbar ranges
        height_cmap_min_v = self.interpolated_data[index]['data'].min()
        height_cmap_max_v = self.interpolated_data[index]['data'].max()
        PSD_cmap_min_v    = self.mean_PSD.min()
        PSD_cmap_max_v    = self.mean_PSD.max()

        display_height_cmap = ax[0].pcolormesh(
            self.meshes['x_mesh'],
            self.meshes['y_mesh'],
            self.interpolated_data[index]['data'],
            cmap = self.cmap_deform,
            vmin = height_cmap_min_v,
            vmax = height_cmap_max_v,
            shading = 'gouraud'
        )

        cb = fig.colorbar(
            display_height_cmap,
            ax       = ax[0],
            location = 'bottom',
            orientation = 'horizontal'
        )

        cb.set_label(
            label    = r'$z$ Displacement $(\mu m)$',
            fontsize = 8
        )

        cb.ax.tick_params(labelsize=6)

        display_PSD_cmap = ax[1].pcolormesh(
            self.meshes['x_frequency_mesh'],
            self.meshes['y_frequency_mesh'],
            self.mean_PSD,
            cmap = self.cmap_psd,
            vmin = PSD_cmap_min_v,
            vmax = PSD_cmap_max_v,
            shading = 'gouraud'
        )

        cb_2 = fig.colorbar(
            display_PSD_cmap,
            ax       = ax[1],
            location = 'bottom',
            orientation = 'horizontal'
        )

        cb_2.set_label(
            label    = r'Power $(\mu m^4)$',
            fontsize = 8
        )

        cb_2.ax.tick_params(labelsize=6)

        ax[2].scatter(
            x     = self.ranged_wavelengths,
            y     = self.mean_flat_PSD[self.masks['remove_zero']][self.masks['range']],
            color = 'k',
            s     = 0.75
        )

        wavelength = self.calculate_theory_wavelength()

        ax[2].axvline(
            x      = wavelength,
            ymin   = 0,
            ymax   = 1,
            color  = 'r',
            lw     = 1,
            alpha  = 0.55,
            ls     = '--',
            label  = r'{:.1f}$\mu m$'.format(wavelength)
        )

        ax[2].legend(
            loc            = 'lower right',
            title          = 'Theory Wavelength',
            fancybox       = True,
            fontsize       = 'small',
            title_fontsize = 'small'
        )

        if display:
            plt.show()

        if file_name:
            fig.savefig(file_name, dpi = 1200)
            plt.close(fig)


    def output_data(self, fpath, n_timesteps): # TODO

        N_MAX_VALUES = 10

        mean_PSD_data = self.mean_flat_PSD[self.masks['remove_zero']][self.masks['range']]

        max_value_index = np.argpartition(mean_PSD_data, -1 * N_MAX_VALUES)[-1 * N_MAX_VALUES:]

        max_values = mean_PSD_data[max_value_index]
        wavelengths = self.ranged_wavelengths[max_value_index]

        pd.DataFrame({
            'max_PSD_values' : max_values,
            'wavelengths'    : wavelengths
        }).to_csv(os.path.join(fpath, 'spatial_wavelength_data.csv'))

        pd.DataFrame({
            'time'                : self.meshes['time'],
            'interpolated_values' : self.temporal_data['interpolated'],
            'windowed_values'     : self.temporal_data['windowed'],
            'frequency'           : self.meshes['frequency'],
            'fft'                 : self.temporal_data['fft'],
            'PSD'                 : self.temporal_data['PSD']
        }).to_csv(os.path.join(fpath, 'temporal_frequency_data.csv'))