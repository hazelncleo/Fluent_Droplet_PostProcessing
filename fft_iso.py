import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fftfreq, fftshift, fft2
from scipy.interpolate import LinearNDInterpolator
from skimage.filters import window
from seaborn import color_palette, dark_palette
from HazelsAwesomeTheme import red_text,green_text,blue_text,yellow_text
plt.style.use('ggplot')


class FFT_ISO:
    def __init__(self, parameters = {'frequency' : 1.63e6,
                                     'amplitude' : 1e-6,
                                     'n_cycles' : 60,
                                     'n_levels_refinement' : 4,
                                     'channel_width' : 50,
                                     'grid_size' : 500}):
                 
        '''
        Docstring for __init__
        
        :param self: Description
        :param data: Description
        :param file: Description
        '''

        self.cmap_psd = color_palette('rocket', as_cmap=True)
        self.cmap_deform = dark_palette("#69d", as_cmap=True)
        
        # Calculate sampling data from provided parameters
        self.sampling_data = {'frequency' : parameters['frequency'],
                              'amplitude' : 1e-6,
                              'n_cycles' : 60,
                              'channel_width' : parameters['channel_width'],
                              'grid_size' : parameters['grid_size'],
                              'elements_along_width' : 2 * np.power(2, parameters['n_levels_refinement']),
                              'elements_along_length' : 20 * np.power(2, parameters['n_levels_refinement'])}
        
        self.sampling_data.update({'sample_spacing' : parameters['grid_size']/self.sampling_data['elements_along_length'],
                                   'n_samples_width' : self.sampling_data['elements_along_width'] + 1,
                                   'n_samples_length' : self.sampling_data['elements_along_length'] + 1})
        
        self.data_loaded = False
        self.frequency_multiplier = np.array([0.5,1,2,3])
        
        self.create_meshes()
        
        self.create_window_functions()
        
        
    def create_meshes(self):
        
        self.meshes = {}
        
        # Get 1D spatial and frequency domain meshes
        self.meshes['channel_width'] = np.linspace(start = -int(self.sampling_data['channel_width']/2),
                                                   stop = int(self.sampling_data['channel_width']/2),
                                                   num = self.sampling_data['n_samples_width'])
        
        self.meshes['channel_length'] = np.linspace(start = -int(self.sampling_data['grid_size']/2),
                                                    stop = int(self.sampling_data['grid_size']/2),
                                                    num = self.sampling_data['n_samples_length'])
        
        self.meshes['channel_width_frequency'] = fftshift(fftfreq(self.sampling_data['n_samples_width'],
                                                                  self.sampling_data['sample_spacing']))
        
        self.meshes['channel_length_frequency'] = fftshift(fftfreq(self.sampling_data['n_samples_length'],
                                                                   self.sampling_data['sample_spacing']))
        
        # Get 2D spatial and frequency domain meshes
        self.meshes['vertical_x_mesh'], self.meshes['vertical_y_mesh'] = np.meshgrid(self.meshes['channel_width'],
                                                                                     self.meshes['channel_length'])
        
        self.meshes['horizontal_x_mesh'], self.meshes['horizontal_y_mesh'] = np.meshgrid(self.meshes['channel_length'],
                                                                                         self.meshes['channel_width'])
        
        self.meshes['vertical_x_frequency_mesh'], self.meshes['vertical_y_frequency_mesh'] = np.meshgrid(self.meshes['channel_width_frequency'],
                                                                                                         self.meshes['channel_length_frequency'])
        
        self.meshes['horizontal_x_frequency_mesh'], self.meshes['horizontal_y_frequency_mesh'] = np.meshgrid(self.meshes['channel_length_frequency'],
                                                                                                             self.meshes['channel_width_frequency'])
        
        
    def create_window_functions(self):
        '''
        Docstring for create_window_functions
        
        :param self: Description
        '''
        self.windows = {'vertical' : window(('tukey', 0.75), self.meshes['vertical_x_mesh'].shape),
                        'horizontal' : window(('tukey', 0.75), self.meshes['horizontal_x_mesh'].shape)} 
    

    def solve(self, data = None, file = None, time_data = [0,0]):
        '''
        Docstring for solve_for_data
        
        :param self: Description
        :param data: Description
        :param file: Description
        '''
        
        if data is not None:
            if file:
                print(yellow_text('Both Data and a File were provided, using the data for subsequent calculations.'))
        elif file:
            # Read coordinate data
            with open(file, 'rb') as f:
                data = np.load(f)
        else:
            print(red_text('No Data or File was provided to perform the FFT on.'))
            raise FileExistsError

        # Convert data to micron and save
        self.time_data = time_data
        self.raw_data = data*1e6
        del data
        
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
        self.raw_data[:,2] = self.raw_data[:,2] - np.mean(self.raw_data[:,2])
    
    
    def interpolate_data(self):
        
        # Create interpolator object
        interpolator = LinearNDInterpolator(self.raw_data[:,:2], self.raw_data[:,2])
        
        # Interpolate onto ordered grid 
        self.data = {'vertical' : interpolator(self.meshes['vertical_x_mesh'],self.meshes['vertical_y_mesh']),
                     'horizontal' : interpolator(self.meshes['horizontal_x_mesh'],self.meshes['horizontal_y_mesh'])}
        
        self.data.update({'vertical_windowed' : self.data['vertical']*self.windows['vertical'],
                          'horizontal_windowed' : self.data['horizontal']*self.windows['horizontal']})

        
    def calculate_ffts(self):
        
        self.fft_data = {'vertical' : fftshift(fft2(self.data['vertical'])),
                         'horizontal' : fftshift(fft2(self.data['horizontal'])),
                         'vertical_windowed' : fftshift(fft2(self.data['vertical_windowed'])),
                         'horizontal_windowed' : fftshift(fft2(self.data['horizontal_windowed']))}
        
        
    def calculate_PSDs(self):
        
        self.PSD = {'vertical' : np.log(np.abs(self.fft_data['vertical'])**2),
                    'horizontal' : np.log(np.abs(self.fft_data['horizontal'])**2),
                    'vertical_windowed' : np.log(np.abs(self.fft_data['vertical_windowed'])**2),
                    'horizontal_windowed' : np.log(np.abs(self.fft_data['horizontal_windowed'])**2)}
        
        self.PSD.update({'vertical_flat' : self.PSD['vertical'].flatten(),
                         'horizontal_flat' : self.PSD['horizontal'].flatten(),
                         'vertical_windowed_flat' : self.PSD['vertical_windowed'].flatten(),
                         'horizontal_windowed_flat' : self.PSD['horizontal_windowed'].flatten()})
        

    def calculate_masks(self):

        self.masks = {'vertical_zero' : self.norms['vertical_frequency_flat'] != 0,
                      'horizontal_zero' : self.norms['horizontal_frequency_flat'] != 0,
                      'vertical_range' : self.norms['vertical_frequency_flat'] > 0.02,
                      'horizontal_range' : self.norms['horizontal_frequency_flat'] > 0.02}
        
        self.masks.update({'vertical_zero-range' : self.masks['vertical_range'][self.masks['vertical_zero']],
                           'horizontal_zero-range' : self.masks['horizontal_range'][self.masks['horizontal_zero']]})

    
    def calculate_normed_wavelengths(self):
        
        self.norms = {
            'vertical_frequency' : np.sqrt(self.meshes['vertical_x_frequency_mesh']**2 + self.meshes['vertical_y_frequency_mesh']**2),
            'horizontal_frequency' : np.sqrt(self.meshes['horizontal_x_frequency_mesh']**2 + self.meshes['horizontal_y_frequency_mesh']**2)
            }
        
        self.norms.update({
            'vertical_frequency_flat' : self.norms['vertical_frequency'].flatten(),
            'horizontal_frequency_flat' : self.norms['horizontal_frequency'].flatten()
            })
        
        self.calculate_masks()
                          
        self.norms.update({
            'vertical_wavelength_flat' : np.reciprocal(self.norms['vertical_frequency_flat'][self.masks['vertical_zero']]),
            'horizontal_wavelength_flat' : np.reciprocal(self.norms['horizontal_frequency_flat'][self.masks['horizontal_zero']])
            })
        
        
    def angles_from_frequency_data(self):
        pass
    
    '''
    ------------------
    PLOTTING
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
        f,ax = plt.subplots(2,4, width_ratios=[1,1,1,1], height_ratios=[1,1], layout='constrained', figsize=(22,13))
        f.get_layout_engine().set(w_pad=0.2, h_pad=0.1, hspace=0.05, wspace=0)

        if title:
            f.suptitle(title, size=50)
        
        # Set Styling
        ax[0,0].set(xlim = [-int(self.sampling_data['grid_size']/2),int(self.sampling_data['grid_size']/2)],
                    ylim = [-int(self.sampling_data['grid_size']/2),int(self.sampling_data['grid_size']/2)],
                    aspect = 'equal')
        
        ax[0,1].set(xlim = [-int(self.sampling_data['grid_size']/2),int(self.sampling_data['grid_size']/2)],
                    ylim = [-int(self.sampling_data['grid_size']/2),int(self.sampling_data['grid_size']/2)],
                    aspect = 'equal')
        
        ax[0,2].set(xlim = [-int(self.sampling_data['grid_size']/2),int(self.sampling_data['grid_size']/2)],
                    ylim = [-int(self.sampling_data['grid_size']/2),int(self.sampling_data['grid_size']/2)],
                    aspect = 'equal')
        
        ax[0,3].set(xlim = [-int(self.sampling_data['grid_size']/2),int(self.sampling_data['grid_size']/2)],
                    ylim = [-int(self.sampling_data['grid_size']/2),int(self.sampling_data['grid_size']/2)],
                    aspect = 'equal')
        
        ax[1,0].set(xlim = [np.min(self.meshes['vertical_y_frequency_mesh']),np.max(self.meshes['vertical_y_frequency_mesh'])],
                    ylim = [np.min(self.meshes['vertical_y_frequency_mesh']),np.max(self.meshes['vertical_y_frequency_mesh'])],
                    aspect = 'equal')
        
        ax[1,1].set(xlim = [np.min(self.meshes['vertical_y_frequency_mesh']),np.max(self.meshes['vertical_y_frequency_mesh'])],
                    ylim = [np.min(self.meshes['vertical_y_frequency_mesh']),np.max(self.meshes['vertical_y_frequency_mesh'])],
                    aspect = 'equal')
        
        ax[1,2].set_aspect(1.75)
        ax[1,3].set_aspect(1.75)

        ax[0,0].set_title('Raw Data for Vertical Section', fontsize=12)
        ax[0,1].set_title('Raw Data for Horizontal Section', fontsize=12)
        ax[0,2].set_title('Windowed Vertical Data', fontsize=12)
        ax[0,3].set_title('Windowed Horizontal Data', fontsize=12)
        ax[1,0].set_title('Logged PSD, Vertical', fontsize=12)
        ax[1,1].set_title('Logged PSD, Horizontal', fontsize=12)
        ax[1,2].set_title('Normed Wavelength Powers, Vertical', fontsize=12)
        ax[1,3].set_title('Normed Wavelength Powers, Horizontal', fontsize=12)

        ax[0,0].set_xlabel(r'Position $(\mu m)$', fontsize=9)
        ax[0,1].set_xlabel(r'Position $(\mu m)$', fontsize=9)
        ax[0,2].set_xlabel(r'Position $(\mu m)$', fontsize=9)
        ax[0,3].set_xlabel(r'Position $(\mu m)$', fontsize=9)
        ax[1,0].set_xlabel(r'Spatial Frequency X $(\frac{1}{\mu m})$', fontsize=9)
        ax[1,1].set_xlabel(r'Spatial Frequency X $(\frac{1}{\mu m})$', fontsize=9)
        ax[1,2].set_xlabel(r'Wavelength $(\mu m)$', fontsize=9)
        ax[1,3].set_xlabel(r'Wavelength $(\mu m)$', fontsize=9)

        ax[0,0].set_ylabel(r'Position $(\mu m)$', fontsize=9)
        ax[0,1].set_ylabel(r'Position $(\mu m)$', fontsize=9)
        ax[0,2].set_ylabel(r'Position $(\mu m)$', fontsize=9)
        ax[0,3].set_ylabel(r'Position $(\mu m)$', fontsize=9)
        ax[1,0].set_ylabel(r'Spatial Frequency Y $(\frac{1}{\mu m})$', fontsize=9)
        ax[1,1].set_ylabel(r'Spatial Frequency Y $(\frac{1}{\mu m})$', fontsize=9)
        ax[1,2].set_ylabel(r'Logged PSD', fontsize=9)
        ax[1,3].set_ylabel(r'Logged PSD', fontsize=9)
        
        # Calculate colorbar ranges
        vmin_1 = min(np.min(self.data['vertical']),np.min(self.data['horizontal']))
        vmax_1 = max(np.max(self.data['vertical']),np.max(self.data['horizontal']))
        vmin_2 = min(np.min(self.PSD['vertical_windowed']),np.min(self.PSD['horizontal_windowed']))
        vmax_2 = max(np.max(self.PSD['vertical_windowed']),np.max(self.PSD['horizontal_windowed']))

        # Plot data
        ax[0,0].pcolormesh(self.meshes['vertical_x_mesh'],
                           self.meshes['vertical_y_mesh'],
                           self.data['vertical'],
                           cmap=self.cmap_deform,
                           vmin=vmin_1,
                           vmax=vmax_1)
        
        ax[0,1].pcolormesh(self.meshes['horizontal_x_mesh'],
                           self.meshes['horizontal_y_mesh'],
                           self.data['horizontal'],
                           cmap=self.cmap_deform,
                           vmin=vmin_1,
                           vmax=vmax_1)
        
        ax[0,2].pcolormesh(self.meshes['vertical_x_mesh'],
                           self.meshes['vertical_y_mesh'],
                           self.data['vertical_windowed'],
                           cmap=self.cmap_deform,
                           vmin=vmin_1,
                           vmax=vmax_1)
        
        pcm_1 = ax[0,3].pcolormesh(self.meshes['horizontal_x_mesh'],
                                   self.meshes['horizontal_y_mesh'],
                                   self.data['horizontal_windowed'],
                                   cmap=self.cmap_deform,
                                   vmin=vmin_1,
                                   vmax=vmax_1)
        
        ax[1,0].pcolormesh(self.meshes['vertical_x_frequency_mesh'],
                           self.meshes['vertical_y_frequency_mesh'],
                           self.PSD['vertical_windowed'],
                           cmap=self.cmap_psd,
                           vmin=vmin_2,
                           vmax=vmax_2)
        
        pcm_2 = ax[1,1].pcolormesh(self.meshes['horizontal_x_frequency_mesh'],
                                   self.meshes['horizontal_y_frequency_mesh'],
                                   self.PSD['horizontal_windowed'],
                                   cmap=self.cmap_psd,
                                   vmin=vmin_2,
                                   vmax=vmax_2)
        
        ax[1,2].scatter(self.norms['vertical_wavelength_flat'][self.masks['vertical_zero-range']],
                        self.PSD['vertical_windowed_flat'][self.masks['vertical_range']],
                        color='k',
                        s=0.75)
        
        ax[1,3].scatter(self.norms['horizontal_wavelength_flat'][self.masks['horizontal_zero-range']],
                        self.PSD['horizontal_windowed_flat'][self.masks['horizontal_range']],
                        color='k',
                        s=0.75)
        
        # Finalizing Styling
        for mult in self.frequency_multiplier:
            ax[1,2].axvline(x = 8.8*mult,
                            ymin = 0,
                            ymax = 1, 
                            color='r',
                            lw=1,
                            alpha=0.75,
                            ls='--',
                            dashes=(4,mult*3), 
                            label = r'{:.1f}$\mu m$'.format(8.8*mult))
            
            ax[1,3].axvline(x = 8.8*mult,
                            ymin = 0,
                            ymax = 1, 
                            color='r',
                            lw=1,
                            alpha=0.75,
                            ls='--',
                            dashes=(4,mult*3), 
                            label = r'{:.1f}$\mu m$'.format(8.8*mult))

        f.colorbar(pcm_1, ax=ax[0,1:3], orientation='horizontal', shrink=0.6, pad=0.06, label=r'$z$ Displacement $(\mu m)$', ticklocation='top')
        f.colorbar(pcm_2, ax=ax[1,:2], orientation='horizontal', shrink=0.6, pad=0.06, label='Power', ticklocation='top')

        ax[1,2].legend(loc='lower right', title='Theory Wavelengths', fancybox=True, fontsize='small', title_fontsize='small')
        ax[1,3].legend(loc='lower right', title='Theory Wavelengths', fancybox=True, fontsize='small', title_fontsize='small')
        
        # Display or save to image
        if display:
            plt.show()
            
        if file_name:
            f.savefig(file_name, dpi = 750)