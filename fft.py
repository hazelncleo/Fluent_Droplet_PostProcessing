import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fftfreq, fftshift, fft2
from scipy.interpolate import LinearNDInterpolator
from skimage.filters import window
from glob import glob
from HazelsAwesomeTheme import red_text,green_text,blue_text,yellow_text
from blessed import Terminal

class FFT_ISO:
    def __init__(self, parameters = {'frequency' : 1.63e6,
                                     'n_levels_refinement' : 4,
                                     'channel_width' : 50,
                                     'grid_size' : 500}):
                 
        '''
        Docstring for __init__
        
        :param self: Description
        :param data: Description
        :param file: Description
        '''
        
        # Calculate sampling data from provided parameters
        self.sampling_data = {'frequency' : parameters['frequency'],
                              'channel_width' : parameters['channel_width'],
                              'grid_size' : parameters['grid_size'],
                              'elements_along_width' : 2 * np.power(2, parameters['n_levels_refinement']),
                              'elements_along_length' : 20 * np.power(2, parameters['n_levels_refinement'])}
        
        self.sampling_data.update({'sample_spacing' : parameters['grid_size']/self.sampling_data['elements_along_length'],
                                   'n_samples_width' : self.sampling_data['elements_along_width'] + 1,
                                   'n_samples_length' : self.sampling_data['elements_along_length'] + 1})
        
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
    

    def solve(self, data = None, file = None):
        '''
        Docstring for solve_for_data
        
        :param self: Description
        :param data: Description
        :param file: Description
        '''
        
        if file:
            if data:
                print(yellow_text('Both Data and a File were provided, using the data for subsequent calculations.'))
            else:
                # Read coordinate data
                with open(file, 'rb') as f:
                    data = np.load(f)
                    
        else:
            print(red_text('No Data or File was provided to perform the FFT on.'))
            raise FileExistsError

        # Convert data to micron and save
        self.raw_data = data*1e6
        del data
        
        self.centre_data()
        
        self.interpolate_data()
        
        self.calculate_ffts()
        
        self.calculate_PSDs()
        
        self.calculate_normed_wavelengths()


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
        
    
    def calculate_normed_wavelengths(self):
        
        self.norms = {
            'vertical_frequency' : np.sqrt(self.meshes['vertical_x_frequency_mesh']**2 + self.meshes['vertical_y_frequency_mesh']**2),
            'horizontal_frequency' : np.sqrt(self.meshes['horizontal_x_frequency_mesh']**2 + self.meshes['horizontal_y_frequency_mesh']**2)
            }
        
        self.norms.update({
            'vertical_frequency_flat' : self.norms['vertical_frequency'].flatten(),
            'horizontal_frequency_flat' : self.norms['horizontal_frequency'].flatten()
            })
        
                          
        self.norms.update({
            'vertical_wavelength_flat' : np.reciprocal(self.norms['vertical_frequency_flat'][self.norms['vertical_frequency_flat'] != 0]),
            'horizontal_wavelength_flat' : np.reciprocal(self.norms['horizontal_frequency_flat'][self.norms['horizontal_frequency_flat'] != 0])
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
        
        f,ax = plt.subplots(2,4, width_ratios=[1,1,1,1], height_ratios=[1,1], figsize=(20,10))
        
        if title:
            f.suptitle(title)
        
        # Set Styling
        ax[0,0].set(xlim = [-int(self.sampling_data['grid_size']/2),int(self.sampling_data['grid_size']/2)],
                    ylim = [-int(self.sampling_data['grid_size']/2),int(self.sampling_data['grid_size']/2)],
                    title = 'Raw Data for Vertical Section',
                    xlabel = r'Position $(\mu m)$',
                    ylabel = r'Position $(\mu m)$')
        
        ax[0,1].set(xlim = [-int(self.sampling_data['grid_size']/2),int(self.sampling_data['grid_size']/2)],
                    ylim = [-int(self.sampling_data['grid_size']/2),int(self.sampling_data['grid_size']/2)],
                    title = 'Raw Data for Horizontal Section',
                    xlabel = r'Position $(\mu m)$',
                    ylabel = r'Position $(\mu m)$')
        
        ax[0,2].set(xlim = [-int(self.sampling_data['grid_size']/2),int(self.sampling_data['grid_size']/2)],
                    ylim = [-int(self.sampling_data['grid_size']/2),int(self.sampling_data['grid_size']/2)],
                    title = 'Windowed Vertical Data',
                    xlabel = r'Position $(\mu m)$',
                    ylabel = r'Position $(\mu m)$')
        
        ax[0,3].set(xlim = [-int(self.sampling_data['grid_size']/2),int(self.sampling_data['grid_size']/2)],
                    ylim = [-int(self.sampling_data['grid_size']/2),int(self.sampling_data['grid_size']/2)],
                    title = 'Windowed Horizontal Data',
                    xlabel = r'Position $(\mu m)$',
                    ylabel = r'Position $(\mu m)$')
        
        ax[1,0].set(xlim = [np.min(self.meshes['vertical_y_frequency_mesh']),np.max(self.meshes['vertical_y_frequency_mesh'])],
                    ylim = [np.min(self.meshes['vertical_y_frequency_mesh']),np.max(self.meshes['vertical_y_frequency_mesh'])],
                    title = 'Logged PSD, Vertical',
                    xlabel = r'Spatial Frequency $(\frac{1}{\mu m})$',
                    ylabel = r'Spatial Frequency $(\frac{1}{\mu m})$')
        
        ax[1,1].set(xlim = [np.min(self.meshes['vertical_y_frequency_mesh']),np.max(self.meshes['vertical_y_frequency_mesh'])],
                    ylim = [np.min(self.meshes['vertical_y_frequency_mesh']),np.max(self.meshes['vertical_y_frequency_mesh'])],
                    title = 'Logged PSD, Horizontal',
                    xlabel = r'Spatial Frequency $(\frac{1}{\mu m})$',
                    ylabel = r'Spatial Frequency $(\frac{1}{\mu m})$')
        
        ax[1,2].set(title = 'Normed Wavelength Powers, Vertical',
                    xlabel = r'Wavelength $(\mu m)$',
                    ylabel = r'Logged PSD')
        
        ax[1,3].set(title = 'Normed Wavelength Powers, Horizontal',
                    xlabel = r'Wavelength $(\mu m)$',
                    ylabel = r'Logged PSD')
        
        # Plot data
        ax[0,0].pcolormesh(self.meshes['vertical_x_mesh'],
                           self.meshes['vertical_y_mesh'],
                           self.data['vertical'])
        
        ax[0,1].pcolormesh(self.meshes['horizontal_x_mesh'],
                           self.meshes['horizontal_y_mesh'],
                           self.data['horizontal'])
        
        ax[0,2].pcolormesh(self.meshes['vertical_x_mesh'],
                           self.meshes['vertical_y_mesh'],
                           self.data['vertical_windowed'])
        
        ax[0,3].pcolormesh(self.meshes['horizontal_x_mesh'],
                           self.meshes['horizontal_y_mesh'],
                           self.data['horizontal_windowed'])
        
        ax[1,0].pcolormesh(self.meshes['vertical_x_frequency_mesh'],
                           self.meshes['vertical_y_frequency_mesh'],
                           self.PSD['vertical_windowed'])
        
        ax[1,1].pcolormesh(self.meshes['horizontal_x_frequency_mesh'],
                           self.meshes['horizontal_y_frequency_mesh'],
                           self.PSD['horizontal_windowed'])
        
        t = self.norms['vertical_wavelength_flat'] < 50
        t_2 = self.norms['vertical_frequency_flat'] > 0.02
        ax[1,2].scatter(self.norms['vertical_wavelength_flat'][t],self.PSD['vertical'].flatten()[t_2],color='k',s=0.75)
        
        if display:
            plt.show()
            
        if file_name:
            plt.savefig(file_name, dpi = 750)
        
         


def main():
    
    a = FFT_ISO()
    a.solve(file='.\\files\\iso_1.npy')
    a.full_plot(display=True)
    # Set sampling data
    print(bruh)
    elements_along_width = 32
    elements_along_length = 320
    sample_spacing = 500/elements_along_length
    N = elements_along_length + 1
    M = elements_along_width + 1

    files = glob('.\\files\\*.npy')

    for i,file in enumerate(files):
        # Read isosurface coordinates
        with open(file, 'rb') as f:
            coords = np.load(f)*1e6

        # Centre Z coordinates around 0
        coords[:,2] = coords[:,2] - np.mean(coords[:,2])
        
        # Create plot
        f,ax = plt.subplots(2,4, width_ratios=[1,1,1,1], height_ratios=[1,1], figsize=(20,10))

        # Create Meshes and other helpers
        channel_width = np.linspace(-25,25,M)
        channel_length = np.linspace(-250,250,N)
        X_vert,Y_vert = np.meshgrid(channel_width,channel_length)
        X_hor,Y_hor = np.meshgrid(channel_length,channel_width)
        windows = window(('tukey', 0.75), X_vert.shape)
        interp = LinearNDInterpolator(coords[:,:2], coords[:,2])

        # Interpolate data onto ordered grid
        Z_vert = interp(X_vert,Y_vert)
        Z_hor = interp(X_hor,Y_hor)

        # Get frequencies
        xf = fftfreq(M, sample_spacing)
        yf = fftfreq(N, sample_spacing)
        Xf,Yf = np.meshgrid(xf,yf)

        # 2D FFT with and without window functions
        Zf_vert_w = fft2(Z_vert*windows)
        Zf_hor_w = fft2(Z_hor*np.transpose(windows))

        Xf, Yf, Zf_vert_w, Zf_hor_w = fftshift(Xf), fftshift(Yf), fftshift(Zf_vert_w), fftshift(Zf_hor_w)

        # Get PSD
        Zf_vert_w = np.log(np.abs(Zf_vert_w)**2)
        Zf_hor_w = np.log(np.abs(Zf_hor_w)**2)

        truth_array = (coords[:,0] > -25.00001) & (coords[:,0] < 25.00001)
        ax[0,0].pcolormesh(X_vert,Y_vert,Z_vert,shading='gouraud')
        ax[0,0].set_xlim(-250,250)
        ax[0,0].set_ylim(-250,250)
        
        ax[0,1].pcolormesh(X_hor,Y_hor,Z_hor)
        ax[0,1].set_xlim(-250,250)
        ax[0,1].set_ylim(-250,250)

        ax[0,2].pcolormesh(X_vert,Y_vert,Z_vert*windows)
        ax[0,2].set_xlim(-250,250)
        ax[0,2].set_ylim(-250,250)

        ax[0,3].pcolormesh(X_hor,Y_hor,Z_hor*np.transpose(windows))
        ax[0,3].set_xlim(-250,250)
        ax[0,3].set_ylim(-250,250)

        ax[1,0].pcolormesh(Xf,Yf,Zf_vert_w)
        ax[1,0].set_xlim(-10*np.max(Xf),10*np.max(Xf))

        ax[1,1].pcolormesh(np.transpose(Yf),np.transpose(Xf),Zf_hor_w,shading='gouraud')
        ax[1,1].set_ylim(-10*np.max(Xf),10*np.max(Xf))

        frequency= np.sqrt(np.power(Xf,2) + np.power(Yf,2)).flatten()
        truth_array = np.logical_and(frequency > 0.022 , frequency <= 0.4) 
        Zf_vert_w = Zf_vert_w.flatten()
        ax[1,2].scatter(1/frequency[truth_array], Zf_vert_w[truth_array],s=0.75,c='k')
        for mult in [0.5,1,2,3]:
            ax[1,2].axvline(8.8*mult,0,1, color='r',lw=1,alpha=0.5,ls='--')
            ax[1,3].axvline(8.8*mult,0,1, color='r',lw=1,alpha=0.5,ls='--')

        frequency = np.transpose(np.sqrt(np.power(Xf,2) + np.power(Yf,2))).flatten()
        truth_array = np.logical_and(frequency > 0.022 , frequency <= 0.4) 
        Zf_hor_w = Zf_hor_w.flatten()
        ax[1,3].scatter(1/frequency[truth_array], Zf_hor_w[truth_array],s=0.75,c='k')
        
        ax[0,0].set_title('Raw Data for Vertical Section')
        ax[0,1].set_title('Raw Data for Horizontal Section')
        ax[0,2].set_title('Windowed Vertical Data')
        ax[0,3].set_title('Windowed Horizontal Data')
        ax[1,0].set_title('Logged PSD, Vertical')
        ax[1,1].set_title('Logged PSD, Horizontal')
        ax[1,2].set_title('Normed Wavelength Powers, Vertical')
        ax[1,3].set_title('Normed Wavelength Powers, Horizontal')
        
        ax[0,0].set_xlabel(r'Position $(\mu m)$')
        ax[0,0].set_ylabel(r'Position $(\mu m)$')
        ax[0,1].set_xlabel(r'Position $(\mu m)$')
        ax[0,1].set_ylabel(r'Position $(\mu m)$')
        ax[0,2].set_xlabel(r'Position $(\mu m)$')
        ax[0,2].set_ylabel(r'Position $(\mu m)$')
        ax[0,3].set_xlabel(r'Position $(\mu m)$')
        ax[0,3].set_ylabel(r'Position $(\mu m)$')
        ax[1,0].set_xlabel(r'Spatial Frequency $(\frac{1}{\mu m})$')
        ax[1,0].set_ylabel(r'Spatial Frequency $(\frac{1}{\mu m})$')
        ax[1,1].set_xlabel(r'Spatial Frequency $(\frac{1}{\mu m})$')
        ax[1,1].set_ylabel(r'Spatial Frequency $(\frac{1}{\mu m})$')
        ax[1,2].set_xlabel(r'Wavelength $(\mu m)$')
        ax[1,2].set_ylabel(r'Logged PSD')
        ax[1,3].set_xlabel(r'Wavelength $(\mu m)$')
        ax[1,3].set_ylabel(r'Logged PSD')
        
        

        plt.savefig('.\\images\\image_{}.png'.format(i),dpi=750)
    




main()





'''
'vertical_frequency_angles' : np.arccos(
np.where(
np.isclose(self.norms['vertical_frequency'],0), 
0, 
np.where(self.meshes['vertical_y_frequency_mesh'] < 0, 
-self.meshes['vertical_x_frequency_mesh'], 
self.meshes['vertical_x_frequency_mesh'])/self.norms['vertical_frequency']
)),

'horizontal_frequency_angles' : np.arccos(
np.where(
np.isclose(self.norms['horizontal_frequency'],0), 
0, 
np.where(self.meshes['horizontal_y_frequency_mesh'] < 0, 
-self.meshes['horizontal_x_frequency_mesh'], 
self.meshes['horizontal_x_frequency_mesh'])/self.norms['horizontal_frequency']
))
})'''