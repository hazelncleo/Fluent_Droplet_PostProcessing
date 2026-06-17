import pandas as pd
import glob
import os
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as ticker
from HazelsAwesomeTheme import red_text,green_text,blue_text,yellow_text

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

        METRECUBED_TO_MICRON = 1e18
        VOLUME_COEFF         = 3 / (4 * np.pi)
        MAX_DIAM             = 150
        FPF_LIMIT            = 5
        QUANTILES_TO_EXTRACT = [0.1,0.5,0.9]

        csv_files = glob.glob(os.path.join(self.folder, 'droplets_data', '*.csv'))
        csv_files.sort(key = alphanum_key)

        self.time_info = np.linspace(5, self.parameters['n_cycles'], len(csv_files))

        quantiles              = np.zeros((len(csv_files), len(QUANTILES_TO_EXTRACT)))
        cumulative_quantiles   = np.zeros(quantiles.shape)
        fine_particle_fraction = np.zeros(len(csv_files))
        times                  = np.zeros(len(csv_files))

        for i,droplet_file in enumerate(csv_files):

            droplet_data = pd.read_csv(droplet_file, skiprows=10)

            with open(droplet_file) as f:
                for j,line in enumerate(f):
                    if j == 1:
                        times[i] = np.float64(line.split(' ')[-1])
                        break

            # Diameter in micron
            droplet_data['diameter'] = 2 * np.cbrt(VOLUME_COEFF * METRECUBED_TO_MICRON * droplet_data.volume)

            # Filter large diameters
            no_large_diameters = droplet_data[droplet_data['diameter'] < MAX_DIAM]

            # Calculate volume quantiles for current step
            if no_large_diameters.shape[0] >= 1:
                quantiles[i,:] = np.quantile(no_large_diameters['diameter'], QUANTILES_TO_EXTRACT, weights=no_large_diameters['volume'], method='inverted_cdf')

            # Append current step diameters & volumes
            if i == 0:
                total_diameters = no_large_diameters['diameter']
                total_volumes = no_large_diameters['volume']
            else:
                total_diameters = np.hstack((total_diameters, no_large_diameters['diameter']))
                total_volumes = np.hstack((total_volumes, no_large_diameters['volume']))

            if total_diameters.shape[0] >= 1:

                # Calculate cumulative quantiles
                cumulative_quantiles[i,:] = np.quantile(total_diameters, QUANTILES_TO_EXTRACT, weights=total_volumes, method='inverted_cdf')

                # Calculate fine particle fraction
                fine_particle_fraction[i] = np.sum(total_volumes, where = (total_diameters <= FPF_LIMIT)) / np.sum(total_volumes)

        self.results_data = pd.DataFrame({
            'timestep'               : np.arange(51, len(csv_files)+51),
            'times'                  : times,
            'cycle_time'             : times * self.parameters['vibration_frequency'],
            'Dv10'                   : quantiles[:,0],
            'Dv50'                   : quantiles[:,1],
            'Dv90'                   : quantiles[:,2],
            'Dv10_cumulative'        : cumulative_quantiles[:,0],
            'Dv50_cumulative'        : cumulative_quantiles[:,1],
            'Dv90_cumulative'        : cumulative_quantiles[:,2],
            'fine_particle_fraction' : fine_particle_fraction
            }
        )

        self.temp_data = {
            'total_diameters'        : total_diameters,
            'total_volumes'          : total_volumes
        }

        if plot_results:
            self.droplet_sizing_plot()

            print(green_text('Droplet sizing plotting completed'))

        if animate_results:
            self.droplet_sizing_animation()



    def droplet_sizing_plot(self):

        N_BINS = 20

        plt.style.use('ggplot')

        f,ax = plt.subplots(1,3,figsize=(16,5))
        f.tight_layout(pad=3.75, h_pad=0.3, w_pad=3.5)

        ax[0].set(
            xlabel = r'Droplet diameter $(\mu m)$',
            ylabel = 'Density',
            title  = 'Volume weighted size distribution',
            xscale = 'log'
        )

        ax[1].set(
            xlabel = 'Time (number of cycles)',
            ylabel = r'Droplet diameter $(\mu m)$',
            title  = 'Quantiles over time'
        )

        ax[2].set(
            xlabel = 'Time (number of cycles)',
            ylabel = r'Droplet diameter $(\mu m)$',
            title  = 'Cumulative quantiles over time'
        )

        sns.histplot(
            x         = self.temp_data['total_diameters'],
            weights   = self.temp_data['total_volumes'],
            stat      = 'density',
            bins      = N_BINS,
            kde       = True,
            kde_kws   = {'bw_adjust' : 3},
            color     = '#9e0d00',
            edgecolor = '#000000',
            ax        = ax[0],
            log_scale = True
        )

        ax[0].xaxis.set_major_locator(ticker.MultipleLocator(1))
        ax[0].xaxis.set_major_formatter(ticker.ScalarFormatter())


        ax[1].plot(self.results_data.cycle_time, self.results_data.Dv90, '-r', label = r'$D_v 90$', lw=0.85)
        ax[1].plot(self.results_data.cycle_time, self.results_data.Dv50, '-g', label = r'$D_v 50$', lw=0.85)
        ax[1].plot(self.results_data.cycle_time, self.results_data.Dv10, '-b', label = r'$D_v 10$', lw=0.85)
        ax[1].legend(loc='upper left', fancybox=True)

        ax[1].sharey(ax[2])


        ax[2].plot(self.results_data.cycle_time, self.results_data.Dv90_cumulative, '-r', label = r'Cumulative $D_v 90$', lw=0.85)
        ax[2].plot(self.results_data.cycle_time, self.results_data.Dv50_cumulative, '-g', label = r'Cumulative $D_v 50$', lw=0.85)
        ax[2].plot(self.results_data.cycle_time, self.results_data.Dv10_cumulative, '-b', label = r'Cumulative $D_v 10$', lw=0.85)
        ax[2].legend(loc='upper left', fancybox=True)

        f.savefig(os.path.join(self.folder, 'output', 'droplet_sizing_plot.png'), dpi=1500)


    def droplet_sizing_animation(self): # TODO
        pass

if __name__ == '__main__':
    sizer = DropletSizingScriptController({'n_cycles' : 40, 'vibration_frequency' : 1.63e6}, 'test_datasets')
    sizer.droplet_sizing_calculation(plot_results=True)