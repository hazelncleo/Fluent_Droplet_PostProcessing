import ansys.pyensight.core as ens

from droplet_sizing_script_controller import DropletSizingScriptController
from post_processor import PostProcessor
from ensight_controller import EnsightController
from fft_iso import FFT_ISO




# TODO, formalise this part of the script
parameters = {
    'frequency'           : 1.63e6,
    'amplitude'           : 1e-6,
    'n_cycles'            : 60,
    'n_levels_refinement' : 4,
    'channel_width'       : 50,
    'grid_size'           : 500
}


if __name__ == '__main__':
    session = ens.LocalLauncher(batch              = True, 
                                ansys_installation = 'C:\\Program Files\\ANSYS Inc\\v251').start()#, use_egl=True, additional_command_line_options=['-v 5']).start()
    
    ensight_pp = EnsightController(session    = session, 
                               parameters = parameters, 
                               fpath      = '')

    ensight_pp.post_process()
else:

    # TODO BETTER WAY TO CHECK IF ENSIGHT EXISTS
    try:
        ensight_pp = EnsightController(ensight = ensight, parameters = parameters)

        ensight_pp.post_process()
    except:
        pass