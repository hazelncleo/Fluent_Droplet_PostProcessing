import warnings

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    import ansys.fluent.core as pyfluent
    warnings.simplefilter("default")

class DropletSizingScriptController:
    def __init__(self, parameters, folder):
        self.parameters = parameters
        self.folder = folder

    def droplet_sizing_calculation(self, plot_results = False, animate_results = False): # TODO
        '''
        Docstring for droplet_sizing
        
        :param self: Description
        '''
        self.droplets_calculated = True


    def droplet_sizing_plot(self): # TODO
        pass


    def droplet_sizing_animation(self): # TODO
        pass

