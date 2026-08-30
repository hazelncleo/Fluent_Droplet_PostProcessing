import warnings

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    import ansys.pyensight.core as ens
    warnings.simplefilter("default")

import numpy as np
import os
import matplotlib.pyplot as plt
import glob
from HazelsAwesomeTheme import red_text,green_text,blue_text,yellow_text
import pandas as pd
import time

from fft_iso import FFT_ISO

class EnsightController:
    def __init__(self, parameters, folder):
        self.folder     = folder
        self.parameters = parameters

    def start_ensight(self):

        # TODO: ansys installation will need to change
        # HPC: ansys_installation = '/apps/ansys/25r1/v251', use_egl = True, use_mpi = 'openmpi', use_sos = 10, interconnect = 'ethernet'
        # MPI IS A MAYBE?? DOESNT SEEM TO WORK
        # Workstation: ansys_installation = 'C:\\Program Files\\ANSYS Inc\\v251'
        session = ens.LocalLauncher(
            batch              = True,
            ansys_installation = 'C:\\Program Files\\ANSYS Inc\\v251',
            additional_command_line_options = ['-glconfig', '-v','3'],
            use_sos = 3,
            use_mpi = 'intel2021'
        ).start()


        self.session = session
        self.ensight = self.session.ensight

        self.eocore = self.ensight.objs.core
        self.eonums = self.ensight.objs.enums
        self.eoutil = self.ensight.utils
        self.parts  = self.ensight.utils.parts
        self.views  = self.ensight.utils.views
        self.query  = self.ensight.utils.query
        self.cwd    = os.getcwd()

        self.load_data()

        self.ensight.file.animation_format("mpeg4")
        self.ensight.solution_time.show_as("step")

        self.t0                  = self.eocore.TIMEVALUES[0]
        self.tn                  = self.eocore.TIMEVALUES[-1]
        self.main_vp             = self.eocore.VPORTS['Main Viewport'][0]
        self.droplets_calculated = False

        self.main_vp.setattrs(
            {
                'BACKGROUNDTYPE' : self.eonums.VPORT_CONS,
                'CONSTANTRGB'    : [0.2, 0.2, 0.2]
            }
        )

        self.init_parts()

        self.init_variables()

        print(green_text('Ensight Instance initialised'))

    '''
    ------------------------
        Helper functions
    ------------------------
    '''

    def set_default_view(self):
        self.views.set_view_direction(1, 1, 1, name="isometric", up_axis=(0,0,1))

        for part in self.eocore.PARTS:
            part.visible = False


    def set_iso_view(self):

        # Set isometric view angle
        self.set_default_view()
        self.ensight.view.perspective("OFF")
        self.ensight.scene.ground_plane_visible("OFF")
        self.ensight.view_transf.zoom(0.6)
        self.ensight.solution_time.time_annotation("ON")
        self.ensight.text.select_begin(0)
        self.ensight.text.change_text(fr"""Time = <\\cnst C1 "%.1f" cycle_time\\>/{self.parameters['n_cycles']} cycles """)

        # Create current vibration position display
        self.top_arrow = self.shape_default.createannot('top_arrow')
        self.top_arrow.setattrs(
            {
                'TYPE'            : self.eonums.ANNOT_SHAPE_ARROW,
                'HEIGHT'          : 0.002,
                'LENGTH'          : 0.150,
                'ARROWTIPLENGTH'  : 0.08,
                'ARROWTIPSIZE'    : 1,
                'LOCATIONX'       : 0.024,
                'LOCATIONY'       : 0.6,
                'ROTATIONALANGLE' : 270,
                'RGB'             : [1, 1, 1],
                'FILL'            : True
            }
        )

        self.bottom_arrow = self.shape_default.createannot('bottom_arrow')
        self.bottom_arrow.setattrs(
            {
                'TYPE'            : self.eonums.ANNOT_SHAPE_ARROW,
                'LENGTH'          : 0.150,
                'HEIGHT'          : 0.002,
                'LOCATIONX'       : 0.024,
                'LOCATIONY'       : 0.3,
                'ROTATIONALANGLE' : 90,
                'ARROWTIPLENGTH'  : 0.08,
                'ARROWTIPSIZE'    : 1,
                'RGB'             : [1, 1, 1],
                'FILL'            : True
            }
        )

        self.top_gauge = self.gauge_default.createannot(self.vibration_state_pos.DESCRIPTION)
        self.top_gauge.setattrs(
            {
                'BACKGROUND'    : True,
                'BACKGROUNDRGB' : [0.2, 0.2, 0.2],
                'LEVELRGB'      : [0.65, 0.07, 0.08],
                'BORDER'        : False,
                'HEIGHT'        : 0.13,
                'LOCATIONX'     : 0.025,
                'LOCATIONY'     : 0.452,
                'MAXIMUM'       : 1,
                'MINIMUM'       : 0,
                'VALUE'         : False,
                'WIDTH'         : 0.02
            }
        )

        self.bottom_gauge = self.gauge_default.createannot(self.vibration_state_neg.DESCRIPTION)
        self.bottom_gauge.setattrs(
            {
                'BACKGROUND'    : True,
                'BACKGROUNDRGB' : [0.65, 0.07, 0.08],
                'BORDER'        : False,
                'HEIGHT'        : 0.13,
                'LOCATIONX'     : 0.025,
                'LOCATIONY'     : 0.318,
                'MAXIMUM'       : 0,
                'MINIMUM'       : -1,
                'LEVELRGB'      : [0.2, 0.2, 0.2],
                'VALUE'         : False,
                'WIDTH'         : 0.02
            }
        )

        self.text_1 = self.text_default.createannot('l')
        self.text_1.setattrs(
            {
                'JUSTIFICATION' : self.eonums.TS_CENTER,
                'SIZE'          : 100,
                'LOCATIONX'     : 0.046,
                'LOCATIONY'     : 0.44,
                'RGB'           : [0.2, 0.2, 0.2]
            }
        )
        self.text_2 = self.text_default.createannot('l')
        self.text_2.setattrs(
            {
                'JUSTIFICATION' : self.eonums.TS_CENTER,
                'SIZE'          : 100,
                'LOCATIONX'     : 0.046,
                'LOCATIONY'     : 0.212,
                'RGB'           : [0.2, 0.2, 0.2]
            }
        )

        self.centre_line = self.line_default.createannot()
        self.centre_line.setattrs(
            {
                'WIDTH'      : 4,
                'LOCATIONX1' : 0.024,
                'LOCATIONY1' : 0.45,
                'LOCATIONX2' : 0.036,
                'LOCATIONY2' : 0.45,
                'RGB'        : [1, 1, 1]
            }
        )


    def create_animation(self, fname = 'bruh'):

        print(f'Starting animation: {fname}')
        self.ensight.solution_time.update_to_first()
        self.ensight.file.animation_rend_offscreen("ON")
        self.ensight.file.animation_numpasses(16)
        self.ensight.file.animation_stereo("current")
        self.ensight.file.animation_screen_tiling(1,1)
        self.ensight.file.animation_file(os.path.join(self.folder, 'output', fname))
        self.ensight.file.animation_window_size("user_defined")
        self.ensight.file.animation_window_xy(2560,1440)
        self.ensight.solution_time.increment(1)
        self.ensight.file.animation_frames(self.tn[0])
        self.ensight.file.animation_start_number(0)
        self.ensight.file.animation_multiple_images("OFF")
        self.ensight.file.animation_raytrace_it("OFF")
        self.ensight.file.animation_play_flipbook("OFF")
        self.ensight.file.animation_play_time("ON")
        self.ensight.file.animation_play_keyframe("OFF")
        self.ensight.file.animation_reset_flipbook("OFF")
        self.ensight.file.animation_reset_traces("OFF")
        self.ensight.file.animation_reset_time("ON")
        self.ensight.file.animation_reset_keyframe("OFF")
        self.ensight.file.save_animation()
        print(f'Animation: {fname}, completed.')


    def load_data(self):
        '''
        Docstring for load_data

        :param self: Description
        '''

        # Ensight case data import
        if True:
            print('Loading data into ensight')
            start_time = time.time()
            self.ensight.data.sos_pass_wildcards("NO")
            self.ensight.data.sos_decompose_type('Temporal')
            self.ensight.data.sos_auto_distrib('dont')
            self.session.load_data(os.path.abspath(os.path.join(self.folder,'data/fluent_model.encas')))
            end_time = time.time()
            print(f'Data loaded into ensight after: {(end_time-start_time):.2f}s')

        # Legacy model fluent data import
        else:
            files = self.get_files()

            print('Loading data into ensight')
            start_time = time.time()

            # Seems to work ok on linux without these???

            self.session.load_data(files + '.cas.h5', result_file = files + '.dat.h5')
            end_time = time.time()

            print(f'Data loaded into ensight after: {(end_time-start_time):.2f}s')


    def get_files(self):
        '''
        Finds the file with shortest file name with extension .cas.h5 in the selected data folder.

        A fluent run using the case file
        '''
        return min(glob.glob(os.path.join(self.folder,'*.cas.h5')), key=lambda x: len(os.path.basename(x)))[:-7] + '-1-*'


    def init_variables(self):
        '''

        '''
        print('Initialising Variables.')
        start_time = time.time()

        self.coords = self.eocore.VARIABLES['Coordinates'][0]
        self.velocity = self.eocore.VARIABLES["velocity"][0]
        self.vf_water = self.eocore.VARIABLES["water_vof"][0]
        self.vf_air = self.eocore.VARIABLES["air_vof"][0]
        self.analysis_time = self.eocore.VARIABLES["Analysis_Time"][0]

        # Calculate time variables
        self.cycle_time = self.eocore.create_variable(
            name    = 'cycle_time',
            sources = [self.fluid_part],
            value   = f"Analysis_Time*{str(self.parameters['vibration_frequency'])}"
        )

        cycle_time_query = self.eoutil.query.create_temporal(
            name        = 'cycle_time',
            query_type  = self.eoutil.query.TEMPORAL_XYZ,
            part_list   = [self.fluid_part],
            variable1   = self.cycle_time,
            xyz         = [0, 0, 0],
            new_plotter = False
        )

        time_data = np.array([[row[1] for row in self.eocore.TIMEVALUES],
                              [row[1] for row in cycle_time_query.QUERY_DATA['xydata']]]).transpose()

        self.results_data = pd.DataFrame(
            data    = time_data,
            columns = ['analysis_time', 'cycle_time']
        )

        self.results_data.index.name = 'timestep_number'

        del time_data, cycle_time_query

        self.vibration_state = self.eocore.create_variable(
            name    = 'vibration_state',
            sources = [self.fluid_part],
            value   = f'COS((2*PI*cycle_time)-PI)',
            private = 1
        )

        self.vibration_state_pos = self.eocore.create_variable(
            name    = 'vibration_state_pos',
            sources = [self.fluid_part],
            value   = f'GT(vibration_state,0)'
        )

        self.vibration_state_neg = self.eocore.create_variable(
            name    = 'vibration_state_neg',
            sources = [self.fluid_part],
            value   = f'LT(vibration_state,0)'
        )

        self.max_velocity = self.eocore.create_variable(
            name    = 'max_velocity',
            value   = 'Max(plist,velocity,[],compute_per_case)',
            sources = [self.fluid_part],
            private = 1
        )

        # Calculate shearrate variables
        self.eocore.create_variable(
            name    = 'temp_1',
            value   = '0*air_vof@MLL/TT',
            sources = [self.fluid_part],
            private = 1
        )

        self.eocore.create_variable(
            name    = 'temp_2',
            value   = '1+temp_1',
            sources = [self.fluid_part],
            private = 1
        )

        self.eocore.create_variable(
            name    = 'water_threshold',
            value   = 'IF_GT(water_vof,0.75)',
            sources = [self.fluid_part]
        )

        self.eocore.create_variable(
            name    = 'temp_shear',
            value   = 'FluidShearMax(plist,velocity,1.0,temp_1,temp_2,1.0)',
            sources = [self.fluid_part],
            private = 1
        )

        self.shearrate = self.eocore.create_variable(
            name    = 'shearrate',
            value   = 'water_threshold*temp_shear@/T',
            sources = [self.fluid_part]
        )

        self.shearrate_vf = self.eocore.create_variable(
            name    = 'shearrate_vf',
            value   = 'water_vof*temp_shear@/T',
            sources = [self.fluid_part]
        )

        # Calculate outlet flowrate variables
        water_velocity_list = [self.fluid_part, *self.outlet_list]
        self.water_velocity = self.eocore.create_variable(
            name    = 'water_velocity',
            value   = 'velocity*water_vof',
            sources = water_velocity_list,
            private = 1
        )

        self.outlet_water_flowrate = self.eocore.create_variable(
            name = 'outlet_water_flowrate',
            value = 'Flow(plist, water_velocity)',
            sources = self.outlet_list
        )

        end_time = time.time()
        print(f'Variables initialised after: {(end_time-start_time):.2f}s')


    def init_parts(self):
        '''

        '''
        self.fluid_part          = self.eocore.PARTS['fluid'][0]
        self.symmetry_part       = self.eocore.PARTS['symmetry'][0]
        self.solid_coupling_part = self.eocore.PARTS['solid_coupling'][0]
        self.outlet_part         = self.eocore.PARTS['outlet'][0]

        try:
            self.outlet_top_part = self.eocore.PARTS['outlet_top'][0]
            self.outlet_list     = [self.outlet_part, self.outlet_top_part]
        except:
            self.outlet_list     = [self.outlet_part]

        self.iso_default         = self.eocore.DEFAULTPARTS[self.ensight.PART_ISO_SURFACE]
        self.shape_default       = self.eocore.DEFAULTANNOTS[self.ensight.ANNOT_SHAPE]
        self.line_default        = self.eocore.DEFAULTANNOTS[self.ensight.ANNOT_LINE]
        self.gauge_default       = self.eocore.DEFAULTANNOTS[self.ensight.ANNOT_GAUGE]
        self.text_default        = self.eocore.DEFAULTANNOTS[self.ensight.ANNOT_TEXT]

    '''
    ---------------------------------------
        Postprocessing output functions
    ---------------------------------------
    '''

    def basic_animation(self):

        start_time = time.time()

        self.basic_iso_part = self.iso_default.createpart(
            name       = "basic_iso",
            sources    = self.fluid_part,
            attributes = [
                ['VARIABLE',   self.vf_water],
                ['COLORBYRGB', [0.2, 0.6, 1]],
                ['OPAQUENESS', 0.55]
            ]
        )[0]

        self.symmetry_part.setattrs({
            'COLORBYPALETTE' : self.vf_air,
            'OPAQUENESS'     : 0.55
        })

        self.fluid_part.VISIBLE          = False
        self.symmetry_part.VISIBLE       = True
        self.outlet_part.VISIBLE         = False
        self.solid_coupling_part.VISIBLE = True

        self.solid_coupling_part.setattrs({
            'OPAQUENESS' : 0.75,
            'COLORBYRGB' : [0.57, 0.57, 0.57]
        })

        self.vf_air_palette = self.vf_air.PALETTE['air_vof<\\\\units>'][0]
        self.vf_air.LEGEND['air_vof<\\\\units>'][0].VISIBLE = False

        self.vf_air_palette.setattrs({
            'NLEVELS'      : 2,
            'MINMAX'       : [0, 0.5],
            'LIMIT_FRINGES' : self.eonums.PALETTE_LIMIT_FRINGES_INVISIBLE,
            'LEVELS_AND_COLORS' : [[0, 0.2, 0.6, 1, 1], [0.5, 0.2, 0.6, 1, 1]]
        })

        self.create_animation('general_animation')

        end_time = time.time()
        print(f'Basic animation completed after: {(end_time-start_time):.2f}s')

        # Cleanup
        self.solid_coupling_part.VISIBLE = False
        self.symmetry_part.VISIBLE       = False
        self.basic_iso_part.VISIBLE      = False

    '''
    -----------------------------
    '''

    def velocity_animation(self):

        start_time = time.time()

        self.velocity_iso_part = self.iso_default.createpart(
            name       = "velocity_iso",
            sources    = self.fluid_part,
            attributes = [
                ['VARIABLE',       self.vf_water],
                ['COLORBYPALETTE', self.velocity]
            ]
        )[0]

        self.fluid_part.VISIBLE          = False
        self.symmetry_part.VISIBLE       = False
        self.outlet_part.VISIBLE         = False
        self.solid_coupling_part.VISIBLE = True
        self.velocity_palette            = self.velocity.PALETTE['velocity<\\\\units>'][0]
        self.velocity_legend             = self.velocity.LEGEND['velocity<\\\\units>'][0]

        self.solid_coupling_part.setattrs(
            {
                'OPAQUENESS' : 0.75,
                'COLORBYRGB' : [0.57, 0.57, 0.57]
            }
        )

        if not hasattr(self, 'max_velocity_query'):
            self.max_velocity_query = self.eoutil.query.create_temporal(
                name        = 'Max velocity vs cycle time',
                query_type  = self.eoutil.query.TEMPORAL_MAXIMUM,
                part_list   = [self.fluid_part],
                variable1   = self.max_velocity,
                variable2   = self.cycle_time,
                new_plotter = False
            )

        final_max_velocity = np.max(np.array(self.max_velocity_query.QUERY_DATA['xydata']).transpose()[1])

        self.velocity_palette.setattrs({
            'MINMAX' : [0, int(np.ceil(final_max_velocity))],
            'NLEVELS': 11
        })

        self.velocity_legend.setattrs({
            'TYPE'    : self.eonums.FNC_CONST,
            'FORMAT'  : '%0.1f',
            'VISIBLE' : True
        })

        self.create_animation('velocity_animation')

        # Cleanup
        self.solid_coupling_part.setattrs({
            'OPAQUENESS' : 1,
            'COLORBYRGB' : [0.57, 0.57, 0.57]
        })

        self.velocity_legend.VISIBLE = False

        self.set_default_view()

        end_time = time.time()
        print(f'Velocity animation completed after {(end_time-start_time):.2f}')

    '''
    -----------------------------
    '''

    def shearrate_calculation(self, plot_results = False, animate_results = False):

        if not hasattr(self, 'max_shearrate_query'):
            self.max_shearrate_query = self.eoutil.query.create_temporal(
                name        = 'Max shear rate vs cycle time',
                query_type  = self.eoutil.query.TEMPORAL_MAXIMUM,
                part_list   = [self.fluid_part],
                variable1   = self.shearrate,
                variable2   = self.cycle_time,
                new_plotter = False
            )

            raw_shearrate_values = np.array(self.max_shearrate_query.QUERY_DATA['xydata']).transpose()[1]

            self.results_data['max_shearrate'] = raw_shearrate_values

            print(green_text('Max shearrate calculation completed'))

        if plot_results:
            self.shearrate_plot()

            print(green_text('Max shearrate plotting completed'))

        if animate_results:
            self.shearrate_animation()

            print(green_text('Max shearrate animation completed'))


    def shearrate_plot(self): # TODO
        pass


    def shearrate_animation(self): # TODO

        # The view must be translated to fit the plot into the display window.
        self.ensight.view_transf.translate(-2e-04,0,0)
        self.ensight.view_transf.zoom(1.2)

        self.shearrate_iso_part = self.iso_default.createpart(
            name="shearrate_iso",
            sources=self.fluid_part,
            attributes=[
                ['VARIABLE',   self.vf_water],
                ['TYPE',       self.eonums.ISO_SURF_SOLID],
                ['CONSTRAINT', self.eonums.CLIP_CHOICE_GREATER]
            ]
        )[0]

        self.shearrate_iso_part.COLORBYPALETTE = self.shearrate
        self.fluid_part.VISIBLE                = False
        self.symmetry_part.VISIBLE             = False
        self.solid_coupling_part.VISIBLE       = False
        self.outlet_part.VISIBLE               = False

        self.shearrate_palette = self.shearrate.PALETTE[0]
        self.shearrate_legend  = self.shearrate.LEGEND[0]

        final_max_shearrate = np.max(self.results_data['max_shearrate'])

        self.shearrate_palette.setattrs(
            {
                'SCALE_METHOD' : self.eonums.PALETTE_SCALE_LOG,
                'MINMAX'       : [1, int(np.power(10, np.ceil(np.log10(final_max_shearrate))))],
                'NLEVELS'      : int(np.ceil(1 + np.log10(final_max_shearrate)))
            }
        )

        self.shearrate_legend.setattrs(
            {
                'TYPE'        : self.eonums.FNC_CONST,
                'FORMAT'      : '%0.1e',
                'HEIGHT'      : 0.45,
                'DESCRIPTION' : 'Shear-rate <\\\\units>',
                'VISIBLE'     : True
            }
        )

        self.max_shearrate_query.setattrs(
            {
                'LINESTYLE' : self.eonums.LINE_SOLID,
                'LINETYPE' : self.eonums.CURVE_LINE_CONNECT,
                'MARKER' : self.eonums.CURVE_MARKER_SQUARE,
                'MARKERSCALE' : 1
            }
        )

        self.max_shearrate_line_plot = self.eocore.defaultplot[0].createplotter()
        self.max_shearrate_query.addtoplot(self.max_shearrate_line_plot)
        self.max_shearrate_line_plot.rescale()

        self.max_shearrate_line_plot.setattrs(
            {
                'PLOTTITLE'          : "Max Shear Rate of Liquid",
                'AXISXTITLE'         : "Number of cycles (n)",
                'AXISYTITLE'         : "Shear-rate (1/s)",
                'LEGENDVISIBLE'      : False,
                'AXISXAUTOSCALE'     : True,
                'AXISYSCALE'         : self.eonums.TRUE,
                'AXISXLABELFORMAT'   : "%.1f",
                'AXISXGRIDTYPE'      : 1,
                'AXISYGRIDTYPE'      : 1,
                'TIMEMARKER'         : True,
                'AXISYAUTOSCALE'     : False,
                'AXISYMIN'           : self.shearrate_palette.MINMAX[0],
                'AXISYMAX'           : self.shearrate_palette.MINMAX[1],
                'ORIGINX'            : 0.5,
                'ORIGINY'            : 0.55,
                'WIDTH'              : 0.48,
                'HEIGHT'             : 0.43,
                'TIMEMARKERRGB'      : [1, 0, 0],
                'TIMEMARKERWIDTH'    : 3,
                'AXISXNUMGRID'       : 1 + int(self.parameters['n_cycles'] / 10),
                'AXISYNUMGRIDLOG'    : 1 + int(np.round(np.log10(self.shearrate_palette.MINMAX[1] / self.shearrate_palette.MINMAX[0]))),
                'AXISYLABELFORMAT'   : '%.0e',
                'AXISYSGRIDTYPE'     : self.eonums.XY_GRID_TICK,
                'AXISXSGRIDTYPE'     : self.eonums.XY_GRID_TICK,
                'AXISYNUMSUBGRIDLOG' : 10,
                'AXISXNUMSUBGRID'    : 4
            }
        )

        self.create_animation('shearrate_animation')

        self.shearrate_iso_part.VISIBLE = False

    '''
    -----------------------------
    '''

    def flowrate_calculation(self, plot_results = False, animate_results = False):
        '''
        Docstring for total_flowrate

        :param self: Description
        :param plot_results: Description
        '''

        if not hasattr(self, 'flowrate_query'):
            self.flowrate_query = self.eoutil.query.create_temporal(
                name        = 'outlet_water_flowrate',
                query_type  = self.eoutil.query.TEMPORAL_XYZ,
                part_list   = [self.fluid_part],
                variable1   = self.outlet_water_flowrate,
                variable2   = self.analysis_time,
                xyz         = [0, 0, 0],
                new_plotter = False
            )

            # Calculate flowrates
            raw_flowrate_values    = np.array(self.flowrate_query.QUERY_DATA['xydata']).transpose()
            t                      = np.append(0, raw_flowrate_values[0])
            dVdt                   = np.append(0, -1 * 1e9 * raw_flowrate_values[1]) # Flowrate in microlitres
            Vol                    = (dVdt[1:] + dVdt[:-1]) * (t[1:] - t[:-1]) / 2
            total_volume_delivered = np.cumsum(Vol)
            volumetric_flowrate    = total_volume_delivered / t[1:]

            self.results_data['total_volume_delivered'] = total_volume_delivered
            self.results_data['volumetric_flowrate']    = volumetric_flowrate

            print(green_text('Outlet flowrate calculation completed'))

        if plot_results:
            self.flowrate_plot()

            print(green_text('Outlet flowrate plotting completed'))

        if animate_results:
            self.flowrate_animation()

            print(green_text('Outlet flowrate animation completed'))


    def flowrate_plot(self): # TODO
        '''

        '''
        f,ax = plt.subplots(1,2)

        ax[0].plot(self.results_data['analysis_time'], self.results_data['total_volume_delivered'])
        ax[1].plot(self.results_data['analysis_time'], self.results_data['volumetric_flowrate'])

        plt.show()


    def flowrate_animation(self): # TODO

        pass

    '''
    -----------------------------
    '''

    def fft_of_surface(self, plot_results = False):

        MAX_N_TIMESTEPS_TO_EXTRACT = 50
        n_timesteps = min(MAX_N_TIMESTEPS_TO_EXTRACT, len(self.eocore.TIMEVALUES))

        times = [self.eocore.TIMEVALUES[-1 * n_timesteps][1], self.eocore.TIMEVALUES[-1][1]]

        fft_calculator = FFT_ISO(parameters = self.parameters, n_timesteps = n_timesteps, times = times)

        self.coord_iso = self.iso_default.createpart(name="coord_iso", sources=self.fluid_part, attributes=[['VARIABLE',self.vf_water]])[0]
        self.ensight.solution_time.show_as("step")
        self.ensight.solution_time.increment(1)
        self.ensight.solution_time.update_to_last()

        for i in range(n_timesteps):

            self.iso_surface_coordinates = self.coord_iso.get_values([self.coords], activate=1)[self.coords]

            fft_calculator.send_data(
                data     = self.iso_surface_coordinates,
                time     = self.eocore.SOLUTIONTIME,
                index    = n_timesteps - i
            )

            self.ensight.solution_time.step_backward()

        fft_calculator.solve()

        fft_calculator.output_data(os.path.join(self.folder, 'output'))

        if plot_results:

            # Spatial plot
            fft_calculator.small_plot(title = 'FFT Plot', file_name = os.path.join(self.folder, 'output', 'spatial_plot.png'), index = n_timesteps)


        print(green_text('FFT of fluid surface completed'))



if __name__ == '__main__':
    parameters = {
        'vibration_frequency' : 1.63e6,
        'vibration_amplitude' : 0.5e-6,
        'n_cycles'            : 40,
        'n_elements'          : 35,
        'channel_width'       : 30,
        'grid_size'           : 500
    }

    folders = [
        'D:\\wavelength_droplet-sizing_project\\simulations\\straight_channel\\rigid_vibration\\163MHz_noise\\wavelengths'
    ]

    for i,folder in enumerate(folders):

        ensig = EnsightController(
            parameters = parameters,
            folder = folder
        )

        os.makedirs(os.path.join(folder,'output'), exist_ok=True)

        ensig.start_ensight()

        ensig.fft_of_surface(plot_results = True)

        ensig.session.close()

        del ensig
