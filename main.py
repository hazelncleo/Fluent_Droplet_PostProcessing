import ansys.pyensight.core as ens
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from fft_iso import FFT_ISO
from scipy.interpolate import LinearNDInterpolator
import pandas as pd

class ensight_class:
    def __init__(
        self, 
        parameters, 
        fpath, 
        session = None, 
        ensight = None
    ):


        if session is not None:
            if ensight is not None:
                print('Both Ensight and Session Object Provided, using session object')

            self.session         = session
            self.ensight         = self.session.ensight
            self.sesion_provided = True

        elif ensight is not None:

            self.ensight         = ensight
            self.sesion_provided = False

        else:
            raise ValueError('ERROR: No session or ensight object specified')


        self.eocore     = self.ensight.objs.core
        self.eonums     = self.ensight.objs.enums
        self.eoutil     = self.ensight.utils
        self.parts      = self.ensight.utils.parts
        self.views      = self.ensight.utils.views
        self.query      = self.ensight.utils.query
        self.parameters = parameters
        
        self.load_data()

        self.ensight.file.animation_format("mpeg4")
        self.ensight.solution_time.show_as("step")

        self.t0      = self.eocore.TIMEVALUES[0]
        self.tn      = self.eocore.TIMEVALUES[-1]
        self.cwd     = os.getcwd()
        self.main_vp = self.eocore.VPORTS['Main Viewport'][0]

        self.main_vp.setattrs(
            {
                'BACKGROUNDTYPE' : self.eonums.VPORT_CONS,
                'CONSTANTRGB'    : [0.2, 0.2, 0.2]
            }
        )

        self.init_parts()

        self.init_variables()

    '''
    ------------------------
        Helper functions
    ------------------------
    '''

    def set_iso_view(self):

        # Set isometric view angle
        self.views.set_view_direction(1, 1, 1, name="isometric", up_axis=(0,0,1))
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

        self.ensight.solution_time.update_to_first()
        self.ensight.file.animation_rend_offscreen("ON")
        self.ensight.file.animation_numpasses(4)
        self.ensight.file.animation_stereo("current")
        self.ensight.file.animation_screen_tiling(1,1)
        self.ensight.file.animation_file(fr"""D:\Uni_Projects\PALM_Projects\Testing\Working_Geometry_Testing\postprocessing\{fname} """)
        self.ensight.file.animation_window_size("user_defined")
        self.ensight.file.animation_window_xy(3840,2160)
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

    '''
    ---------------------------------------
        Postprocessing output functions
    ---------------------------------------
    '''

    def basic_animation(self): # TODO

        pass


    def velocity_animation(self):
        
        self.velocity_iso_part = self.iso_default.createpart(
            name       = "velocity_iso", 
            sources    = self.fluid_part, 
            attributes = {
                            'VARIABLE'       : self.vf_water,
                            'COLORBYPALETTE' : self.velocity
                         }
        )[0]

        self.fluid_part.VISIBLE    = False
        self.symmetry_part.VISIBLE = False
        self.outlet_part.VISIBLE   = False
        self.velocity_palette      = self.velocity.PALETTE['Velocity<\\\\units>'][0]
        self.velocity_legend       = self.velocity.LEGEND['Velocity<\\\\units>'][0]

        self.solid_coupling_part.setattrs(
            {
                'OPAQUENESS' : 0.75,
                'COLORBYRGB' : [0.57, 0.57, 0.57]
            }
        )

        self.velocity_palette.set_range_to_over_time_minmax(self.t0[0], self.tn[0])

        self.velocity_palette.setattrs(
            {
                'MINMAX'  : [0, np.ceil(self.velocity_palette.MINMAX[1])],
                'NLEVELS' : 11
            }
        )

        self.velocity_legend.setattrs(
            {
                'TYPE'    : self.eonums.FNC_CONST,
                'FORMAT'  : '%0.1f'
            }
        )

        self.create_animation('vel_anim')

        self.velocity_iso_part.VISIBLE = False


    def shearrate_animation(self):
        
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
                ])[0]
        
        self.shearrate_iso_part.COLORBYPALETTE = self.shearrate
        self.fluid_part.VISIBLE                = False
        self.symmetry_part.VISIBLE             = False
        self.solid_coupling_part.VISIBLE       = False
        self.outlet_part.VISIBLE               = False

        self.shearrate_palette = self.shearrate.PALETTE[0]
        self.shearrate_legend  = self.shearrate.LEGEND[0]
        
        self.shearrate_palette.set_range_to_over_time_minmax(self.t0[0], self.tn[0])
        self.shearrate_palette.setattrs(
            {
                'SCALE_METHOD' : self.eonums.PALETTE_SCALE_LOG,
                'MINMAX'       : [1, np.power(10, np.ceil(np.log10(self.shearrate_palette.MINMAX[1])))],
                'NLEVELS'      : int(1 + np.log10(np.round(self.shearrate_palette.MINMAX[1])))
            }
        )

        self.shearrate_legend.setattrs(
            {
                'TYPE'        : self.eonums.FNC_CONST,
                'FORMAT'      : '%0.1e',
                'HEIGHT'      : 0.45,
                'DESCRIPTION' : 'Shear-rate <\\\\units>'
            }
        )

        self.max_shearrate_query = self.eoutil.query.create_temporal(
            name        = 'Max shear rate vs cycle time',
            query_type  = self.eoutil.query.TEMPORAL_MAXIMUM,
            part_list   = [self.fluid_part],
            variable1   = self.shearrate,
            variable2   = self.cycle_time,
            new_plotter = False
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


    def droplet_sizing(self): # TODO
        '''
        Docstring for droplet_sizing
        
        :param self: Description
        '''
        pass


    def total_flowrate(self, plot_results = True): # TODO
        
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
        raw_flowrate_values = np.array(self.flowrate_query.QUERY_DATA['xydata']).transpose()
        t                   = np.append(0, raw_flowrate_values[0])
        dVdt                = np.append(0, (-1 * 1e9 / 997) * raw_flowrate_values[1]) # Flowrate in microlitres
        V                   = (dVdt[1:] + dVdt[:-1]) * (t[1:] - t[:-1]) / 2
        V_total             = np.cumsum(V)
        dVdt_total          = V_total / t[1:]

        if plot_results:
            
            f,ax = plt.subplots(1,4)
            
            
            ax[0].plot(t[1:], dVdt)
            ax[2].plot(t[1:], V)
            ax[1].plot(t[1:], V_total)
            ax[3].plot(t[1:], dVdt_total)
            plt.show()
            return 0
        else:
            return np.vstack([t,
                              dVdt,
                              V,
                              V_total,
                              dVdt_total])            


    def fft_of_surface(self):

        
        fft_calculator = FFT_ISO(parameters = self.parameters)

        self.coord_iso = self.iso_default.createpart(name="coord_iso", sources=self.fluid_part, attributes=[['VARIABLE',self.vf_water]])[0]
        self.ensight.solution_time.show_as("step")
        self.ensight.solution_time.increment(1)
        self.ensight.solution_time.update_to_first()

        for i in range(self.tn[0]):
            self.iso_surface_coordinates = self.coord_iso.get_values([self.coords], activate=1)[self.coords]
            
            fft_calculator.solve(data=self.iso_surface_coordinates, time_data=[self.eocore.TIMESTEP+1,self.eocore.SOLUTIONTIME])

            fft_calculator.full_plot(f'FFT Plot t={i}', f'.\\images\\iso_plot_{int(self.eocore.TIMESTEP+1)}.png')
            
            self.ensight.solution_time.step_forward()


    def get_output_data(self): # TODO
        '''
        Data columns:
            time,
            cycle_time,
            max_shearrate,
            total_flowrate,
            fpf,
            dv10,
            dv50,
            dv90
        '''
        
        # Times
        analysis_time_query = self.eoutil.query.create_temporal(
            name        = 'Analysis_Time',
            query_type  = self.eoutil.query.TEMPORAL_XYZ,
            part_list   = [self.fluid_part],
            variable1   = self.analysis_time,
            xyz         = [0, 0, 0],
            new_plotter = False
        )

        cycle_time_query = self.eoutil.query.create_temporal(
            name        = 'cycle_time',
            query_type  = self.eoutil.query.TEMPORAL_XYZ,
            part_list   = [self.fluid_part],
            variable1   = self.cycle_time,
            xyz         = [0, 0, 0],
            new_plotter = False
        )

        if not hasattr(self, 'max_shearrate_query'):
            self.max_shearrate_query = self.eoutil.query.create_temporal(
                name        = 'Max shear rate vs cycle time',
                query_type  = self.eoutil.query.TEMPORAL_MAXIMUM,
                part_list   = [self.fluid_part],
                variable1   = self.shearrate,
                variable2   = self.cycle_time,
                new_plotter = False
            )

        analysis_time_values = np.array([row[1] for row in analysis_time_query.QUERY_DATA['xydata']])
        cycle_time_values    = np.array([row[1] for row in cycle_time_query.QUERY_DATA['xydata']])
        max_shearrate_values = np.array([row[1] for row in self.max_shearrate_query.QUERY_DATA['xydata']]) 


    def save_data_to_csv(self): # TODO
        
        pass


    '''
    -------------------------
        Control functions
    -------------------------
    '''

    def post_process(self, options = None): # TODO

        self.set_iso_view()

        self.total_flowrate()

        self.get_output_data()
        
        #self.velocity_animation()

        #self.shearrate_animation()

        #self.fft_of_surface()

    '''
    -----------------------------
        init Helper functions
    -----------------------------
    '''

    def load_data(self): # TODO
        '''
        Docstring for load_data
        
        :param self: Description
        '''
    
        if self.sesion_provided:
            
            #self.session.load_data('D:\\Uni_Projects\\PALM_Projects\\Simulations\\Fluent_Droplet_PostProcessing\\data\\test_simple_vibrate-1-*.cas.h5',
            #                    result_file='D:\\Uni_Projects\\PALM_Projects\\Simulations\\Fluent_Droplet_PostProcessing\\data\\test_simple_vibrate-1-*.dat.h5')
            self.session.load_data('D:\\Uni_Projects\\PALM_Projects\\Simulations\\Fluent_Droplet_PostProcessing\\data_droplets\\test_simple_vibrate-1-*.cas.h5',
                                   result_file='D:\\Uni_Projects\\PALM_Projects\\Simulations\\Fluent_Droplet_PostProcessing\\data_droplets\\test_simple_vibrate-1-*.dat.h5')
            #self.session.load_data('D:\\Uni_Projects\\PALM_Projects\\Testing\\Working_Geometry_Testing\\adaptive_64\\test_simple_vibrate-1-*.cas.h5',
            #                       result_file='D:\\Uni_Projects\\PALM_Projects\\Testing\\Working_Geometry_Testing\\adaptive_64\\test_simple_vibrate-1-*.dat.h5')

        else:
            self.ensight.part.select_default()
            self.ensight.part.modify_begin()
            self.ensight.part.elt_representation("3D_feature_2D_full")
            self.ensight.part.modify_end()
            self.ensight.data.binary_files_are("native")
            self.ensight.data.format("Fluent_HDF5")
            self.ensight.data.reader_option("'Load Internal Parts' OFF")
            self.ensight.data.reader_option("'Load _M1 _M2 vars' OFF")
            self.ensight.data.reader_option("'Load all cell types' OFF")
            self.ensight.data.reader_option("'Poly to Regular Cell' ON")
            self.ensight.data.reader_option("'Poly faced Hex to Poly' OFF")
            self.ensight.data.reader_option("'Fix Hanging Nodes' ON")
            self.ensight.data.reader_option("'Use Zone IDs for Parts' OFF")
            self.ensight.data.reader_option("'Show Fluent Variable Display Name' ON")
            self.ensight.data.reader_option("'Enable Part Grouping' OFF")
            self.ensight.data.reader_option("'Console Output' 'Normal'")
            self.ensight.data.reader_option("'Time Values' 'Read Time Values'")
            self.ensight.data.result(r"""/fs04/scratch2/pa16/submodels/simple_vibration_runs/adaptive_64_750nm/test_simple_vibrate-1-*.dat.h5 """)
            self.ensight.data.shift_time(1.000000, 0.000000, 0.000000)
            self.ensight.solution_time.monitor_for_new_steps("off")
            self.ensight.data.replace(r"""/fs04/scratch2/pa16/submodels/simple_vibration_runs/adaptive_64_750nm/test_simple_vibrate-1-*.cas.h5 """)
            self.ensight.command.delay_refresh("ON")


    def init_variables(self):
        '''
        
        '''

        self.coords = self.eocore.VARIABLES['Coordinates'][0]
        self.velocity = self.eocore.VARIABLES["Velocity"][0]
        self.vf_water = self.eocore.VARIABLES["Volume_fraction_water"][0]
        self.vf_air = self.eocore.VARIABLES["Volume_fraction_air"][0]
        self.analysis_time = self.eocore.VARIABLES["Analysis_Time"][0]

        # Calculate time variables
        self.cycle_time = self.eocore.create_variable(
            name    = 'cycle_time',
            sources = [self.fluid_part],
            value   = f"Analysis_Time*{str(self.parameters['frequency'])}"
        )
        
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
        
        # Calculate shearrate variables
        self.eocore.create_variable(
            name    = 'temp_1',
            value   = '0*Volume_fraction_air@MLL/TT', 
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
            value   = 'IF_GT(Volume_fraction_water,0.75)', 
            sources = [self.fluid_part]
        )
        
        self.eocore.create_variable(
            name    = 'temp_shear',
            value   = 'FluidShearMax(plist,Velocity,1.0,temp_1,temp_2,1.0)', 
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
            value   = 'Volume_fraction_water*temp_shear@/T', 
            sources = [self.fluid_part]
        )

        # Calculate outlet flowrate variables
        self.water_velocity = self.eocore.create_variable(
            name    = 'water_velocity', 
            value   = 'Velocity*Volume_fraction_water', 
            sources = [self.fluid_part, self.outlet_part],
            private = 1
        )

        self.outlet_water_flowrate = self.eocore.create_variable(
            name = 'outlet_water_flowrate',
            value = 'Flow(plist, water_velocity)',
            sources = [self.outlet_part]
        )


    def init_parts(self):
        '''
        
        '''
        self.fluid_part          = self.eocore.PARTS['fluid'][0]
        self.symmetry_part       = self.eocore.PARTS['symmetry'][0]
        self.solid_coupling_part = self.eocore.PARTS['solid_coupling'][0]
        self.outlet_part         = self.eocore.PARTS['outlet'][0]
        self.iso_default         = self.eocore.DEFAULTPARTS[self.ensight.PART_ISO_SURFACE]
        self.shape_default       = self.eocore.DEFAULTANNOTS[self.ensight.ANNOT_SHAPE]
        self.line_default        = self.eocore.DEFAULTANNOTS[self.ensight.ANNOT_LINE]
        self.gauge_default       = self.eocore.DEFAULTANNOTS[self.ensight.ANNOT_GAUGE]
        self.text_default        = self.eocore.DEFAULTANNOTS[self.ensight.ANNOT_TEXT]









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
    
    ensight_pp = ensight_class(session    = session, 
                               parameters = parameters, 
                               fpath      = '')

    ensight_pp.post_process()
else:

    # TODO BETTER WAY TO CHECK IF ENSIGHT EXISTS
    try:
        ensight_pp = ensight_class(ensight = ensight, parameters = parameters)

        ensight_pp.post_process()
    except:
        pass