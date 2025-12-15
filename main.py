import ansys.pyensight.core as ens
#from ansys.pyensight.core import libuserd
import numpy as np
import os
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
import matplotlib.cm as cm
from fft_iso import FFT_ISO
from scipy.interpolate import LinearNDInterpolator
import pandas as pd

class ensight_class:
    def __init__(self, parameters, fpath, session = None, ensight = None):

        if session is not None:
            if ensight is not None:
                print('Both Ensight and Session Object Provided, using session object')
            self.session = session
            self.ensight = self.session.ensight
            self.sesion_provided = True
        elif ensight is not None:
            self.ensight = ensight
            self.sesion_provided = False
        else:
            raise ValueError('ERROR: No session or ensight object specified')
            

        self.eocore = self.ensight.objs.core
        self.eonums = self.ensight.objs.enums
        self.eoutil = self.ensight.utils
        self.parts = self.ensight.utils.parts
        self.views = self.ensight.utils.views
        self.query = self.ensight.utils.query
        self.parameters = parameters

        if self.sesion_provided:
            self.session.load_data('D:\\Uni_Projects\\PALM_Projects\\Testing\\Working_Geometry_Testing\\postprocessing\\data_droplets\\test_simple_vibrate-1-*.cas.h5',
                                result_file='D:\\Uni_Projects\\PALM_Projects\\Testing\\Working_Geometry_Testing\\postprocessing\\data_droplets\\test_simple_vibrate-1-*.dat.h5')
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
            self.ensight.data.shift_time(1.000000,0.000000,0.000000)
            self.ensight.solution_time.monitor_for_new_steps("off")
            self.ensight.data.replace(r"""/fs04/scratch2/pa16/submodels/simple_vibration_runs/adaptive_64_750nm/test_simple_vibrate-1-*.cas.h5 """)
            self.ensight.command.delay_refresh("ON")

        self.ensight.file.animation_format("mpeg4")
        self.ensight.solution_time.show_as("step")
        self.t0 = self.eocore.TIMEVALUES[0]
        self.tn = self.eocore.TIMEVALUES[-1]

        self.main_vp = self.eocore.VPORTS['Main Viewport'][0]
        self.main_vp.setattrs(dict([['BACKGROUNDTYPE',self.eonums.VPORT_CONS],
                                    ['CONSTANTRGB',[0.2,0.2,0.2]]]))

        self.fluid_part = self.eocore.PARTS['fluid'][0]
        self.symmetry_part = self.eocore.PARTS['symmetry'][0]
        self.solid_coupling_part = self.eocore.PARTS['solid_coupling'][0]
        self.outlet_part = self.eocore.PARTS['outlet'][0]
        self.iso_default = self.eocore.DEFAULTPARTS[self.ensight.PART_ISO_SURFACE]
        self.shape_default = self.eocore.DEFAULTANNOTS[self.ensight.ANNOT_SHAPE]
        self.line_default = self.eocore.DEFAULTANNOTS[self.ensight.ANNOT_LINE]
        self.gauge_default = self.eocore.DEFAULTANNOTS[self.ensight.ANNOT_GAUGE]
        self.text_default = self.eocore.DEFAULTANNOTS[self.ensight.ANNOT_TEXT]
        self.coords = self.eocore.VARIABLES['Coordinates'][0]
        self.velocity = self.eocore.VARIABLES["Velocity"][0]
        self.vf_water = self.eocore.VARIABLES["Volume_fraction_water"][0]
        self.vf_air = self.eocore.VARIABLES["Volume_fraction_air"][0]

        self.cycle_time = self.eocore.create_variable(name='cycle_time',
                                                      sources=[self.fluid_part],
                                                      value=f"Analysis_Time*{str(self.parameters['frequency'])}")
        self.vibration_state = self.eocore.create_variable(name='vibration_state',
                                                      sources=[self.fluid_part],
                                                      value=f'COS((2*PI*cycle_time)-PI)',
                                                      private=1)
        self.vibration_state_pos = self.eocore.create_variable(name='vibration_state_pos',
                                                      sources=[self.fluid_part],
                                                      value=f'GT(vibration_state,0)')
        self.vibration_state_neg = self.eocore.create_variable(name='vibration_state_neg',
                                                      sources=[self.fluid_part],
                                                      value=f'LT(vibration_state,0)')
        


    def set_iso_view(self):
        self.views.set_view_direction(1, 1, 1, name="isometric", up_axis=(0,0,1))
        self.ensight.view.perspective("OFF")
        self.ensight.scene.ground_plane_visible("OFF")
        self.ensight.view_transf.zoom(0.6)
        self.ensight.solution_time.time_annotation("ON")
        self.ensight.text.select_begin(0)
        self.ensight.text.change_text(fr"""Time = <\\cnst C1 "%.1f" cycle_time\\>/{self.parameters['n_cycles']} cycles """)

        self.top_arrow = self.shape_default.createannot('top_arrow')
        self.top_arrow.setattrs(dict([['TYPE',self.eonums.ANNOT_SHAPE_ARROW],
                                      ['HEIGHT',0.002],
                                      ['LENGTH',0.150],
                                      ['ARROWTIPLENGTH',0.08],
                                      ['ARROWTIPSIZE',1],
                                      ['LOCATIONX',0.024],
                                      ['LOCATIONY',0.6],
                                      ['ROTATIONALANGLE',270],
                                      ['RGB',[1,1,1]],
                                      ['FILL',True]]))
                                                        
        self.bottom_arrow = self.shape_default.createannot('bottom_arrow')
        self.bottom_arrow.setattrs(dict([['TYPE',self.eonums.ANNOT_SHAPE_ARROW],
                                         ['LENGTH',0.150],
                                         ['HEIGHT',0.002],
                                         ['LOCATIONX',0.024],
                                         ['LOCATIONY',0.3],
                                         ['ROTATIONALANGLE',90],
                                         ['ARROWTIPLENGTH',0.08],
                                         ['ARROWTIPSIZE',1],
                                         ['RGB',[1,1,1]],
                                         ['FILL',True]]))

        self.top_gauge = self.gauge_default.createannot(self.vibration_state_pos.DESCRIPTION)
        self.top_gauge.setattrs(dict([['BACKGROUND',True],
                                      ['BACKGROUNDRGB',[0.2,0.2,0.2]],
                                      ['LEVELRGB',[0.65,0.07,0.08]],
                                      ['BORDER',False],
                                      ['HEIGHT',0.13],
                                      ['LOCATIONX',0.025],
                                      ['LOCATIONY',0.452],
                                      ['MAXIMUM',1],
                                      ['MINIMUM',0],
                                      ['VALUE',False],
                                      ['WIDTH',0.02]]))
        
        self.bottom_gauge = self.gauge_default.createannot(self.vibration_state_neg.DESCRIPTION)
        self.bottom_gauge.setattrs(dict([['BACKGROUND',True],
                                         ['BACKGROUNDRGB',[0.65,0.07,0.08]],
                                         ['BORDER',False],
                                         ['HEIGHT',0.13],
                                         ['LOCATIONX',0.025],
                                         ['LOCATIONY',0.318],
                                         ['MAXIMUM',0],
                                         ['MINIMUM',-1],
                                         ['LEVELRGB',[0.2,0.2,0.2]],
                                         ['VALUE',False],
                                         ['WIDTH',0.02]]))
        
        self.text_1 = self.text_default.createannot('l')
        self.text_1.setattrs(dict([['JUSTIFICATION',self.eonums.TS_CENTER],
                                      ['SIZE',100],
                                      ['LOCATIONX',0.046],
                                      ['LOCATIONY',0.44],
                                      ['RGB',[0.2,0.2,0.2]]]))
        self.text_2 = self.text_default.createannot('l')
        self.text_2.setattrs(dict([['JUSTIFICATION',self.eonums.TS_CENTER],
                                      ['SIZE',100],
                                      ['LOCATIONX',0.046],
                                      ['LOCATIONY',0.212],
                                      ['RGB',[0.2,0.2,0.2]]]))

        self.centre_line = self.line_default.createannot()
        self.centre_line.setattrs(dict([['WIDTH',4],
                                         ['LOCATIONX1',0.024],
                                         ['LOCATIONY1',0.45],
                                         ['LOCATIONX2',0.036],
                                         ['LOCATIONY2',0.45],
                                         ['RGB',[1,1,1]]]))


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


    def basic_animation(self): # TODO

        pass


    def velocity_animation(self):
        
        self.velocity_iso_part = self.iso_default.createpart(name="velocity_iso", sources=self.fluid_part, attributes=[['VARIABLE',self.vf_water]])[0]
        self.velocity_iso_part.COLORBYPALETTE = self.velocity

        self.fluid_part.VISIBLE = False
        self.symmetry_part.VISIBLE = False
        self.solid_coupling_part.OPAQUENESS = 0.75
        self.solid_coupling_part.COLORBYRGB=[0.57,0.57,0.57]
        self.outlet_part.VISIBLE = False

        self.velocity_palette = self.velocity.PALETTE['Velocity<\\\\units>'][0]
        self.velocity_legend = self.velocity.LEGEND['Velocity<\\\\units>'][0]

        self.velocity_palette.set_range_to_over_time_minmax(self.t0[0],self.tn[0])
        self.velocity_palette.MINMAX = [0,np.ceil(self.velocity_palette.MINMAX[1])]
        self.velocity_palette.NLEVELS = 11
        self.velocity_legend.TYPE = self.eonums.FNC_CONST
        self.velocity_legend.FORMAT = '%0.1f'
        
        
        

        self.create_animation('vel_anim')
        self.velocity_iso_part.VISIBLE = False

        
    def shearrate_animation(self):
        
        self.ensight.view_transf.translate(-2e-04,0,0)
        self.ensight.view_transf.zoom(1.2)

        self.eocore.create_variable(name='temp_1',
                                    value='0*Volume_fraction_air@MLL/TT', 
                                    sources=[self.fluid_part], 
                                    private=1)
        self.eocore.create_variable(name='temp_2',
                                    value='1+temp_1', 
                                    sources=[self.fluid_part], 
                                    private=1)
        self.eocore.create_variable(name='water_threshold',
                                    value='IF_GT(Volume_fraction_water,0.75)', 
                                    sources=[self.fluid_part])
        self.eocore.create_variable(name='temp_shear',
                                    value='FluidShearMax(plist,Velocity,1.0,temp_1,temp_2,1.0)', 
                                    sources=[self.fluid_part], 
                                    private=1)

        self.shearrate = self.eocore.create_variable(name='shearrate', 
                                                     value='water_threshold*temp_shear@/T', 
                                                     sources=[self.fluid_part])
        self.shearrate_vf = self.eocore.create_variable(name='shearrate_vf', 
                                                        value='Volume_fraction_water*temp_shear@/T', 
                                                        sources=[self.fluid_part])
        

        self.shearrate_iso_part = self.iso_default.createpart(name="shearrate_iso", sources=self.fluid_part, attributes=[['VARIABLE',self.vf_water],['TYPE',self.eonums.ISO_SURF_SOLID],['CONSTRAINT',self.eonums.CLIP_CHOICE_GREATER]])[0]
        self.shearrate_iso_part.COLORBYPALETTE = self.shearrate

        self.fluid_part.VISIBLE = False
        self.symmetry_part.VISIBLE = False
        self.solid_coupling_part.VISIBLE = False
        self.outlet_part.VISIBLE = False

        self.shearrate_palette = self.shearrate.PALETTE[0]
        self.shearrate_legend = self.shearrate.LEGEND[0]

        self.shearrate_palette.set_range_to_over_time_minmax(self.t0[0],self.tn[0])
        self.shearrate_palette.SCALE_METHOD = self.eonums.PALETTE_SCALE_LOG

        self.shearrate_palette.MINMAX = [1, np.power(10,np.ceil(np.log10(self.shearrate_palette.MINMAX[1])))]
        self.shearrate_palette.NLEVELS = int(1+np.log10(np.round(self.shearrate_palette.MINMAX[1])))

        self.shearrate_legend.TYPE = self.eonums.FNC_CONST
        self.shearrate_legend.FORMAT = '%0.1e'
        self.shearrate_legend.HEIGHT = 0.45
        self.shearrate_legend.DESCRIPTION = 'Shear-rate <\\\\units>'

        self.max_shearrate_query = self.eoutil.query.create_temporal('Max shear rate vs cycle time',
                                                                     query_type=self.eoutil.query.TEMPORAL_MAXIMUM,
                                                                     part_list=[self.fluid_part],
                                                                     variable1=self.shearrate,
                                                                     variable2=self.cycle_time,
                                                                     new_plotter=False)

        self.max_shearrate_query.LINESTYLE = self.eonums.LINE_SOLID
        self.max_shearrate_query.LINETYPE = self.eonums.CURVE_LINE_CONNECT
        self.max_shearrate_query.MARKER = self.eonums.CURVE_MARKER_SQUARE
        self.max_shearrate_query.MARKERSCALE = 1

        self.max_shearrate_line_plot = self.eocore.defaultplot[0].createplotter()
        self.max_shearrate_query.addtoplot(self.max_shearrate_line_plot)
        self.max_shearrate_line_plot.rescale()
        self.max_shearrate_line_plot.PLOTTITLE = "Max Shear Rate of Liquid"
        self.max_shearrate_line_plot.AXISXTITLE = "Number of cycles (n)"
        self.max_shearrate_line_plot.AXISYTITLE = "Shear-rate (1/s)"
        self.max_shearrate_line_plot.LEGENDVISIBLE = False
        self.max_shearrate_line_plot.AXISXAUTOSCALE = True
        self.max_shearrate_line_plot.AXISYSCALE = self.eonums.TRUE
        self.max_shearrate_line_plot.AXISXLABELFORMAT = "%.1f"
        self.max_shearrate_line_plot.AXISXGRIDTYPE = 1
        self.max_shearrate_line_plot.AXISYGRIDTYPE = 1
        self.max_shearrate_line_plot.TIMEMARKER = True
        self.max_shearrate_line_plot.AXISYAUTOSCALE = False
        self.max_shearrate_line_plot.AXISYMIN = self.shearrate_palette.MINMAX[0]
        self.max_shearrate_line_plot.AXISYMAX = self.shearrate_palette.MINMAX[1]
        self.max_shearrate_line_plot.ORIGINX = 0.5
        self.max_shearrate_line_plot.ORIGINY = 0.55
        self.max_shearrate_line_plot.WIDTH = 0.48
        self.max_shearrate_line_plot.HEIGHT = 0.43
        self.max_shearrate_line_plot.TIMEMARKERRGB = [1,0,0]
        self.max_shearrate_line_plot.TIMEMARKERWIDTH = 3
        self.max_shearrate_line_plot.AXISXNUMGRID = 1+int(self.parameters['n_cycles']/10)
        self.max_shearrate_line_plot.AXISYNUMGRIDLOG = 1+int(np.round(np.log10(self.shearrate_palette.MINMAX[1]/self.shearrate_palette.MINMAX[0])))
        self.max_shearrate_line_plot.AXISYLABELFORMAT = '%.0e'
        self.max_shearrate_line_plot.AXISYSGRIDTYPE = self.eonums.XY_GRID_TICK
        self.max_shearrate_line_plot.AXISXSGRIDTYPE = self.eonums.XY_GRID_TICK
        self.max_shearrate_line_plot.AXISYNUMSUBGRIDLOG = 10
        self.max_shearrate_line_plot.AXISXNUMSUBGRID = 4


        self.create_animation('shear_anim')
        self.shearrate_iso_part.VISIBLE = False


    def droplet_sizing(self): # TODO
        
        self.eocore.create_variable(name='temp_coords',
                                    sources=[self.fluid_part],
                                    value="Coordinates",
                                    private=1)

        self.element_coords = self.eocore.create_variable(name='element_coords',
                                                          sources=[self.fluid_part],
                                                          value="NodeToElem(plist, temp_coords)")
        
        self.element_volume = self.eocore.create_variable(name='element_volume',
                                                          sources=[self.fluid_part],
                                                          value="EleSize(plist)")
        
        self.nodal_vf_water = self.eocore.create_variable(name='nodal_vf_water',
                                                          sources=[self.fluid_part],
                                                          value="ElemToNode(plist, Volume_fraction_water)")
        
        


        #for i in range(self.tn[0]):
        self.default_filter = self.eocore.DEFAULTPARTS[self.ensight.PART_FILTER_PART]
        self.filter = self.default_filter.createpart(name='filter', sources=[self.fluid_part], attributes = [['ELTFILTER1ACTIVE',True],
                                                                                                           ['ELTFILTER1VARIABLE',self.vf_water],
                                                                                                           ['ELTFILTER1VARCOMP',self.eonums.VECTOR_MAGNITUDE],
                                                                                                           ['ELTFILTER1TESTVALUE',0.75],
                                                                                                           ['ELTFILTER1TESTOP',self.eonums.ELE_FAILED_LESS]])[0]
        

        self.ensight.solution_time.show_as("step")
        self.ensight.solution_time.increment(1)
        self.ensight.solution_time.update_to_last()

        self.fluid_data_2 = self.fluid_part.get_values([self.coords, self.element_coords, self.element_volume, self.vf_water, self.nodal_vf_water], activate=1)
        self.fluid_data = self.filter.get_values([self.coords, self.element_coords, self.element_volume, self.vf_water, self.nodal_vf_water, self.cycle_time], activate=1)

        #self.ensight.part.select_begin(1)
        #a = self.fluid_data['ELEMENT_IDS'][100]
        #self.ensight.show_info.element(a)

        a = self.ensight.query_parts(parts=[self.fluid_part])
        print(a)
        
        self.bruh = self.fluid_data[self.coords]*1e6
        t = self.bruh[:,2] > 465
        self.bruh = self.bruh[t]
        clustering = DBSCAN(eps=5,min_samples=2).fit(self.bruh)


        f,ax = plt.subplots(1,2,subplot_kw={"projection": "3d"})
        ax[0].scatter(xs=self.bruh[:,0],ys=self.bruh[:,1],zs=self.bruh[:,2],c=self.fluid_data[self.nodal_vf_water][t])

        ax[1].scatter(xs=self.bruh[:,0],ys=self.bruh[:,1],zs=self.bruh[:,2],c=clustering.labels_)

        plt.show()

        print('a')

        #print(np.max(self.fluid_data[self.coords][:,2])*1e6,np.min(self.fluid_data[self.coords][:,2])*1e6)
        

        #self.fluid_data_2 = self.fluid_part.get_values([self.coords, self.vf_water], activate=1)

        
        #print(np.max(self.fluid_data_2[self.coords][:,2]),np.min(self.fluid_data_2[self.coords][:,2]))

        #with open('coords.npy','wb') as f1, open('vof.npy','wb') as f2:
        #    np.save(f1,self.fluid_coordinates)
        '''
        channel_width = np.linspace(start = -25,
                        stop = 25,
                        num = 65)
        
        channel_length = np.linspace(start = -250,
                        stop = 250,
                        num = 641)
        
        z = np.linspace(start = 300,
                        stop = 550,
                        num = 321)
        print('here1')
        X,Y,Z = np.meshgrid(channel_width,channel_length,z)
        print('here2')

        interp = LinearNDInterpolator(self.fluid_data[self.coords]*1e6,self.fluid_data[self.nodal_vf_water])
        print('here3')
        points = interp(X,Y,Z)
        print('here4')
        f,ax = plt.subplots(1,1, subplot_kw={"projection": "3d"})
        print('here5')
        ax.scatter(xs=X,ys=Y,zs=Z,c=points)
        print('here6')
        '''




    def total_flowrate(self): # TODO

        pass


    def fft_of_surface(self):

        
        fft_calculator = FFT_ISO(parameters = self.parameters)

        self.coord_iso = self.iso_default.createpart(name="coord_iso", sources=self.fluid_part, attributes=[['VARIABLE',self.vf_water]])[0]
        self.ensight.solution_time.show_as("step")
        self.ensight.solution_time.increment(1)
        self.ensight.solution_time.update_to_first()

        for i in range(self.tn[0]):
            self.iso_surface_coordinates = self.coord_iso.get_values([self.coords,'Analysis_Time'], activate=1)[self.coords]
            
            fft_calculator.solve(data=self.iso_surface_coordinates, time_data=[self.eocore.TIMESTEP+1,self.eocore.SOLUTIONTIME])

            fft_calculator.full_plot(f'FFT Plot t={i}', f'.\\images\\iso_plot_{int(self.eocore.TIMESTEP+1)}.png')
            
            self.ensight.solution_time.step_forward()
    

def self_distance(x):
    
    arr = np.zeros(x.shape[0])

    for i,row in enumerate(x):
        a = np.linalg.norm(x-row,axis=1)
        arr[i] = np.min(a[np.nonzero(a)])

    return arr

def post_process(ensight_pp, options = None): # TODO

    ensight_pp.set_iso_view()
    
    #ensight_pp.velocity_animation()

    #ensight_pp.shearrate_animation()

    #ensight_pp.fft_of_surface()

    ensight_pp.droplet_sizing()


# TODO, formalise this part of the script
parameters = {'frequency' : 1.63e6,
              'amplitude' : 1e-6,
              'n_cycles' : 60,
              'n_levels_refinement' : 4,
              'channel_width' : 50,
              'grid_size' : 500}


if __name__ == '__main__':
    session = ens.LocalLauncher(batch=True, ansys_installation='C:\\Program Files\\ANSYS Inc\\v251').start()#, use_egl=True, additional_command_line_options=['-v 5']).start()
    
    ensight_pp = ensight_class(session= session, parameters=parameters, fpath='')

    post_process(ensight_pp)
else:

    # TODO BETTER WAY TO CHECK IF ENSIGHT EXISTS
    try:
        ensight_pp = ensight_class(ensight = ensight, parameters = parameters)

        post_process(ensight_pp)
    except:
        pass

'''
a=libuserd.LibUserd()
a.initialize()
b=a.load_data('D:\\Uni_Projects\\PALM_Projects\\Testing\\Working_Geometry_Testing\\postprocessing\\data\\test_simple_vibrate-1-00323.cas.h5',
              result_file='D:\\Uni_Projects\\PALM_Projects\\Testing\\Working_Geometry_Testing\\postprocessing\\data\\test_simple_vibrate-1-00323.dat.h5',
              file_format='Fluent_HDF5')
c=b.parts()
conn=c[0].element_conn(libuserd.ElementType.HEX08)
conn.shape = (len(conn)//8, 8)
'''