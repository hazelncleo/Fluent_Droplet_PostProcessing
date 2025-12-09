import ansys.pyensight.core as ensight
import numpy as np
import os
from scipy.fft import rfft2, rfftfreq, fftn, fftshift
from scipy.interpolate import LinearNDInterpolator
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import inspect


class ensight_class:
    def __init__(self, session, parameters, fpath):
        self.session = session
        self.eocore = session.ensight.objs.core
        self.eonums = session.ensight.objs.enums
        self.eoutil = session.ensight.utils
        self.parts = session.ensight.utils.parts
        self.views = session.ensight.utils.views
        self.query = session.ensight.utils.query
        self.parameters = parameters

        self.session.load_data('D:\\Uni_Projects\\PALM_Projects\\Testing\\Working_Geometry_Testing\\postprocessing\\test_simple_vibrate-1-*.cas.h5',
                               result_file='D:\\Uni_Projects\\PALM_Projects\\Testing\\Working_Geometry_Testing\\postprocessing\\test_simple_vibrate-1-*.dat.h5')
        #self.session.load_data('D:\\Uni_Projects\\PALM_Projects\\Testing\\Working_Geometry_Testing\\adaptive_64\\test_simple_vibrate-1-*.cas.h5',
        #                       result_file='D:\\Uni_Projects\\PALM_Projects\\Testing\\Working_Geometry_Testing\\adaptive_64\\test_simple_vibrate-1-*.dat.h5')

        self.session.ensight.file.animation_format("mpeg4")
        self.session.ensight.solution_time.show_as("time")
        self.t0 = self.eocore.TIMEVALUES[0]
        self.tn = self.eocore.TIMEVALUES[-1]

        self.main_vp = self.eocore.VPORTS['Main Viewport'][0]
        self.main_vp.setattrs(dict([['BACKGROUNDTYPE',self.eonums.VPORT_CONS],
                                    ['CONSTANTRGB',[0.2,0.2,0.2]]]))

        self.fluid_part = self.eocore.PARTS['fluid'][0]
        self.symmetry_part = self.eocore.PARTS['symmetry'][0]
        self.solid_coupling_part = self.eocore.PARTS['solid_coupling'][0]
        self.outlet_part = self.eocore.PARTS['outlet'][0]
        self.iso_default = self.eocore.DEFAULTPARTS[session.ensight.PART_ISO_SURFACE]
        self.shape_default = self.eocore.DEFAULTANNOTS[self.session.ensight.ANNOT_SHAPE]
        self.line_default = self.eocore.DEFAULTANNOTS[self.session.ensight.ANNOT_LINE]
        self.gauge_default = self.eocore.DEFAULTANNOTS[self.session.ensight.ANNOT_GAUGE]
        self.text_default = self.eocore.DEFAULTANNOTS[self.session.ensight.ANNOT_TEXT]
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
        self.session.ensight.view.perspective("OFF")
        self.session.ensight.scene.ground_plane_visible("OFF")
        self.session.ensight.view_transf.zoom(0.6)
        self.session.ensight.solution_time.show_as("time")
        self.session.ensight.solution_time.time_annotation("ON")
        self.session.ensight.text.select_begin(0)
        self.session.ensight.text.change_text(fr"""Time = <\\cnst C1 "%.1f" cycle_time\\>/{self.parameters['n_cycles']} cycles """)

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

        self.session.ensight.solution_time.update_to_first()
        self.session.ensight.file.animation_rend_offscreen("ON")
        self.session.ensight.file.animation_numpasses(4)
        self.session.ensight.file.animation_stereo("current")
        self.session.ensight.file.animation_screen_tiling(1,1)
        self.session.ensight.file.animation_file(fr"""D:\Uni_Projects\PALM_Projects\Testing\Working_Geometry_Testing\postprocessing\{fname} """)
        self.session.ensight.file.animation_window_size("user_defined")
        self.session.ensight.file.animation_window_xy(3840,2160)
        self.session.ensight.solution_time.increment(1)
        self.session.ensight.file.animation_frames(self.tn[0])
        self.session.ensight.file.animation_start_number(0)
        self.session.ensight.file.animation_multiple_images("OFF")
        self.session.ensight.file.animation_raytrace_it("OFF")
        self.session.ensight.file.animation_play_flipbook("OFF")
        self.session.ensight.file.animation_play_time("ON")
        self.session.ensight.file.animation_play_keyframe("OFF")
        self.session.ensight.file.animation_reset_flipbook("OFF")
        self.session.ensight.file.animation_reset_traces("OFF")
        self.session.ensight.file.animation_reset_time("ON")
        self.session.ensight.file.animation_reset_keyframe("OFF")
        self.session.ensight.file.save_animation()


    def basic_animation(self):

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
        #self.session.ensight.view_transf.translate(-8e-05,8e-05,0)
        self.session.ensight.view_transf.translate(-2e-04,0,0)
        self.session.ensight.view_transf.zoom(1.2)

        self.eocore.create_variable(name='temp_1',
                                    value='0*Volume_fraction_air@MLL/TT', 
                                    sources=[self.fluid_part], 
                                    private=1)
        self.eocore.create_variable(name='temp_2',
                                    value='1+temp_1', 
                                    sources=[self.fluid_part], 
                                    private=1)
        self.eocore.create_variable(name='water_threshold',
                                    value='IF_GT(Volume_fraction_water,0.3)', 
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
        self.shearrate_palette.LIMIT_FRINGES = self.eonums.PALETTE_LIMIT_FRINGES_INVISIBLE

        self.shearrate_palette.MINMAX = [1, np.power(10,np.ceil(np.log10(self.shearrate_palette.MINMAX[1])))]
        self.shearrate_palette.NLEVELS = int(1+np.log10(np.round(self.shearrate_palette.MINMAX[1])))

        self.shearrate_legend.TYPE = self.eonums.FNC_CONST
        self.shearrate_legend.FORMAT = '%0.1e'

        self.max_shearrate_query = self.eoutil.query.create_temporal('Max shear rate vs cycle time',
                                                                     query_type=self.eoutil.query.TEMPORAL_MAXIMUM,
                                                                     part_list=[self.fluid_part],
                                                                     variable1=self.shearrate,
                                                                     variable2=self.cycle_time,
                                                                     new_plotter=False)

        self.max_shearrate_query.LINESTYLE = self.eonums.LINE_SOLID
        self.max_shearrate_query.LINETYPE = self.eonums.CURVE_LINE_CONNECT
        self.max_shearrate_query.MARKER = self.eonums.CURVE_MARKER_SQUARE
        self.max_shearrate_query.MARKERSCALE = 4

        self.max_shearrate_line_plot = self.eocore.defaultplot[0].createplotter()
        self.max_shearrate_query.addtoplot(self.max_shearrate_line_plot)
        self.max_shearrate_line_plot.rescale()
        self.max_shearrate_line_plot.PLOTTITLE = "Max Shear Rate of Liquid"
        self.max_shearrate_line_plot.AXISXTITLE = "Number of cycles (n)"
        self.max_shearrate_line_plot.AXISYTITLE = "Shear rate (1/s)"
        self.max_shearrate_line_plot.LEGENDVISIBLE = False
        self.max_shearrate_line_plot.AXISXAUTOSCALE = True
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


        self.create_animation('shear_anim')
        self.shearrate_iso_part.VISIBLE = False


    def droplet_sizing(self):

        pass


    def total_flowrate(self):

        pass


    def fft_of_surface(self):

        self.coord_iso = self.iso_default.createpart(name="coord_iso", sources=self.fluid_part, attributes=[['VARIABLE',self.vf_water]])[0]
        self.session.ensight.solution_time.show_as("step")
        self.session.ensight.solution_time.increment(1)
        self.session.ensight.solution_time.update_to_first()

        for i in range(10):
            self.iso_surface_coordinates = self.coord_iso.get_values([self.coords,'Analysis_Time'], activate=1)[self.coords]
            with open('iso_{}.npy'.format(i),'wb') as f:
                np.save(f,self.iso_surface_coordinates)
            
            self.session.ensight.solution_time.step_forward()
  
        if False:
            N = 641
            M = 65
            xf_vert = np.zeros((N,M))
            xf_hor = np.zeros((M,N))
            f,(ax,ax_f) = plt.subplots(2,3)

            x = np.linspace(-250e-6,250e-6,N)
            y = np.linspace(-250e-6,250e-6,N)
            X,Y = np.meshgrid(x,y)

            interp = LinearNDInterpolator(self.iso_surface_coordinates[:,:2], self.iso_surface_coordinates[:,2])
            xf_hor = interp(X[288:353,:],Y[288:353,:])
            xf_vert = interp(X[:,288:353],Y[:,288:353])

            
            ax[0].scatter(X[288:353,:],Y[288:353,:],s=0.01,c=xf_hor, vmin=475e-6,vmax=485e-6)
            ax[1].scatter(X[:,288:353],Y[:,288:353],s=0.01,c=xf_vert, vmin=475e-6,vmax=485e-6)
            ax[2].scatter(self.iso_surface_coordinates[:,0],self.iso_surface_coordinates[:,1],c=self.iso_surface_coordinates[:,2],s=0.1, vmin=475e-6,vmax=485e-6)
            ax[0].set_xlim(-250e-6,250e-6)
            ax[0].set_ylim(-250e-6,250e-6)
            ax[1].set_xlim(-250e-6,250e-6)
            ax[1].set_ylim(-250e-6,250e-6)
            ax[2].set_xlim(-250e-6,250e-6)
            ax[2].set_ylim(-250e-6,250e-6)
            
            #ax[0].plot(Y[:,150],normalised_distance[:,150],'-r')
            #ax[0].plot(Y[:,150],new_distance[:,150],'-k')
            #ax[0].set_xlim(-25e-6,25e-6)

            Z = fftn(xf_hor)[15:M//2+1,15:N//2+1]
            Z_2 = fftn(xf_vert)[15:N//2+1,15:M//2+1]
            freq_rows = rfftfreq(N, 50/64)[15:]
            freq_cols = rfftfreq(M, 50/64)[15:]

            F_R,F_C = np.meshgrid(freq_rows,freq_cols)

            wavelength_normed = np.sqrt(np.power(F_R,2)+np.power(F_C,2)).flatten()
            Z_flat = np.log(np.power(np.abs(Z),2)).flatten()
            ax_f[1].scatter(wavelength_normed,Z_flat)
            #ax_f[0].scatter(F_R,F_C,s=1,c=np.log(np.power(np.abs(Z),2)))
            #ax_f[0].set_xlim(0,)
            #ax_f[1].contour(F_R,F_C,np.log(np.power(np.abs(Z),2)))
            ax_f[0].imshow(np.log(np.power(np.abs(Z),2)))
            #ax_f[1].contour(F_C,F_R,np.log(np.power(np.abs(Z_2),2)))
            #ax_f[1].plot(1/freq_rows[15:],np.abs(np.real(Z[50,15:])))

            #ax[1].scatter(X,Y,s=10,c=new_distance)
            #ax[1].plot(lines[:,0],lines[:,1],'-r')
            #ax[1].imshow(np.real(Z), cmap=cm.gray)


            plt.show()
    




def main():
    
    session = ensight.LocalLauncher(batch=True).start()#, use_egl=True, additional_command_line_options=['-v 5']).start()

    parameters = {'n_cycles' : 60,
                 'frequency' : 1.63e6,
                 'amplitude' : 1e-6}
    
    ensight_pp = ensight_class(session, parameters, '')

    ensight_pp.set_iso_view()
    
    ensight_pp.velocity_animation()

    #ensight_pp.shearrate_animation()

    ensight_pp.fft_of_surface()
    

    '''
    
    '''






#ax[0].scatter([:,0], iso_surface_coordinates[coords][:,2])
#ax[0].set_xlim(-250e-6,250e-6)
#ax[0].set_ylim(450e-6,500e-6)


#freq = np.fft.rfft2(iso_surface_coordinates[coords], axes=[0,1])
#abs_freq = np.abs(freq)
#ax[1].plot(abs_freq[:,0], abs_freq[:,1])





'''



'''
#fluid_part.VISIBLE = False
#symmetry_part.VISIBLE = False
#solid_coupling_part.OPAQUENESS = 0.75
#outlet_part.VISIBLE = False
'''
session.ensight.function.palette("Velocity")
session.ensight.function.set_palette_to_minmax()
session.ensight.function.range(0,np.ceil(eocore.ANNOTS["Velocity<\\\\units>"][0].RANGE[1]))
session.ensight.function.number_of_levels(11)
session.ensight.function.type("banded")
session.ensight.legend.select_palette_begin("Velocity")
session.ensight.legend.format("%0.1f")
session.ensight.legend.select_palette_end()
'''


#session.ensight.file.animation_rend_offscreen("ON")
#session.ensight.file.animation_numpasses(4)
#session.ensight.file.animation_stereo("current")
#session.ensight.file.animation_screen_tiling(1,1)
#session.ensight.file.animation_file(r"""D:\Uni_Projects\PALM_Projects\Testing\Working_Geometry_Testing\postprocessing\velocity_anim """)
#session.ensight.file.animation_window_size("user_defined")
#session.ensight.file.animation_window_xy(3840,2160)
#session.ensight.file.animation_frames(eocore.TIMESTEP_LIMITS[1]+1)
#session.ensight.file.animation_start_number(0)
#session.ensight.file.animation_multiple_images("OFF")
#session.ensight.file.animation_raytrace_it("OFF")
#session.ensight.file.animation_play_flipbook("OFF")
#session.ensight.file.animation_play_time("ON")
#session.ensight.file.animation_play_keyframe("OFF")
#session.ensight.file.animation_reset_flipbook("OFF")
#session.ensight.file.animation_reset_traces("OFF")
#session.ensight.file.animation_reset_time("ON")
#session.ensight.file.animation_reset_keyframe("OFF")
#session.ensight.file.save_animation()

#iso_part.COLORBYPALETTE = None
#iso_part.COLORBYRGB = [0.0,0.6666666666666666,1.0]
#session.ensight.part.clone(0,2)

#image_data = session.render(1920, 1080, aa=4)
#with open("D:\\Uni_Projects\\PALM_Projects\\Testing\\Working_Geometry_Testing\\adaptive_32\\postprocessing\\image.png", "wb") as f:
#    f.write(image_data)



if __name__ == '__main__':
    main()