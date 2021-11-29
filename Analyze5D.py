#!/usr/local/bin python3
"""addon for analyze5d"""

import os
from gi.repository import Gtk
import numpy as np
import h5py
import hdf5plugin
from matplotlib import pyplot, colors, cm
from matplotlib.figure import Figure
from matplotlib.backends.backend_gtk3agg import FigureCanvasGTK3Agg as FigureCanvas
from mpl_toolkits.axes_grid1 import make_axes_locatable
from misc_dude import * 
from multiprocessing import Pool, cpu_count

def Display2Data(Axe, x, y):
    """display coordinate to data coordinate"""

    return Axe.transData.inverted().transform(np.array([(x,y)]))[0]


class MainWindow:

    def Generate(self, widget, dude):

        ymax = int(dude.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.get_value())
        ymin = int(dude.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.get_value())
        xmax = int(dude.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.get_value())
        xmin = int(dude.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.get_value())   
        shifts = np.loadtxt(os.path.join(self.analysis_folder, "shifts.txt"), dtype=int)
        full = self.Full_ToggleButton.get_active()
        if float(self.TwoTheta0_Entry.get_text())<0:
            self.gamma = np.radians(np.genfromtxt(os.path.join(self.analysis_folder, "gamma.csv"), delimiter=',')[ymin:ymax+1, xmin:xmax+1])
            self.twotheta = np.radians(np.genfromtxt(os.path.join(self.analysis_folder, "twotheta.csv"), delimiter=',')[ymin:ymax+1, xmin:xmax+1])
        else:
            if dude.dimY == 516 and dude.dimX == 516:
                pixelsize = 55e-6
            elif dude.dimY == 195 and dude.dimX == 487:
                pixelsize = 125e-6
            elif dude.dimY == 1062 and dude.dimX == 1028:
                pixelsize = 75e-6
            twotheta0 = float(self.TwoTheta0_Entry.get_text())
            gamma0 = float(self.Gamma0_Entry.get_text())
            distance = float(self.Distance_Entry.get_text())
            self.twotheta, self.gamma = np.meshgrid(np.arange(xmin, xmax+1),np.arange(ymin, ymax+1))
            self.twotheta = ((self.twotheta-(xmin+xmax)/2.)*pixelsize/distance+np.radians(twotheta0))
            self.gamma = (((ymin+ymax)/2.-self.gamma)*pixelsize/distance+np.radians(gamma0))
        self.theta = np.radians(90-shifts[:,1]/1000.)[:,np.newaxis,np.newaxis]
        ystep = shifts[:,8].mean()
        print("ystep: ", ystep)
        xstep = shifts[:,9].mean()
        print("xstep:", xstep)
        aspect = np.sin(np.radians(float(self.Theta0_Entry.get_text())))*np.abs(ystep/xstep)
        print("aspect:", aspect)
        scan_ny = shifts[0,3]-shifts[0,2]+1
        scan_nx = shifts[0,5]-shifts[0,4]+1
        if full:
            print("using the full extend of the data")
            scan_ny += shifts[:,2].max()+(shifts[:,6]-1-shifts[:,3]).max()
            scan_nx += shifts[:,4].max()+(shifts[:,7]-1-shifts[:,5]).max()
        det_ny = ymax-ymin+1
        det_nx = xmax-xmin+1
        det_bin = int(self.Detector_Bin_Entry.get_text())
        self.data5d = np.zeros((shifts.shape[0], scan_ny, scan_nx, int(det_ny/det_bin), int(det_nx/det_bin)), dtype=np.float64)
        self.wavelength = 12.398/float(self.Energy_Entry.get_text())
        K = 2*np.pi/self.wavelength
        magic_angle = np.arccos(np.cos(self.twotheta)/np.cos(self.gamma))
        print(magic_angle.mean()/3.14*180)
        self.qy = K*np.sin(self.gamma)*np.ones(self.theta.shape)
        self.qx = K*np.cos(self.gamma)*np.sin(magic_angle)
        self.qz = K*np.cos(self.gamma)*np.cos(magic_angle)-K
        self.qz, self.qx = np.cos(self.theta)*self.qz-np.sin(self.theta)*self.qx, np.sin(self.theta)*self.qz+np.cos(self.theta)*self.qx

        if det_bin > 1:
            self.qy = self.qy[:,:int(det_ny/det_bin)*det_bin, :int(det_nx/det_bin)*det_bin].reshape(-1,int(det_ny/det_bin), det_bin, int(det_nx/det_bin), det_bin).mean(4).mean(2)
            self.qz = self.qz[:,:int(det_ny/det_bin)*det_bin, :int(det_nx/det_bin)*det_bin].reshape(-1,int(det_ny/det_bin), det_bin, int(det_nx/det_bin), det_bin).mean(4).mean(2)
            self.qx = self.qx[:,:int(det_ny/det_bin)*det_bin, :int(det_nx/det_bin)*det_bin].reshape(-1,int(det_ny/det_bin), det_bin, int(det_nx/det_bin), det_bin).mean(4).mean(2)
            self.twotheta = self.twotheta[:int(det_ny/det_bin)*det_bin, :int(det_nx/det_bin)*det_bin].reshape(int(det_ny/det_bin), det_bin, int(det_nx/det_bin), det_bin).mean(3).mean(1)
            self.gamma = self.gamma[:int(det_ny/det_bin)*det_bin, :int(det_nx/det_bin)*det_bin].reshape(int(det_ny/det_bin), det_bin, int(det_nx/det_bin), det_bin).mean(3).mean(1)

        for i in range(shifts.shape[0]):
            print("grabbing data from scan {0}".format(shifts[i,0]))
            if dude.eiger_enabled:
                h5_filename = [f.name for f in os.scandir(dude.h5_folder) if "scan_{0}".format(shifts[i,0]) in f.name]
                h5 = h5py.File(os.path.join(dude.h5_folder, h5_filename[0]), 'r')
                datatmp= h5['/entry/data/data'][:,ymin:ymax+1,xmin:xmax+1]
                datatmp[datatmp>2**31]=0
                h5.close()
            else:
                image_path = os.path.join(dude.Image_folder, str(shifts[i,0]))
                image_list = [imagefile.name for imagefile in os.scandir(image_path) if imagefile.name.endswith('.tif')]
                tif_index = np.array([int(filename.split(".")[0].split("_")[-1]) for filename in image_list])
                image_list = np.take(image_list , tif_index.argsort())
                image_list = list(map(lambda x:os.path.join(image_path,x), image_list))
                datatmp = np.zeros((len(image_list), det_ny, det_nx), dtype="float64")
                pool = Pool(processes=cpu_count()) 
                for j, tmparr in enumerate(pool.imap(image_loader, image_list)):
                    datatmp[j] = tmparr[ymin:ymax+1, xmin:xmax+1]
                pool.close()
                datatmp[datatmp<0] = 0
                
            datatmp = datatmp.reshape(shifts[i,6], shifts[i,7], det_ny, det_nx)
            if full:
                datatmp = np.pad(datatmp,((shifts[:,2].max()-shifts[i,2],shifts[i,3]-shifts[:,3].min()), (shifts[:,4].max()-shifts[i,4],shifts[i,5]-shifts[:,5].min()), (0,0), (0,0)))
            else:
                datatmp = datatmp[shifts[i][2]:shifts[i][3]+1,shifts[i][4]:shifts[i][5]+1]
            if det_bin > 1:
                datatmp = datatmp[:,:,:int(det_ny/det_bin)*det_bin, :int(det_nx/det_bin)*det_bin].reshape(scan_ny,scan_nx,int(det_ny/det_bin), det_bin, int(det_nx/det_bin), det_bin).sum(5).sum(3)
            self.data5d[i] = datatmp

        if self.SwapXY_ToggleButton.get_active():
            print("swapping scanning directions")
            self.data5d = np.swapaxes(self.data5d, 1, 2)
        if self.FlipX_ToggleButton.get_active():
            print("flipping image horizontally")
            self.data5d = np.flip(self.data5d, 2)
        if self.FlipY_ToggleButton.get_active():
            print("flipping image vertically")
            self.data5d = np.flip(self.data5d, 1)

        print("writing to data5D.npz, this might take a while")
        np.savez(os.path.join(self.analysis_folder, "data5D.npz"), data=self.data5d, qx=self.qx, qy=self.qy, qz=self.qz, twotheta=self.twotheta, gamma=self.gamma, theta=self.theta)
        self.data_12 = self.data5d.sum(4).sum(3).sum(0)
        self.Result_Image = self.Result_Axe.imshow(self.data_12, cmap="viridis", aspect=aspect, origin='lower', interpolation = 'nearest')
        self.Result_Canvas.draw()
        self.Threshold1_HScale_Adjustment.handler_block(self.Threshold1_Changed_Handler)
        self.Threshold1_HScale_Adjustment.set_lower(self.data_12.min())
        self.Threshold1_HScale_Adjustment.set_upper(self.data_12.max())
        self.Threshold1_HScale_Adjustment.set_value(self.data_12.min())
        self.Threshold1_HScale_Adjustment.handler_unblock(self.Threshold1_Changed_Handler)
        self.mask1 = np.zeros(self.data_12.shape)
        self.Threshold2_HScale_Adjustment.handler_block(self.Threshold2_Changed_Handler)
        self.Threshold2_HScale_Adjustment.set_lower(self.data_12.min())
        self.Threshold2_HScale_Adjustment.set_upper(self.data_12.max())
        self.Threshold2_HScale_Adjustment.set_value(self.data_12.min())
        self.Threshold2_HScale_Adjustment.handler_unblock(self.Threshold2_Changed_Handler)
        self.mask2 = np.zeros(self.data_12.shape)


    def Save(self, widget):
        
        np.savez(os.path.join(self.analysis_folder, "Result.npz"), I=self.data_12, mask1=self.mask1, d=self.d_spacing.data, \
                            tilt_lr=self.tilt_lr.data, tilt_ud=self.tilt_ud.data, mask2=self.mask2)
        self.Result_Figure.savefig(os.path.join(self.analysis_folder, "Result.png"), dpi=300)
        self.Result_Figure.savefig(os.path.join(self.analysis_folder, "Result.svg"), dpi=300, transparent=True)
    

    def Load(self, widget, dude):

        d5 = np.load(os.path.join(self.analysis_folder, "data5D.npz"))
        shifts = np.loadtxt(os.path.join(self.analysis_folder, "shifts.txt"), dtype=int)
        ystep = shifts[:,8].mean()
        print("ystep: ", ystep)
        xstep = shifts[:,9].mean()
        print("xstep:", xstep)
        aspect = np.sin(np.radians(float(self.Theta0_Entry.get_text())))*np.abs(ystep/xstep)
        print("aspect:", aspect)
        self.wavelength = 12.398/float(self.Energy_Entry.get_text())
        self.data5d = d5["data"]
        self.qx = d5["qx"]
        self.qy = d5["qy"]
        self.qz = d5["qz"]
        self.twotheta = d5["twotheta"]
        self.gamma = d5["gamma"]
        self.theta = d5["theta"]
        self.data_12 = self.data5d.sum(4).sum(3).sum(0)
        self.Result_Axe.cla()
        self.Result_Axe.set_axis_off()
        self.cax.set_axis_off()
        self.Result_Image = self.Result_Axe.imshow(self.data_12, cmap="viridis", aspect=aspect, origin='lower', interpolation = 'nearest')
        self.Threshold1_HScale_Adjustment.handler_block(self.Threshold1_Changed_Handler)
        self.Threshold1_HScale_Adjustment.set_lower(self.data_12.min())
        self.Threshold1_HScale_Adjustment.set_upper(self.data_12.max())
        self.Threshold1_HScale_Adjustment.set_value(self.data_12.min())
        self.Threshold1_HScale_Adjustment.handler_unblock(self.Threshold1_Changed_Handler)
        self.mask1 = np.zeros(self.data_12.shape)
        self.Threshold2_HScale_Adjustment.handler_block(self.Threshold2_Changed_Handler)
        self.Threshold2_HScale_Adjustment.set_lower(self.data_12.min())
        self.Threshold2_HScale_Adjustment.set_upper(self.data_12.max())
        self.Threshold2_HScale_Adjustment.set_value(self.data_12.min())
        self.Threshold2_HScale_Adjustment.handler_unblock(self.Threshold2_Changed_Handler)
        self.mask2 = np.zeros(self.data_12.shape)
        self.Result_Canvas.draw()
        try:
            self.Result_Canvas.mpl_disconnect(self.Result_Canvas_Mouse_Hover_Event)
            self.Result_Canvas.mpl_disconnect(self.Result_Canvas_Button_Pressed_Event)
        except:
            pass
        
        
    def Calculate(self, widget):

        data5d_sum = self.data5d.sum(axis=(0,3,4))
        self.qx_com = (self.data5d * self.qx[:,np.newaxis,np.newaxis,:,:]).sum(axis=(0,3,4))/data5d_sum
        self.qy_com = (self.data5d * self.qy[:,np.newaxis,np.newaxis,:,:]).sum(axis=(0,3,4))/data5d_sum
        self.qz_com = (self.data5d * self.qz[:,np.newaxis,np.newaxis,:,:]).sum(axis=(0,3,4))/data5d_sum

        # check for nan(s)
        self.qx_com[np.isnan(self.qx_com)] = self.qx_com[~np.isnan(self.qx_com)].mean()
        self.qy_com[np.isnan(self.qy_com)] = self.qy_com[~np.isnan(self.qy_com)].mean()
        self.qz_com[np.isnan(self.qz_com)] = self.qz_com[~np.isnan(self.qz_com)].mean()

        q_com = np.sqrt(self.qx_com**2+self.qy_com**2+self.qz_com**2)


        self.d_spacing = np.ma.MaskedArray(data = 2*np.pi/q_com, mask = self.mask1)
        #self.twotheta_com = (self.data5d.sum(0) * self.twotheta).sum(axis=(2,3))/self.data5d.sum(axis=(0,3,4))
        #self.d_spacing = np.ma.MaskedArray(data = self.wavelength/2/np.sin(self.twotheta_com/2), mask = self.mask1)
 
        self.dmin_HScale_Adjustment.handler_block(self.dmin_Changed_Handler)
        self.dmax_HScale_Adjustment.handler_block(self.dmax_Changed_Handler)
        self.dmin_HScale_Adjustment.set_lower(self.d_spacing.data.min()-1e-4)
        self.dmin_HScale_Adjustment.set_upper(self.d_spacing.data.max())
        self.dmin_HScale_Adjustment.set_value(self.d_spacing.min())
        self.dmax_HScale_Adjustment.set_lower(self.d_spacing.data.min())
        self.dmax_HScale_Adjustment.set_upper(self.d_spacing.data.max()+1e-4)
        self.dmax_HScale_Adjustment.set_value(self.d_spacing.max())
        self.dmin_HScale_Adjustment.handler_unblock(self.dmin_Changed_Handler)
        self.dmax_HScale_Adjustment.handler_unblock(self.dmax_Changed_Handler)
        self.Result_Image.set_array(self.d_spacing)
        self.Result_Image.set_clim(self.d_spacing.min(), self.d_spacing.max())
        self.Result_Image.set_cmap(self.cmap)
        self.cb = self.Result_Figure.colorbar(self.Result_Image, cax=self.cax, orientation='vertical', format="%.4f")
        self.cax.set_axis_on()

        #tilt_ud = (self.data5d.sum(0) * self.gamma).sum(axis=(2,3))/self.data5d.sum(axis=(0,3,4)) # angular version of the calculation
        tilt_ud = np.arcsin(self.qy_com/q_com)
        tilt_ud -= tilt_ud.flatten()[self.data_12.argmax()]
        self.tilt_ud = np.ma.MaskedArray(data = tilt_ud, mask = self.mask2)
        
        #tilt_lr = (self.data5d.sum(axis=(3,4))*self.theta).sum(0)/self.data5d.sum(axis=(0,3,4)) # angular version of the calculation
        tilt_lr = np.arctan2(-self.qx_com,-self.qz_com)
        tilt_lr -= tilt_lr.flatten()[self.data_12.argmax()]
        self.tilt_lr = np.ma.MaskedArray(data = tilt_lr, mask = self.mask2)

        self.quiver_changed(None)

        self.Result_Canvas_Mouse_Hover_Event = self.Result_Canvas.mpl_connect('motion_notify_event', self.Result_Canvas_Mouse_Hover)
        self.Result_Canvas_Button_Press_Event = self.Result_Canvas.mpl_connect('button_press_event', self.Result_Canvas_Button_Pressed)


    def quiver_changed(self, widget):

        try:
            self.Result_Quiver.remove()
            del self.Result_Quiver
        except Exception as e:
            pass
        x, y = np.meshgrid(np.arange(self.d_spacing.shape[1]), np.arange(self.d_spacing.shape[0]))
        bins = int(self.Tilt_Bin_Entry.get_text())
        ny = int(y.shape[0]/bins)
        nx = int(y.shape[1]/bins)

        self.Result_Quiver = self.Result_Axe.quiver(x[:ny*bins,:nx*bins].reshape(ny,bins,nx,bins).mean(3).mean(1), \
                                                    y[:ny*bins,:nx*bins].reshape(ny,bins,nx,bins).mean(3).mean(1), \
                                                    (self.tilt_lr)[:ny*bins,:nx*bins].reshape(ny,bins,nx,bins).mean(3).mean(1)+float(self.Uoffset_Entry.get_text()),\
                                                    (self.tilt_ud)[:ny*bins,:nx*bins].reshape(ny,bins,nx,bins).mean(3).mean(1)+float(self.Voffset_Entry.get_text()),\
                                                    scale = float(self.Tilt_Scale_Entry.get_text()),
                                                    width = float(self.Tilt_Width_Entry.get_text()), scale_units = "inches")
        self.Result_Canvas.draw()
        
 
    def d_limit_changed(self, widget):

        self.Result_Image.set_clim(self.dmin_HScale_Adjustment.get_value(), self.dmax_HScale_Adjustment.get_value())
        self.Result_Canvas.draw()


    def set_threshold_otsu(self, widget):

        self.Threshold1_HScale_Adjustment.set_value(threshold_otsu(self.data_12))
        self.Threshold2_HScale_Adjustment.set_value(threshold_otsu(self.data_12))
    

    def threshold1_changed(self, widget):

        self.mask1 = self.data_12 < self.Threshold1_HScale_Adjustment.get_value()
        self.d_spacing.mask = self.mask1
        self.Result_Image.set_array(self.d_spacing)
        self.Result_Canvas.draw()

        
    def threshold2_changed(self, widget):

        self.mask2 = self.data_12 < self.Threshold2_HScale_Adjustment.get_value()
        self.tilt_ud.mask = self.mask2
        self.tilt_lr.mask = self.mask2
        self.quiver_changed(None)

    
    def Result_Canvas_Figure_Entered(self, event):

        #Otherwise button_press_event.key produces none
        self.Result_Canvas.grab_focus()


    def Result_Canvas_Button_Pressed(self, event):

        x, y = list(map(lambda x: int(round(x, 0)), (Display2Data(self.Result_Axe, event.x, event.y))))
        dim_y, dim_x = self.data_12.shape

        if (y < dim_y-0.5) and (y > -.5) and (x > -0.5) and (x < dim_x-0.5):
            if event.button == 3 and event.key != "shift":
                self.mask1[y,x] = ~self.mask1[y,x]
                self.d_spacing.mask = self.mask1
                self.Result_Image.set_array(self.d_spacing)
                self.Result_Canvas.draw()
            elif event.button == 3 and event.key == "shift":
                self.mask2[y,x] = True
                self.tilt_ud.mask = self.mask2
                self.tilt_lr.mask = self.mask2
                self.quiver_changed(None)

    def Result_Canvas_Mouse_Hover(self, event):

        x, y = list(map(lambda x: int(round(x, 0)), (Display2Data(self.Result_Axe, event.x, event.y))))
        dim_y, dim_x = self.data_12.shape

        if (y < dim_y-0.5) and (y > -.5) and (x > -0.5) and (x < dim_x-0.5):
            self.label.set_text('[Intensity={0:.3f}]  [d-spacing={1:.3f}]  [tilt-lr={2:.3f}]  [tilt-ud={3:.3f}]'.format(self.data_12[y,x], self.d_spacing.data[y,x], self.tilt_lr.data[y,x], self.tilt_ud.data[y,x]))
            

    def MainWindow_Destroy(self, widget, dude): 

        self.data5D = None
        dude.Analyze5D = None
        self.win.destroy()

            
    def __init__(self, dude):

        self.cmap = cm.coolwarm.copy()
        self.cmap.set_under("w")
        self.cmap.set_over("w")

        self.label = Gtk.Label()

        self.Result_Figure = Figure()
        self.Result_Axe = self.Result_Figure.add_axes([0, 0.02, 0.88, .96])
        self.Result_Axe.set_axis_off()
        divider = make_axes_locatable(self.Result_Axe)
        self.cax = divider.append_axes("right", size="5%", pad=0.05)
        self.cax.set_axis_off()
        self.Result_Canvas = FigureCanvas(self.Result_Figure)
        self.Result_ScrolledWindow = Gtk.ScrolledWindow()
        self.Result_ScrolledWindow.set_policy(Gtk.PolicyType.AUTOMATIC, Gtk.PolicyType.AUTOMATIC)
        self.Result_ScrolledWindow.add(self.Result_Canvas)
        self.Result_ScrolledWindow_EventBox = Gtk.EventBox()
        self.Result_ScrolledWindow_EventBox.add(self.Result_ScrolledWindow)
        self.Result_ScrolledWindow.set_size_request(600,600)

        self.Result_Canvas_Figure_Enter_Event = self.Result_Canvas.mpl_connect('figure_enter_event', self.Result_Canvas_Figure_Entered)
        
        
        Energy_Label = Gtk.Label(" E (keV) ")
        self.Energy_Entry = Gtk.Entry()
        self.Energy_Entry.set_text(str(10))
        self.Energy_Entry.set_width_chars(5)

        self.Full_ToggleButton = Gtk.ToggleButton("Full")
        self.Full_ToggleButton.set_tooltip_text("active: full extend of the 2D data will be used, original size of data should be equal. \ninactive: useable area with complete overlap, original size of data does not need to be equal.")
        
        self.SwapXY_ToggleButton = Gtk.ToggleButton("X-Y")
        self.SwapXY_ToggleButton.set_tooltip_text("toggle this when you are scanning samy or hybridy as the inner loop motor.")
        self.FlipX_ToggleButton = Gtk.ToggleButton("-X")
        self.FlipX_ToggleButton.set_tooltip_text("toggle this when you are scanning hybridx from negative to positive, or attoz from positive to negative.")
        self.FlipY_ToggleButton = Gtk.ToggleButton("-Y")
        self.FlipY_ToggleButton.set_tooltip_text("toggle this when you are scanning hybridy from positive to negative, or samy from negative to positive.")

        Generate_Button = Gtk.Button("Gen")
        Generate_Button.set_tooltip_text("Take the shift correction file and angular values of each pixel, to produce the 5D space.")
        Generate_Button.connect("clicked", self.Generate, dude)

        Load_Button = Gtk.Button("Load")
        Load_Button.set_tooltip_text("Load data5D.npz file")
        Load_Button.connect("clicked", self.Load, dude)

        Calculate_Button = Gtk.Button("Calc")
        Calculate_Button.connect("clicked", self.Calculate)

        Save_Button = Gtk.Button("Save")
        Save_Button.set_tooltip_text("save results into Results.npz and figure into Result.png")
        Save_Button.connect("clicked", self.Save)

        #Threshold_Otsu_Button = Gtk.Button("Otsu")
        #Threshold_Otsu_Button.connect("clicked", self.set_threshold_otsu)
        self.Threshold1_HScale_Adjustment = Gtk.Adjustment(value = 1, lower = 100000, upper = 200000, step_increment = 1, page_increment = 100, page_size = 0)
        self.Threshold1_Changed_Handler = self.Threshold1_HScale_Adjustment.connect("value_changed", self.threshold1_changed)
        Threshold1_HScale = Gtk.HScale()
        Threshold1_HScale.set_adjustment(self.Threshold1_HScale_Adjustment)
        Threshold1_HScale.set_size_request(150, 20)
        Threshold1_HScale.set_value_pos(Gtk.PositionType.LEFT)
        Threshold1_HScale.set_digits(0)
        Threshold1_HScale.set_tooltip_text("thresholding for d-spacing plot")
        self.Threshold2_HScale_Adjustment = Gtk.Adjustment(value = 1, lower = 100000, upper = 200000, step_increment = 1, page_increment = 100, page_size = 0)
        self.Threshold2_Changed_Handler = self.Threshold2_HScale_Adjustment.connect("value_changed", self.threshold2_changed)
        Threshold2_HScale = Gtk.HScale()
        Threshold2_HScale.set_adjustment(self.Threshold2_HScale_Adjustment)
        Threshold2_HScale.set_size_request(150, 20)
        Threshold2_HScale.set_value_pos(Gtk.PositionType.LEFT)
        Threshold2_HScale.set_digits(0)
        Threshold2_HScale.set_tooltip_text("thresholding for lattice tilt plot")


        dmin_Label = Gtk.Label(" dmin ")
        self.dmin_HScale_Adjustment = Gtk.Adjustment(value = 1, lower = 1, upper = 20, step_increment = .0001, page_increment = .01, page_size = 0)
        self.dmin_Changed_Handler = self.dmin_HScale_Adjustment.connect("value_changed", self.d_limit_changed)
        dmin_HScale = Gtk.HScale()
        dmin_HScale.set_adjustment(self.dmin_HScale_Adjustment)
        dmin_HScale.set_size_request(150, 20)
        dmin_HScale.set_value_pos(Gtk.PositionType.LEFT)
        dmin_HScale.set_digits(4)

        dmax_Label = Gtk.Label(" dmax ")
        self.dmax_HScale_Adjustment = Gtk.Adjustment(value = 1, lower = 1, upper = 20, step_increment = .0001, page_increment = .01, page_size = 0)
        self.dmax_Changed_Handler = self.dmax_HScale_Adjustment.connect("value_changed", self.d_limit_changed)
        dmax_HScale = Gtk.HScale()
        dmax_HScale.set_adjustment(self.dmax_HScale_Adjustment)
        dmax_HScale.set_size_request(150, 20)
        dmax_HScale.set_value_pos(Gtk.PositionType.LEFT)
        dmax_HScale.set_digits(4)

        Tilt_Bin_Label = Gtk.Label(" bin ")
        self.Tilt_Bin_Entry = Gtk.Entry()
        self.Tilt_Bin_Entry.set_width_chars(5)
        self.Tilt_Bin_Entry.set_text('2')
        self.Tilt_Bin_Entry.connect("changed", self.quiver_changed)
        Tilt_Scale_Label = Gtk.Label(" scale ")
        self.Tilt_Scale_Entry = Gtk.Entry()
        self.Tilt_Scale_Entry.set_tooltip_text("0.1 means a tilt of 1 mrad will be 0.01 inches long, keep this consistent when plotting multiple data")
        self.Tilt_Scale_Entry.set_width_chars(5)
        self.Tilt_Scale_Entry.set_text('0.1')
        self.Tilt_Scale_Entry.connect("changed", self.quiver_changed)
        Tilt_Width_Label = Gtk.Label(" width ")
        self.Tilt_Width_Entry = Gtk.Entry()
        self.Tilt_Width_Entry.set_width_chars(7)
        self.Tilt_Width_Entry.set_text('0.005')
        self.Tilt_Width_Entry.connect("changed", self.quiver_changed)
        Uoffset_Label = Gtk.Label(" dU ")
        self.Uoffset_Entry = Gtk.Entry()
        self.Uoffset_Entry.set_width_chars(7)
        self.Uoffset_Entry.set_text('0.000')
        self.Uoffset_Entry.connect("changed", self.quiver_changed)
        Voffset_Label = Gtk.Label(" dV ")
        self.Voffset_Entry = Gtk.Entry()
        self.Voffset_Entry.set_width_chars(7)
        self.Voffset_Entry.set_text('0.000')
        self.Voffset_Entry.connect("changed", self.quiver_changed)
        
        Theta0_Label = Gtk.Label(" Theta0 ")
        self.Theta0_Entry = Gtk.Entry()
        self.Theta0_Entry.set_tooltip_text("aspect ratio will be stretched accordingly, 90 means no stretching")
        self.Theta0_Entry.set_width_chars(5)
        self.Theta0_Entry.set_text('90')
        TwoTheta0_Label = Gtk.Label(" TwoTheta0 ")
        self.TwoTheta0_Entry = Gtk.Entry()
        self.TwoTheta0_Entry.set_tooltip_text("negative number means loading from twotheta.txt, if not this defines the twotheta of the center of the RoI")
        self.TwoTheta0_Entry.set_width_chars(5)
        self.TwoTheta0_Entry.set_text('-1')
        Gamma0_Label = Gtk.Label(" Gamma0 ")
        self.Gamma0_Entry = Gtk.Entry()
        self.Gamma0_Entry.set_tooltip_text("if twotheta0 is negative means loading from gamma.txt, if not this defines the gamma of the center of the RoI")
        self.Gamma0_Entry.set_width_chars(5)
        self.Gamma0_Entry.set_text('-1')
        Distance_Label = Gtk.Label(" distance ")
        self.Distance_Entry = Gtk.Entry()
        self.Distance_Entry.set_tooltip_text("if twotheta0 is negative  means ignore this value, if not this defines the distance of the center of the RoI in meters")
        self.Distance_Entry.set_width_chars(5)
        self.Distance_Entry.set_text('-1')
        Detector_Bin_Label = Gtk.Label(" bin ")
        self.Detector_Bin_Entry = Gtk.Entry()
        self.Detector_Bin_Entry.set_tooltip_text("detector binning to reduce RAM usage")
        self.Detector_Bin_Entry.set_width_chars(5)
        self.Detector_Bin_Entry.set_text('1')

        HBox1 = Gtk.HBox(homogeneous = False, spacing = 3)
        HBox1.pack_start(Energy_Label, False, False, 3)
        HBox1.pack_start(self.Energy_Entry, False, False, 3)
        HBox1.pack_start(self.Full_ToggleButton, False, False, 3)
        HBox1.pack_start(self.SwapXY_ToggleButton, False, False, 3)
        HBox1.pack_start(self.FlipX_ToggleButton, False, False, 3)
        HBox1.pack_start(self.FlipY_ToggleButton, False, False, 3)
        vsep = Gtk.VSeparator()
        HBox1.pack_start(vsep, True, True, 3)
        HBox1.pack_start(Generate_Button, False, False, 3)
        HBox1.pack_start(Load_Button, False, False, 3)
        vsep = Gtk.VSeparator()
        HBox1.pack_start(vsep, True, True, 3)
        HBox1.pack_start(Calculate_Button, False, False, 3)
        vsep = Gtk.VSeparator()
    
        HBox2 = Gtk.HBox(homogeneous = False, spacing = 3)
        HBox2.pack_start(dmin_Label, False, False, 3)
        HBox2.pack_start(dmin_HScale, True, True, 3)
        HBox2.pack_start(dmax_Label, False, False, 3)
        HBox2.pack_start(dmax_HScale, True, True, 3)

        HBox3 = Gtk.HBox(homogeneous = False, spacing = 3)
        HBox3.pack_start(Tilt_Bin_Label, False, False, 3)
        HBox3.pack_start(self.Tilt_Bin_Entry, False, False, 3)
        HBox3.pack_start(Tilt_Scale_Label, False, False, 3)
        HBox3.pack_start(self.Tilt_Scale_Entry, False, False, 3)
        HBox3.pack_start(Tilt_Width_Label, False, False, 3)
        HBox3.pack_start(self.Tilt_Width_Entry, False, False, 3)
        HBox3.pack_start(Uoffset_Label, False, False, 3)
        HBox3.pack_start(self.Uoffset_Entry, False, False, 3)
        HBox3.pack_start(Voffset_Label, False, False, 3)
        HBox3.pack_start(self.Voffset_Entry, False, False, 3)

        HBox4 = Gtk.HBox(homogeneous = False, spacing = 3)
        #HBox4.pack_start(Threshold_Otsu_Button, False, False, 3)
        HBox4.pack_start(Threshold1_HScale, True, True, 3)
        HBox4.pack_start(Threshold2_HScale, True, True, 3)
        #HBox4.pack_start(vsep, True, True, 3)
        HBox4.pack_start(Save_Button, False, False, 3)


        HBox5 = Gtk.HBox(homogeneous = False, spacing = 3)
        HBox5.pack_start(Theta0_Label, False, False, 3)
        HBox5.pack_start(self.Theta0_Entry, False, False, 3)
        HBox5.pack_start(TwoTheta0_Label, False, False, 3)
        HBox5.pack_start(self.TwoTheta0_Entry, False, False, 3)
        HBox5.pack_start(Gamma0_Label, False, False, 3)
        HBox5.pack_start(self.Gamma0_Entry, False, False, 3)
        HBox5.pack_start(Distance_Label, False, False, 3)
        HBox5.pack_start(self.Distance_Entry, False, False, 3)
        HBox5.pack_start(Detector_Bin_Label, False, False, 3)
        HBox5.pack_start(self.Detector_Bin_Entry, False, False, 3)

        VBox1 = Gtk.VBox(homogeneous = False, spacing = 3)
        VBox1.pack_start(self.label, False, False, 3)
        VBox1.pack_start(self.Result_ScrolledWindow_EventBox, False, False, 3)
        VBox1.pack_start(HBox1, False, False, 3)
        VBox1.pack_start(HBox5, False, False, 3)
        VBox1.pack_start(HBox2, False, False, 3)
        VBox1.pack_start(HBox3, False, False, 3)
        VBox1.pack_start(HBox4, False, False, 3)

        self.analysis_folder = os.path.join(os.path.abspath(os.path.join(dude.MDA_folder, os.pardir)),"Analysis")

        self.win = Gtk.Window()
        self.win.connect("destroy", self.MainWindow_Destroy, dude)
        self.win.add(VBox1)
        self.win.set_title("Analyze 5D")
        self.win.show_all()
