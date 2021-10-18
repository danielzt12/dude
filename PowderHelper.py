#!/usr/local/bin python3
"""addon for powder helper"""

from gi.repository import Gtk
import numpy as np
from matplotlib import pyplot, colors
from matplotlib.figure import Figure
from matplotlib.backends.backend_gtk3agg import FigureCanvasGTK3Agg as FigureCanvas
from scipy.ndimage import morphology
from math import *
from misc_dude import *
import os
import time

class MainWindow:

    def save(self, widget, dude):

        if self.parent.sparse_enabled:
            dat = self.parent.image[:, self.indices.flatten()==index+1]
        else:
            dat = self.parent.image[:, self.indices==index+1]
        gam = self.gamma[self.indices==index+1]
        nbin = 100
        gamrange = np.linspace(gam.min(), gam.max(), nbin)
        indi = np.digitize(gam, gamrange)
        binc = np.bincount(indi.flatten())
        I = np.zeros((dim0*dim1, nbin-1))
        for i in range(nbin-1):
            I[:,i] = dat[:,indi==i+1].sum(1)/binc[i+1]
        np.save("Q{0:.0f}.npz".format(self.Q[index]*1000), data=I.reshape(dim0,dim1,nbin-1), tilt=gamrange)        
        


    def Integrate(self, widget, dude):
    
        t0 = time.time()
        K = 2*np.pi/(12.398418746/float(self.Energy_Entry.get_text()))
        Q = np.fabs(K*np.sin(np.radians(self.twotheta)))
        nbin = int(self.Nbin_Entry.get_text())
        Qrange = np.linspace(Q.min(), Q.max(), nbin)
        self.I = np.zeros(nbin-1)
        self.indices = np.digitize(Q, Qrange)
        self.bincount = np.bincount(self.indices.flatten())
        for i in range(nbin-1):
            if dude.sparse_enabled:
                self.I[i] = dude.image[:,self.indices.flatten()==i+1].sum()/self.bincount[i+1]
            else:
                self.I[i] = dude.image[:,self.indices==i+1].sum()/self.bincount[i+1]
        self.Image_Axe.cla()
        self.Q = (Qrange[1:]+Qrange[:-1])/2.
        print (time.time()-t0)
        #self.Image_Axe.plot(self.Q, self.I, "k-.")
        self.Image_Axe.semilogy(self.Q, self.I, "k-.")
        self.Image_Axe.set_xlabel(r"Q ($\AA^{ -1}$)")
        self.Image_Canvas.draw()

        np.savetxt(os.path.join(os.path.abspath(os.path.join(dude.MDA_folder, os.pardir, 'Analysis', "{0}.csv".format(dude.MDA_File_ListStore[dude.mda_selection_path[0]][0])))), np.vstack((self.Q, self.I)))


    def Image_Canvas_Mouse_Hover(self, event):

        xpointer, ypointer = Display2Data(self.Image_Axe, event.x, event.y)
        try:
            index = np.fabs(self.Q - xpointer).argmin()
        except:
            return
        if self.parent.sparse_enabled:
            self.parent.Image_Image.set_array(np.multiply(self.image,(self.indices!=index+1)))
        else:
            self.parent.Image_Image.set_array(self.image*(self.indices!=index+1))
        self.parent.Image_Image.set_norm(colors.LogNorm())
        self.parent.Image_Canvas.draw()
        dims = self.parent.data[0]["dimensions"]
        ndim = len(dims)
        if dims[-1] == 2048:
            ndim -= 1

        if event.button == 3:
            if ndim == 2:
                if self.parent.sparse_enabled:
                    self.parent.Plot2D_Image.set_array(self.parent.image[:,self.indices.flatten()==index+1].sum(1).reshape(dims[0], dims[1]))
                else:
                    self.parent.Plot2D_Image.set_array(self.parent.image[:,self.indices==index+1].sum(1).reshape(dims[0], dims[1]))
                self.parent.Plot2D_Image.set_norm(colors.LogNorm())
                self.parent.Plot2D_Canvas.draw()
            elif ndim == 1:
                if not self.parent.Scan_ToolBox_Y_ToggleButton.get_active():
                    self.parent.Plot1D_Axe.cla()
                if self.parent.sparse_enabled:
                    self.parent.Plot1D_Axe.semilogy(np.array(self.parent.data[1].p[0].data), self.parent.image[:,self.indices.flatten()==index+1].sum(1))
                else:
                    self.parent.Plot1D_Axe.semilogy(np.array(self.parent.data[1].p[0].data), self.parent.image[:,self.indices==index+1].sum(1))
                self.parent.Plot1D_Canvas.draw()


    def MainWindow_Destroy(self, widget, dude): 

        dude.PowderHelper = None
        self.win.destroy()


    def __init__(self, dude):

        if not dude.AllImagesLoaded:
            dude.LoadAllImages(None)
        self.parent = dude
        if dude.sparse_enabled:
            self.image = dude.image.sum(0).reshape(1062,1028)
        else:
            self.image = dude.image.mean(0)
        dude.Image_Image.set_array(self.image)
        dude.Image_Image.set_norm(colors.LogNorm())
        dude.Image_Canvas.draw()
        self.Preset_indices = np.zeros((4,dude.dimY,dude.dimX), dtype=np.bool)
        self.Preset_bincount = np.ones(4)
        analysis_folder = os.path.join(dude.MDA_folder[:-4], "Analysis")
        self.twotheta = np.loadtxt(os.path.join(analysis_folder, "twotheta.csv"), delimiter=",")
        self.gamma = np.loadtxt(os.path.join(analysis_folder, "gamma.csv"), delimiter=",")
        self.Image_Figure = Figure()
        self.Image_Axe = self.Image_Figure.add_axes([0.1, 0.1, .85, .85])
        self.Image_Axe.set_xlabel(r"Q ($\AA^{ -1}$)")
        self.Image_Canvas = FigureCanvas(self.Image_Figure)
        self.Image_ScrolledWindow = Gtk.ScrolledWindow()
        self.Image_ScrolledWindow.set_policy(Gtk.PolicyType.AUTOMATIC, Gtk.PolicyType.AUTOMATIC)
        self.Image_ScrolledWindow.add(self.Image_Canvas)
        self.Image_ScrolledWindow_EventBox = Gtk.EventBox()
        self.Image_ScrolledWindow_EventBox.add(self.Image_ScrolledWindow)
        self.Image_ScrolledWindow.set_size_request(700,500)
        self.Image_Canvas_Mouse_Hover_Event = self.Image_Canvas.mpl_connect('motion_notify_event', self.Image_Canvas_Mouse_Hover)
        self.Image_Canvas_Button_Press_Event = self.Image_Canvas.mpl_connect('button_press_event', self.Image_Canvas_Mouse_Hover) 

        Energy_Label = Gtk.Label(" E(keV) ")
        self.Energy_Entry = Gtk.Entry()
        self.Energy_Entry.set_text("10.0")
        self.Energy_Entry.set_width_chars(5)
        Nbin_Label = Gtk.Label(" Nbin ")
        self.Nbin_Entry = Gtk.Entry()
        self.Nbin_Entry.set_text("100")
        self.Nbin_Entry.set_width_chars(5)

        #self.Threshold_Spin_Adjustment = Gtk.Adjustment(value = 0.01, lower = 0, upper = 10, step_increment = 0.001, page_increment = 0.01, page_size = 0)
        #Threshold_SpinButton = Gtk.SpinButton()
        #Threshold_SpinButton.set_adjustment(self.Threshold_Spin_Adjustment)
        #Threshold_SpinButton.set_digits(3)

        self.Preset0_RadioButton = Gtk.RadioButton(label = " 0 ")
        self.Preset1_RadioButton = Gtk.RadioButton(label = " 1 ", group = self.Preset0_RadioButton)
        self.Preset2_RadioButton = Gtk.RadioButton(label = " 2 ", group = self.Preset0_RadioButton)
        self.Preset3_RadioButton = Gtk.RadioButton(label = " 3 ", group = self.Preset0_RadioButton)
        self.Preset4_RadioButton = Gtk.RadioButton(label = " 4 ", group = self.Preset0_RadioButton)
        
        Integrate_Button = Gtk.Button(" Integrate ")
        Integrate_Button.connect("clicked", self.Integrate, dude)

        Save_Button = Gtk.Button(" Save ")
        Save_Button.connect("clicked", self.save, dude)

        HBox1 = Gtk.HBox(homogeneous = False, spacing = 3)
        HBox1.pack_start(Energy_Label, False, False, 0)
        HBox1.pack_start(self.Energy_Entry, False, False, 0)
        HBox1.pack_start(Nbin_Label, False, False, 0)
        HBox1.pack_start(self.Nbin_Entry, False, False, 0)
        #HBox1.pack_start(Threshold_SpinButton, False, False, 0)
        #HBox1.pack_start(self.Preset0_RadioButton, False, False, 0)
        #HBox1.pack_start(self.Preset1_RadioButton, False, False, 0)
        #HBox1.pack_start(self.Preset2_RadioButton, False, False, 0)
        #HBox1.pack_start(self.Preset3_RadioButton, False, False, 0)
        #HBox1.pack_start(self.Preset4_RadioButton, False, False, 0)
        HBox1.pack_start(Integrate_Button, False, False, 0)
        #HBox1.pack_start(Save_Button, False, False, 0)
        #HBox1.pack_start(Reload_Button, False, False, 0)

        VBox1 = Gtk.VBox(homogeneous = False, spacing = 3)
        VBox1.pack_start(self.Image_ScrolledWindow_EventBox, False, False, 0)
        VBox1.pack_start(HBox1, False, False, 0)
        
        self.win = Gtk.Window()
        self.win.connect("destroy", self.MainWindow_Destroy, dude)
        self.win.add(VBox1)
        self.win.set_title("PowderHelper")
        self.win.show_all()


"""
 else:
                if self.Preset1_RadioButton.get_active():
                    self.Preset_indices[0] = self.indices==index+1
                    self.Preset_bincount[0] = self.bincount[index+1]
                    try:
                        self.Preset1_Vline.remove()
                    except:
                        pass
                    self.Preset1_Vline = self.Image_Axe.axvline(self.Q[index], color="#FF0000")
                elif self.Preset2_RadioButton.get_active():
                    self.Preset_indices[1] = self.indices==index+1
                    self.Preset_bincount[1] = self.bincount[index+1]
                    try:
                        self.Preset2_Vline.remove()
                    except:
                        pass
                    self.Preset2_Vline = self.Image_Axe.axvline(self.Q[index], color="#80FF00")
                elif self.Preset3_RadioButton.get_active():
                    self.Preset_indices[2] = self.indices==index+1
                    self.Preset_bincount[2] = self.bincount[index+1]
                    try:
                        self.Preset3_Vline.remove()
                    except:
                        pass
                    self.Preset3_Vline = self.Image_Axe.axvline(self.Q[index], color="#00FFFF")
                elif self.Preset4_RadioButton.get_active():
                    self.Preset_indices[3] = self.indices==index+1
                    self.Preset_bincount[3] = self.bincount[index+1]
                    try:
                        self.Preset4_Vline.remove()
                    except:
                        pass
                    self.Preset4_Vline = self.Image_Axe.axvline(self.Q[index], color="#8000FF")
                data = np.zeros((dim0, dim1, 4))
                hsv = np.ones((dim0, dim1, 3))
                for i in range(4):
                    if self.parent.sparse_enabled:
                        data[:,:,i] = self.parent.image[:,self.Preset_indices[i].flatten()].sum(1).reshape(dim0, dim1)/self.Preset_bincount[i]
                    else:
                        data[:,:,i] = self.parent.image[:,self.Preset_indices[i]].sum(1).reshape(dim0, dim1)//self.Preset_bincount[i]
                hsv[:,:,1][data.std(2)<self.Threshold_Spin_Adjustment.get_value()] = 0
                datamax = data.max(2)
                hsv[:,:,0] = data.argmax(2) * .25
                #hsv[:,:,1] = np.log(datamax-datamax.min()+1)
                #hsv[:,:,1] = hsv[:,:,1]/hsv[:,:,1].max()
                rgb = colors.hsv_to_rgb(hsv)
                self.parent.Plot2D_Image.set_array(rgb)
                if self.parent.Plot2D_Log_ToggleButton.get_active():
                    self.parent.Plot2D_Image.set_norm(colors.LogNorm())
                else:
                    self.parent.Plot2D_Image.set_norm(colors.Normalize())
                self.Image_Canvas.draw()
"""
