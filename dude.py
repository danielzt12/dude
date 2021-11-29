#!/usr/local/bin python3
"""dude.py: Just Another Data Viewer"""

import os
import sys
import time
from readMDA import *
from math import *
import numpy as np
from multiprocessing import Pool, cpu_count
import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk
from matplotlib import pyplot, colors, patches, font_manager
from matplotlib.figure import Figure
from matplotlib.backends.backend_gtk3agg import FigureCanvasGTK3Agg as FigureCanvas
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import fabio
import h5py
import hdf5plugin
from scipy import sparse
from misc_dude import *
import Shortcuts
import PowderHelper
import ShowMetadata
import ShiftCorrection
import Analyze5D
try:
    import PtychoLib
except:
    pass
else:
    import ProbeLibrary
    

class MyMainWindow:

    #------------------------Add-ons-------------------------------#

    def Addon_PtychoLib(self, widget):
        
        try:
            self.PtychoLib.win.present()
        except:
            self.PtychoLib = PtychoLib.MainWindow(self)

    def Addon_Shortcuts(self, widget):
        
        try:
            self.Shortcuts.win.present()
        except:
            self.Shortcuts = Shortcuts.MainWindow(self)

    def Addon_ProbeLibrary(self, widget):
        
        try:
            self.ProbeLibrary.win.present()
        except:
            self.ProbeLibrary = ProbeLibrary.MainWindow(self)

    def Addon_PowderHelper(self, widget):
        
        try:
            self.PowderHelper.win.present()
        except:
            self.PowderHelper = PowderHelper.MainWindow(self)

    def Addon_ShiftCorrection(self, widget):
        
        try:
            self.ShiftCorrection.win.present()
        except:
            self.ShiftCorrection = ShiftCorrection.MainWindow(self)

    def Addon_Analyze5D(self, widget):
        
        try:
            self.Analyze5D.win.present()
        except:
            self.Analyze5D = Analyze5D.MainWindow(self)

    def Addon_ShowMetadata(self, widget):
        
        try:
            self.ShowMetadata.win.present()
        except:
            self.ShowMetadata = ShowMetadata.MainWindow(self)

    #-------------------------miscellaneous------------------------#

    def Console_Command_Sent(self, widget):

        if self.Console_Entry.get_text() == 'clear':
            widget.get_completion().get_model().clear()
        elif self.Console_Entry.get_text() == 'history':
            for command in widget.get_completion().get_model():
                print(command[0])
        else:
            try:
                exec(self.Console_Entry.get_text())
            except Exception as e:
                print(e)
            else:
                widget.get_completion().get_model().append([self.Console_Entry.get_text()])
                self.console_cmd_index = len(widget.get_completion().get_model())
                self.Console_Entry.set_text('')


    def Image_Canvas_Draw_Rectangle(self, xmin, xmax, ymin, ymax, color, animated = False):

        Rectangle =  patches.Rectangle((xmin, ymin), width = xmax-xmin, height = ymax - ymin,\
        alpha = 1, edgecolor = color, fill = False, linewidth= 3, animated = animated)
        self.Image_Axe.add_patch(Rectangle)
        return Rectangle

    #------------detector image related callbacks--------------#

    def LoadAllImages(self, widget):
        
        t0 = time.time()
        if self.eiger_enabled:
            dd = self.h5["/entry/instrument/detector/data"]
            nd0 = dd.maxshape[0]-dd.shape[0] # maxshape is the intended shape, before interruption of scans
            if self.pump_probe:
                nd0 = int(nd0/2)
            if self.sparse_enabled:
                self.image = None
                nd1 = dd.shape[0]
                nd3 = 0
                nblock = 1000
                for i in range(int(nd1/nblock)+(nd1%nblock>0)):
                    nd2 = nblock if nd1>nblock else nd1
                    nd1 -= nblock
                    if not self.pump_probe:
                        ddd = dd[nd3:nd3+nd2][()].reshape(nd2, 1062*1028)
                        ddd[ddd>2**31]=0
                    else:
                        ddd = np.zeros((int(nd2/2), 1062, 1028))
                        tmp = dd[nd3:nd3+nd2,-512:,:][()].reshape(int(nd2/2), 2, 512, 1028)
                        ddd[:,:512,:] = tmp[:,0]
                        ddd[:,-512:,:] = tmp[:,1]
                        ddd = ddd.reshape(int(nd2/2), 1062*1028)
                        ddd[ddd>65530]=0
                    if i == 0:
                        self.image = sparse.csr_matrix(ddd)
                    else:
                        self.image = sparse.vstack((self.image, sparse.csr_matrix(ddd)))
                    nd3 += nblock
                if nd0:
                    self.image = sparse.vstack((self.image, sparse.csr_matrix(np.zeros((nd0,1062*1028)))))
                #self.image = None
                #tmpdata = self.h5["/entry/instrument/detector/data"][()]
                #tmpdata[tmpdata>2**13]=0
                #if self.Scan_ToolBox_ImageData_MaskAbove_ToggleButton.get_active():
                #    tmpdata *= self.custom_mask
                #self.image = sparse.csr_matrix(tmpdata.reshape(self.h5["/entry/instrument/detector/data"].shape[0], 1062*1028))
                #tmpdata = None
            else:
                self.image = dd[()]
                self.image[self.image>2**31] = 0
                if nd0:
                    self.image = np.vstack((self.image, np.zeros((nd0, 1062, 1028))))
                if self.Scan_ToolBox_ImageData_MaskAbove_ToggleButton.get_active():
                    self.image *= self.custom_mask
        else:
            image_list = list(map(lambda x:os.path.join(self.image_path,x), self.image_list))
            self.image = np.zeros((self.ntiff,self.dimY, self.dimX), dtype="int32") # edist 20211031 from int16
            pool = Pool(processes=cpu_count()) 
            for i, tmparr in enumerate(pool.imap(image_loader, image_list)):
                self.image[i] = tmparr
            pool.close()
            self.image[self.image<0] = 0
            if self.Scan_ToolBox_ImageData_MaskAbove_ToggleButton.get_active():
                self.image *= self.custom_mask
        print (time.time()-t0)
        self.AllImagesLoaded = True
        self.Scan_ToolBox_U_ToggleButton.set_sensitive(True)
        self.Scan_ToolBox_CustomROI_Entry.set_sensitive(True)
        self.Scan_ToolBox_CustomROI_Add_Button.set_sensitive(True)
        self.Scan_ToolBox_CustomROI_Sum_Button.set_sensitive(True)


    def Image_Canvas_Button_Released(self, event):

        self.Image_xend, self.Image_yend = list(map(lambda x: round(x, 0), (Display2Data(self.Image_Axe, event.x, event.y))))
       
        self.Image_Canvas.mpl_disconnect(self.Image_Canvas_Button_Release_Event)
        self.Image_Canvas_Button_Press_Event = self.Image_Canvas.mpl_connect('button_press_event', self.Image_Canvas_Button_Pressed)

        if self.Image_Rectangle_Drawing:
            self.Image_Rectangle_Drawing = False
            self.non_animated_background = None
            try:
                self.Image_ShowROI_Rectangle.remove()
            except (AttributeError, ValueError):
                pass
            xmin = int(min(self.Image_xend, self.Image_xstart))
            xmax = int(max(self.Image_xend, self.Image_xstart))
            ymin = int(min(self.Image_yend, self.Image_ystart))
            ymax = int(max(self.Image_yend, self.Image_ystart))
            
            self.Image_ShowROI_Rectangle = self.Image_Canvas_Draw_Rectangle(xmin, xmax, ymin, ymax, 'b')
            self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.handler_block(self.Scan_ToolBox_CustomROI_YMin_SpinButton_Changed)
            self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.handler_block(self.Scan_ToolBox_CustomROI_XMin_SpinButton_Changed)
            self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.handler_block(self.Scan_ToolBox_CustomROI_XMax_SpinButton_Changed)
            self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.set_upper(self.dimX-1)
            self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.set_value(xmin)
            self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.set_lower(0)
            self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.set_value(xmax)
            self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.set_upper(self.dimY-1)
            self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.set_value(ymin)
            self.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.set_lower(0)
            self.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.set_value(ymax)
            self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.handler_unblock(self.Scan_ToolBox_CustomROI_YMin_SpinButton_Changed)
            self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.handler_unblock(self.Scan_ToolBox_CustomROI_XMin_SpinButton_Changed)
            self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.handler_unblock(self.Scan_ToolBox_CustomROI_XMax_SpinButton_Changed) 
        elif self.Image_Zooming:
            self.Image_Zooming = False
            self.Image_Zoomed = True
            self.Image_ZoomOut_Button.set_sensitive(True)
            self.non_animated_background = None
            try:
                self.Image_ShowROI_Rectangle.remove()
            except (AttributeError, ValueError):
                pass
            xmin = int(min(self.Image_xend, self.Image_xstart))
            xmax = int(max(self.Image_xend, self.Image_xstart))
            ymin = int(min(self.Image_yend, self.Image_ystart))
            ymax = int(max(self.Image_yend, self.Image_ystart))
            dim_y, dim_x = self.Image_Image.get_array().shape
            dx = xmax - xmin
            dy = ymax - ymin
            if 1.*dim_y/dim_x > dy/dx:
                self.Image_Axe.set_xlim(xmin-0.5, xmax-0.5)
                self.Image_Axe.set_ylim((ymin+ymax)/2.+.5*dx/dim_x*dim_y-0.5, (ymin+ymax)/2.-.5*dx/dim_x*dim_y-0.5)
            else:
                self.Image_Axe.set_ylim(ymax-0.5, ymin-0.5)
                self.Image_Axe.set_xlim((xmin+xmax)/2.-.5*dy/dim_y*dim_x-0.5, (xmin+xmax)/2.+.5*dy/dim_y*dim_x-0.5)
            self.Image_Canvas.draw()
        elif self.Image_Panning:
            self.Image_Panning = False
            self.Image_ScrolledWindow_EventBox.get_window().set_cursor(Gdk.Cursor(Gdk.CursorType.TCROSS))
       

    def Image_Canvas_Button_Scrolled(self, event):

        xdata, ydata = list(map(lambda x: round(x, 0), (Display2Data(self.Image_Axe, event.x, event.y))))
        old_x_range = self.Image_Axe.get_xlim()[1]-self.Image_Axe.get_xlim()[0]
        old_y_range = self.Image_Axe.get_ylim()[1]-self.Image_Axe.get_ylim()[0]
        range_x_new = old_x_range*(1-0.1*event.step)
        range_y_new = old_y_range*(1-0.1*event.step)

        self.Image_Axe.set_xlim(xdata-range_x_new*(xdata-self.Image_Axe.get_xlim()[0])/old_x_range, xdata+range_x_new*(self.Image_Axe.get_xlim()[1]-xdata)/old_x_range)
        self.Image_Axe.set_ylim(ydata-range_y_new*(ydata-self.Image_Axe.get_ylim()[0])/old_y_range, ydata+range_y_new*(self.Image_Axe.get_ylim()[1]-ydata)/old_y_range)
        self.Image_Canvas.draw()
        self.Image_ZoomOut_Button.set_sensitive(True)
        self.Image_Zoomed = True
                

    def Image_Canvas_Figure_Entered(self, event):

        #Otherwise button_press_event.key produces none
        self.Image_Canvas.grab_focus()


    def Image_Canvas_Mouse_Hover(self, event):

        self.Image_xend, self.Image_yend = list(map(lambda x: int(round(x, 0)), (Display2Data(self.Image_Axe, event.x, event.y))))
        try:
            dim_y, dim_x = self.Image_Image.get_array().shape
        except:
            return

        if (self.Image_yend < dim_y-0.5) and (self.Image_yend > -.5) and (self.Image_xend > -0.5) and (self.Image_xend < dim_x-0.5):
            if self.show_angle:
                self.Image_Value_Label.set_text('[2The={0:.2f}]  [Gam={1:.2f}]  [Z={2}]'.format(self.twotheta[self.Image_yend, self.Image_xend], self.gamma[self.Image_yend, self.Image_xend], int(self.Image_Image.get_array().data[self.Image_yend, self.Image_xend])))
                # since i don't know which version of matplotlib, somehow get_array() returns a masked array hence the new .data
            else:
                self.Image_Value_Label.set_text('[X = {0:d}]  [Y = {1:d}]  [Z = {2}]'.format(self.Image_xend, self.Image_yend, int(self.Image_Image.get_array().data[self.Image_yend, self.Image_xend])))

            if self.Image_Rectangle_Drawing or self.Image_Zooming:
                if self.non_animated_background != None:
                    # restore the clean slate background
                    self.Image_Canvas.restore_region(self.non_animated_background)
                    if self.Image_xstart > self.Image_xend: #modify the starting point
                        self.Image_Rectangle.set_x(self.Image_xend)
                    self.Image_Rectangle.set_width(abs(self.Image_xend-self.Image_xstart))
                    if self.Image_ystart > self.Image_yend: #modify the starting point
                        self.Image_Rectangle.set_y(self.Image_yend)
                    self.Image_Rectangle.set_height(abs(self.Image_yend-self.Image_ystart))
                    self.Image_Axe.draw_artist(self.Image_Rectangle)
                    self.Image_Canvas.blit(self.Image_Axe.bbox)
                else:
                    xmax = int(max(self.Image_xend, self.Image_xstart))
                    xmin = int(min(self.Image_xend, self.Image_xstart))
                    ymax = int(max(self.Image_yend, self.Image_ystart))
                    ymin = int(min(self.Image_yend, self.Image_ystart))
                    self.Image_Rectangle = self.Image_Canvas_Draw_Rectangle(xmin, xmax, ymin, ymax, 'b', animated = True)
                    self.Image_Canvas.draw()
                    self.non_animated_background = self.Image_Canvas.copy_from_bbox(self.Image_Axe.bbox)
            elif self.Image_Panning:
                self.Image_Axe.set_xlim(self.Image_Axe.get_xlim()[0]+self.Image_xstart-self.Image_xend, \
                                      self.Image_Axe.get_xlim()[1]+self.Image_xstart-self.Image_xend)
                self.Image_Axe.set_ylim(self.Image_Axe.get_ylim()[0]+self.Image_ystart-self.Image_yend, \
                                      self.Image_Axe.get_ylim()[1]+self.Image_ystart-self.Image_yend)
                self.Image_xstart, self.Image_ystart = list(map(lambda x: int(round(x, 0)), (Display2Data(self.Image_Axe, event.x, event.y))))
                self.Image_Canvas.draw()
            
    def Image_ZoomOut(self, widget):

        dim_y, dim_x = self.Image_Image.get_array().shape
        self.Image_Axe.set_xlim(-0.5, dim_x-0.5)
        self.Image_Axe.set_ylim(dim_y-0.5, -0.5)
        self.Image_Canvas.draw()
        self.Image_Zoomed = False
        self.Image_ZoomOut_Button.set_sensitive(False)

    def Image_Canvas_Button_Pressed(self, event):
        
        self.Image_xstart, self.Image_ystart = list(map(lambda x: round(x, 0), (Display2Data(self.Image_Axe, event.x, event.y))))

        if event.button == 2:
            if self.Image_Zoomed:
                dim_y, dim_x = self.Image_Image.get_array().shape
                self.Image_Axe.set_xlim(-0.5, dim_x-0.5)
                self.Image_Axe.set_ylim(dim_y-0.5, -0.5)
                self.Image_Canvas.draw()
                self.Image_Zoomed = False
        else:
            if event.button == 1:
                if self.Image_Zoomed == False:
                    self.Image_Zooming = True
                else:
                    self.Image_Panning = True
                    self.Image_ScrolledWindow_EventBox.get_window().set_cursor(Gdk.Cursor(Gdk.CursorType.FLEUR))
            elif event.button == 3:
                self.Image_Rectangle_Drawing = True
                try:
                    self.Image_ShowROI_Rectangle.remove()
                except (AttributeError, ValueError):
                    pass
            self.Image_Canvas.mpl_disconnect(self.Image_Canvas_Button_Press_Event)
            self.Image_Canvas_Button_Release_Event = self.Image_Canvas.mpl_connect('button_release_event', self.Image_Canvas_Button_Released)


    def Image_AutoScale_toggled(self, widget):

        if widget.get_active():
            data = self.Image_Image.get_array()
            vmin = data.min()
            vmax = data.max()
            self.Image_Vmin_HScale.set_sensitive(False)
            self.Image_Vmax_HScale.set_sensitive(False)
            self.Image_Vmin_HScale_Adjustment.handler_block(self.Image_Vmin_Changed_Handler)
            self.Image_Vmax_HScale_Adjustment.handler_block(self.Image_Vmax_Changed_Handler)
            self.Image_Vmin_HScale_Adjustment.set_upper(vmax-1)
            self.Image_Vmin_HScale_Adjustment.set_value(vmin)
            self.Image_Vmax_HScale_Adjustment.set_lower(vmin+1)
            self.Image_Vmax_HScale_Adjustment.set_value(vmax)
            self.Image_Vmin_HScale_Adjustment.handler_unblock(self.Image_Vmin_Changed_Handler)
            self.Image_Vmax_HScale_Adjustment.handler_unblock(self.Image_Vmax_Changed_Handler)
            if self.Image_Log_ToggleButton.get_active():
                self.Image_Image.set_norm(colors.LogNorm(vmin = vmin+.1, vmax = vmax+.1))
            else:
                self.Image_Image.set_norm(colors.Normalize(vmin = vmin, vmax = vmax))
            self.Image_Canvas.draw()
        else:
            self.Image_Vmin_HScale.set_sensitive(True)
            self.Image_Vmax_HScale.set_sensitive(True)


    def Image_Vscale_Changed(self, widget):
  
        self.Image_Vmax_HScale_Adjustment.set_lower(self.Image_Vmin_HScale_Adjustment.get_value()+1)
        self.Image_Vmin_HScale_Adjustment.set_upper(self.Image_Vmax_HScale_Adjustment.get_value()-1)
        if self.Image_Log_ToggleButton.get_active():
            self.Image_Image.norm = colors.LogNorm(vmin = self.Image_Vmin_HScale_Adjustment.get_value()+.1, vmax = self.Image_Vmax_HScale_Adjustment.get_value()+.1)
        else:
            self.Image_Image.set_norm(colors.Normalize(vmin = self.Image_Vmin_HScale_Adjustment.get_value(), vmax = self.Image_Vmax_HScale_Adjustment.get_value()))
        self.Image_Canvas.draw()


    def Image_Plot_HScale_Changed(self, widget):

        index = int((widget.get_value()-self.image_index_min)/self.nbin)
        if self.eiger_enabled:
            if self.image_index_max > 0:
                if self.pump_probe:
                    image_data = np.zeros((1062,1028), dtype=np.uint16)
                    image_data[:512] = self.h5["/entry/instrument/detector/data"][int(index*2),-512:]
                    image_data[-512:] = self.h5["/entry/instrument/detector/data"][int(index*2+1), -512:]
                else:
                    image_data = self.h5["/entry/instrument/detector/data"][index].astype(np.int32)
            else:
                return # h5 not available
        else:
            image_data = fabio.open(os.path.join(self.image_path,self.image_list[index])).data
        image_data[image_data<0]=0
        if self.Scan_ToolBox_ImageData_MaskAbove_ToggleButton.get_active():
            image_data *= self.custom_mask
        try:
            self.Image_Image.set_array(image_data)
        except:
            self.Image_Image = self.Image_Axe.imshow(image_data, aspect = 'equal', interpolation = 'nearest', cmap = self.cm)
        else:
            pass

        vmin = image_data.min()
        vmax = image_data.max()
        self.Image_Vmin_HScale_Adjustment.set_lower(vmin)
        self.Image_Vmax_HScale_Adjustment.set_upper(vmax)
        if self.Image_AutoScale_ToggleButton.get_active():
            self.Image_Vmin_HScale_Adjustment.handler_block(self.Image_Vmin_Changed_Handler)
            self.Image_Vmax_HScale_Adjustment.handler_block(self.Image_Vmax_Changed_Handler)
            self.Image_Vmin_HScale_Adjustment.set_upper(vmax-1)
            self.Image_Vmin_HScale_Adjustment.set_value(vmin)
            self.Image_Vmax_HScale_Adjustment.set_lower(vmin+1)
            self.Image_Vmax_HScale_Adjustment.set_value(vmax)
            self.Image_Vmin_HScale_Adjustment.handler_unblock(self.Image_Vmin_Changed_Handler)
            self.Image_Vmax_HScale_Adjustment.handler_unblock(self.Image_Vmax_Changed_Handler)
        if self.Image_Log_ToggleButton.get_active():
            self.Image_Image.set_norm(colors.LogNorm(vmin = self.Image_Vmin_HScale_Adjustment.get_value()+.1, vmax = self.Image_Vmax_HScale_Adjustment.get_value()+.1))
        else:
            self.Image_Image.set_norm(colors.Normalize(vmin = self.Image_Vmin_HScale_Adjustment.get_value(), vmax = self.Image_Vmax_HScale_Adjustment.get_value()))
        self.Image_Canvas.draw()

    #----------------1D Plot related callbacks----------------#

    def Plot1D_Canvas_Mouse_Hover(self, event):

        xpointer, ypointer = Display2Data(self.Plot1D_Axe, event.x, event.y)
        try:
            index = np.fabs(self.Plot1D_xdata - xpointer).argmin()
        except:
            return
        if index != self.Image_Plot_HScale_Adjustment.get_value()+self.image_index_min:
            if event.button == 3: 
                try:
                    for hl in self.Plot1D_Hightlight:
                        hl.remove()
                except (AttributeError, ValueError):
                    pass
                x, y = self.Plot1D_xdata[index], self.Plot1D_ydata[index]
                self.Plot1D_Hightlight = self.Plot1D_Axe.plot(x, y, "ro", markersize=10)
                self.Plot1D_Canvas.draw()
                self.Plot1D_Value_Label.set_text('[X = {0:.3f}] [Y = {1:.3f}]'.format(x, y))
                index = self.image_index_min + index * self.nbin # added *self.nbin
                self.Image_Plot_HScale_Adjustment.set_value(index)


    def Plot1D_Vscale_Changed(self, widget):

        self.Plot1D_Axe.set_xscale('log' if self.Plot1D_LogX_ToggleButton.get_active() else 'linear')
        self.Plot1D_Axe.set_yscale('log' if self.Plot1D_LogY_ToggleButton.get_active() else 'linear')
        self.Plot1D_Canvas.draw()


    def Scan_Plot1D(self, xdata, ydata, hold=False):

        if not self.Scan_ToolBox_Y_ToggleButton.get_active() and not hold:
            self.Plot1D_Axe.cla()
        self.Plot1D_Axe.plot(xdata, ydata, "-o")
        self.Plot1D_Axe.set_xscale('log' if self.Plot1D_LogX_ToggleButton.get_active() else 'linear')
        self.Plot1D_Axe.set_yscale('log' if self.Plot1D_LogY_ToggleButton.get_active() else 'linear')
        self.Plot1D_xdata = xdata
        self.Plot1D_ydata = ydata

        
    #----------------2D Plot related callbacks----------------#
    def Plot2D_Canvas_Mouse_Pressed(self, event):

        global _xstart, _ystart
        _xstart, _ystart = map(lambda x: int(round(x, 0)), (Display2Data(self.Plot2D_Axe, event.x, event.y)))
        self.Plot2D_Canvas.mpl_disconnect(self.Plot2D_Canvas_Mouse_Pressed)
        self.Plot2D_Canvas_Button_Release_Event = self.Plot2D_Canvas.mpl_connect('button_release_event', self.Plot2D_Canvas_Mouse_Released)

    def Plot2D_Canvas_Mouse_Released(self, event):

        _xend, _yend = map(lambda x: int(round(x, 0)), (Display2Data(self.Plot2D_Axe, event.x, event.y)))
        self.ttlist += [[(_xstart+event.xdata)/2, (_ystart+event.ydata)/2]]
        self.Plot2D_Canvas.mpl_disconnect(self.Plot2D_Canvas_Mouse_Released)
        self.Plot2D_Canvas_Button_Press_Event = self.Plot2D_Canvas.mpl_connect('button_press_event', self.Plot2D_Canvas_Mouse_Pressed)
        self.Plot2D_Axe.plot(self.ttlist[-1][0], self.ttlist[-1][1], 'r.')
        self.Plot2D_Canvas.draw()
    
    def Plot2D_Canvas_Mouse_Hover(self, event):

        x1, x2 = map(lambda x: int(round(x, 0)), (Display2Data(self.Plot2D_Axe, event.x, event.y)))
        try:
            dim_y, dim_x = self.Plot2D_Image.get_array().shape
        except:
            return
        if (x2 < dim_y-0.5) and (x2 > -.5) and (x1 > -0.5) and (x1 < dim_x-0.5):
            index = int(self.image_index_min+(x1+x2*self.Plot2D_ydata.shape[1])*self.nbin)
            if index != self.Image_Plot_HScale_Adjustment.get_value():
                if event.button == 3:
                    self.Plot2D_P0_Label.set_text('P0 [X = {0:.3f}] [Y = {1:.3f}] [Z = {2:.3f}]'.format(self.Plot2D_xdata2[x2,x1], self.Plot2D_xdata1[x2,x1], self.Plot2D_ydata[x2,x1]))
                    self.Image_Plot_HScale_Adjustment.set_value(index)
                elif event.button == 1:
                    lelements = self.Plot2D_P0_Label.get_text().split()
                    x0, y0 = float(lelements[3][:-1]), float(lelements[6][:-1])
                    self.Plot2D_P1_Label.set_text('P1 [X = {0:.3f}] [Y = {1:.3f}] [dX = {2:.3f}] [dY = {3:.3f}]'.format(self.Plot2D_xdata2[x2,x1], self.Plot2D_xdata1[x2,x1], self.Plot2D_xdata2[x2,x1]-x0, self.Plot2D_xdata1[x2,x1]-y0))

        #elif event.xdata != None and event.ydata != None:
        #    self.Plot2D_Value_Label.set_text('[X = {0:.3f}]  [Y = {1:.3f}]'.format(event.xdata, event.ydata))
    

    def Plot2D_AutoScale_toggled(self, widget):

        if widget.get_active():
            data = self.Plot2D_Image.get_array()
            vmin = data.min()
            vmax = data.max()
            self.Plot2D_Vmin_HScale.set_sensitive(False)
            self.Plot2D_Vmax_HScale.set_sensitive(False)
            self.Plot2D_Vmin_HScale_Adjustment.handler_block(self.Plot2D_Vmin_Changed_Handler)
            self.Plot2D_Vmax_HScale_Adjustment.handler_block(self.Plot2D_Vmax_Changed_Handler)
            self.Plot2D_Vmin_HScale_Adjustment.set_upper(vmax-1)
            self.Plot2D_Vmin_HScale_Adjustment.set_value(vmin)
            self.Plot2D_Vmax_HScale_Adjustment.set_lower(vmin+1)
            self.Plot2D_Vmax_HScale_Adjustment.set_value(vmax)
            self.Plot2D_Vmin_HScale_Adjustment.handler_unblock(self.Plot2D_Vmin_Changed_Handler)
            self.Plot2D_Vmax_HScale_Adjustment.handler_unblock(self.Plot2D_Vmax_Changed_Handler)
            if self.Plot2D_Log_ToggleButton.get_active():
                self.Plot2D_Image.set_norm(colors.LogNorm(vmin = vmin+.1, vmax = vmax+.1))
            else:
                self.Plot2D_Image.set_norm(colors.Normalize(vmin = vmin, vmax = vmax))
            self.Plot2D_Canvas.draw()
        else:
            self.Plot2D_Vmin_HScale.set_sensitive(True)
            self.Plot2D_Vmax_HScale.set_sensitive(True)

    def Plot2D_ScaleBar_toggled(self, widget):
        
        try:
            self.Plot2D_scalebar.remove()
        except:
            pass
        if widget.get_active() and self.data[1].curr_pt > 2:
            um_per_pixel = np.abs(self.data[1].p[0].data[1]-self.data[1].p[0].data[0])
            um_choice = np.array([0.1,0.2,0.5,1,2,5,10,20,30,40,50,100,200,500])
            x_um = um_choice[np.argmin(np.abs(um_choice / (um_per_pixel * self.data[2].npts) - 1/4))]
            x_pixel = int(x_um / um_per_pixel)
            self.Plot2D_scalebar = AnchoredSizeBar(self.Plot2D_Axe.transData,\
                                                   x_pixel, u'{0} μm'.format(x_um), 'lower right', \
                                                   pad=1, label_top = True,\
                                                   color='white',\
                                                   frameon=False,\
                                                   size_vertical=self.Plot2D_xdata2.shape[0]/80., \
                                                   fontproperties = font_manager.FontProperties(size=20))
            self.Plot2D_Axe.add_artist(self.Plot2D_scalebar)
        else:
            pass
        self.Plot2D_Canvas.draw()


    def Plot2D_Vscale_Changed(self, widget):
  
        self.Plot2D_Vmax_HScale_Adjustment.set_lower(self.Plot2D_Vmin_HScale_Adjustment.get_value()+1)
        self.Plot2D_Vmin_HScale_Adjustment.set_upper(self.Plot2D_Vmax_HScale_Adjustment.get_value()-1)
        if self.Plot2D_Log_ToggleButton.get_active():
            self.Plot2D_Image.set_norm(colors.LogNorm(vmin = self.Plot2D_Vmin_HScale_Adjustment.get_value()+.1, vmax = self.Plot2D_Vmax_HScale_Adjustment.get_value()+.1))
        else:
            self.Plot2D_Image.set_norm(colors.Normalize(vmin = self.Plot2D_Vmin_HScale_Adjustment.get_value(), vmax = self.Plot2D_Vmax_HScale_Adjustment.get_value()))
        self.Plot2D_Canvas.draw()


    def Scan_PlotSpiral(self, xdata1, xdata2, ydata):

        if not self.Scan_ToolBox_Y_ToggleButton.get_active():
            self.Plot2D_Axe.cla()
        self.Plot2D_Axe.set_axis_off()
        self.Plot2D_Image = self.Plot2D_Axe.scatter(xdata2, xdata1, ydata, s=2, edgecolor='none', cmap=self.cm)
        if self.Plot2D_Log_ToggleButton.get_active():
            self.Plot2D_Image.set_norm(colors.LogNorm(vmin = self.Plot2D_Vmin_HScale_Adjustment.get_value()+.1, vmax = self.Plot2D_Vmax_HScale_Adjustment.get_value()+.1))
        else:
            self.Plot2D_Image.set_norm(colors.Normalize(vmin = self.Plot2D_Vmin_HScale_Adjustment.get_value(), vmax = self.Plot2D_Vmax_HScale_Adjustment.get_value()))
        self.Plot2D_xdata1 = xdata1
        self.Plot2D_xdata2 = xdata2
        self.Plot2D_ydata = ydata


    def Scan_Plot2D(self, xdata1, xdata2, ydata):

        if not self.Scan_ToolBox_Y_ToggleButton.get_active():
            self.Plot2D_Axe.cla()
        self.Plot2D_Axe.set_axis_off()
        self.Plot2D_Image = self.Plot2D_Axe.imshow(ydata,  aspect = 'equal', interpolation = 'nearest', cmap=self.cm, origin="lower")
        if self.Plot2D_Log_ToggleButton.get_active():
            self.Plot2D_Image.set_norm(colors.LogNorm(vmin = max(0.1,self.Plot2D_Vmin_HScale_Adjustment.get_value()), vmax = self.Plot2D_Vmax_HScale_Adjustment.get_value()))
        else:
            self.Plot2D_Image.set_norm(colors.Normalize(vmin = self.Plot2D_Vmin_HScale_Adjustment.get_value(), vmax = self.Plot2D_Vmax_HScale_Adjustment.get_value()))
        self.Plot2D_xdata1 = xdata1
        self.Plot2D_xdata2 = xdata2
        self.Plot2D_ydata = ydata
        if self.Plot2D_ScaleBar_ToggleButton.get_active():
            self.Plot2D_ScaleBar_toggled(self.Plot2D_ScaleBar_ToggleButton)

    #-------------------main menu callbacks--------------------#

    def FileDialog_Construction(self, widget, flag = 0):

        if flag == 0: #MDA
            FileDialog = Gtk.FileChooserDialog(title = "Choose MDA Folder", action = Gtk.FileChooserAction.SELECT_FOLDER)
            FileDialog.set_transient_for(self.Main_Window)
            FileDialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_OPEN, Gtk.ResponseType.OK)
            response = FileDialog.run()
            if response == Gtk.ResponseType.OK:
                self.MDA_folder = FileDialog.get_filename()
                if "mda" in self.MDA_folder:
                    self.Folder_Open()
                    self.Image_folder = os.path.join(os.path.abspath(os.path.join(self.MDA_folder, os.pardir)),"Images")
                    self.h5_folder = os.path.join(os.path.abspath(os.path.join(self.MDA_folder, os.pardir)),"h5")                
                FileDialog.destroy()
            else:
                FileDialog.destroy()
        elif flag == 1: # previously for selecting Images folder is not found, not useless
            FileDialog = Gtk.FileChooserDialog(title = "Choose Image Folder", action = Gtk.FileChooserAction.SELECT_FOLDER)
            FileDialog.set_transient_for(self.Main_Window)
            FileDialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_OPEN, Gtk.ResponseType.OK)
            response = FileDialog.run()
            if response == Gtk.ResponseType.OK:
                self.Image_folder = FileDialog.get_filename()
            FileDialog.destroy()
        elif flag == 2:
            FileDialog = Gtk.FileChooserDialog(title = "Choose a file", action = Gtk.FileChooserAction.OPEN)
            FileDialog.set_transient_for(self.Main_Window)
            FileDialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_OPEN, Gtk.ResponseType.OK)
            response = FileDialog.run()
            if response == Gtk.ResponseType.OK:
                filename = FileDialog.get_filename()
            else:
                filename = None
            FileDialog.destroy()
            return filename
        elif flag == 3:
            FileDialog = Gtk.FileChooserDialog(title = "Choose a file", action = Gtk.FileChooserAction.SAVE)
            FileDialog.set_transient_for(self.Main_Window)
            FileDialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_SAVE, Gtk.ResponseType.OK)
            response = FileDialog.run()
            if response == Gtk.ResponseType.OK:
                filename = FileDialog.get_filename()
            else:
                filename = None
            FileDialog.destroy()
            return filename


    def Detector_Changed(self, widget):

        if not widget.get_active():
            return
        text = widget.get_label()
        self.pilatus_enabled = False
        self.eiger_enabled = False
        if text == "1028 x 1062":
            self.dimY = 1062
            self.dimX = 1028
            self.eiger_enabled = True
        elif text == "516 x 516":
            self.dimY = 516
            self.dimX = 516
        elif text == "515 x 515":
            self.dimY = 515
            self.dimX = 515
        elif text == "256 x 256":
            self.dimY = 256
            self.dimX = 256
        elif text == "487 x 195":
            self.dimY = 195
            self.dimX = 487
            self.pilatus_enabled = True
        elif text == "1360 x 1024":
            self.dimY = 1024
            self.dimX = 1360
        self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.set_upper(self.dimX-1)
        self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.set_upper(self.dimY-1)
        self.Image_Axe.cla()
        self.Image_Axe.set_axis_off()
        del self.Image_Image
        #self.Image_Axe.xaxis.set_ticklabels([])
        #self.Image_Axe.yaxis.set_ticklabels([])


    def Change_Colormap(self, widget):

        if widget.get_active(): #the callback is invoked two times, deselecting one and selecing another
            self.cm = widget.get_label()
            self.Image_Image.set_cmap(self.cm)
            self.Image_Canvas.draw()
            self.Plot2D_Image.set_cmap(self.cm)
            self.Plot2D_Canvas.draw()

    def SparseToggled(self, widget):

        if widget.get_active():
            self.sparse_enabled = True
        else:
            self.sparse_enabled = False

    def ShowAngleToggled(self, widget):

        if widget.get_active():
            analysis_folder = os.path.join(os.path.abspath(os.path.join(self.MDA_folder, os.pardir, 'Analysis')))    
            if os.path.exists(os.path.join(analysis_folder, "gamma.csv")):
                self.gamma = np.genfromtxt(os.path.join(analysis_folder, "gamma.csv"), delimiter=',')
                self.twotheta = np.genfromtxt(os.path.join(analysis_folder, "twotheta.csv"), delimiter=',')
            else:
                self.gamma = np.zeros((1062,1028))
                self.twotheta = np.copy(self.gamma)
            self.show_angle = True
        else:
            self.show_angle = False

    def DirtyFixToggled(self, widget):
        
        if widget.get_active():
            self.dirty_fix = True
        else:
            self.dirty_fix = False

    def PumpProbeToggled(self, widget):
        
        if widget.get_active():
            self.pump_probe = True
        else:
            self.pump_probe = False


    def MainWindow_Destroy(self, widget): 

        self.image = None
        try:
            self.h5.close()
        except:
            pass
        try:
            self.Shortcuts.win.destroy()
        except:
            pass
        try:
            self.ProbeLibrary.win.destroy()
        except:
            pass
        try:
            self.PowderHelper.win.destroy()
        except:
            pass
        try:
            self.ShiftCorrection.win.destroy()
        except:
            pass
        try:
            self.Analyze5D.win.destroy()
        except:
            pass
        try:
            self.PtychoLib.win.destroy()
        except:
            pass
        try:
            self.ShowMetadata.win.destroy()
        except:
            pass
        self.Main_Window.destroy()
        Gtk.main_quit()


    def AboutThisProgram(self, widget):

        dialog = Gtk.AboutDialog()
        dialog.set_program_name("dude")
        dialog.set_comments("Diffraction User Data Explorer")
        dialog.set_version("2.3.5")
        dialog.set_website("https://hub.docker.com/r/danielzt12/dude")
        dialog.set_website_label("https://hub.docker.com/r/danielzt12/dude")
        dialog.get_child().get_children()[1].hide()
        dialog.get_child().get_children()[0].get_children()[0].hide()
        dialog.run()
        

    #---------------scan toolbox callbacks--------------------#   

    def Scan_ToolBox_ImageData_MaskAbove(self, widget):

        if widget.get_active():
            threshold = int(self.Scan_ToolBox_ImageData_HotPixel_Adjustment.get_value())
            self.custom_mask = self.Image_Image.get_array()<threshold
        else:
            self.custom_mask = None


    def Scan_ToolBox_ImageData_MaskBelow(self, widget):

        threshold = int(self.Scan_ToolBox_ImageData_HotPixel_Adjustment.get_value())
        if self.sparse_enabled:
            self.image.data[self.image.data<threshold]=0
            self.image.eliminate_zeros()
        else:
            self.image[self.image<threshold] = 0


    def Scan_ToolBox_U_Toggled(self, widget):

        if not widget.get_active():
            try:
                self.Image_ShowROI_Rectangle.remove()
            except (AttributeError, ValueError):
                pass
            self.Image_Canvas.draw()
        self.Scan_ToolBox_Plot_Changed(None)


    def Scan_ToolBox_CustomROI_Summed(self, widget):
        
        ndim = len(self.data)-1 if self.data[0]["dimensions"][-1] != 2048 else len(self.data)-2
        xmin = int(self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.get_value())
        xmax = int(self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.get_value())
        ymin = int(self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.get_value())
        ymax = int(self.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.get_value())
        if ndim == 1:
            xdata =  np.array(self.data[1].p[0].data)
            if self.sparse_enabled:
                indice0 = np.zeros((1062,1028))
                indice0[ymin:ymax+1,xmin:xmax+1]=1
                indice0 = indice0.flatten().nonzero()[0]
                ydata = np.array(self.image[:,indice0].sum(1))[:,0]
            else:
                ydata = self.image[:, ymin:ymax+1,xmin:xmax+1].sum(1).sum(1)
            self.Scan_Plot1D(xdata[(xdata!=0)*(ydata!=0)], ydata[(xdata!=0)*(ydata!=0)])
            if self.pump_probe:
                if ymin>550:
                    ymin -= 550
                    ymax -= 550
                else:
                    ymin += 550
                    ymax += 550
                indice0 = np.zeros((1062,1028))
                indice0[ymin:ymax+1,xmin:xmax+1]=1
                indice0 = indice0.flatten().nonzero()[0]
                ydata2 = np.array(self.image[:,indice0].sum(1))[:,0]
                self.Scan_Plot1D(xdata[(xdata!=0)*(ydata!=0)], ydata2[(xdata!=0)*(ydata!=0)], hold=True) 
            self.Plot1D_Canvas.draw()
        elif ndim == 2:
            xdata2 = np.zeros((self.data[0]["dimensions"][0],self.data[0]["dimensions"][1]))
            xdata1 = np.copy(xdata2)+1
            ydata = np.copy(xdata2)
            datatmp = np.array(self.data[2].p[0].data)
            xdata2[:datatmp.shape[0]] = datatmp
            xdata1 *= np.array(self.data[1].p[0].data)[:,np.newaxis]
            if self.sparse_enabled:
                indice0 = np.zeros((1062,1028))
                indice0[ymin:ymax+1,xmin:xmax+1]=1
                indice0 = indice0.flatten().nonzero()[0]
                ydata = self.image[:,indice0].sum(1).reshape(xdata1.shape)
            else:
                ydata = self.image[:, ymin:ymax+1,xmin:xmax+1].sum(1).sum(1).reshape(xdata1.shape)
            self.Plot2D_Axe.cla()
            self.Plot2D_Axe.set_axis_off()
            if self.pump_probe:
                if ymin>550:
                    ymin -= 550
                    ymax -= 550
                else:
                    ymin += 550
                    ymax += 550
                indice0 = np.zeros((1062,1028))
                indice0[ymin:ymax+1,xmin:xmax+1]=1
                indice0 = indice0.flatten().nonzero()[0]
                # for some reasons the images aren't right, and needs to be shifted
                ydata2 = self.image[:,indice0].sum(1).reshape(xdata1.shape)
                ydata_flat = ydata.flatten()
                ydata_flat[0,:-1] = ydata_flat[0,1:]
                ydata = ydata_flat.reshape(ydata.shape)
                ydata_flat = ydata2.flatten()
                ydata_flat[0,:-1] = ydata_flat[0,1:]
                ydata2 = ydata_flat.reshape(ydata2.shape)
                ydata = np.hstack((ydata,ydata2))
                self.Plot2D_xdata1 = np.hstack((self.Plot2D_xdata1,self.Plot2D_xdata1))
                self.Plot2D_xdata2 = np.hstack((self.Plot2D_xdata2,self.Plot2D_xdata2))
            self.Plot2D_Image = self.Plot2D_Axe.imshow(ydata,  aspect = 'equal', interpolation = 'nearest', cmap=self.cm, origin="lower")
            self.Plot2D_ydata = ydata
            if self.Plot2D_Log_ToggleButton.get_active():
                self.Plot2D_Image.set_norm(colors.LogNorm())
            else:
                self.Plot2D_Image.set_norm(colors.Normalize())
            self.Plot2D_Canvas.draw()
        
    def Scan_ToolBox_CustomROI_Added(self, widget):
        
        roi_name = self.Scan_ToolBox_CustomROI_Entry.get_text()
        roi_info = [roi_name, \
                        int(self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.get_value()),\
                        int(self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.get_value()),\
                        int(self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.get_value()),\
                        int(self.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.get_value())]
        for row in self.CustomROI_store:
            if roi_info[0] == row[0]:
                roi_info[0]+="*"
        self.CustomROI_store.append(roi_info)
        self.Scan_ToolBox_U_ComboBox.set_active(len(self.Scan_ToolBox_U_ComboBox.get_model())-1)
        self.Scan_ToolBox_U_ToggleButton.set_active(True)


    def Scan_ToolBox_CustomROI_Changed(self, widget):
        
        xmin = int(self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.get_value())
        xmax = int(self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.get_value())
        ymin = int(self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.get_value())
        ymax = int(self.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.get_value())
        self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.set_upper(xmax-1)
        self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.set_lower(xmin+1)
        self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.set_upper(ymax-1)
        self.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.set_lower(ymin+1)
        
        try:
            self.Image_ShowROI_Rectangle.remove()
        except (AttributeError, ValueError):
            pass
        self.Image_ShowROI_Rectangle = self.Image_Canvas_Draw_Rectangle(xmin, xmax, ymin, ymax, 'b')
        self.Image_Canvas.draw()
        
     
    def Scan_ToolBox_Plot_Changed(self, widget):

        ndim = len(self.data)-1 if self.data[0]["dimensions"][-1] != 2048 else len(self.data)-2
        if ndim == 1:
            xdata =  np.array(self.data[1].p[0].data)

            if self.Scan_ToolBox_U_ToggleButton.get_active():
                roi_info = self.CustomROI_store[self.Scan_ToolBox_U_ComboBox.get_active()]
                if self.sparse_enabled:   # disabling sparse for 1d
                    indice0 = np.zeros((1062,1028))
                    indice0[roi_info[3]:roi_info[4]+1,roi_info[1]:roi_info[2]+1]=1
                    indice0 = indice0.flatten().nonzero()[0]
                    ydata = np.array(self.image[:,indice0].sum(1))[:,0]
                    #print(ydata.shape)
                else:
                    ydata = self.image[:, roi_info[3]:roi_info[4]+1,roi_info[1]:roi_info[2]+1].sum(1).sum(1)
                self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.handler_block(self.Scan_ToolBox_CustomROI_YMin_SpinButton_Changed)
                self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.handler_block(self.Scan_ToolBox_CustomROI_XMin_SpinButton_Changed)
                self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.handler_block(self.Scan_ToolBox_CustomROI_XMax_SpinButton_Changed)
                self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.set_upper(self.dimX-1)
                self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.set_value(roi_info[1])
                self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.set_lower(0)
                self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.set_value(roi_info[2])
                self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.set_upper(self.dimY-1)
                self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.set_value(roi_info[3])
                self.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.set_lower(0)
                self.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.set_value(roi_info[4])
                self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.handler_unblock(self.Scan_ToolBox_CustomROI_YMin_SpinButton_Changed)
                self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.handler_unblock(self.Scan_ToolBox_CustomROI_XMin_SpinButton_Changed)
                self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.handler_unblock(self.Scan_ToolBox_CustomROI_XMax_SpinButton_Changed) 
            else:
                ydata =  np.array(self.data[1].d[self.Scan_ToolBox_Y_ComboBox.get_active()].data)
                # dirty fix first point bug
                if self.dirty_fix and \
                   ("eiger" in self.data[1].d[self.Scan_ToolBox_Y_ComboBox.get_active()].name or \
                   "QMPX3" in self.data[1].d[self.Scan_ToolBox_Y_ComboBox.get_active()].name ):
                    ydata[:-1] = ydata[1:]
                    # print("dirty fixing first point bug")
            if self.Scan_ToolBox_M_ToggleButton.get_active():
                mdata = np.array(self.data[1].d[self.Scan_ToolBox_M_ComboBox.get_active()].data)
                mdata /= mdata.mean()
                ydata /= mdata
            self.Scan_Plot1D(xdata[(xdata!=0)*(ydata!=0)], ydata[(xdata!=0)*(ydata!=0)]) 
            self.Plot1D_Canvas.draw()
        else:
            # this is to avoid crashing on unfinished data
            xdata2 = np.zeros((self.data[0]["dimensions"][0],self.data[0]["dimensions"][1]))
            xdata1 = np.copy(xdata2)+1
            ydata = np.copy(xdata2)
            mdata = np.copy(xdata2)
            datatmp = np.array(self.data[2].p[0].data)
            xdata2[:datatmp.shape[0]] = datatmp
            xdata1 *= np.array(self.data[1].p[0].data)[:,np.newaxis]

            if self.Scan_ToolBox_U_ToggleButton.get_active():
                roi_info = self.CustomROI_store[self.Scan_ToolBox_U_ComboBox.get_active()]
                if self.sparse_enabled:
                    indice0 = np.zeros((1062,1028))
                    indice0[roi_info[3]:roi_info[4]+1,roi_info[1]:roi_info[2]+1]=1
                    indice0 = indice0.flatten().nonzero()[0]
                    ydata = self.image[:,indice0].sum(1).reshape(xdata1.shape)
                else:
                    ydata = self.image[:, roi_info[3]:roi_info[4]+1,roi_info[1]:roi_info[2]+1].sum(1).sum(1).reshape(xdata1.shape)
                self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.handler_block(self.Scan_ToolBox_CustomROI_YMin_SpinButton_Changed)
                self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.handler_block(self.Scan_ToolBox_CustomROI_XMin_SpinButton_Changed)
                self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.handler_block(self.Scan_ToolBox_CustomROI_XMax_SpinButton_Changed)
                self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.set_upper(self.dimX-1)
                self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.set_value(roi_info[1])
                self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.set_lower(0)
                self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.set_value(roi_info[2])
                self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.set_upper(self.dimY-1)
                self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.set_value(roi_info[3])
                self.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.set_lower(0)
                self.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.set_value(roi_info[4])
                self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.handler_unblock(self.Scan_ToolBox_CustomROI_YMin_SpinButton_Changed)
                self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.handler_unblock(self.Scan_ToolBox_CustomROI_XMin_SpinButton_Changed)
                self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.handler_unblock(self.Scan_ToolBox_CustomROI_XMax_SpinButton_Changed) 
            else:
                datatmp =  np.array(self.data[2].d[self.Scan_ToolBox_Y_ComboBox.get_active()].data)
                ydata[:datatmp.shape[0]] = datatmp
                if self.dirty_fix and \
                   ("eiger" in self.data[2].d[self.Scan_ToolBox_Y_ComboBox.get_active()].name or \
                    "QMPX3" in self.data[2].d[self.Scan_ToolBox_Y_ComboBox.get_active()].name ) and \
                   self.data[1].time != "whatever":
                    print("dirty fixing first point bug")
                    ydata_flat = ydata.flatten()
                    ydata_flat[:-1] = ydata_flat[1:]
                    ydata = ydata_flat.reshape(ydata.shape)
            if self.Scan_ToolBox_M_ToggleButton.get_active():
                datatmp = np.array(self.data[2].d[self.Scan_ToolBox_M_ComboBox.get_active()].data)
                mdata[:datatmp.shape[0]] = datatmp
                mdata /= mdata[mdata>0].mean()
                mdata[mdata==0] = 1
                ydata /= mdata

            ydata_min = ydata[ydata!=0].min()
            self.Plot2D_Vmin_HScale_Adjustment.set_lower(ydata_min)
            self.Plot2D_Vmax_HScale_Adjustment.set_upper(ydata.max())
            if self.Plot2D_AutoScale_ToggleButton.get_active():
                self.Plot2D_Vmin_HScale_Adjustment.handler_block(self.Plot2D_Vmin_Changed_Handler)
                self.Plot2D_Vmax_HScale_Adjustment.handler_block(self.Plot2D_Vmax_Changed_Handler)
                self.Plot2D_Vmin_HScale_Adjustment.set_upper(ydata.max()-.01)
                self.Plot2D_Vmin_HScale_Adjustment.set_value(ydata_min)
                self.Plot2D_Vmax_HScale_Adjustment.set_lower(ydata_min+.01)
                self.Plot2D_Vmax_HScale_Adjustment.set_value(ydata.max())
                self.Plot2D_Vmin_HScale_Adjustment.handler_unblock(self.Plot2D_Vmin_Changed_Handler)
                self.Plot2D_Vmax_HScale_Adjustment.handler_unblock(self.Plot2D_Vmax_Changed_Handler)
            if self.MDA_File_ListStore[self.mda_selection_path[0]][1] != "ao01":
                self.Scan_Plot2D(xdata1, xdata2, ydata)
            else:
                self.Scan_PlotSpiral(xdata1, xdata2, ydata)
            self.Plot2D_Canvas.draw()


    def Scan_TreeView_Selection_Changed(self, treeselection):

        treemodelfilter, filterpathlist = treeselection.get_selected_rows()
        self.mda_selection_path = list(map(treemodelfilter.convert_path_to_child_path, filterpathlist))
        mdapath = self.MDA_File_ListStore[self.mda_selection_path[0]][-2]
        self.Scan_ToolBox_U_ToggleButton.set_active(False)
        self.Scan_ToolBox_U_ToggleButton.set_sensitive(False)
        self.AllImagesLoaded = False
        self.Scan_Load(mdapath)


    def Scan_TreeView_Select_Hijack(self, treeselection, treefilter, irow, selected):
        
        irow, icolumn = self.MDA_File_TreeView.get_cursor()
        if icolumn == self.MDA_File_TreeView.get_column(10): 
            self.MDA_File_ListStore[irow][-1] = not self.MDA_File_ListStore[irow][-1]
            return False
        else:
            return True


    def Scan_Load(self, mdapath):

        Scan_ToolBox_Y_ComboBox_Select = self.Scan_ToolBox_Y_ComboBox.get_active()
        Scan_ToolBox_M_ComboBox_Select = self.Scan_ToolBox_M_ComboBox.get_active()
        self.Scan_ToolBox_Y_ComboBox.handler_block(self.Scan_ToolBox_Y_Changed_Handler)
        self.Scan_ToolBox_M_ComboBox.handler_block(self.Scan_ToolBox_M_Changed_Handler)
        self.MDA_Det_store.clear()
        self.data = readMDA(mdapath, verbose=0)
        ndim = len(self.data)-1 if self.data[0]["dimensions"][-1] != 2048 else len(self.data)-2
        for i in range(self.data[ndim].nd):
            dname = self.data[ndim].d[i].name
            self.MDA_Det_store.append([dname.replace("s26_eiger_cnm","eiger"), 1])
        self.Scan_ToolBox_M_ComboBox.set_active(Scan_ToolBox_M_ComboBox_Select)
        self.Scan_ToolBox_Y_ComboBox.handler_unblock(self.Scan_ToolBox_Y_Changed_Handler)
        self.Scan_ToolBox_M_ComboBox.handler_unblock(self.Scan_ToolBox_M_Changed_Handler)
        self.Plot_Notebook.set_current_page(ndim-1)        
        self.Scan_ToolBox_Y_ComboBox.set_active(Scan_ToolBox_Y_ComboBox_Select)

        if self.eiger_enabled:
            h5_filename = [f.name for f in os.scandir(self.h5_folder) if "scan_{0}".format(self.MDA_File_ListStore[self.mda_selection_path[0]][0]) in f.name]
            if len(h5_filename) == 1:
                h5_filename = h5_filename[0]
            else:
                h5_index = np.array([int(f.split(".")[0].split("_")[-1]) for f in h5_filename])
                h5_filename = np.take(h5_filename, h5_index.argsort())[0]
            try:
                self.h5.close()
            except Exception as e:
                pass
            try:
                self.h5 = h5py.File(os.path.join(self.h5_folder, h5_filename), 'r')
            except Exception as e:
                self.image_index_max = -1
            else:
                self.image_index_max = self.h5['entry/data/data'].shape[0]-1
                self.image_index_min = 0
                self.nbin = 1
                if hasattr(self, "ShowMetadata"):
                    self.ShowMetadata.h5 = self.h5
                    self.ShowMetadata.update_keystore()
                    print("updating metadata")
                
        else:
            self.image_path = os.path.join(self.Image_folder, mdapath.split(".")[0].split("_")[-1].lstrip("0"))
            # os.chmod(self.image_path, 0o777)
            if self.pilatus_enabled == True:
                self.image_list = [imagefile.name for imagefile in os.scandir(self.image_path) if ('.tif' in imagefile.name and 'il' in imagefile.name)]
            else:
                self.image_list = [imagefile.name for imagefile in os.scandir(self.image_path) if ('.tif' in imagefile.name and 'il' not in imagefile.name)]
            self.tif_index = np.array([int(filename.split(".")[0].split("_")[-1]) for filename in self.image_list])
            self.image_list = np.take(self.image_list , self.tif_index.argsort())
            self.nbin = 1
            if ndim == 2:
                self.ntiff = self.data[0]['dimensions'][0] * self.data[0]['dimensions'][1]
                for i in range(self.data[ndim].nd):
                    dname = self.data[ndim].d[i].name
                    if "FileNumber" in dname:
                        mda_index = np.array(self.data[2].d[i].data).flatten()
                #print("a,",mda_index)
                print (mda_index.max(), self.tif_index.max(), mda_index.min(), self.tif_index.min())
            else:
                self.ntiff = self.data[0]['dimensions'][0]
            if len(self.image_list) != self.ntiff:
                print ("I only print this because the file saving here has been messy!")
                print ("Expecting {0} files but got {1} instead".format(self.ntiff, len(self.image_list)))
                if ndim == 1:
                    if self.ntiff > len(self.image_list):
                        print ("I suppose the scan is interrupted and will not attempt to correct that")
                    else:
                        if len(self.image_list)%self.ntiff == 0:
                            self.nbin = max(1,int(len(self.image_list)/self.ntiff))
                            print ("I suppose the scan has {0} repeated measurement for each point".format(self.nbin))
                            self.image_list = self.image_list[0::self.nbin]
                        else:
                            print ("There is a chance that files are being saved into the wrong folder. Please correct manually 1!")
                else:
                    if len(self.image_list) > self.ntiff:
                        self.nbin = max(1,int(len(self.image_list)/self.ntiff))
                        if len(self.image_list)%self.ntiff:
                            self.nbin+=1
                        if self.nbin>1:
                            if not len(self.image_list)%self.ntiff:
                                print ("I suppose the scan has {0} repeated measurement for each point".format(self.nbin))
                                self.image_list = self.image_list[0::self.nbin]
                            else:
                                if mda_index.min() == self.tif_index.min():
                                    inc_index = mda_index[1:] - mda_index[:-1]
                                    hiccup_index = (inc_index == 2).nonzero()[0]
                                    if inc_index.max() == 2:
                                        print ("I suppose the scan has {0} repeated measurement for each point, but there are {1} images missing here!".format(self.nbin, self.ntiff*self.nbin-len(self.image_list)))
                                        if inc_index[hiccup_index+1].sum() == 0 : # because a hiccup is always followed by a 0                                 
                                            print("Experiencing {0} hiccups. Use only tif index.".format((inc_index==2).sum()))
                                            mda_index[hiccup_index+1] -= 1
                                            self.image_list = self.image_list[(mda_index-mda_index.min()).astype(int)]
                                            self.image_list = self.image_list[0::self.nbin]
                                        else:
                                            print("Experiencing {0} hiccups, unable to correct.".format((inc_index==2).sum()))
                                    elif inc_index.min() < 0 or inc_index.max() > 2:
                                        print ("This should not happen. Please correct manually 5!")
                                    elif inc_index.min() == 0 and (inc_index==0).sum() == mda_index.shape[0] - self.tif_index.shape[0] + hiccup_index.shape[0]:
                                        print("Lost {0} images\nDeciding to trust mda index.".format((inc_index==0).sum()- hiccup_index.shape[0]))
                                        self.image_list = self.image_list[0::self.nbin]
                                    else:
                                        print("Dunno what to do 7!")
                                else:
                                    print("Dunno what to do 6!")
                        else:
                            if mda_index.max() < self.tif_index.max():
                                if mda_index.min() == self.tif_index.min():
                                    print ("Ignoring {0} extra images in the folder.\nThey might belong to the next scan.".format(self.tif_index.max() - mda_index.max()))
                                    self.image_list = self.image_list[self.tif_index<=mda_index.max()]
                                    self.tif_index = self.tif_index[self.tif_index<=mda_index.max()]
                                else:
                                    print ("There is a chance that files are being saved into the wrong folder. Please correct manually 2!")
                            else:
                                if mda_index.max() == self.tif_index.max() and mda_index.min() == self.tif_index.min() +1:
                                    print ("This bug is new in 2019, the mda index is shifted by 1 with regard to the tif index, but strangely for these long scans there is one more point in the beginning")
                                    self.image_list = self.image_list[:-1]
                                else:
                                    print ("This should not happen. Please correct manually 3!")
                    else:
                        if mda_index.min() == self.tif_index.min():
                            inc_index = mda_index[1:] - mda_index[:-1]
                            hiccup_index = (inc_index == 2).nonzero()[0]
                            if inc_index.min() == 0 and (inc_index==0).sum() == mda_index.shape[0] - self.tif_index.shape[0] + hiccup_index.shape[0]:
                                print("Lost {0} images\nDeciding to trust mda index.".format((inc_index==0).sum()- hiccup_index.shape[0]))
                                self.image_list = self.image_list[(mda_index-mda_index.min()).astype(int)]
                            elif inc_index.max() == 2:
                                if inc_index[hiccup_index+1].sum() == 0 : # because a hiccup is always followed by a 0
                                    print("Experiencing {0} hiccups. Use only tif index.".format((inc_index==2).sum()))
                                    mda_index[hiccup_index+1] -= 1
                                    self.image_list = self.image_list[mda_index]
                                else:
                                    print("Experiencing {0} hiccups, unable to correct.".format((inc_index==2).sum()))
                            elif inc_index.min() < 0 or inc_index.max() > 2:
                                print ("This should not happen. Please correct manually 4!")
                            else:
                                print("Dunno what to do 8!")
                        else:
                            if mda_index.min() == self.tif_index.min() +1:
                                print ("This bug is new in 2019, the mda index is shifted by 1 with regard to the tif index, and I suppose you interrupted the scan. Other than that all is well")
                            else:
                                print ("There is a chance that files are being saved into the wrong folder. Please correct manually!")
            if self.tif_index.shape[0]:
                self.image_index_min = self.tif_index.min()
                self.image_index_max = self.tif_index.max()
                ###self.image_index_min = mda_index.min()
                ###self.image_index_max = mda_index.max()

        self.Image_Plot_HScale_Adjustment.set_upper(self.image_index_max)
        self.Image_Plot_HScale_Adjustment.set_lower(self.image_index_min)
        if self.Image_Plot_HScale_Adjustment.get_value() == 0:
            self.Image_Plot_HScale_Changed(self.Image_Plot_HScale_Adjustment)
        else:
            self.Image_Plot_HScale_Adjustment.set_value(self.image_index_min)
        self.Image_Plot_HScale_Adjustment.set_step_increment(self.nbin)
        self.Scan_ToolBox_CustomROI_Entry.set_sensitive(False)
        self.Scan_ToolBox_CustomROI_Add_Button.set_sensitive(False)
        self.Scan_ToolBox_CustomROI_Sum_Button.set_sensitive(False)
        

    def Folder_Refresh(self, widget):

        self.MDA_cursor_current = self.MDA_File_TreeView.get_cursor()[0] #unstable
        self.MDA_File_TreeView.get_selection().handler_block(self.MDA_File_TreeView_Selection_Changed_Handler) #unstable
        MDAfile_list = [mdafile for mdafile in os.listdir(self.MDA_folder) if mdafile.endswith('.mda')]
        index_list = np.array([int(filename.split(".")[0].split("_")[-1]) for filename in MDAfile_list])
        last_index = int(self.MDA_File_ListStore[-1][0])
        MDAfile_list = np.array(MDAfile_list)[index_list>=last_index]
        index_list = index_list[index_list>=last_index]
        MDAfile_list = np.take(MDAfile_list , index_list.argsort())
        self.Folder_Scan(MDAfile_list) 


    def Folder_Open(self):

        self.MDA_File_TreeView.get_selection().handler_block(self.MDA_File_TreeView_Selection_Changed_Handler)
        self.MDA_File_ListStore.clear()
        self.MDA_File_TreeView.get_selection().handler_unblock(self.MDA_File_TreeView_Selection_Changed_Handler)
        MDAfile_list = [mdafile for mdafile in os.listdir(self.MDA_folder) if mdafile.endswith('mda')]
        index_list = np.array([int(filename.split(".")[0].split("_")[-1]) for filename in MDAfile_list])
        MDAfile_list = np.take(MDAfile_list , index_list.argsort())
        self.Folder_Scan(MDAfile_list)


    def Folder_Scan(self, MDAfile_list):

        if len(self.MDA_File_ListStore):
            self.MDA_File_ListStore.remove(self.MDA_File_ListStore[-1].iter)
        mda_list = list(map(lambda x:os.path.join(self.MDA_folder,x), MDAfile_list))
        pool = Pool(processes=cpu_count()) 
        for i, tmplist in enumerate(pool.imap(mda_loader, mda_list)):
            self.MDA_File_ListStore.append(tmplist)
        pool.close()
        self.MDA_File_TreeView.get_selection().handler_unblock(self.MDA_File_TreeView_Selection_Changed_Handler)#unstable
        if hasattr(self, "MDA_cursor_current"):
            self.MDA_File_TreeView.set_cursor(self.MDA_cursor_current)#unstable
        
               
    def Console_KeyPressed(self, widget, event):

        keyname = Gdk.keyval_name(event.keyval)
        if keyname == 'Up' and (event.state & Gdk.ModifierType.SHIFT_MASK):
            if self.console_cmd_index > 0:
                self.console_cmd_index -= 1
            else:
                self.console_cmd_index = len(widget.get_completion().get_model())-1
            widget.set_text(widget.get_completion().get_model()[self.console_cmd_index][0])
        elif keyname == 'Down' and (event.state & Gdk.ModifierType.SHIFT_MASK):
            if self.console_cmd_index < len(widget.get_completion().get_model()) -1:
                self.console_cmd_index += 1
            else:
                self.console_cmd_index = 0
            widget.set_text(widget.get_completion().get_model()[self.console_cmd_index][0])
 

    def __init__(self):

        self.cm = "viridis"
        self.dimY = 1062
        self.dimX = 1028
        self.Image_Rectangle_Drawing = False
        self.Image_Zoomed = False
        self.Image_Zooming = False
        self.Image_Panning = False
        self.non_animated_background = None
        self.custom_mask = None

        #------------------------MDA File Window------------------------------#    
        MDA_File_ScrolledWindow = Gtk.ScrolledWindow()
        MDA_File_ScrolledWindow.set_policy(Gtk.PolicyType.ALWAYS, Gtk.PolicyType.ALWAYS)
        MDA_File_ScrolledWindow.set_overlay_scrolling(False)
        self.MDA_File_ListStore = Gtk.ListStore(str, str, str, str, str, str, str, str, str, str, str, bool)
        self.MDA_File_TreeView_Filter = self.MDA_File_ListStore.filter_new()
        self.MDA_File_TreeView = Gtk.TreeView.new_with_model(self.MDA_File_TreeView_Filter)
        self.MDA_File_TreeView.get_selection().set_mode(Gtk.SelectionMode.SINGLE)
        self.MDA_File_TreeView_Selection_Changed_Handler = self.MDA_File_TreeView.get_selection().connect("changed", self.Scan_TreeView_Selection_Changed)
        # hijack the default function to allow multiselect without loading all the scans
        self.MDA_File_TreeView.get_selection().set_select_function(self.Scan_TreeView_Select_Hijack)

        i = 0      
        for name in ["N#","M1","MIN","MAX","NP","M2","MIN","MAX","NP","CT"]:
            MDA_File_CellRendererText = Gtk.CellRendererText()     
            MDA_File_CellRendererText.set_property("xalign", 1)
            MDA_File_TreeViewColumn = Gtk.TreeViewColumn(name, MDA_File_CellRendererText, text = i)
            i += 1
            self.MDA_File_TreeView.append_column(MDA_File_TreeViewColumn)
        MDA_File_CellRendererToggle = Gtk.CellRendererToggle()     
        MDA_File_CellRendererToggle.set_property("xalign", .5)
        MDA_File_TreeViewColumn = Gtk.TreeViewColumn("", MDA_File_CellRendererToggle, active = 11)
        self.MDA_File_TreeView.append_column(MDA_File_TreeViewColumn)
        self.MDA_File_TreeView.set_enable_search(False)
        MDA_File_ScrolledWindow.add(self.MDA_File_TreeView)
        MDA_File_ScrolledWindow.set_size_request(450, 400)

        
        # ROI Plot : X, Y, Mon
        self.MDA_Det_store = Gtk.ListStore(str, int)
        renderer_text = Gtk.CellRendererText()    
        self.CustomROI_store = Gtk.ListStore(str, int, int, int, int)
        renderer_text2 = Gtk.CellRendererText()    

        self.Scan_ToolBox_M_ToggleButton = Gtk.ToggleButton(label = " MON ")
        self.Scan_ToolBox_M_ComboBox = Gtk.ComboBox.new_with_model(self.MDA_Det_store)
        self.Scan_ToolBox_M_ComboBox.pack_start(renderer_text, True)
        self.Scan_ToolBox_M_ComboBox.add_attribute(renderer_text, "text", 0)
        self.Scan_ToolBox_M_ComboBox.add_attribute(renderer_text, "visible", 1)
        self.Scan_ToolBox_M_ToggleButton.connect("toggled", self.Scan_ToolBox_Plot_Changed)
        self.Scan_ToolBox_Y_ToggleButton = Gtk.ToggleButton(label = " DET ")
        self.Scan_ToolBox_Y_ComboBox = Gtk.ComboBox.new_with_model(self.MDA_Det_store)
        self.Scan_ToolBox_Y_ComboBox.pack_start(renderer_text, True)
        self.Scan_ToolBox_Y_ComboBox.add_attribute(renderer_text, "text", 0)
        self.Scan_ToolBox_Y_ComboBox.add_attribute(renderer_text, "visible", 1)
        self.Scan_ToolBox_U_ToggleButton = Gtk.ToggleButton(label = " ROI ")
        self.Scan_ToolBox_U_ComboBox = Gtk.ComboBox.new_with_model(self.CustomROI_store)
        self.Scan_ToolBox_U_ComboBox.pack_start(renderer_text2, True)
        self.Scan_ToolBox_U_ComboBox.add_attribute(renderer_text2, "text", 0)
        self.Scan_ToolBox_U_ToggleButton.connect("toggled", self.Scan_ToolBox_U_Toggled)
        self.Scan_ToolBox_U_ToggleButton.set_active(False)
        self.Scan_ToolBox_U_ToggleButton.set_sensitive(False)

        Scan_Det_HBox1 = Gtk.HBox(homogeneous = False, spacing = 0)
        Scan_Det_HBox1.set_border_width(3)
        Scan_Det_HBox1.pack_start(self.Scan_ToolBox_Y_ToggleButton, False, False, 3)
        Scan_Det_HBox1.pack_start(self.Scan_ToolBox_Y_ComboBox, True, True, 3)
        Scan_Det_HBox2 = Gtk.HBox(homogeneous = False, spacing = 0)
        Scan_Det_HBox2.set_border_width(3)
        Scan_Det_HBox2.pack_start(self.Scan_ToolBox_U_ToggleButton, False, False, 3)
        Scan_Det_HBox2.pack_start(self.Scan_ToolBox_U_ComboBox, True, True, 3)
        Scan_Det_HBox3 = Gtk.HBox(homogeneous = False, spacing = 0)
        Scan_Det_HBox3.set_border_width(3)
        Scan_Det_HBox3.pack_start(self.Scan_ToolBox_M_ToggleButton, False, False, 3)
        Scan_Det_HBox3.pack_start(self.Scan_ToolBox_M_ComboBox, True, True, 3)
        
        self.Scan_ToolBox_Y_Changed_Handler = self.Scan_ToolBox_Y_ComboBox.connect("changed", self.Scan_ToolBox_Plot_Changed)
        self.Scan_ToolBox_M_Changed_Handler = self.Scan_ToolBox_M_ComboBox.connect("changed", self.Scan_ToolBox_Plot_Changed)
        self.Scan_ToolBox_U_Changed_Handler = self.Scan_ToolBox_U_ComboBox.connect("changed", self.Scan_ToolBox_Plot_Changed)
        
         # ROI Management
        self.Scan_ToolBox_CustomROI_Entry = Gtk.Entry()
        self.Scan_ToolBox_CustomROI_Entry.set_size_request(40, 20)
        self.Scan_ToolBox_CustomROI_Sum_Button = Gtk.Button(label = " Sum ")
        self.Scan_ToolBox_CustomROI_Sum_Button.connect("clicked", self.Scan_ToolBox_CustomROI_Summed)
        self.Scan_ToolBox_CustomROI_Add_Button = Gtk.Button(label = " Add ")
        self.Scan_ToolBox_CustomROI_Add_Button.connect("clicked", self.Scan_ToolBox_CustomROI_Added)

        Scan_ToolBox_ImageData_LoadAll_Button = Gtk.Button(label = 'Load All')
        Scan_ToolBox_ImageData_LoadAll_Button.connect("clicked", self.LoadAllImages)

        Scan_ToolBox_CustomROI_HBox = Gtk.HBox(homogeneous = False, spacing = 3)
        Scan_ToolBox_CustomROI_HBox.set_border_width(3)
        Scan_ToolBox_CustomROI_HBox.pack_start(self.Scan_ToolBox_CustomROI_Entry, True, True, 0)
        Scan_ToolBox_CustomROI_HBox.pack_start(self.Scan_ToolBox_CustomROI_Add_Button, False, True, 0)
        Scan_ToolBox_CustomROI_HBox.pack_start(self.Scan_ToolBox_CustomROI_Sum_Button, False, True, 0)
        Scan_ToolBox_CustomROI_HBox.pack_start(Scan_ToolBox_ImageData_LoadAll_Button, False, True, 0)
        self.Scan_ToolBox_CustomROI_Entry.set_sensitive(False)
        self.Scan_ToolBox_CustomROI_Add_Button.set_sensitive(False)
        self.Scan_ToolBox_CustomROI_Sum_Button.set_sensitive(False)
        
        # ROI Definition
        self.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment = Gtk.Adjustment(value = 200, lower = 0, upper = self.dimY-1, step_increment = 1, page_increment = 2, page_size = 0)
        self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment = Gtk.Adjustment(value = 100, lower = 0, upper = self.dimY-1, step_increment = 1, page_increment = 2, page_size = 0)
        self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment = Gtk.Adjustment(value = 200, lower = 0, upper = self.dimX-1, step_increment = 1, page_increment = 2, page_size = 0)
        self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment = Gtk.Adjustment(value = 100, lower = 0, upper = self.dimX-1, step_increment = 1, page_increment = 2, page_size = 0)
        Scan_ToolBox_CustomROI_YMax_SpinButton = Gtk.SpinButton()
        Scan_ToolBox_CustomROI_YMax_SpinButton.set_adjustment(self.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment)
        Scan_ToolBox_CustomROI_YMax_SpinButton.set_digits(0)
        Scan_ToolBox_CustomROI_YMin_SpinButton = Gtk.SpinButton()
        Scan_ToolBox_CustomROI_YMin_SpinButton.set_adjustment(self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment)
        Scan_ToolBox_CustomROI_YMin_SpinButton.set_digits(0)
        Scan_ToolBox_CustomROI_XMax_SpinButton = Gtk.SpinButton()
        Scan_ToolBox_CustomROI_XMax_SpinButton.set_adjustment(self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment)
        Scan_ToolBox_CustomROI_XMax_SpinButton.set_digits(0)
        Scan_ToolBox_CustomROI_XMin_SpinButton = Gtk.SpinButton()
        Scan_ToolBox_CustomROI_XMin_SpinButton.set_adjustment(self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment)
        Scan_ToolBox_CustomROI_XMin_SpinButton.set_digits(0)
        Scan_ToolBox_CustomROI_Cut_Button = Gtk.Button(label = ' Cut ')

        self.Scan_ToolBox_CustomROI_YMax_SpinButton_Changed = self.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.connect('value-changed', self.Scan_ToolBox_CustomROI_Changed)
        self.Scan_ToolBox_CustomROI_YMin_SpinButton_Changed = self.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.connect('value-changed', self.Scan_ToolBox_CustomROI_Changed)
        self.Scan_ToolBox_CustomROI_XMax_SpinButton_Changed = self.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.connect('value-changed', self.Scan_ToolBox_CustomROI_Changed)
        self.Scan_ToolBox_CustomROI_XMin_SpinButton_Changed = self.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.connect('value-changed', self.Scan_ToolBox_CustomROI_Changed)

        Scan_ToolBox_CustomROI_Table = Gtk.Table(n_rows = 3, n_columns = 3, homogeneous = False)
        Scan_ToolBox_CustomROI_Table.attach(Scan_ToolBox_CustomROI_YMax_SpinButton, left_attach = 1, right_attach = 2, top_attach = 2, bottom_attach = 3)
        Scan_ToolBox_CustomROI_Table.attach(Scan_ToolBox_CustomROI_YMin_SpinButton, left_attach = 1, right_attach = 2, top_attach = 0, bottom_attach = 1)  
        Scan_ToolBox_CustomROI_Table.attach(Scan_ToolBox_CustomROI_XMin_SpinButton, left_attach = 0, right_attach = 1, top_attach = 1, bottom_attach = 2)  
        Scan_ToolBox_CustomROI_Table.attach(Scan_ToolBox_CustomROI_XMax_SpinButton, left_attach = 2, right_attach = 3, top_attach = 1, bottom_attach = 2)

        # Customized HotPixel Mask
        self.Scan_ToolBox_ImageData_HotPixel_Adjustment = Gtk.Adjustment(value = 0, lower = 0, upper = 65536, step_increment = 1, page_increment = 10, page_size = 0)
        Scan_ToolBox_ImageData_HotPixel_SpinButton = Gtk.SpinButton()
        Scan_ToolBox_ImageData_HotPixel_SpinButton.set_adjustment(self.Scan_ToolBox_ImageData_HotPixel_Adjustment)
        Scan_ToolBox_ImageData_HotPixel_SpinButton.set_digits(0)
        self.Scan_ToolBox_ImageData_MaskAbove_ToggleButton = Gtk.ToggleButton(label = 'Mask Above')
        self.Scan_ToolBox_ImageData_MaskAbove_ToggleButton.connect('clicked', self.Scan_ToolBox_ImageData_MaskAbove)
        self.Scan_ToolBox_ImageData_MaskAbove_ToggleButton.set_tooltip_text('When this button is pressed down, all pixels on the current image which are hotter than the defined threshold will be masked, the mask will then be applied to all images as long as this button is kept pressed down.\n This is independent of the static hot pixel correction in Menu->Settings') 

        self.Scan_ToolBox_ImageData_MaskBelow_Button = Gtk.Button(label = 'Mask Below')
        self.Scan_ToolBox_ImageData_MaskBelow_Button.connect('clicked', self.Scan_ToolBox_ImageData_MaskBelow)
        self.Scan_ToolBox_ImageData_MaskBelow_Button.set_tooltip_text('When this button is pressed, all loaded data with value below the set value is set to 0') 

        Scan_ToolBox_ImageData_HBox = Gtk.HBox(homogeneous = False, spacing = 3)
        Scan_ToolBox_ImageData_HBox.set_border_width(3)
        Scan_ToolBox_ImageData_HBox.pack_start(Scan_ToolBox_ImageData_HotPixel_SpinButton, True, True, 0)
        Scan_ToolBox_ImageData_HBox.pack_start(self.Scan_ToolBox_ImageData_MaskAbove_ToggleButton, True, False, 0)
        Scan_ToolBox_ImageData_HBox.pack_start(self.Scan_ToolBox_ImageData_MaskBelow_Button, True, False, 0)
        #Scan_ToolBox_ImageData_HBox.pack_start(Scan_ToolBox_ImageData_LoadAll_Button, True, False, 0)

        Spec_ToolBox_VBox = Gtk.VBox(homogeneous = False, spacing = 0)
        Spec_ToolBox_VBox.set_border_width(3)
        Spec_ToolBox_VBox.pack_start(Scan_ToolBox_ImageData_HBox, False, False, 5)
        Spec_ToolBox_VBox.pack_start(Scan_ToolBox_CustomROI_HBox, False, False, 5)
        Spec_ToolBox_VBox.pack_start(Scan_ToolBox_CustomROI_Table, False, False, 5)
        Spec_ToolBox_VBox.pack_start(Scan_Det_HBox1, False, False, 5)
        Spec_ToolBox_VBox.pack_start(Scan_Det_HBox2, False, False, 5)
        Spec_ToolBox_VBox.pack_start(Scan_Det_HBox3, False, False, 5)
        Spec_ToolBox_VBox.set_size_request(350, 400)

        ############################  Scan Plot 1D  ###############################
        
        self.Plot1D_Figure = Figure()
        self.Plot1D_Axe = self.Plot1D_Figure.add_axes([0.08, 0.08, 0.87, 0.87])
        self.Plot1D_Axe.set_autoscale_on(True)
        self.Plot1D_Canvas = FigureCanvas(self.Plot1D_Figure)
        Plot1D_ScrolledWindow = Gtk.ScrolledWindow()
        Plot1D_ScrolledWindow.set_policy(Gtk.PolicyType.AUTOMATIC, Gtk.PolicyType.AUTOMATIC)
        Plot1D_ScrolledWindow.add(self.Plot1D_Canvas)

        self.Plot1D_Canvas_Mouse_Hover_Event = self.Plot1D_Canvas.mpl_connect('motion_notify_event', self.Plot1D_Canvas_Mouse_Hover)
        self.Plot1D_Canvas_Button_Press_Event = self.Plot1D_Canvas.mpl_connect('button_press_event', self.Plot1D_Canvas_Mouse_Hover) # so that a single right click also changes the images

        self.Plot1D_LogX_ToggleButton = Gtk.ToggleButton(label = " LogX ")
        self.Plot1D_LogX_ToggleButton.connect("toggled", self.Plot1D_Vscale_Changed)
        self.Plot1D_LogY_ToggleButton = Gtk.ToggleButton(label = " LogY ")
        self.Plot1D_LogY_ToggleButton.connect("toggled", self.Plot1D_Vscale_Changed)
        Plot1D_Canvas_Toolbar_VSeparator1 = Gtk.VSeparator()
        Plot1D_Canvas_Toolbar_VSeparator2 = Gtk.VSeparator()

        self.Plot1D_Value_Label = Gtk.Label()     

        Plot1D_Canvas_Toolbar_HBox = Gtk.HBox(homogeneous = False, spacing = 3)
        Plot1D_Canvas_Toolbar_HBox.set_border_width(3)
        Plot1D_Canvas_Toolbar_HBox.pack_start(Plot1D_Canvas_Toolbar_VSeparator1, False, False, 3)
        Plot1D_Canvas_Toolbar_HBox.pack_start(self.Plot1D_LogX_ToggleButton, False, False, 0)
        Plot1D_Canvas_Toolbar_HBox.pack_start(self.Plot1D_LogY_ToggleButton, False, False, 0)
        Plot1D_Canvas_Toolbar_HBox.pack_start(Plot1D_Canvas_Toolbar_VSeparator2, False, False, 3)
        Plot1D_Canvas_Toolbar_HBox.pack_end(self.Plot1D_Value_Label, False, False, 0)

        Plot1D_VBox = Gtk.VBox(homogeneous = False, spacing = 3)
        Plot1D_VBox.set_border_width(3)
        Plot1D_VBox.pack_start(Plot1D_Canvas_Toolbar_HBox, False, False, 0)
        Plot1D_VBox.pack_start(Plot1D_ScrolledWindow, True, True, 0)
        
        Plot1D_VBox.set_size_request(800,600)

        #------------------------Plot 2D-------------------------------#

        self.Plot2D_Figure = Figure()
        self.Plot2D_Axe = self.Plot2D_Figure.add_axes([0, 0, 1, 1])
        self.Plot2D_Canvas = FigureCanvas(self.Plot2D_Figure)
        Plot2D_ScrolledWindow = Gtk.ScrolledWindow()
        Plot2D_ScrolledWindow.set_policy(Gtk.PolicyType.AUTOMATIC, Gtk.PolicyType.AUTOMATIC)
        Plot2D_ScrolledWindow.add(self.Plot2D_Canvas)

        Plot2D_ScrolledWindow_EventBox = Gtk.EventBox()
        Plot2D_ScrolledWindow_EventBox.add(Plot2D_ScrolledWindow)
        
        self.Plot2D_Canvas_Mouse_Hover_Event = self.Plot2D_Canvas.mpl_connect('motion_notify_event', self.Plot2D_Canvas_Mouse_Hover)
        self.Plot2D_Canvas_Button_Press_Event = self.Plot2D_Canvas.mpl_connect('button_press_event', self.Plot2D_Canvas_Mouse_Hover) # so that a single right click also changes the images
        #self.Plot2D_Canvas_Button_Press_Event = self.Plot2D_Canvas.mpl_connect('button_press_event', self.Plot2D_Canvas_Mouse_Pressed)
        
        self.Plot2D_Log_ToggleButton = Gtk.ToggleButton(label = "Log")
        self.Plot2D_Log_ToggleButton.set_active(True)
        self.Plot2D_Log_ToggleButton.connect("toggled", self.Plot2D_Vscale_Changed)
        
        self.Plot2D_AutoScale_ToggleButton = Gtk.ToggleButton(label = "Auto")
        self.Plot2D_AutoScale_ToggleButton.set_active(True)
        self.Plot2D_AutoScale_ToggleButton.connect("toggled", self.Plot2D_AutoScale_toggled)

        self.Plot2D_ScaleBar_ToggleButton = Gtk.ToggleButton(label = "Scale")
        self.Plot2D_ScaleBar_ToggleButton.set_active(False)
        self.Plot2D_ScaleBar_ToggleButton.connect("toggled", self.Plot2D_ScaleBar_toggled)

        Plot2D_VSeperator = Gtk.VSeparator()

        self.Plot2D_Vmin_HScale_Adjustment = Gtk.Adjustment(value = 0, lower = 0, upper = 10000, step_increment = .1, page_increment = 1, page_size = 0)
        self.Plot2D_Vmin_Changed_Handler = self.Plot2D_Vmin_HScale_Adjustment.connect("value_changed", self.Plot2D_Vscale_Changed)
        self.Plot2D_Vmin_HScale = Gtk.HScale()
        self.Plot2D_Vmin_HScale.set_adjustment(self.Plot2D_Vmin_HScale_Adjustment)
        self.Plot2D_Vmin_HScale.set_size_request(150, 20)
        self.Plot2D_Vmin_HScale.set_value_pos(Gtk.PositionType.LEFT)
        self.Plot2D_Vmin_HScale.set_digits(1)

        self.Plot2D_Vmax_HScale_Adjustment = Gtk.Adjustment(value = 10000, lower = 10000, upper = 20000, step_increment = .1, page_increment = 1, page_size = 0)
        self.Plot2D_Vmax_Changed_Handler = self.Plot2D_Vmax_HScale_Adjustment.connect("value_changed", self.Plot2D_Vscale_Changed)
        self.Plot2D_Vmax_HScale = Gtk.HScale()
        self.Plot2D_Vmax_HScale.set_adjustment(self.Plot2D_Vmax_HScale_Adjustment)
        self.Plot2D_Vmax_HScale.set_size_request(150, 20)
        self.Plot2D_Vmax_HScale.set_value_pos(Gtk.PositionType.LEFT)
        self.Plot2D_Vmax_HScale.set_digits(1)
        
        self.Plot2D_Vmin_HScale.set_sensitive(False)
        self.Plot2D_Vmax_HScale.set_sensitive(False)

        self.Plot2D_P0_Label = Gtk.Label()
        self.Plot2D_P1_Label = Gtk.Label()
        
        Plot2DToolbar_HBox1 = Gtk.HBox(homogeneous = False, spacing = 3)
        Plot2DToolbar_HBox1.set_border_width(3)
        Plot2DToolbar_HBox1.pack_start(self.Plot2D_Log_ToggleButton, False, False, 0)
        Plot2DToolbar_HBox1.pack_start(self.Plot2D_AutoScale_ToggleButton, False, False, 0)
        Plot2DToolbar_HBox1.pack_start(self.Plot2D_ScaleBar_ToggleButton, False, False, 0)
        Plot2DToolbar_HBox1.pack_start(Plot2D_VSeperator, False, False, 3)
        Plot2DToolbar_HBox1.pack_start(self.Plot2D_Vmin_HScale, True, True, 0)
        Plot2DToolbar_HBox1.pack_start(self.Plot2D_Vmax_HScale, True, True, 0)

        Plot2DToolbar_HBox2 = Gtk.HBox(homogeneous = False, spacing = 3)
        Plot2DToolbar_HBox2.set_border_width(3)
        Plot2DToolbar_HBox2.pack_end(self.Plot2D_P0_Label, False, False, 0)
        Plot2DToolbar_HBox2.pack_start(self.Plot2D_P1_Label, False, False, 0)
 
        Plot2D_VBox = Gtk.VBox(homogeneous = False, spacing = 3)
        Plot2D_VBox.set_border_width(3)
        Plot2D_VBox.pack_start(Plot2DToolbar_HBox1, False, False, 0)
        Plot2D_VBox.pack_start(Plot2DToolbar_HBox2, False, False, 0)
        Plot2D_VBox.pack_start(Plot2D_ScrolledWindow_EventBox, True, True, 0)

        self.PlotSpare_Figure = Figure()
        self.PlotSpare_Axe = self.PlotSpare_Figure.add_axes([0, 0.05, 1, .9])
        self.PlotSpare_Canvas = FigureCanvas(self.PlotSpare_Figure)

        self.Plot_Notebook = Gtk.Notebook()
        self.Plot_Notebook.set_show_tabs(False)
        self.Plot_Notebook.append_page(Plot1D_VBox)
        self.Plot_Notebook.append_page(Plot2D_VBox)
        self.Plot_Notebook.append_page(self.PlotSpare_Canvas)

        self.Image_Figure = Figure()
        self.Image_Axe = self.Image_Figure.add_axes([0, 0, 1, 1])
        self.Image_Axe.set_axis_off()
        #self.Image_Axe.xaxis.set_ticklabels([])
        #self.Image_Axe.yaxis.set_ticklabels([])
        self.Image_Canvas = FigureCanvas(self.Image_Figure)
        
        self.Image_Canvas_Figure_Enter_Event = self.Image_Canvas.mpl_connect('figure_enter_event', self.Image_Canvas_Figure_Entered)
        self.Image_Canvas_Button_Press_Event = self.Image_Canvas.mpl_connect('button_press_event', self.Image_Canvas_Button_Pressed)
        self.Image_Canvas_Button_Scroll_Event = self.Image_Canvas.mpl_connect('scroll_event', self.Image_Canvas_Button_Scrolled)
        self.Image_Canvas_Mouse_Hover_Event = self.Image_Canvas.mpl_connect('motion_notify_event', self.Image_Canvas_Mouse_Hover)

        self.Image_ScrolledWindow = Gtk.ScrolledWindow()
        self.Image_ScrolledWindow.set_policy(Gtk.PolicyType.AUTOMATIC, Gtk.PolicyType.AUTOMATIC)
        self.Image_ScrolledWindow.add(self.Image_Canvas)
        self.Image_ScrolledWindow_EventBox = Gtk.EventBox()
        self.Image_ScrolledWindow_EventBox.add(self.Image_ScrolledWindow)

        self.Image_Value_Label = Gtk.Label()

        self.Image_Plot_HScale_Adjustment = Gtk.Adjustment(value = 0, lower = 0, upper = 100000, step_increment = 1, page_increment = 1, page_size = 0)
        self.Image_Plot_HScale_Adjustment.connect("value_changed", self.Image_Plot_HScale_Changed)
        self.Image_Plot_HScale = Gtk.HScale()
        self.Image_Plot_HScale.set_adjustment(self.Image_Plot_HScale_Adjustment)
        self.Image_Plot_HScale.set_size_request(350, 20)
        self.Image_Plot_HScale.set_value_pos(Gtk.PositionType.LEFT)
        self.Image_Plot_HScale.set_digits(0)

        self.Image_ZoomOut_Button = Gtk.Button(label = "Full")
        self.Image_ZoomOut_Button.connect("clicked", self.Image_ZoomOut)
        self.Image_ZoomOut_Button.set_sensitive(False)
        self.Image_AutoScale_ToggleButton = Gtk.ToggleButton(label = "Auto")
        self.Image_AutoScale_ToggleButton.set_active(True)
        self.Image_AutoScale_ToggleButton.connect("toggled", self.Image_AutoScale_toggled)

        self.Image_Vmin_HScale_Adjustment = Gtk.Adjustment(value = 0, lower = 0, upper = 10000, step_increment = 1, page_increment = 10, page_size = 0)
        self.Image_Vmin_Changed_Handler = self.Image_Vmin_HScale_Adjustment.connect("value_changed", self.Image_Vscale_Changed)
        self.Image_Vmin_HScale = Gtk.Scale(orientation = Gtk.Orientation.HORIZONTAL)
        self.Image_Vmin_HScale.set_adjustment(self.Image_Vmin_HScale_Adjustment)
        self.Image_Vmin_HScale.set_size_request(20, 20)
        self.Image_Vmin_HScale.set_value_pos(Gtk.PositionType.LEFT)
        self.Image_Vmin_HScale.set_digits(0)
        Image_Vmin_Label = Gtk.Label()
        Image_Vmin_Label.set_text('Vmin:')

        self.Image_Vmax_HScale_Adjustment = Gtk.Adjustment(value = 10000, lower = 10000, upper = 20000, step_increment = 1, page_increment = 10, page_size = 0)
        self.Image_Vmax_Changed_Handler = self.Image_Vmax_HScale_Adjustment.connect("value_changed", self.Image_Vscale_Changed)
        self.Image_Vmax_HScale = Gtk.HScale()
        self.Image_Vmax_HScale.set_adjustment(self.Image_Vmax_HScale_Adjustment)
        self.Image_Vmax_HScale.set_size_request(20, 20)
        self.Image_Vmax_HScale.set_value_pos(Gtk.PositionType.LEFT)
        self.Image_Vmax_HScale.set_digits(0)
        Image_Vmax_Label = Gtk.Label()
        Image_Vmax_Label.set_text('Vmax:')

        self.Image_Log_ToggleButton = Gtk.ToggleButton(label = "Log")
        self.Image_Log_ToggleButton.set_active(True)
        self.Image_Log_ToggleButton.connect("toggled", self.Image_Vscale_Changed)

        self.Image_Vmin_HScale.set_sensitive(False)
        self.Image_Vmax_HScale.set_sensitive(False)
        
        Image_Toolbar_HBox1 = Gtk.HBox(homogeneous = False, spacing = 3)
        Image_Toolbar_HBox1.set_border_width(3)
        Image_Toolbar_HBox1.pack_start(Image_Vmin_Label, False, False, 3)
        Image_Toolbar_HBox1.pack_start(self.Image_Vmin_HScale, True, True, 3)
        Image_Toolbar_HBox1.pack_start(Image_Vmax_Label, False, False, 3)
        Image_Toolbar_HBox1.pack_start(self.Image_Vmax_HScale, True, True, 3)

        Image_Toolbar_HBox2 = Gtk.HBox(homogeneous = False, spacing = 3)
        Image_Toolbar_HBox2.set_border_width(3)
        Image_Toolbar_HBox2.pack_start(self.Image_ZoomOut_Button, False, False, 3)
        Image_Toolbar_HBox2.pack_start(self.Image_Log_ToggleButton, False, False, 3)
        Image_Toolbar_HBox2.pack_start(self.Image_AutoScale_ToggleButton, False, False, 3)
        Image_Toolbar_HBox2.pack_start(self.Image_Plot_HScale, False, False, 0)
        Image_Toolbar_HBox2.pack_end(self.Image_Value_Label, False, False, 0)

        Image_Toolbar_VBox = Gtk.VBox()
        Image_Toolbar_VBox.set_border_width(3)
        Image_Toolbar_VBox.pack_start(Image_Toolbar_HBox1, False, False, 0)
        Image_Toolbar_VBox.pack_start(Image_Toolbar_HBox2, False, False, 0)

        Image_VBox = Gtk.VBox(homogeneous = False, spacing = 3)
        Image_VBox.set_border_width(3)
        Image_VBox.pack_start(self.Image_ScrolledWindow_EventBox, True, True, 0)
        Image_VBox.pack_start(Image_Toolbar_VBox, False, False, 0)

        Image_VBox.set_size_request(920,1000)

        self.Console_Entry = Gtk.Entry()
        self.Console_Entry.connect('activate', self.Console_Command_Sent)
        self.Console_Entry.connect('key-press-event', self.Console_KeyPressed)
        Console_EntryCompletion = Gtk.EntryCompletion()
        Console_EntryCompletion.set_model(Gtk.ListStore(str))
        Console_EntryCompletion.set_text_column(0)
        self.Console_Entry.set_completion(Console_EntryCompletion)
        Console_EntryCompletion.set_inline_completion(True)
        Console_EntryCompletion.set_popup_completion(False)

        AccelGrp = Gtk.AccelGroup.new() 
        FileMenuItem = Gtk.MenuItem(label = "File")
        FileMenu = Gtk.Menu()
        FileMenuItem.set_submenu(FileMenu)
        FileOpenMenuItem = Gtk.MenuItem(label = "Open")
        FileOpenMenuItem.add_accelerator("activate", AccelGrp, Gdk.KEY_F2, 0, Gtk.AccelFlags.VISIBLE)
        FileOpenMenuItem.connect("activate", self.FileDialog_Construction, 0)
        FileRefreshMenuItem = Gtk.MenuItem(label = "Refresh")
        FileRefreshMenuItem.add_accelerator("activate", AccelGrp, Gdk.KEY_F5, 0, Gtk.AccelFlags.VISIBLE)
        FileRefreshMenuItem.connect("activate", self.Folder_Refresh)
        FileQuitMenuItem = Gtk.MenuItem(label = "Quit")
        FileQuitMenuItem.add_accelerator("activate", AccelGrp, Gdk.KEY_Q, Gdk.ModifierType.CONTROL_MASK, Gtk.AccelFlags.VISIBLE)
        FileQuitMenuItem.connect("activate", self.MainWindow_Destroy)
        FileMenu.add(FileOpenMenuItem)
        FileMenu.add(FileRefreshMenuItem)
        FileMenu.add(FileQuitMenuItem)
        
        SettingMenuItem = Gtk.MenuItem(label = "Setting")
        SettingMenu = Gtk.Menu()
        SettingMenuItem.set_submenu(SettingMenu)
        SettingDetectorMenuItem = Gtk.MenuItem(label = "Detector")
        SettingDetectorMenu = Gtk.Menu()
        SettingDetectorMenuItem.set_submenu(SettingDetectorMenu)
        SettingDetector0MenuItem = Gtk.RadioMenuItem(label="1028 x 1062")
        SettingDetector0MenuItem.set_active(True)
        SettingDetector0MenuItem.connect("activate", self.Detector_Changed)
        SettingDetector1MenuItem = Gtk.RadioMenuItem.new_from_widget(SettingDetector0MenuItem)
        SettingDetector1MenuItem.set_label("516 x 516")
        SettingDetector1MenuItem.connect("activate", self.Detector_Changed)
        SettingDetector2MenuItem = Gtk.RadioMenuItem.new_from_widget(SettingDetector0MenuItem)
        SettingDetector2MenuItem.set_label("515 x 515")
        SettingDetector2MenuItem.connect("activate", self.Detector_Changed)
        SettingDetector3MenuItem = Gtk.RadioMenuItem.new_from_widget(SettingDetector0MenuItem)
        SettingDetector3MenuItem.set_label("256 x 256")
        SettingDetector3MenuItem.connect("activate", self.Detector_Changed)
        SettingDetector4MenuItem = Gtk.RadioMenuItem.new_from_widget(SettingDetector0MenuItem)
        SettingDetector4MenuItem.set_label("487 x 195")
        SettingDetector4MenuItem.connect("activate", self.Detector_Changed)
        SettingDetector5MenuItem = Gtk.RadioMenuItem.new_from_widget(SettingDetector0MenuItem)
        SettingDetector5MenuItem.set_label("1360 x 1024")
        SettingDetector5MenuItem.connect("activate", self.Detector_Changed)
        SettingDetectorMenu.add(SettingDetector0MenuItem)
        SettingDetectorMenu.add(SettingDetector1MenuItem)
        SettingDetectorMenu.add(SettingDetector2MenuItem)
        SettingDetectorMenu.add(SettingDetector3MenuItem)
        SettingDetectorMenu.add(SettingDetector4MenuItem)
        SettingDetectorMenu.add(SettingDetector5MenuItem)
        SettingColormapMenuItem = Gtk.MenuItem(label = "Colormap")
        SettingColormapMenu = Gtk.Menu()
        SettingColormapMenuItem.set_submenu(SettingColormapMenu)
        SettingColormap0MenuItem = Gtk.RadioMenuItem(label="jet")
        SettingColormap0MenuItem.connect("activate", self.Change_Colormap)
        SettingColormap1MenuItem = Gtk.RadioMenuItem.new_from_widget(SettingColormap0MenuItem)
        SettingColormap1MenuItem.set_label("binary")
        SettingColormap1MenuItem.connect("activate", self.Change_Colormap)
        SettingColormap2MenuItem = Gtk.RadioMenuItem.new_from_widget(SettingColormap0MenuItem)
        SettingColormap2MenuItem.set_label("gray")
        SettingColormap2MenuItem.connect("activate", self.Change_Colormap)
        SettingColormap3MenuItem = Gtk.RadioMenuItem.new_from_widget(SettingColormap0MenuItem)
        SettingColormap3MenuItem.set_label("hot")
        SettingColormap3MenuItem.connect("activate", self.Change_Colormap)
        SettingColormap4MenuItem = Gtk.RadioMenuItem.new_from_widget(SettingColormap0MenuItem)
        SettingColormap4MenuItem.set_label("coolwarm")
        SettingColormap4MenuItem.connect("activate", self.Change_Colormap)
        SettingColormap5MenuItem = Gtk.RadioMenuItem.new_from_widget(SettingColormap0MenuItem)
        SettingColormap5MenuItem.set_label("viridis")
        SettingColormap5MenuItem.set_active(True)
        SettingColormap5MenuItem.connect("activate", self.Change_Colormap)
        SettingColormap6MenuItem = Gtk.RadioMenuItem.new_from_widget(SettingColormap0MenuItem)
        SettingColormap6MenuItem.set_label("plasma")
        SettingColormap6MenuItem.connect("activate", self.Change_Colormap)
        SettingColormap7MenuItem = Gtk.RadioMenuItem.new_from_widget(SettingColormap0MenuItem)
        SettingColormap7MenuItem.set_label("inferno")
        SettingColormap7MenuItem.connect("activate", self.Change_Colormap)
        SettingColormap8MenuItem = Gtk.RadioMenuItem.new_from_widget(SettingColormap0MenuItem)
        SettingColormap8MenuItem.set_label("magma")
        SettingColormap8MenuItem.connect("activate", self.Change_Colormap)
        SettingColormap9MenuItem = Gtk.RadioMenuItem.new_from_widget(SettingColormap0MenuItem)
        SettingColormap9MenuItem.set_label("cividis")
        SettingColormap9MenuItem.connect("activate", self.Change_Colormap)
        SettingColormapMenu.add(SettingColormap0MenuItem)
        SettingColormapMenu.add(SettingColormap1MenuItem)
        SettingColormapMenu.add(SettingColormap2MenuItem)
        SettingColormapMenu.add(SettingColormap4MenuItem)
        SettingColormapMenu.add(SettingColormap5MenuItem)
        SettingColormapMenu.add(SettingColormap6MenuItem)
        SettingColormapMenu.add(SettingColormap7MenuItem)
        SettingColormapMenu.add(SettingColormap8MenuItem)
        SettingColormapMenu.add(SettingColormap9MenuItem)
        SettingSparseItem = Gtk.CheckMenuItem(label = "Sparse Matrix")
        SettingSparseItem.set_active(True)
        SettingSparseItem.connect("activate", self.SparseToggled)
        SettingShowAngleItem = Gtk.CheckMenuItem(label = "Show Angle")
        SettingShowAngleItem.set_active(False)
        SettingShowAngleItem.connect("activate", self.ShowAngleToggled)
        SettingDirtyFixItem = Gtk.CheckMenuItem(label = "Dirty Fix")
        SettingDirtyFixItem.set_active(True)
        SettingDirtyFixItem.connect("activate", self.DirtyFixToggled)
        SettingPumpProbeItem = Gtk.CheckMenuItem(label = "Pump Probe")
        SettingPumpProbeItem.set_active(False)
        SettingPumpProbeItem.connect("activate", self.PumpProbeToggled)

        AddonsMenuItem = Gtk.MenuItem(label = "Add-ons")
        AddonsMenu = Gtk.Menu()
        AddonsMenuItem.set_submenu(AddonsMenu)
        AddonMenuItem = Gtk.MenuItem(label = "User Shortcuts")
        AddonMenuItem.connect("activate", self.Addon_Shortcuts)
        AddonsMenu.add(AddonMenuItem)
        AddonMenuItem = Gtk.MenuItem(label = "Show Metadata")
        AddonMenuItem.connect("activate", self.Addon_ShowMetadata)
        AddonsMenu.add(AddonMenuItem)
        AddonMenuItem = Gtk.MenuItem(label = "Run PtychoLib")
        AddonMenuItem.connect("activate", self.Addon_PtychoLib)
        AddonsMenu.add(AddonMenuItem)
        AddonMenuItem = Gtk.MenuItem(label = "Probe Library")
        AddonMenuItem.connect("activate", self.Addon_ProbeLibrary)
        AddonsMenu.add(AddonMenuItem)
        AddonMenuItem = Gtk.MenuItem(label = "Powder Helper")
        AddonMenuItem.connect("activate", self.Addon_PowderHelper)
        AddonsMenu.add(AddonMenuItem)
        AddonMenuItem = Gtk.MenuItem(label = "Shift Correction")
        AddonMenuItem.connect("activate", self.Addon_ShiftCorrection)
        AddonsMenu.add(AddonMenuItem)
        AddonMenuItem = Gtk.MenuItem(label = "Analyze 5D")
        AddonMenuItem.connect("activate", self.Addon_Analyze5D)
        AddonsMenu.add(AddonMenuItem)
            
        HelpMenuItem = Gtk.MenuItem(label = "Help")
        HelpMenu = Gtk.Menu()
        HelpMenuItem.set_submenu(HelpMenu)
        HelpAboutMenuItem = Gtk.MenuItem(label = "About")
        HelpAboutMenuItem.connect("activate", self.AboutThisProgram)
        HelpMenu.add(HelpAboutMenuItem)
        
        SettingMenu.add(SettingDetectorMenuItem)
        SettingMenu.add(SettingColormapMenuItem)
        SettingMenu.add(SettingSparseItem)
        SettingMenu.add(SettingShowAngleItem)
        SettingMenu.add(SettingDirtyFixItem)
        SettingMenu.add(SettingPumpProbeItem)
        
        MenuBar = Gtk.MenuBar()
        MenuBar.add(FileMenuItem)
        MenuBar.add(SettingMenuItem)
        MenuBar.add(AddonsMenuItem)
        MenuBar.add(HelpMenuItem)
        
        self.Main_Window = Gtk.Window()
        self.Main_Window.add_accel_group(AccelGrp)

        HBox1 = Gtk.HBox(homogeneous = False, spacing = 3)
        HBox1.pack_start(MDA_File_ScrolledWindow, False, False, 0)
        HBox1.pack_start(Spec_ToolBox_VBox, False, False, 0)
        VBox1 = Gtk.VBox(homogeneous = False, spacing = 3)
        VBox1.pack_start(HBox1, False, False, 0)
        VBox1.pack_start(self.Plot_Notebook, False, False, 0)
        HBox2 = Gtk.HBox(homogeneous = False, spacing = 3)
        HBox2.pack_start(VBox1, False, False, 0)
        HBox2.pack_start(Image_VBox, False, False, 0)

        # LV1
        Main_VBox = Gtk.VBox(homogeneous = False, spacing = 3)
        Main_VBox.set_border_width(3)
        Main_VBox.pack_start(MenuBar, False, False, 0)
        Main_VBox.pack_start(HBox2, True, True, 0)
        Main_VBox.pack_start(self.Console_Entry, False, False, 0)

        # Main Window
        self.Main_Window.connect("destroy", self.MainWindow_Destroy)
        self.Main_Window.set_title("Diffraction User Data Explorer")
        self.Main_Window.add(Main_VBox)
        self.Main_Window.set_resizable(False)
        self.Main_Window.show_all()

        self.Image_ScrolledWindow_EventBox.get_window().set_cursor(Gdk.Cursor(Gdk.CursorType.TCROSS))
        self.pilatus_enabled = False
        self.eiger_enabled = True
        self.sparse_enabled= True
        self.show_angle = False
        self.dirty_fix = True
        self.pump_probe = False

#------------------------Main------------------------------#          

def main():
    
    Gtk.main()
    return 0

if __name__ == "__main__":

    MyMainWindow()
    main()
