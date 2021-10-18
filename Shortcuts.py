#!/usr/local/bin python3
"""addon for user shortcuts"""

from gi.repository import Gtk
import numpy as np
from scipy import ndimage, io, sparse
from matplotlib import colors
from math import *
from readMDA import *
import time
import h5py

try:
    import epics
except:
    pass


class MainWindow:

    def SumImages(self, widget, dude):

        if not dude.AllImagesLoaded:
            dude.LoadAllImages(None)
        if dude.sparse_enabled:
            dude.Image_Image.set_array(dude.image.mean(0).reshape(1062,1028))
        else:
            dude.Image_Image.set_array(dude.image.mean(0))
        dude.Image_Image.set_norm(colors.LogNorm())
        dude.Image_Canvas.draw()


    def MaxMap(self, widget, dude):

        if not dude.AllImagesLoaded:
            dude.LoadAllImages(None)
        dude.Plot2D_Image.set_array(dude.image.max(1).max(1).reshape(dude.Plot2D_ydata.shape))
        dude.Plot2D_Image.set_norm(colors.Normalize())
        dude.Plot2D_Canvas.draw()
        dude.Plot2D_ydata = data
        dude.Plot_Notebook.set_current_page(1)

        
    def OverviewMap(self, widget, dude):

        dude.Plot2D_Axe.cla()
        dude.Plot2D_Axe.set_axis_off()
        dude.Plot2D_Axe.set_aspect("equal")
        ymin = 1e10
        ymax = 1
        lims = (9e9,-9e9,9e9,-9e9)
        IMs = []
        for row in dude.MDA_File_ListStore:
            if row[11] == True: # row[11] is the multiselect toggle button
                if dude.Scan_ToolBox_U_ToggleButton.get_active(): # if custom ROI is selected
                    pass # do nothing = not implemented
                else:
                    mda = readMDA(row[10], verbose=0) # row[10] contains the path of the mda file
                    if dude.eiger_enabled:
                        h5_filename = [f.name for f in os.scandir(dude.h5_folder) if "scan_{0}".format(row[0]) in f.name][0]
                        #print(h5_filename)
                        h5 = h5py.File(os.path.join(dude.h5_folder, h5_filename), 'r')
                        xdata2_0 = h5["/entry/instrument/26-ID-C/ATTO SAM Z"][()] 
                        theta = h5["/entry/instrument/26-ID-C/SAMPLE THETA"][()]
                        #print(np.sin(np.radians(theta)))
                    else:
                        for d in mda[2].d:
                            if d.name == "atto2:m3.RBV":
                                xdata2_0 = np.array(d.data)
                            if d.name == "atto2:PIC867:1:m1.RBV": 
                                theta = np.array(d.data)#-90
                                print (theta.mean())
                    datatmp = np.array(mda[2].p[0].data)
                    dimy, dimx = datatmp.shape[0], datatmp.shape[1]
                    xdata2 = np.ones((dimy,dimx))
                    xdata1 = np.copy(xdata2)#+1
                    #ydata = np.copy(xdata2)
                    xdata2 *= datatmp
                    xdata1 *= np.array(mda[1].p[0].data)[:dimy,np.newaxis]
                    if mda[1].p[0].name == "26idcnpi:Y_HYBRID_SP.VAL":
                        xdata1_0 = xdata1
                    elif mda[1].p[0].name == "26idcnpi:X_HYBRID_SP.VAL":
                        #xdata2_0 = xdata1 - xdata2_0 * np.cos(np.radians(theta))
                        xdata2_0 =  -xdata1 /  np.sin(np.radians(theta)) + xdata2_0
                    if mda[2].p[0].name == "26idcnpi:Y_HYBRID_SP.VAL":
                        xdata1_0 = xdata2
                    elif mda[2].p[0].name == "26idcnpi:X_HYBRID_SP.VAL":
                        xdata2_0 =  -xdata2 /  np.sin(np.radians(theta)) + xdata2_0
                    if mda[2].p[0].name == "26idbATTO:m3.VAL":
                        for d in mda[2].d:
                            if d.name == "26idcnpi:m35.RBV":
                                xdata1_0 = np.array(d.data)
                            if d.name == "26idcnpi:m34.RBV":
                                xdata2_0 = np.array(d.data)
                        xdata2_0 = -xdata2_0 / np.sin(np.radians(theta)) + xdata2
                    ydata = np.array(mda[2].d[dude.Scan_ToolBox_Y_ComboBox.get_active()].data)
                    if dude.dirty_fix and "eiger" in mda[2].d[dude.Scan_ToolBox_Y_ComboBox.get_active()].name and mda[1].time != "whatever":
                        #print("dirty fixing first point bug")
                        ydata_flat = ydata.flatten()
                        ydata_flat[:-1] = ydata_flat[1:]
                        ydata = ydata_flat.reshape(ydata.shape)
                    #dude.Plot2D_Axe.pcolormesh(xdata2, xdata1, ydata, cmap=dude.cm, alpha=0.5)
                    IMs += [dude.Plot2D_Axe.imshow(ydata,cmap=dude.cm, alpha=0.5, extent=(xdata2_0[0,0],xdata2_0[-1,-1],xdata1_0[0,0],xdata1_0[-1,-1]), origin="lower")]
                    lims = (min(lims[0],min(xdata2_0[0,0],xdata2_0[-1,-1])),\
                            max(lims[1],max(xdata2_0[0,0],xdata2_0[-1,-1])),\
                            min(lims[2],min(xdata1_0[0,0],xdata1_0[-1,-1])),\
                            max(lims[3],max(xdata1_0[0,0],xdata1_0[-1,-1])))
                    ymin = max(1,min(ymin, ydata.min()))
                    ymax = max(ymax, ydata.max())
                    
        for IM in IMs:
            if dude.Plot2D_Log_ToggleButton.get_active():
                IM.set_norm(colors.LogNorm(ymin, ymax))
            else:
                IM.set_norm(colors.Normalize(ymin, ymax))
        dude.Plot2D_Axe.set_xlim(lims[0],lims[1])
        dude.Plot2D_Axe.set_ylim(lims[2],lims[3])
        dude.Plot2D_Axe.set_aspect(1)
        dude.Plot2D_Canvas.draw()
        dude.Plot2D_Figure.savefig("overview.png", dpi=300)


    def CoMX(self, widget, dude):

        if not dude.AllImagesLoaded:
            dude.LoadAllImages(None)
        ymax = int(dude.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.get_value())
        ymin = int(dude.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.get_value())
        xmax = int(dude.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.get_value())
        xmin = int(dude.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.get_value())
        if not dude.sparse_enabled:
            data = (((dude.image[:,ymin:ymax+1,xmin:xmax+1].sum(1)*np.arange(xmin,xmax+1)).sum(1)) /\
                    (dude.image[:,ymin:ymax+1,xmin:xmax+1].sum(1).sum(1))).reshape(dude.Plot2D_ydata.shape)
        else:
            m0 = np.arange(xmin,xmax+1)*np.ones((ymax-ymin+1,1))
            m0 = sparse.csr_matrix(m0.flatten()*np.ones((dude.image.shape[0],1)))
            m1 = np.zeros((1062,1028))
            m1[ymin:ymax+1,xmin:xmax+1] = 1
            m1 = m1.flatten().nonzero()[0]
            data = ((m0.multiply(dude.image[:,m1])).sum(1) / dude.image[:,m1].sum(1)).reshape(dude.Plot2D_ydata.shape)
        dude.Plot2D_Image.set_array(data)
        dude.Plot2D_Image.set_norm(colors.Normalize())
        dude.Plot2D_Image.set_cmap("coolwarm")
        dude.Plot2D_Canvas.draw()
        dude.Plot2D_ydata = data
        dude.Plot_Notebook.set_current_page(1)


    def CoMY(self, widget, dude):

        if not dude.AllImagesLoaded:
            dude.LoadAllImages(None)
        ymax = int(dude.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.get_value())
        ymin = int(dude.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.get_value())
        xmax = int(dude.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.get_value())
        xmin = int(dude.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.get_value())
        if not dude.sparse_enabled:
            data = (((dude.image[:,ymin:ymax+1,xmin:xmax+1].sum(2)*np.arange(ymin,ymax+1)).sum(1)) /\
                    (dude.image[:,ymin:ymax+1,xmin:xmax+1].sum(1).sum(1))).reshape(dude.Plot2D_ydata.shape)
        else:
            m0 = np.arange(ymin,ymax+1)[:,np.newaxis]*np.ones(xmax-xmin+1)
            m0 = sparse.csr_matrix(m0.flatten()*np.ones((dude.image.shape[0],1)))
            m1 = np.zeros((1062,1028))
            m1[ymin:ymax+1,xmin:xmax+1] = 1
            m1 = m1.flatten().nonzero()[0]
            data = ((m0.multiply(dude.image[:,m1])).sum(1) / dude.image[:,m1].sum(1)).reshape(dude.Plot2D_ydata.shape)
        dude.Plot2D_Image.set_array(data)
        dude.Plot2D_Image.set_norm(colors.Normalize())
        dude.Plot2D_Canvas.draw()
        dude.Plot2D_ydata = data
        dude.Plot_Notebook.set_current_page(1)


    def SumX(self, widget, dude):

        dude.Plot1D_Axe.cla()
        ymax = int(dude.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.get_value())
        ymin = int(dude.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.get_value())
        xmax = int(dude.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.get_value())
        xmin = int(dude.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.get_value())
        ydata = np.copy(dude.Image_Image.get_array()[ymin:ymax+1,xmin:xmax+1]).sum(1)
        xdata = np.arange(ymin,ymax+1)
        dude.Plot1D_Axe.plot(xdata, ydata)
        dude.Plot1D_Axe.set_xscale('log' if dude.Plot1D_LogX_ToggleButton.get_active() else 'linear')
        dude.Plot1D_Axe.set_yscale('log' if dude.Plot1D_LogY_ToggleButton.get_active() else 'linear')
        dude.Plot1D_xdata = xdata
        dude.Plot1D_ydata = ydata
        dude.Plot1D_Canvas.draw()
        dude.Plot_Notebook.set_current_page(0)

        
    def SumY(self, widget, dude):

        dude.Plot1D_Axe.cla()
        ymax = int(dude.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.get_value())
        ymin = int(dude.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.get_value())
        xmax = int(dude.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.get_value())
        xmin = int(dude.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.get_value())
        ydata = np.copy(dude.Image_Image.get_array().data[ymin:ymax+1,xmin:xmax+1]).sum(0)
        xdata = np.arange(xmin,xmax+1)
        dude.Plot1D_Axe.plot(xdata, ydata)
        dude.Plot1D_Axe.set_xscale('log' if dude.Plot1D_LogX_ToggleButton.get_active() else 'linear')
        dude.Plot1D_Axe.set_yscale('log' if dude.Plot1D_LogY_ToggleButton.get_active() else 'linear')
        dude.Plot1D_xdata = xdata
        dude.Plot1D_ydata = ydata
        dude.Plot1D_Canvas.draw()
        dude.Plot_Notebook.set_current_page(0)

        
    def Qxz(self, widget, dude):

        theta = -np.radians(np.array(dude.data[1].p[0].data))
        twotheta = np.radians(dude.h5["/entry/instrument/26-ID-C/Det Two Theta"][()])
        wavelength = 12398./dude.h5["/entry/instrument/26-ID-B/DCM Energy"][()]

        if not dude.AllImagesLoaded:
            dude.LoadAllImages(None)
        ymax = int(dude.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.get_value())
        ymin = int(dude.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.get_value())
        xmax = int(dude.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.get_value())
        xmin = int(dude.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.get_value())
        if not dude.sparse_enabled:
            data = dude.image[:,ymin:ymax+1,xmin:xmax+1].sum(1).T
        else:
            data = (dude.image.toarray().reshape(theta.shape[0], 1062, 1028))[:,ymin:ymax+1, xmin:xmax+1].sum(1).T

        ###
        d0 = 0.18
        distance = np.c_[ np.sqrt(d0**2 + ((np.arange(xmin, xmax+1)-514)*75e-6+0.085)**2) ]
        ###
        print(distance.min(), distance.max())

        #twotheta = twotheta + np.c_[(np.arange(xmin, xmax+1)-514)]*75e-6/distance
        twotheta = np.c_[twotheta + np.arctan(((np.arange(xmin, xmax+1)-514)*75e-6+0.085)/d0)]
        gamma = 0
        #gamma = gamma + np.c_[258-(np.arange(ymin, ymax+1))]*55e-6/distance
        print(twotheta.min()/3.14*180, twotheta.max()/3.14*180)

        K = 2*np.pi/wavelength
        qx = K*np.sin(twotheta)*np.cos(gamma)
        qz = K*np.cos(twotheta)*np.cos(gamma)-K

        qz, qx = np.cos(theta)*qz-np.sin(theta)*qx, np.sin(theta)*qz+np.cos(theta)*qx
        print(qz.max(), qz.min(), qx.max(), qx.min())
        np.savez_compressed("Qxy.npz", data=data, qz=qz, qx=qx, theta=theta, twotheta = twotheta)

        dude.Plot2D_Axe.cla()
        dude.Plot2D_Axe.set_axis_off()
        #dude.Plot2D_Image = dude.Plot2D_Axe.scatter(qx, qz, c=data, s=1, edgecolor='none', cmap=dude.cm, norm=colors.LogNorm())
        dude.Plot2D_Image = dude.Plot2D_Axe.pcolormesh(qz, qx, data, cmap=dude.cm, norm=colors.LogNorm(), shading='auto')
        dude.Plot2D_Axe.set_aspect(1)
        dude.Plot2D_Canvas.draw()
        dude.Plot_Notebook.set_current_page(1)
    

    def defRoI(self, widget, dude, flag):

        ymax = int(dude.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.get_value())
        ymin = int(dude.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.get_value())
        xmax = int(dude.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.get_value())
        xmin = int(dude.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.get_value())
        if dude.dimY == 516 and dude.dimX == 516:
            epics.caput("QMPX3:ROI{0}:MinX".format(flag), xmin, wait=True)
            epics.caput("QMPX3:ROI{0}:SizeX".format(flag), xmax-xmin, wait=True)
            epics.caput("QMPX3:ROI{0}:MinY".format(flag), ymin, wait=True)
            epics.caput("QMPX3:ROI{0}:SizeY".format(flag), ymax-ymin, wait=True)
        elif dude.dimY == 195 and dude.dimX == 487:
            epics.caput("dp_pilatus4:ROI{0}:MinX".format(flag), xmin, wait=True)
            epics.caput("dp_pilatus4:ROI{0}:SizeX".format(flag), xmax-xmin, wait=True)
            epics.caput("dp_pilatus4:ROI{0}:MinY".format(flag), ymin, wait=True)
            epics.caput("dp_pilatus4:ROI{0}:SizeY".format(flag), ymax-ymin, wait=True)
        elif dude.dimY == 1062 and dude.dimX == 1028:
            epics.caput("s26_eiger_cnm:ROI{0}:MinX".format(flag), xmin, wait=True)
            epics.caput("s26_eiger_cnm:ROI{0}:SizeX".format(flag), xmax-xmin, wait=True)
            epics.caput("s26_eiger_cnm:ROI{0}:MinY".format(flag), ymin, wait=True)
            epics.caput("s26_eiger_cnm:ROI{0}:SizeY".format(flag), ymax-ymin, wait=True)

            
    def ShowLatest(self, widget, dude):

        if dude.dimY == 516 and dude.dimX == 516:
            img = epics.caget("QMPX3:image1:ArrayData", as_numpy=1, count=516*516).reshape(516,516)
        elif dude.dimY == 195 and dude.dimX == 487:
            img = epics.caget("dp_pilatus4:image1:ArrayData", as_numpy=1, count=195*487).reshape(195,487)
        elif dude.dimY == 1062 and dude.dimX == 1028:
            img = epics.caget("s26_eiger_cnm:image1:ArrayData", as_numpy=1, count=1062*1028).reshape(1062,1028)
        dude.Image_Image.set_array(img.astype(np.int32))
        dude.Image_Image.set_norm(colors.LogNorm())
        dude.Image_Canvas.draw()


    def GotoP0(self, widget, dude):

        if dude.MDA_File_ListStore[dude.mda_selection_path[0]][1] == "hybridy" and dude.MDA_File_ListStore[dude.mda_selection_path[0]][5] == "hybridx":
            lelements = dude.Plot2D_P0_Label.get_text().split()
            x0, y0 = float(lelements[3][:-1]), float(lelements[6][:-1])
            dialog = Gtk.MessageDialog(dude.Main_Window, 0, Gtk.MessageType.QUESTION, Gtk.ButtonsType.YES_NO, "Moving hybridy to {0:.3f} and hybridx to {1:.3f}?".format(y0, x0))
            response = dialog.run()
            if response == Gtk.ResponseType.YES:
                epics.caput("26idcnpi:X_HYBRID_SP.VAL", x0)
                time.sleep(.1)
                epics.caput("26idcnpi:Y_HYBRID_SP.VAL", y0)
            dialog.destroy()
        elif dude.MDA_File_ListStore[dude.mda_selection_path[0]][1] == "samy":
            if epics.caget("26idc:1:userCalc7.VAL") or epics.caget("26idc:1:userCalc5.VAL"):
                print("unlock hybrid before doing this")
                return
            lelements = dude.Plot2D_P0_Label.get_text().split()
            x0, y0 = float(lelements[3][:-1]), float(lelements[6][:-1])
            if dude.MDA_File_ListStore[dude.mda_selection_path[0]][5][:4] == "atto":
                z_or_x = dude.MDA_File_ListStore[dude.mda_selection_path[0]][5][4]
                dialog = Gtk.MessageDialog(dude.Main_Window, 0, Gtk.MessageType.QUESTION, Gtk.ButtonsType.YES_NO, "Moving samy to {0:.3f} and atto{1} to {2:.3f}?".format(y0, z_or_x, x0))
                response = dialog.run()
                if response == Gtk.ResponseType.YES:
                    epics.caput("26idcnpi:m17.VAL", y0)
                    time.sleep(.1)
                    if z_or_x == "z":
                        epics.caput("atto2:m3.VAL", x0)
                    elif z_or_x == "x":
                        epics.caput("atto2:m4.VAL", x0)
                    else:
                        print("impossible!", z_or_x)
                dialog.destroy()
        else:
            print("unrecognized motors!")

 
    def ShowSpiral(self, widget, dude):

        for i in range(dude.data[1].nd):
            dname = dude.data[1].d[i].name
            if "m34" in dname:
                xx = dude.data[1].d[i].data
            elif "m35" in dname:
                yy = dude.data[1].d[i].data
        dude.PlotSpare_Axe.cla()
        dude.PlotSpare_Axe.set_aspect(1)
        dude.PlotSpare_Axe.scatter(xx,yy,c=dude.Plot1D_ydata, s=1500/sqrt(len(yy)), cmap=dude.cm, norm=colors.LogNorm())
        dude.PlotSpare_Canvas.draw()
        dude.Plot_Notebook.set_current_page(2)


    def SaveMat(self, widget, dude):

        if not dude.AllImagesLoaded:
            dude.LoadAllImages(None)
        mdict = {}
        ymax = int(dude.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.get_value())
        ymin = int(dude.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.get_value())
        xmax = int(dude.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.get_value())
        xmin = int(dude.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.get_value())
        if dude.sparse_enabled:
            m1 = np.zeros((1062,1028))
            m1[ymin:ymax+1,xmin:xmax+1] = 1
            m1 = m1.flatten().nonzero()[0]
            t0=time.time()
            mdict["data"] = dude.image[:,m1].toarray().reshape(dude.image.shape[0],ymax-ymin+1,xmax-xmin+1)
            if "/entry/scan/Positioner" in dude.h5:
                pos_keys = list(dude.h5["/entry/scan/Positioner"].keys())
                for key in pos_keys:
                    mdict[key] = np.array(dude.h5["/entry/scan/Positioner"].get(key))
            mdict["theta"] = np.array(dude.h5["/entry/instrument/26-ID-C/SAMPLE THETA"])
        else:
            mdict["data"] = dude.image[:,ymin:ymax+1,xmin:xmax+1]
        mdict["roi"] = [ymin, ymax, xmin, xmax]

        mat_folder = os.path.join(dude.h5_folder, "mat")
        if not os.path.exists(mat_folder):
            os.mkdir(mat_folder)
            os.chmod(mat_folder, 0o777)

        io.savemat(os.path.join(mat_folder,"scan_{0:04d}.mat".format(int(dude.MDA_File_ListStore[dude.mda_selection_path[0]][0]))), mdict)
        print("done")

        
    def Flip(self, widget, dude):

        if dude.MDA_File_ListStore[dude.mda_selection_path[0]][1] == "hybridy":
            xdata = np.array(dude.data[2].p[0].data)
            if xdata[0,1]>xdata[0,0] and not dude.Plot2D_Axe.xaxis_inverted():
                dude.Plot2D_Axe.invert_xaxis()
            ydata = np.array(dude.data[1].p[0].data)
            if ydata[1]<ydata[0] and not dude.Plot2D_Axe.yaxis_inverted():
                dude.Plot2D_Axe.invert_yaxis()
        elif dude.MDA_File_ListStore[dude.mda_selection_path[0]][1] == "samy":
            xdata = np.array(dude.data[2].p[0].data)
            if xdata[0,1]<xdata[0,0]  and not dude.Plot2D_Axe.xaxis_inverted():
                dude.Plot2D_Axe.invert_xaxis()
            ydata = np.array(dude.data[1].p[0].data)
            if ydata[1]>ydata[0] and not dude.Plot2D_Axe.yaxis_inverted():
                dude.Plot2D_Axe.invert_yaxis()
        dude.Plot2D_Canvas.draw()


    def Stretch(self, widget, dude):

        if dude.eiger_enabled:
            theta = dude.h5["/entry/instrument/26-ID-C/SAMPLE THETA"][()]
        else:
            for d in dude.data[2].d:
                if d.name == "atto2:PIC867:1:m1.RBV":
                    theta = np.array(d.data).mean()
        dude.Plot2D_Axe.set_aspect(sin(radians(theta)))
        dude.Plot2D_Canvas.draw()


    def MainWindow_Destroy(self, widget, dude): 

        dude.Shortcuts = None
        self.win.destroy()

            
    def __init__(self, dude):

        button_list = []

        SumImages_Button = Gtk.Button("Sum Images")
        SumImages_Button.set_tooltip_text("Average over all the images in the scan, and show it in log scale")
        SumImages_Button.connect("clicked", self.SumImages, dude)
        button_list += [SumImages_Button]

        MaxMap_Button = Gtk.Button("Max Map")
        MaxMap_Button.set_tooltip_text("Show the max intensity of each point on the map")
        MaxMap_Button.connect("clicked", self.MaxMap, dude)
        button_list += [MaxMap_Button]

        OverviewMap_Button = Gtk.Button("Overview Map")
        OverviewMap_Button.set_tooltip_text("Show overview map on multiselected scans.")
        OverviewMap_Button.connect("clicked", self.OverviewMap, dude)
        button_list += [OverviewMap_Button]

        ShowSpiral_Button = Gtk.Button("ShowSpiral")
        ShowSpiral_Button.set_tooltip_text("1")
        ShowSpiral_Button.connect("clicked", self.ShowSpiral, dude)
        button_list += [ShowSpiral_Button]
        
        CoMX_Button = Gtk.Button("X CoM in RoI")
        CoMX_Button.set_tooltip_text("Calculate the CoM in the X direction in the RoI defined by your right mouse button, the data in the Y direction is directly integrated")
        CoMX_Button.connect("clicked", self.CoMX, dude)
        button_list += [CoMX_Button]

        CoMY_Button = Gtk.Button("Y CoM in RoI")
        CoMY_Button.set_tooltip_text("Calculate the CoM in the Y direction in the RoI defined by your right mouse button, the data in the X direction is directly integrated")
        CoMY_Button.connect("clicked", self.CoMY, dude)
        button_list += [CoMY_Button]

        SumX_Button = Gtk.Button("X Sum in RoI")
        SumX_Button.set_tooltip_text("Sum the data in the X direction in the RoI defined by your right mouse button, and show the result along the Y direction")
        SumX_Button.connect("clicked", self.SumX, dude)
        button_list += [SumX_Button]

        SumY_Button = Gtk.Button("Y Sum in RoI")
        SumY_Button.set_tooltip_text("Sum the data in the Y direction in the RoI defined by your right mouse button, and show the result along the X direction")
        SumY_Button.connect("clicked", self.SumY, dude)
        button_list += [SumY_Button]

        for i in range(1,5):
            RoI_Button = Gtk.Button("Redefine RoI{0}".format(i))
            RoI_Button.set_tooltip_text("Set the RoI you scan with as the RoI you just drew")
            RoI_Button.connect("clicked", self.defRoI, dude, i)
            button_list += [RoI_Button]

        Qxz_Button = Gtk.Button("Generate Qxz")
        Qxz_Button.set_tooltip_text("")
        Qxz_Button.connect("clicked", self.Qxz, dude)
        button_list += [Qxz_Button]

        ShowLatest_Button = Gtk.Button("Latest Image")
        ShowLatest_Button.set_tooltip_text("Show the latest image to help you redefine your RoIs.")
        ShowLatest_Button.connect("clicked", self.ShowLatest, dude)
        button_list += [ShowLatest_Button]

        GotoP0_Button = Gtk.Button("Goto P0")
        GotoP0_Button.set_tooltip_text("Move Piezo X Y to position indicated by P0")
        GotoP0_Button.connect("clicked", self.GotoP0, dude)
        button_list += [GotoP0_Button]

        SaveMat_Button = Gtk.Button("RoI to mat")
        SaveMat_Button.set_tooltip_text("save roi data to matlab file under the h5 folder")
        SaveMat_Button.connect("clicked", self.SaveMat, dude)
        button_list += [SaveMat_Button]

        Flip_Button = Gtk.Button("Flip Plot2D")
        Flip_Button.set_tooltip_text("automatically flipping, the cursor values are still correct after flipping")
        Flip_Button.connect("clicked", self.Flip, dude)
        button_list += [Flip_Button]

        Stretch_Button = Gtk.Button("Stretch Plot2D")
        Stretch_Button.set_tooltip_text("automatically stretching, the cursor values are still correct after stretching")
        Stretch_Button.connect("clicked", self.Stretch, dude)
        button_list += [Stretch_Button]

        grid = Gtk.Grid()
        nrow = int(sqrt(len(button_list)))

        for i in range(len(button_list)):
            grid.attach(button_list[i], i%nrow, i/nrow, 1, 1)
        
        self.win = Gtk.Window()
        self.win.connect("destroy", self.MainWindow_Destroy, dude)
        self.win.add(grid)
        self.win.set_title("User Shortcuts")
        self.win.show_all()

"""
 def SumSubsets(self, widget, dude):

        istart = int(dude.Image_Plot_HScale_Adjustment.get_value())
        image_file = dude.image_list[0]
        prefix = image_file[:-10]
        for i in range(istart, istart+dude.nbin):
            image_file2 = prefix+"{0:06d}.tif".format(i)
            try:
                data += fabio.open(os.path.join(dude.image_path,image_file2)).data
            except:
                data = fabio.open(os.path.join(dude.image_path,image_file2)).data
        dude.Image_Image.set_array(data)
        dude.Image_Image.set_norm(colors.LogNorm())
        dude.Image_Canvas.draw()

    def flipfly(self, widget, dude):

        tt = dude.Plot2D_Image.get_array();
        tt[1::2] = tt[1::2][:,::-1]
        dude.Plot2D_Image.set_array(tt);
        dude.Plot2D_Canvas.draw()


    def Joe1(self, widget, dude):

        if not dude.AllImagesLoaded:
            dude.LoadAllImages(None)
        mdict = {}
        mdict["data"] = {}
        ymax = int(dude.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.get_value())
        ymin = int(dude.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.get_value())
        xmax = int(dude.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.get_value())
        xmin = int(dude.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.get_value())
        mdict["data"]["xcen_mean"] = (((dude.image[:,ymin:ymax+1,xmin:xmax+1].sum(1)*np.arange(xmin,xmax+1)).sum(1)) /\
                                     (dude.image[:,ymin:ymax+1,xmin:xmax+1].sum(1).sum(1))).reshape(dude.Plot2D_ydata.shape) 
        mdict["data"]["ycen_mean"] = (((dude.image[:,ymin:ymax+1,xmin:xmax+1].sum(2)*np.arange(ymin,ymax+1)).sum(1)) /\
                                      (dude.image[:,ymin:ymax+1,xmin:xmax+1].sum(1).sum(1))).reshape(dude.Plot2D_ydata.shape)
        mdict["data"]["Int_mean"] = dude.image.mean((1,2)).reshape(dude.Plot2D_ydata.shape)
        mdict["data"]["Int_var"] = dude.image.var((1,2)).reshape(dude.Plot2D_ydata.shape)
        io.savemat("scan_{0}.mat".format(dude.MDA_File_ListStore[dude.mda_selection_path[0]][0]), mdict)
        print("done")


    def Wen1(self, widget, dude):

        dude.Plot1D_Axe.cla()
        ydata = 0
        cnt = 0
        for row in dude.MDA_File_ListStore:
            if row[11] == True: # row[11] is the multiselect toggle button
                if dude.Scan_ToolBox_U_ToggleButton.get_active(): # if custom ROI is selected
                    pass # do nothing = not implemented
                else:
                    mda = readMDA(row[10], verbose=0) # row[10] contains the path of the mda file
                    xdata = np.array(mda[1].p[0].data)
                    ydata += np.array(mda[1].d[dude.Scan_ToolBox_Y_ComboBox.get_active()].data)
                    cnt += 1

        dude.Plot1D_Axe.plot(xdata, ydata/cnt)
        dude.Plot1D_Canvas.draw()



   def Ahn1(self, widget, dude):

        i = 0
        for row in dude.MDA_File_ListStore:
            if row[11] == True: # row[11] is the multiselect toggle button
                if dude.Scan_ToolBox_U_ToggleButton.get_active(): # if custom ROI is selected
                    pass # do nothing = not implemented
                else:
                    if i==0:
                        mda = np.array((readMDA(row[10], verbose=0))[1].d[dude.Scan_ToolBox_Y_ComboBox.get_active()].data)
                    else:
                        mda += np.array((readMDA(row[10], verbose=0))[1].d[dude.Scan_ToolBox_Y_ComboBox.get_active()].data)
                    i += 1 
        mda /= i
        dude.Plot1D_Axe.plot(dude.Plot1D_xdata, mda)
        dude.Plot1D_Canvas.draw()

    def Tao3(self, widget, dude):

        scannum = dude.MDA_File_ListStore[dude.mda_selection_path[0]][0]
        data = dude.image.reshape(200,200,129,4,129,4).sum(5).sum(3)
        shiftdir = "/CNMshare/savedata/pythonscripts/Tao/IH-MA-73/shiftcorr3"
        shift = np.load(os.path.join(shiftdir, "shift{0}.npz".format(scannum)))

        data = np.roll(data, shift["y"], axis=0)
        for i in range(200):
            data[i] = np.roll(data[i], int(shift["x"][i]), axis=0)
    
        #data = data[53:193,30:170]
        homox = 195
        homoy = 303
        xx, yy = np.meshgrid(np.arange(516).reshape(129,4).mean(1),np.arange(516).reshape(129,4).mean(1))
        K = 2*np.pi/1.5498
        eta0 = np.radians(34.8018)
        qz = K*np.sin(np.arctan((homoy-yy)*55e-6/0.98)+eta0)+K*np.sin(eta0)
        qy = -K*np.cos(np.arctan((homoy-yy)*55e-6/0.98)+eta0)+K*np.cos(eta0)
        qx = K*np.sin(np.arctan((xx-homox)*55e-6/0.98))
        qy += (K-K*np.cos(np.arctan((xx-homox)*55e-6/0.98)))
        #d_eta = (int(scannum)-35)*0.04/180*np.pi
        d_eta = (int(scannum)-8)*0.04/180*np.pi
        #d_eta = (int(scannum)-105)*0.05/180*np.pi
        qy, qz = qz*sin(-d_eta)+qy*cos(-d_eta), qz*cos(-d_eta)-qy*sin(-d_eta)
        np.savez_compressed(shiftdir+"/cor{0}.npz".format(scannum), data=data, qx=qx, qy=qy, qz=qz)
        dude.Plot2D_Image.set_array(data.sum(2).sum(2))
        dude.Plot2D_Image.set_norm(colors.LogNorm())
        dude.Plot2D_Canvas.draw()

    def Tao1(self, widget, dude):

        scannum = dude.MDA_File_ListStore[dude.mda_selection_path[0]][0]
        data = dude.image.reshape(200,200,129,4,129,4).sum(5).sum(3)
        shiftdir = "/CNMshare/savedata/pythonscripts/Tao/IH-MA-73/shiftcorr2"
        shift = np.load(os.path.join(shiftdir, "shift{0}.npz".format(scannum)))

        data = np.roll(data, shift["y"], axis=0)
        for i in range(200):
            data[i] = np.roll(data[i], int(shift["x"][i]), axis=0)
    
        #data = data[53:193,30:170]
        homox = 195
        homoy = 303
        xx, yy = np.meshgrid(np.arange(516).reshape(129,4).mean(1),np.arange(516).reshape(129,4).mean(1))
        K = 2*np.pi/1.5498
        eta0 = np.radians(34.8018)
        qz = K*np.sin(np.arctan((homoy-yy)*55e-6/0.98)+eta0)+K*np.sin(eta0)
        qy = -K*np.cos(np.arctan((homoy-yy)*55e-6/0.98)+eta0)+K*np.cos(eta0)
        qx = K*np.sin(np.arctan((xx-homox)*55e-6/0.98))
        qy += (K-K*np.cos(np.arctan((xx-homox)*55e-6/0.98)))
        #d_eta = (int(scannum)-35)*0.04/180*np.pi
        d_eta = (int(scannum)-105)*0.05/180*np.pi
        qy, qz = qz*sin(-d_eta)+qy*cos(-d_eta), qz*cos(-d_eta)-qy*sin(-d_eta)
        np.savez_compressed("cor{0}.npz".format(scannum), data=data, qx=qx, qy=qy, qz=qz)
        dude.Plot2D_Image.set_array(data.sum(2).sum(2))
        dude.Plot2D_Image.set_norm(colors.LogNorm())
        dude.Plot2D_Canvas.draw()


#from skimage.measure import compare_ssim as ssim
    def SSIM_1st(self, widget, dude):

        if not dude.AllImagesLoaded:
            dude.LoadAllImages(None)
        ymax = int(dude.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.get_value())
        ymin = int(dude.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.get_value())
        xmax = int(dude.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.get_value())
        xmin = int(dude.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.get_value())
        d1 = np.zeros(dude.Plot2D_ydata.shape)
        data = dude.image[:,ymin:ymax+1,xmin:xmax+1].reshape(dude.Plot2D_ydata.shape[0], dude.Plot2D_ydata.shape[1],\
                                                             ymax-ymin+1, xmax-xmin+1)
        for i in range(1,dude.Plot2D_ydata.shape[0]-1):
            for j in range(1,dude.Plot2D_ydata.shape[1]-1):
                d1[i,j] = ssim(data[i,j], data[i+1,j])+\
                          ssim(data[i,j], data[i-1,j])+\
                          ssim(data[i,j], data[i,j+1])+\
                          ssim(data[i,j], data[i,j-1])
        dude.Plot2D_Image.set_array(d1)
        dude.Plot2D_Image.set_norm(colors.Normalize())
        dude.Plot2D_Canvas.draw()
        dude.Plot2D_ydata = d1



   def RotImages(self, widget, dude):

        if not dude.AllImagesLoaded:
            dude.LoadAllImages(None)
        gamma = epics.caget("26idcSOFT:sm4.RBV")
        twotheta = epics.caget("26idcSOFT:sm3.RBV")
        rotangle = -degrees(asin(sin(radians(gamma))/tan(radians(twotheta))))
        dude.image = ndimage.rotate(dude.image, rotangle, axes=(1,2), reshape=False)


    def OverviewFly(self, widget, dude):

        dude.Plot2D_Axe.cla()
        dude.Plot2D_Axe.set_axis_off()
        dude.Plot2D_Axe.set_aspect("equal")
        ymin = 1e10
        ymax = 1
        lims = (9e9,-9e9,9e9,-9e9)
        IMs = []
        for row in dude.MDA_File_ListStore:
            if row[11] == True: # row[11] is the multiselect toggle button
                mda = readMDA(row[10], verbose=0)
                datatmp = np.array(mda[2].p[0].data)
                dimy, dimx = datatmp.shape[0], datatmp.shape[1]
                xdata2 = np.ones((dimy,dimx))
                xdata1 = np.copy(xdata2)#+1
                #ydata = np.copy(xdata2)
                xdata2 *= datatmp
                xdata1 *= np.array(mda[1].p[0].data)[:dimy,np.newaxis]
                ydata = np.array(mda[2].d[dude.Scan_ToolBox_Y_ComboBox.get_active()].data)
                IMs += [dude.Plot2D_Axe.imshow(ydata.T,cmap=dude.cm, extent=(xdata1[0,0],xdata1[-1,-1],xdata2[0,0],xdata2[-1,-1]), origin="lower")]
                lims = (min(lims[0],min(xdata1[0,0],xdata1[-1,-1])),\
                        max(lims[1],max(xdata1[0,0],xdata1[-1,-1])),\
                        min(lims[2],min(xdata2[0,0],xdata2[-1,-1])),\
                        max(lims[3],max(xdata2[0,0],xdata2[-1,-1])))
                ymin = max(1,min(ymin, ydata[ydata!=0].min()))
                ymax = max(ymax, ydata.max())
                    
        for IM in IMs:
            if dude.Plot2D_Log_ToggleButton.get_active():
                IM.set_norm(colors.LogNorm(ymin, ymax))
            else:
                IM.set_norm(colors.Normalize(ymin, ymax))
        dude.Plot2D_Axe.set_xlim(lims[0],lims[1])
        dude.Plot2D_Axe.set_ylim(lims[3],lims[2])
        dude.Plot2D_Axe.set_aspect(1)
        dude.Plot2D_Canvas.draw()


    def AngleRoI(self, widget, dude):

        if not dude.AllImagesLoaded:
            dude.LoadAllImages(None)
        ymax = int(dude.Scan_ToolBox_CustomROI_YMax_Spin_Adjustment.get_value())
        ymin = int(dude.Scan_ToolBox_CustomROI_YMin_Spin_Adjustment.get_value())
        xmax = int(dude.Scan_ToolBox_CustomROI_XMax_Spin_Adjustment.get_value())
        xmin = int(dude.Scan_ToolBox_CustomROI_XMin_Spin_Adjustment.get_value())
        if not dude.sparse_enabled:
            comx = (((dude.image[:,ymin:ymax+1,xmin:xmax+1].sum(1)*np.arange(xmin,xmax+1)).sum(1)) /\
                    (dude.image[:,ymin:ymax+1,xmin:xmax+1].sum(1).sum(1))).reshape(dude.Plot2D_ydata.shape)
            comy = (((dude.image[:,ymin:ymax+1,xmin:xmax+1].sum(2)*np.arange(ymin,ymax+1)).sum(1)) /\
                    (dude.image[:,ymin:ymax+1,xmin:xmax+1].sum(1).sum(1))).reshape(dude.Plot2D_ydata.shape)
        else:
            m0 = np.arange(ymin,ymax+1)[:,np.newaxis]*np.ones(xmax-xmin+1)
            m0 = sparse.csr_matrix(m0.flatten()*np.ones((dude.image.shape[0],1)))
            m1 = np.zeros((1062,1028))
            m1[ymin:ymax+1,xmin:xmax+1] = 1
            m1 = m1.flatten().nonzero()[0]
            comy = ((m0.multiply(dude.image[:,m1])).sum(1) / dude.image[:,m1].sum(1)).reshape(dude.Plot2D_ydata.shape)

            m0 = np.arange(xmin,xmax+1)*np.ones((ymax-ymin+1,1))
            m0 = sparse.csr_matrix(m0.flatten()*np.ones((dude.image.shape[0],1)))
            comx = ((m0.multiply(dude.image[:,m1])).sum(1) / dude.image[:,m1].sum(1)).reshape(dude.Plot2D_ydata.shape)
        comx -= comx.mean()
        comy -= comy.mean()
        hsv = np.ones((dude.Plot2D_ydata.shape[0], dude.Plot2D_ydata.shape[1], 3))
        hsv[:,:,0] = (np.arctan2(comy, comx)/np.pi+1.)/2.
        print (hsv[:,:,0].min(),hsv[:,:,0].max() )
        rgb = colors.hsv_to_rgb(hsv)
        dude.Plot2D_Image.set_array(rgb)
        dude.Plot2D_Canvas.draw()
"""
