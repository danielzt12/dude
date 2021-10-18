#!/usr/local/bin python3
"""addon for user shortcuts"""

from gi.repository import Gtk
import h5py

class MainWindow:

    def Open_MetadataFile(self, widget, dude):

        self.h5 = h5py.File(dude.FileDialog_Construction(None, 2), "r")
        self.update_keystore()

    def update_keystore(self):

        list_keys = []

        keys_lv1 = self.h5.keys()
        for key1 in keys_lv1:
            item1 = self.h5.get(key1)
            keys_lv2 = item1.keys()
            for key2 in keys_lv2:
                #print(key2)
                item2 = item1.get(key2)
                if isinstance(item2, h5py.Dataset):
                    if not len(item2.shape):
                        list_keys += [key1]
                        break
                else:
                    keys_lv3 = item2.keys()
                    for key3 in keys_lv3:
                        #print(key3)
                        item3 = item2.get(key3)
                        if isinstance(item3, h5py.Dataset):
                            if not len(item3.shape):
                                list_keys += ["/".join([key1, key2])]
                                break
                        else:
                            keys_lv4 = item3.keys()
                            for key4 in keys_lv4:
                                #print(key4)
                                item4 = item3.get(key4)
                                if isinstance(item4, h5py.Dataset):
                                    if not len(item4.shape):
                                        list_keys += ["/".join([key1, key2, key3])]
                                        break
        tmp_i = self.Metadata_Key_ComboBox.get_active()
        self.Metadata_Key_ComboBox.handler_block(self.Key_Changed_Handler)
        self.Metadata_Key_ComboBox.get_model().clear()
        for key in list_keys:
            self.Metadata_Key_ComboBox.get_model().append([key.split("/")[-1], key])
        self.Metadata_Key_ComboBox.handler_unblock(self.Key_Changed_Handler)
        if tmp_i < len(self.Metadata_Key_ComboBox.get_model()):
            self.Metadata_Key_ComboBox.set_active(tmp_i)
        else:
            self.Metadata_Key_ComboBox.set_active(0)

    def update(self, widget):
                                    
        self.Metadata_ListStore.clear()
        key = widget.get_model()[widget.get_active()][1]
        for item in self.h5.get(key).keys():
            if isinstance(self.h5.get(key).get(item), h5py.Dataset):
                if not len(self.h5.get(key).get(item).shape):
                    self.Metadata_ListStore.append([item, "{0:.3f}".format(self.h5.get(key).get(item)[()])])

    def MainWindow_Destroy(self, widget, dude): 

        #dude.ShowMetadata = None
        delattr(dude, "ShowMetadata")
        self.win.destroy()

    def __init__(self, dude):

        self.Metadata_keystore = Gtk.ListStore(str, str)
        renderer_text = Gtk.CellRendererText()    
        self.Metadata_Key_ComboBox = Gtk.ComboBox.new_with_model(self.Metadata_keystore)
        self.Metadata_Key_ComboBox.pack_start(renderer_text, True)
        self.Metadata_Key_ComboBox.add_attribute(renderer_text, "text", 0)

        Metadata_ScrolledWindow = Gtk.ScrolledWindow()
        Metadata_ScrolledWindow.set_policy(Gtk.PolicyType.ALWAYS, Gtk.PolicyType.ALWAYS)
        Metadata_ScrolledWindow.set_overlay_scrolling(False)
        self.Metadata_ListStore = Gtk.ListStore(str, str)
        self.Metadata_TreeView_Filter = self.Metadata_ListStore.filter_new()
        self.Metadata_TreeView = Gtk.TreeView.new_with_model(self.Metadata_TreeView_Filter)
        self.Metadata_TreeView.get_selection().set_mode(Gtk.SelectionMode.NONE)

        Metadata_CellRendererText = Gtk.CellRendererText()     
        #Metadata_CellRendererText.set_property("xalign", 1)
        Metadata_TreeViewColumn = Gtk.TreeViewColumn("key", Metadata_CellRendererText, text = 0)
        self.Metadata_TreeView.append_column(Metadata_TreeViewColumn)

        Metadata_CellRendererText = Gtk.CellRendererText()     
        #Metadata_CellRendererText.set_property("xalign", 1)
        Metadata_TreeViewColumn = Gtk.TreeViewColumn("value", Metadata_CellRendererText, text = 1)
        self.Metadata_TreeView.append_column(Metadata_TreeViewColumn)

        self.Key_Changed_Handler = self.Metadata_Key_ComboBox.connect("changed", self.update)

        self.Metadata_TreeView.set_headers_visible(False)
        self.Metadata_TreeView.set_enable_search(False)
        Metadata_ScrolledWindow.add(self.Metadata_TreeView)
        Metadata_ScrolledWindow.set_size_request(300, 400)

        Metadata_FileOpen_Button = Gtk.Button(" Open ")
        Metadata_FileOpen_Button.connect("clicked", self.Open_MetadataFile, dude)
        
        Main_HBox = Gtk.HBox(homogeneous = False, spacing = 3)
        Main_HBox.set_border_width(3)
        Main_HBox.pack_start(Metadata_FileOpen_Button, False, False, 0)
        Main_HBox.pack_start(self.Metadata_Key_ComboBox, True, True, 0)

        Main_VBox = Gtk.VBox(homogeneous = False, spacing = 3)
        Main_VBox.set_border_width(3)
        Main_VBox.pack_start(Main_HBox, False, False, 0)
        Main_VBox.pack_start(Metadata_ScrolledWindow, False, False, 0)

        self.win = Gtk.Window()
        self.win.connect("destroy", self.MainWindow_Destroy, dude)
        self.win.add(Main_VBox)
        self.win.set_title("Show Metadata")
        self.win.show_all()


def main():
    
    Gtk.main()
    return 0

if __name__ == "__main__":

    MyMainWindow()
    main()
