# style.py module
from tkinter import ttk
from misc import all_children

# consider making a class style constructor for each widget class used in the main app

#class MainHeaderStyle:
   # def __init__(self):
               
class LayoutManager: 
    def __init__(self, widget, children = None): # widget will probably be mostly root window and toplevel windows
        
        self.widget = widget
        
        if children:
            self.children_widgets = all_children(self.widget)
        else:
            self.children_widgets = self.widget
    
    def apply_default_layout(self):
        
        for child in self.children_widgets:        
            self.generic_padding(child)

    def generic_padding(self, widget):
        widget.grid(ipadx = 2, ipady = 2, padx = 5, pady = 5, sticky = 'NESW')

#class RootLayout:
    
#class TopLevelLayout:
    
        
    
def window_style():
    s = ttk.Style()#('winnative', 'clam', 'alt', 'default', 'classic', 'vista', 'xpnative')
    s.theme_use("winnative") #will not work on any widgets that are not from the ttk library
            
def check_size(wid):
    print('checking window size')
    # make sure the window isn't larger than the size of the screen
    wid.update_idletasks()
    max_width = int(wid.winfo_screenwidth()*0.9)
    max_height = int(wid.winfo_screenheight()*0.9)

    width = int(wid.winfo_reqwidth())
    height = int(wid.winfo_reqheight())

    if (wid.winfo_reqwidth() > max_width) or (wid.winfo_reqheight() > max_height):

        print('window size ({}, {}) too large'.format(wid.winfo_reqwidth(), wid.winfo_reqheight()))    

        if (wid.winfo_reqwidth() > max_width) and (wid.winfo_reqheight() > max_height):
            new_width = max_width
            new_height = max_height

        elif wid.winfo_reqwidth() > max_width:
            new_width = max_width
            new_height = height

        elif wid.winfo_reqheight() > max_height:
            new_width = width
            new_height = height

        x = int((wid.winfo_screenwidth() // 2) - (new_width // 2))
        y = int((wid.winfo_screenheight() // 2) - (new_height // 2))
        wid.geometry('{}x{}+{}+{}'.format(new_width, new_height, x, y))  

        print('window resized..., new size ({}, {})'.format(wid.winfo_reqwidth(), wid.winfo_reqheight()))     

#class HeaderStyle: