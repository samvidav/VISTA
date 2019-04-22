# main.py module
import tkinter as tk
from tkinter import font, ttk, messagebox

from style import window_style, LayoutManager
import template

class InitiateProgram(tk.Tk):
    '''Create the main architecture of the main screen of the UI... this is where the user inputs the structure of their data'''
    
    def __init__(self, master):
        self.master = master
               
        master.title('Generate Files for Interactome')
        master.wm_attributes('-topmost', True)
        master.wm_attributes('-topmost', False)
        master.protocol("WM_DELETE_WINDOW", self.on_close)
        self.create_UI()
    

    def create_UI(self):
        '''Create the actual layout of the main-page UI'''
        
        master = self.master     
        
        #variables
        self.master.conditions = tk.StringVar()
        self.master.weightcutoff = tk.DoubleVar()
        self.master.abundance = tk.BooleanVar()

        self.master.conditions.set('e.g. "24, 48, 72"')
        self.master.weightcutoff.set(0.5)
        self.master.abundance.set(False)
        
        # Frames/LabelFrames (Frames with a border around them)
        master.topframe = ttk.LabelFrame(master, text = 'Intructions for Use', labelanchor = tk.N)
        master.bottomframe = ttk.LabelFrame(master, text ='Tell us about your data', labelanchor = tk.N)
        
        # Labels
        master.instructions = ttk.Label(master.topframe, text = 'Use this tool to generate nodes and edge files for exploratory analysis and\nvisualization of spatial and temporal interaction data.', wraplength = 1000 )
        master.label = ttk.Label(master.bottomframe, text = 'Enter your time-points separated by commas.',justify = tk.LEFT)
        master.label0 = ttk.Label(master.bottomframe, text = 'Enter an edge confidence threshold between 0 and 1\n(edges below this value will not be visualized)',justify = tk.LEFT)
        master.label1 = ttk.Label(master.bottomframe, text = 'Do you have abundance data for your prey?', wraplength = 500,justify = tk.LEFT)
        
        # Entries
        master.entry = ttk.Entry(master.bottomframe, textvariable = self.master.conditions, justify = tk.CENTER) # conditions separated by commas
        master.entry0 = ttk.Entry(master.bottomframe, textvariable = self.master.weightcutoff, justify = tk.CENTER)

        # radiobuttons
        radiobutton1 = tk.Radiobutton(master.bottomframe, text = 'Yes' , variable = self.master.abundance, value = True, indicatoron=False)
        radiobutton2 = tk.Radiobutton(master.bottomframe, text = 'No' , variable = self.master.abundance, value = False, indicatoron=False)
        
        # Buttons
        master.gobutton = ttk.Button(master.bottomframe, text = "Go", command = self.go_button)
        
        # layout
        master.topframe.grid(row =0, column = 0, columnspan = 3, ipadx=5, ipady=5, padx=10, pady=10, sticky = tk.N+tk.E+tk.S+tk.W)
        master.bottomframe.grid(row =1, column = 0, columnspan = 3, rowspan = 5, ipadx=5, ipady=5, padx=10, pady=10, sticky = tk.N+tk.E+tk.S+tk.W)
        
        master.instructions.grid(row=0,column=0, columnspan=3, ipadx=5, ipady=5, padx=10, pady=10)
        master.label.grid(row=0, column = 0, sticky = tk.W, padx=5, pady=5)
        master.label0.grid(row=1, column = 0, sticky = tk.W, padx=5, pady=5)
        master.label1.grid(row=2, column=0, sticky = tk.W, padx=5, pady=5)
        
        master.entry.grid(row = 0, column=1, columnspan=2, sticky = tk.N+tk.E+tk.S+tk.W)
        master.entry0.grid(row = 1, column=1, columnspan=2, sticky = tk.N+tk.E+tk.S+tk.W)
        radiobutton1.grid(row=2, column=1,ipadx=5, ipady=5, sticky = tk.N+tk.E+tk.S+tk.W)
        radiobutton2.grid(row=2, column=2,ipadx=5, ipady=5, sticky = tk.N+tk.E+tk.S+tk.W)
        
        master.gobutton.grid(row = 3, column = 1, columnspan=2, sticky=tk.N+tk.E+tk.S+tk.W)
        
        #styling
        window_style()
        LM = LayoutManager(master, children = True)
        LM.apply_default_layout()

        # make the window and widgets scalable
        for row in range(0, 4):
            for col in range(0, 3):
                for wid in [master, master.topframe, master.bottomframe]:
                    wid.rowconfigure(row, weight=1)
                    wid.columnconfigure(col, weight=1)

    def go_button(self):
          
        if self.master.conditions.get() == 'e.g. "24, 48, 72"':
            messagebox.showerror('ERROR', 'Please enter your conditions')
        else:
            conditions = self.master.conditions.get()
            # split by commas and remove any spurious leading or trailing spaces
            conditions = [x.strip() for x in conditions.split(',') if len(x.strip())>0]
            
            if len([x for x in conditions if conditions.count(x) > 1])>1: # force condition names to be unique
                messagebox.showerror('ERROR', 'Please enter unique names for each of your conditions')

            else:
                self.input_template = template.Template(self.master, conditions, abundance = self.master.abundance.get(), weightcutoff=self.master.weightcutoff.get())
            
    def on_close(self):
        
        self.master.destroy()
        print('exit')
    
if __name__ == '__main__':
    
    root = tk.Tk()
    InitiateProgram(root)
    root.mainloop()