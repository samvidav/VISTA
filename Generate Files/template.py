#template.py module

import tkinter as tk
from tkinter import ttk
import tkinter.filedialog
from tkinter import messagebox
from pandas import concat, read_csv

from style import window_style, LayoutManager, check_size
from dragdrop import SourceDragManager, TargetDragManager
import misc

class Template(tk.Toplevel):

    def __init__(self, master, conditions, abundance=None, localization=None, n=False, weightcutoff=0.5):
               
        self.master = master
        self.n = n
        self.abundance = abundance
        self.localization = localization
        self.conditions = conditions
        self.weightcutoff = weightcutoff

        master.iconify()
        
        tk.Toplevel.__init__(self, master)
        self.populate_template()
        
        window_style()
        LM = LayoutManager(self, children = 1)
        LM.apply_default_layout()
    
    def populate_template(self):
        
        self.files = tk.StringVar()
        self.files.set('Choose file(s) to import')
        
        #frames
        self.frame = ttk.LabelFrame(self, text = 'Instructions', labelanchor = tk.N)
        
        # label
        self.instructions = ttk.Label(self.frame, text = '1. Click "Open" to choose the .csv or .txt file(s) with edge data you want to import.\n'
                                    '2. Click "Expand All" to view the imported column headers.\n'
                                    '3. Drag and drop the appropriate columns into the correct headers in the window to the right.\n'
                                    '4. Click "Continue" once done to check that the information has been imported correctly.\n'
                                    '\n'
                                    'Available shortcuts and mouse combos for drag-and-drop: \n'
                                    'Shift-LeftClick will select consecutive table items \n'
                                    'Control-LeftClick will add items one by one to the selection \n'
                                    'Backspace or delete will clear items from the tables \n'
                                    'LeftClick and drag will attempt to copy a single item to the data table \n'
                                    'RightClick or MiddleClick and drag will copy a group of items to consecutive cells on the other table \n'
                                    'Items can only be copied unidirectionally from the import table to the export table \n'
                                    'If you accidentally delete an item from the data table, re-copy it from the import table \n'
                                    'If you accidentally delete an item from the import table, re-import the file and only the missing node(s) will be re-added to the table \n',
                                     wraplength = 1000, justify = tk.LEFT)
        
        # Entries
        self.file_entry = ttk.Entry(self, textvariable = self.files)
        
        # Treeview tables
        
        # source file tree
        self.file_tree = TreeviewMaker(self)
        
        self.file_tree.config(columns = ['Column Number', 'Header'])
        self.file_tree.column('#0', width = 150)
        self.file_tree.heading('#0', text = 'click to sort', command = self.file_tree.sort_AZ)
        self.file_tree.column('Column Number', anchor = tk.CENTER)
        
        for col in self.file_tree['columns']:
            self.file_tree.heading(col, text = col, command = lambda col=col: self.file_tree.hide_show_cols(col))
            self.file_tree.column(col, width = 150)

        self.file_tree.df_dict = {}
        self.file_tree.moveable_nodes = []
        
        # target file tree
        self.target_tree = TreeviewMaker(self)
        self.target_tree.config(columns = ['Source', 'File Name', 'Header', 'Column Number'])        
        self.target_tree.populate_target_tree(self.conditions, self.abundance, self.localization, self.n)
        
        for col in self.target_tree['columns']:
            self.target_tree.heading(col, text = col, command = lambda col=col: self.target_tree.hide_show_cols(col))
            self.target_tree.column(col, width = 150)
        
        self.target_tree['displaycolumns'] = ['File Name', 'Header']
        self.target_tree.expand_all(self.target_tree.get_children()[0])
        self.target_tree.tag_configure('curritem', background = 'sky blue')
        
        # Buttons
        self.open_button = ttk.Button(self, text = 'Open', command = lambda: self.file_tree.populate_file_tree(self.files))
        
        self.file_tree_expand = ttk.Button(self, text = 'Expand All', command = self.file_tree.expand_all)
        self.file_tree_collapse = ttk.Button(self, text = 'Collapse All', command = self.file_tree.collapse_all)
        
        self.target_tree_expand = ttk.Button(self, text = 'Expand All', command = self.target_tree.expand_all)
        self.target_tree_collapse = ttk.Button(self, text = 'Collapse All', command = self.target_tree.collapse_all)
        
        self.GO_button = ttk.Button(self, text = 'Continue', command = lambda: self.ConfirmInput())
        
        # layout
        self.frame.grid(row=0, column=0, columnspan = 4, sticky = 'NSEW')
        self.instructions.pack(side = tk.LEFT)
        
        self.file_entry.grid(row = 1, column = 0, sticky = 'NSEW')
        self.open_button.grid(row = 1, column = 1, sticky = 'NSEW')
        self.file_tree.grid(row = 3, column = 0, columnspan = 2, rowspan =2, sticky = 'NSEW')
        
        self.GO_button.grid(row = 1, column = 2, columnspan =2, sticky ='NSEW')
        
        self.file_tree_expand.grid(row = 2, column = 0, sticky = 'NSEW')
        self.file_tree_collapse.grid(row = 2, column = 1, sticky = 'NSEW')
        self.target_tree_expand.grid(row = 2, column = 2, sticky = 'NSEW')
        self.target_tree_collapse.grid(row = 2, column = 3, sticky = 'NSEW')
        
        self.target_tree.grid(row=3, column=2, rowspan = 2, columnspan=2, sticky = 'NSEW')

        # make the window and widgets scalable
        for row in range(0, 4):
            for col in range(0, 4):
                self.rowconfigure(row, weight=1)
                self.columnconfigure(col, weight=1)
        
        self.frame.rowconfigure(0, weight=1)
        for col in range(0, 4):
            self.frame.columnconfigure(col, weight=1)

        # make sure window size is not larger than the screen
        check_size(self)
        
        # configure dnd
        sourceDND = SourceDragManager()
        sourceDND.make_dragable(self.file_tree)
        
        targetDND = TargetDragManager()
        targetDND.make_dragable(self.target_tree)

    ############################################################################# 
    def ConfirmInput(self):
        '''Purpose is to show the user how the program is going to input the data and make sure that it looks correct'''
        
        self.top = tk.Toplevel(self.master)
        self.top.iconify()
        
        if self.n:
            self.nodes = self.get_tree_data()

            # some sort of check here to make sure that all the nodes are represented in the nodes file
            if self.nodes is not None:
                # visualize nodes data in GUI   
                self.tree_from_df(self.nodes)
                self.populate_confirm()
                self.top.deiconify()
            else:
                self.top.destroy()

        else:
            self.edges = self.get_tree_data()

            if self.edges is not None:
                # visualize edges data in GUI
                self.tree_from_df(self.edges)
                self.populate_confirm()
                self.top.deiconify()
            
            else:
                self.top.destroy()
    
    def get_nodes(self, conditions):

        # clear existing widgets
        for wid in misc.all_children(self):
            wid.destroy()

        # variables
        localization = tk.BooleanVar() 

        # labels
        label2 = ttk.Label(self, text = 'Do you have localization data for your bait and preys?\nIf not, we will assign them using UniProt localization annotations', justify = tk.LEFT, wraplength=500)

        # RadioButtons
        radiobutton3 = tk.Radiobutton(self, text = 'Yes' , variable = localization, value = True, indicatoron=False)
        radiobutton4 = tk.Radiobutton(self, text = 'No' , variable = localization, value = False, indicatoron=False)
        
        # buttons
        GO = ttk.Button(self, text = "Submit", command = lambda: self.repopulate_template(localization.get()))
        
        # layout
        label2.grid(row=1, column=0, sticky = tk.W, padx=5, pady=5)

        radiobutton3.grid(row=1, column=1,ipadx=5, ipady=5, sticky = tk.N+tk.E+tk.S+tk.W)
        radiobutton4.grid(row=1, column=2,ipadx=5, ipady=5, sticky = tk.N+tk.E+tk.S+tk.W)
        
        GO.grid(row=2, column=1, columnspan=2, sticky=tk.N+tk.E+tk.S+tk.W)
    
    def repopulate_template(self, loc):
        self.localization = loc
        self.populate_template()

    def populate_confirm(self):

        # Buttons
        continue_button = ttk.Button(self.top, text = 'Continue', command = lambda: self.cont(self.top))
        go_back_button = ttk.Button(self.top, text = 'Return to input template', command = self.destroy)

        # Labels
        l1 = ttk.Label(self.top, text = 'This is the data that the program will analyze, does it look correct?\n If so, click continue. \n Otherwise, click Return to input template to return to the previous screen.')
        
        #scrollbar
        scroll = ttk.Scrollbar(self.top, orient = tk.VERTICAL)
        self.tree.config(yscrollcommand = scroll.set)
        scroll.config(command = self.tree.yview)

        x_scroll = ttk.Scrollbar(self.top, orient = tk.HORIZONTAL)
        self.tree.config(xscrollcommand = x_scroll.set)
        x_scroll.config(command = self.tree.xview)

        # layout
        l1.grid(row = 0, column = 0, columnspan = 2)

        continue_button.grid(row=1, column = 1)
        go_back_button.grid(row = 1, column = 0)

        self.tree.grid(row = 2, column = 0, columnspan = 2)
        scroll.grid(row = 2, column = 2)
        x_scroll.grid(row=3, column=0, columnspan=2, sticky=tk.N+tk.E+tk.S+tk.W)

        # styling
        window_style()
        LM = LayoutManager(self.top, children = True)
        LM.apply_default_layout() 

        # make the window and widgets scalable
        for row in range(0, 3):
            for col in range(0, 2):
                self.top.rowconfigure(row, weight=1)
                self.top.columnconfigure(col, weight=1)
        
        # make sure window size is not larger than the screen
        check_size(self.top)
        
    
    def get_tree_data(self):
            
        dfs = []
        ERR = None
        ERR_node = []
        
        for C in self.target_tree.get_children():
            dfs_ = []

            for node in self.target_tree.get_children(C):
                x = self.target_tree.item(node)
                vals = x['values']
                if vals:
                    source = vals[0]
                    col = vals[3]-1
                    series = self.file_tree.df_dict[source].iloc[:, col]
                    series.name = node
                    dfs_.append(series)
                    
                else:
                    ERR = 1
                    ERR_node.append(node)

            if len(dfs_)>0:
                df = concat(dfs_, axis=1)

                if self.n:
                    df = df.set_index(['Uniprot Accession #'])
                
                else:
                    df = df.set_index(['Bait Uniprot Accession #_{}'.format(C), 'Prey Uniprot Accession #_{}'.format(C)])
                    df.index.names = ['Bait Uniprot Accession #', 'Prey Uniprot Accession #']
                
                dfs.append(df)
        
        if ERR:
            messagebox.showerror("Error", "Please enter data for {}".format(str(ERR_node)))
            return None

        elif len(dfs)>0:
            return concat(dfs, axis=1)
        
    def tree_from_df(self, df):

        df = df.reset_index()
        tree = TreeviewMaker(self.top)
        tree.config(columns = list(df.columns.values))
        tree.column('#0', width=0)

        for col in df.columns.values:
            tree.heading(col, text = col)
            tree.column(col, width = 200)

        for row in range(0, df.shape[0]):
            iid = tree.insert('', 'end')
            for i, col in enumerate(df.columns.values):
                tree.set(iid, col, value = df.iloc[row, i])

        self.tree = tree                 
        
    def cont(self, top):
        if self.n:
            self.nodes.to_csv('nodes.csv')
            self.master.destroy()
        else:
            self.edges.to_csv('edges_{}.csv'.format(self.weightcutoff))
            self.n = True
            self.get_nodes(self.target_tree.get_children())
            
            top.destroy()
        
    
class TreeviewMaker(ttk.Treeview): ## call the method from outside the class structure
    '''Creates Treeview widgets for the data insertion template window
    Has methods: populate_target_tree \n get_files \n populate_file_tree \n tree_to_df \n
    expand_all \n collapse_all \n all_nodes \n hide_show_cols \n sort_AZ'''
    
    def __init__(self, master):
        ttk.Treeview.__init__(self, master, height = 10)
        self.master = master
        
    def populate_target_tree(self, conditions, abundance, localization, n):
        '''Creates an empty treeview defined by the user parameters provided for conditions, abundance (if available), localization (if available), and if proteins are from HCMV. 
        Also defines the attribute Treeview.clearable_nodes to be used to define what data in the tree can be deleted by the user'''
        
        self.clearable_nodes = []
        
        if not n:
            l = ['Bait Uniprot Accession #', 'Prey Uniprot Accession #', 'Edge weight']
                
            if abundance:
                l = l + ['Prey abundance']
            
            for C in conditions:
            
                self.insert('', 'end', text = C, iid = C, open=True)

                for text in l:
                    self.insert(C, 'end', text = text, iid = text+'_'+C)
                    self.clearable_nodes.append(text+'_'+C)
        else:
            iid = self.insert('', 'end', open=True)
            l = ['Uniprot Accession #','Gene Name', 'Species Id\n(e.g. 9606 for human)']
            if localization:
                l = l + ['Localization']
            
            for text in l:
                self.insert(iid, 'end', text = text, iid = text)
                self.clearable_nodes.append(text)

    def get_files(self, stringvar):
        '''opens a tkinter filedialog window for the user to select one or more files to import. Sets the Entry textvariable
        to the string of filenames'''
        
        self.filenames =  tk.filedialog.askopenfilenames(parent = self.master, initialdir = "/",title = "Select file",filetypes = (("CSV files","*.csv"),("Text files","*.txt"),("all files","*.*")))
        stringvar.set(list(self.filenames))
        
    def populate_file_tree(self, stringvar):
        
        '''populates the Treeivew with headers for each CSV column generated from conversion to a pandas dataframe by
        the read_CSV method'''
        
        self.get_files(stringvar)
        files = self.filenames
        
        for file in files:
            
            df = read_csv(file, index_col = None, header = 0, sep=None, engine = 'python')
            self.df_dict[file] = df
                    
            if not self.exists(file):
                self.insert('', 'end', text = file.split('/')[-1], iid=file)

                for i in range(0, df.shape[1]):
                    myiid = '{}_col{}'.format(file, i)
                    
                    self.moveable_nodes.append(myiid)
                    
                    self.insert(file, i, iid = myiid)
                    self.set(myiid, 'Column Number', value = i+1)
                    self.set(myiid, 'Header', value = df.columns.values[i]) 
            else:
                 for i in range(0, df.shape[1]):
                    myiid = '{}_col{}'.format(file, i)
                                        
                    if not self.exists(myiid):
                        self.moveable_nodes.append(myiid)
                        self.insert(file, i, iid = myiid)
                        self.set(myiid, 'Column Number', value = i+1)
                        self.set(myiid, 'Header', value = df.columns.values[i])
        
    def expand_all(self, node = None):
        '''expands all nodes by recursively expanding each child of the treeview'''
        if node == None:
            node = ''
            
        iids = self.all_nodes(node)
        for iid in iids:
            self.item(iid, open = True)
        
    def collapse_all(self, node = None):
        '''collapses all nodes by recursively collapsing each child of the treeview'''
        if node == None:
            node = ''
            
        for node in self.get_children(node):
            self.item(node, open = False)
            
    def all_nodes(self, node):
        '''recursively generates a list of all children of a treeview or of a particular node'''
        
        node_list = [item for item in self.get_children(node)]

        for item in node_list :
            try:
                node_list.extend(self.get_children(item))
            except AttributeError:
                pass
        
        return node_list
    
    def hide_show_cols(self, col):
        '''hides or shows a column by adjusting its width''' #modify to also expand the other column?
        
        width = self.column(col)['width']
        
        if width > 100:
            self.column(col, width = 20)
        else:
            self.column(col, width = 150)
    
    def sort_AZ(self):
        '''sorts source node children by filename'''
        iids = self.get_children('')
        files = {}
        for iid in iids:
            files[iid] = iid.split('/')[-1]
        
        sort_files = sorted(files, key = files.get)
        
        for i, iid in enumerate(sort_files):
            self.move(iid, '', i)


    
        
        
        
   
        
        