# dragdrop.py
import template

# configure errors for when single value cells are already filled, etc.

class SourceDragManager(): #to be used on source widget
    '''This is basically a binding manager for the file-import Treeview. It defines a quasi "drag-and-drop" protocol
    that allows the user to copy data headers from the file tree to the target tree. It also defines bindings for 
    delete and backspace so that nodes values can be deleted, but does not allow deletion of entire nodes'''
    
    def make_dragable(self, widget):
        '''This is the method that should be called on an instance of SourceDragManager to apply the below bindings to a widget
        ex: dnd = SourceDragManager()
            dnd.make_draggable(widget)'''
        
        widget.bind("<ButtonPress-1>", self.on_start)
        widget.bind("<B1-Motion>", self.on_drag)
        widget.bind("<ButtonRelease-1>", self.on_drop)
        
        widget.bind("<ButtonPress-3>", self.B2_drag_start)
        widget.bind("<B3-Motion>", self.on_drag)
        widget.bind("<ButtonRelease-3>",self.on_drop)
        
        widget.bind("<ButtonPress-2>", self.B2_drag_start)
        widget.bind("<B2-Motion>", self.on_drag)
        widget.bind("<ButtonRelease-2>",self.on_drop)
        
        widget.bind('<Delete>', self.delete)
        widget.bind('Backspace', self.delete)
    
        self.source = widget

    def on_start(self, event):
        # you could use this method to create a floating window
        # that represents what is being dragged.
        self.source.configure(cursor = 'hand1')
        self.last_tags = None
        self.last_row = None
        self.target = None
    
    def B2_drag_start(self, event):
        '''Unbinds the default bindings of Buttons 2 and 3 so that they can be used, instead, to move groups of items'''
        
        #self.source.configure(cursor = 'hand1')
        pass
            
        
    def on_drag(self, event):
        '''Highlights the target node of the target widget, basically binds Button1 + hovering to a function that
        highlights the node that the mouse is over'''
        
        # you could use this method to move a floating window that
        # represents what you're dragging
        
        #self.source.configure(cursor="hand1")
        
        x,y = event.widget.winfo_pointerxy()
        target = event.widget.winfo_containing(x,y)
       
        if not target == self.source:
            if isinstance(target, template.TreeviewMaker):
                target_row = target.identify_row(event.y)
                self.target = target
                
                if self.last_row:
                    if not target_row == self.last_row:
                        target.item(self.last_row, tags = self.last_tags) # reset the old tags
                        self.last_row = target_row
                        self.last_tags = target.item(target_row, 'tags')
                        target.item(target_row, tags=['curritem'])
                else:
                    self.last_row = target_row
                    self.last_tags = target.item(target_row, 'tags')
                    target.item(target_row, tags=['curritem'])
                       
            else:
                if self.target:
                    self.target.item(self.last_row, tags = self.last_tags)
        else:
            if self.target:
                self.target.item(self.last_row, tags = self.last_tags)
      

    def on_drop(self, event):
        '''configures what happends when a user attempts to drag and drop data onto the target tree'''
        
        #self.source.configure(cursor="arrow")
        
        # find the widget under the cursor
        x,y = event.widget.winfo_pointerxy()
        target = event.widget.winfo_containing(x,y)
        
        stagnate = 0
        
        if isinstance(target, template.TreeviewMaker):
            if not target == self.source:
                target_row = target.identify_row(event.y)
                selection = self.source.selection()
                
                #resets the row color
                target.item(self.last_row, tags=self.last_tags)
                
                if target_row in target.clearable_nodes:
                    if len(selection) == 1:
                        if selection[0] in self.source.moveable_nodes: #the simplest case
                            target.set(target_row, 'Source', value = self.source.parent(selection[0]))
                            target.set(target_row, 'Header', value = self.source.item(selection[0])['values'][1])
                            target.set(target_row, 'Column Number', value = self.source.item(selection[0])['values'][0])
                            target.set(target_row, 'File Name', value =self.source.parent(selection[0]).split('/')[-1])
                        
                    else: 
                        for n, item in enumerate(selection):
                            if item in self.source.moveable_nodes:
                                if n == 0:
                                    target.set(target_row, 'Source', value = self.source.parent(item))
                                    target.set(target_row, 'Header', value = self.source.item(item)['values'][1])
                                    target.set(target_row, 'Column Number', value = self.source.item(item)['values'][0])
                                    target.set(target_row, 'File Name', value = self.source.parent(item).split('/')[-1])
                                else:
                                    if not stagnate == 1:
                                        target_row = target.next(target_row)
                                    if target_row in target.clearable_nodes:
                                        stagnate = 0
                                        target.set(target_row, 'Source', value = self.source.parent(item))
                                        target.set(target_row, 'Header', value = self.source.item(item)['values'][1])
                                        target.set(target_row, 'Column Number', value = self.source.item(item)['values'][0])
                                        target.set(target_row, 'File Name', value = self.source.parent(item).split('/')[-1])
                                    
                            else:
                                stagnate = 1
        
        self.last_target_row = None
            
    def delete(self, event): # delete from file_tree
        selection = self.source.selection()
        
        for iid in selection:
            self.source.delete(iid)
            
            if iid in self.source.moveable_nodes:
                self.source.moveable_nodes.remove(iid)
            
class TargetDragManager():
    def make_dragable(self, widget):
        self.source = widget
        
        widget.bind('<Delete>', self.delete)
        widget.bind('Backspace', self.delete)
        
    def delete(self, event): # delete from the data_tree
       
        selection = self.source.selection()
        
        for iid in selection:
            
            if iid in self.source.clearable_nodes:
                self.source.item(iid, values = [])

        self.source.selection_set()  
            
            