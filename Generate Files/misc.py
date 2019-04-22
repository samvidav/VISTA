# misc.py module

def all_children(wid):
    wid_list = wid.winfo_children()

    for item in wid_list :
        if item.winfo_children():
            wid_list.extend(item.winfo_children())
        
    return wid_list