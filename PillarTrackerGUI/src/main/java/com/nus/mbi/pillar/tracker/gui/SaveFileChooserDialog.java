package com.nus.mbi.pillar.tracker.gui;

import fiji.util.gui.GenericDialogPlus;

/**
 *
 * @author xiaochun
 */
public class SaveFileChooserDialog {
    public String save_fname = "";
    private String defaultPath = "";

    public SaveFileChooserDialog(String defaultPath) {
        this.defaultPath  = defaultPath;
    }
    
    public boolean showDialog(){
        GenericDialogPlus gd = new GenericDialogPlus("Save a file:");
        gd.addFileField("save file to:", defaultPath, 50);
        gd.showDialog();
        if(gd.wasCanceled()) return false;        
        save_fname = gd.getNextString();
        return true;
    }
            
}
