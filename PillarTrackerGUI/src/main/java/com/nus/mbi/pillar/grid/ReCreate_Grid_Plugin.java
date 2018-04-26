package com.nus.mbi.pillar.grid;

import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.plugin.PlugIn;

/**
 *
 * @author xiaochun
 */
public class ReCreate_Grid_Plugin implements PlugIn{
    public double threshold = 0;
    public boolean process_all_frames = false;
    public GridProperty grid_property = new GridProperty();
    public String grid_fname;
    
    @Override
    public void run(String arg) {
        
        
        
    }
    
    public boolean showDialog(boolean set_grid){
        double threshold_deflection = threshold;
        GenericDialogPlus gd = new GenericDialogPlus("Settings for Grid Recreation");
        
        if(set_grid && grid_property==null) set_grid = false;
        double s = grid_property.spacing;
        double a = grid_property.oblique;
        double b = grid_property.grid_angle;
        if(set_grid){
            gd.addNumericField("        spacing:", s, 2, 10, "pixel");
            gd.addNumericField("        grid_oblique:", a, 2, 10, "degree");		
            gd.addNumericField("        grid_angle:", b, 2, 10, "degree");	
        }
        
        gd.addNumericField("       deflection threshold:", threshold_deflection, 2, 10, "nm");
        gd.addCheckbox(" process all frames?", false);
        String new_grid_fname = grid_fname;
        if(grid_fname.endsWith(".dxy")) new_grid_fname = grid_fname.replaceAll(".dxy", ".grid");
        
        gd.addFileField("save to:", new_grid_fname, 50);
        gd.showDialog();
        if (gd.wasCanceled()) return false;
        
        if(set_grid){
            s = gd.getNextNumber();
            a = gd.getNextNumber();
            b = gd.getNextNumber();
            if(s<0){
                IJ.showMessage("the spacing must be larger than 0");
                return false;
            }            
        }
        
        threshold_deflection = gd.getNextNumber();
        if(threshold_deflection<=0){
            IJ.showMessage("the threshold must be larger than 0");
            return false;
        }
        
        boolean is_all_frames = gd.getNextBoolean();
        new_grid_fname = gd.getNextString();
        
        if(set_grid){            
            grid_property.spacing = s;
            grid_property.oblique = a;
            grid_property.grid_angle = b;
        }
        
        threshold = threshold_deflection;
        process_all_frames = is_all_frames;
        grid_fname = new_grid_fname;
        return true;
    }
}
