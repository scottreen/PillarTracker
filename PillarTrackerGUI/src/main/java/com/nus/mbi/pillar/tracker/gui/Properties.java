package com.nus.mbi.pillar.tracker.gui;

/**
 *
 * @author xiaochun
 */
public class Properties {
    public double spacing;
    public double grid_oblique;
    public double grid_angle = 90;
    public double diameter;
    public double gauss_sigma;
    public double pixel_size = Double.NaN;
    public boolean dark_pillar;
    
    public static double convertNanoToPixel(double nm, double pixel_size){
        return nm/pixel_size;
    }
    
    public static double convertPixelToNano(double pixel, double pixel_size){
        return pixel*pixel_size;
    }
    
    public double convertNanoToPixel(double nm){
        return nm/pixel_size;
    }
    
    public double convertPixelToNano(double pixel){
        return pixel*pixel_size;
    }
}
