package com.nus.mbi.pillar.detection;

import ij.ImagePlus;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.process.FloatPolygon;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import java.awt.AWTEvent;
import java.awt.Label;

/**
 *
 * @author xiaochun
 */
public class Threshold_Multi_Points implements DialogListener {
    public FloatPolygon points = null;
    private boolean previewing;
    public ImagePlus imp;
    private Label messageArea;
    private double threshold_high;
    private double threshold_low;
    private double maximum;
    private double minimum;
    
    @Override
    public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {        
        if(previewing){            
            threshold_low=gd.getNextNumber();               
            threshold_high=gd.getNextNumber();
            
//            if(threshold_low<minimum || threshold_high>maximum){                
//                return false;
//            }
            FloatPolygon new_points = new FloatPolygon();   
            int npoints = points.npoints;            
            ImageProcessor img = imp.getProcessor();
            for (int i=0; i<npoints; i++) {                
                double fx = points.xpoints[i];
                double fy = points.ypoints[i];
                int x = (int)Math.round(fx);
                int y = (int)Math.round(fy);
                float p = img.getf(x, y);
                if(p>=threshold_low && p<=threshold_high) 
                    new_points.addPoint(fx, fy);                
            }        
                 
            if(imp!=null && imp.isVisible()){
                imp.setRoi(new PointRoi(new_points));
            }
            
            messageArea.setText("the number of points: " + new_points.npoints);
        }

       return (!gd.invalidNumber());
    }
    
    public boolean showDialog(ImagePlus ip) {        
            imp = ip;
            Roi roi = ip.getRoi();
            if(roi==null) return false;
            if(PointRoi.class.isInstance(roi)) points = roi.getFloatPolygon();            
            else return false;
            
            ImageStatistics stat = ip.getProcessor().getStatistics();
            maximum = stat.max;
            minimum = stat.min;//ip.getProcessor().getMin();
            GenericDialog gd = new GenericDialog("Threshold for Multi-Point Selections");	
            
            gd.addNumericField("threshold_low:", minimum, 2, 10, ">"+minimum);	      
            gd.addNumericField("threshold_high:", maximum, 2, 10, "<"+maximum);	        
              
            
            gd.addMessage("the number of points: " + points.npoints);    
            messageArea = (Label)gd.getMessage();                    
            gd.addDialogListener(this);
            previewing = true;

            gd.showDialog();
            if (gd.wasCanceled()){
                ip.setRoi(roi);
                return false;
            }
            previewing = false;
            
            threshold_low=gd.getNextNumber();      
            threshold_high=gd.getNextNumber();
                     

            return true;
    }
}
