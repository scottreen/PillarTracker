
package com.nus.mbi.pillar.detection;

import com.nus.mbi.pillar.tracker.pillar_detector_match_filter;
import ij.IJ;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.process.ImageProcessor;
import java.awt.AWTEvent;

/**
 *
 * @author xiaochun
 */
public class LocalMaximaThreashold implements DialogListener{
    //private ImagePlus imp;
    private ImageProcessor image;
    
    public double threshold;
    private boolean previewing;
    private pillar_detector_match_filter detector;

    public LocalMaximaThreashold(ImageProcessor image, double threshold, pillar_detector_match_filter detector) {
        //this.imp = imp;
        this.image = image;
        this.threshold = threshold;
        this.detector = detector;
    }    
          
    public boolean showDialog(){
        GenericDialog gd = new GenericDialog("Settings for Threshold");            
        
        gd.addNumericField("threshold:",  threshold, 0, 10, "");        
        gd.addDialogListener(this);
        previewing = true;
        gd.showDialog();
        if (gd.wasCanceled()) return false;        
        
        threshold = gd.getNextNumber();    
        if(threshold<0){
            IJ.log("the threshold must be larger than zero!");
            return false;
        }        
        previewing = false;        
        draw_mask(threshold);        
        return true;
    }
    
    private void draw_mask(double t){        
        detector.drawLocalMaxima(image, t);
        //imp.updateAndDraw();
    }
    
    /** Read the parameters (during preview or after showing the dialog) */
    public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
        if(previewing){
            double t = gd.getNextNumber();                        
            draw_mask(t);
        }
        
       return (!gd.invalidNumber());
    }  
}
