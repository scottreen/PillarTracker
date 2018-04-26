package com.nus.mbi.pillar.tracker;

import ij.ImagePlus;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.process.FloatPolygon;
import ij.process.ImageProcessor;
import java.awt.AWTEvent;
import java.awt.Color;
import java.awt.Label;

/**
 *
 * @author xiaochun
 */
public class SetMaskFilterRadius_Plugin implements DialogListener {
    private boolean previewing = false;
    private Label messageArea;                  // reference to the textmessage area for displaying the status of mask filter
    public ImagePlus imp;    
    public FloatPolygon points;
    public int low_radius = 5;
    public int larger_mask_radius = 0;
    public int step_radius = 1;
    public boolean apply_mask_filter = false;
    
        private void draw_mask(int max_radius, int step, Color color_off_center){            
            //draw_mask(imp, points, off_center_radius, color_off_center);
            int xc = imp.getWidth()/2;
            int yc = imp.getHeight()/2;            
            Overlay fft_overlay = new Overlay();
            for(int r=low_radius; r<=max_radius; r += step){
                Overlay overlay = get_mask_off_center(points, xc, yc, r, color_off_center);
                for(int i=0; i<overlay.size(); i++) fft_overlay.add(overlay.get(i));
            }
            imp.setOverlay(fft_overlay);
            imp.updateAndDraw();
        }
    
        public static void draw_mask(ImagePlus imp, FloatPolygon points, int off_center_radius, Color color_off_center){

            int xc = imp.getWidth()/2;
            int yc = imp.getHeight()/2;

            Overlay fft_overlay = get_mask_off_center(points, xc, yc, off_center_radius, color_off_center);

            imp.setOverlay(fft_overlay);
            imp.updateAndDraw();
        }
        
        public static void draw_mask(ImageProcessor imp, FloatPolygon points, int off_center_radius, Color color_off_center){

            int xc = imp.getWidth()/2;
            int yc = imp.getHeight()/2;

            Overlay fft_overlay = get_mask_off_center(points, xc, yc, off_center_radius, color_off_center);
            imp.drawOverlay(fft_overlay);
            //imp.setOverlay(fft_overlay);
        }
        
        public static Overlay get_mask_off_center(FloatPolygon points, int xc, int yc, int off_center_radius, Color color_off_center){
            Overlay fft_overlay = new Overlay();
            int w = off_center_radius*2+1;
            if(off_center_radius>0 && points!=null){
                for (int i=0; i<points.npoints; i++) {
                    int x = Math.round(points.xpoints[i]);
                    int y = Math.round(points.ypoints[i]);
                    int dx = xc - x;
                    int dy = yc - y;
                    if(dx*dx+dy*dy>9){
                        OvalRoi oval = new OvalRoi(x-off_center_radius, y-off_center_radius, w, w);
                        //oval.setFillColor(color_off_center);
                        oval.setStrokeColor(color_off_center);
                        fft_overlay.add(oval);
                    }
                }
            }
            return fft_overlay;
        }

    /** Read the parameters (during preview or after showing the dialog) */
        public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
            if(previewing){            
                int radius = (int)gd.getNextNumber(); 
                int step = (int)gd.getNextNumber();  
                if(radius<0 || step<0){                
                    return false;
                }      
                if(imp!=null && imp.isVisible()){
                    draw_mask(radius, step, Color.WHITE);
                }

                //String msg = radius<=low_radius? " no mask will be applied, use raw image" : " will apply mask filter instead of raw image";
                messageArea.setText(message(low_radius, radius, step));
            }

           return (!gd.invalidNumber());
        }
    
        private String message(int start, int end, int step){
            int times = 0;
            for(int i=start; i<=end; i+=step) times++;
            String msg = ""+times + " iterations for mask radius will be applied";
            return msg;
        }
        
    public boolean showFFTDialog() {        
            GenericDialog gd = new GenericDialog("Mask Setting for Pillar Localization");	
            int ref_mask_radius = low_radius;
            gd.addMessage("the mask radius currently used for pillar reconstruction: " + ref_mask_radius);    
            gd.addMessage("the number of masks: " + points.npoints);    
            //gd.addCheckbox("process mask filtered image instead of raw image?", false);
            gd.addNumericField("Maximum_Radius:", larger_mask_radius, 0, 10, ">"+ref_mask_radius);	        
            gd.addNumericField("Radius_Step:", step_radius, 0, 10, ">0");	        
            
            gd.addMessage(message(low_radius, larger_mask_radius, step_radius));
            messageArea = (Label)gd.getMessage();        
            //gd.addCheckbox("apply mask filter with maximum radius?", apply_mask_filter);
            gd.addDialogListener(this);
            previewing = true;

            gd.showDialog();
            if (gd.wasCanceled()) return false;
            previewing = false;

            larger_mask_radius=(int)gd.getNextNumber();
            step_radius=(int)gd.getNextNumber();
            //apply_mask_filter = gd.getNextBoolean();        
            
            return true;
    }

}
