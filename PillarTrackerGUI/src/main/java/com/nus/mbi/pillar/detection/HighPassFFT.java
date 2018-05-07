package com.nus.mbi.pillar.detection;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Undo;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.measure.Measurements;
import ij.plugin.PlugIn;
import ij.process.ColorProcessor;
import ij.process.FHT;
import ij.process.FloatPolygon;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import java.awt.AWTEvent;
import java.awt.Color;
import java.util.ArrayList;
import java.util.List;
import com.nus.mbi.pillar.stat.IntPoint;

/**
 *  rewrite the FFT class in ij.plugin.FFT to implement the Fourier mask filter.
 * @author xiaochun
 */
public class HighPassFFT implements  PlugIn, Measurements, DialogListener  {
    private ImagePlus imp;
    private boolean padded;
    private int originalWidth;
    private int originalHeight;

    //private Overlay FFT_Overlay = null;
    public FloatPolygon points = null;
    public int off_center_radius = 5;
    public int center_radius = 5;
    private boolean previewing = false;
    //private boolean white_mask_center = true;
    //private boolean white_mask_off_center = true;
    private boolean white_mask = true;
    private boolean process_all_frames = true;
    
    public int step_radius_off_center = 1;
    public int end_radius_off_center = 0;
    
    private boolean use_phase_correlation = false;
    private FHT first_frame_fht = null;
    private boolean show_phase_difference = false;

    public HighPassFFT() {
    }

    public void setup(ImagePlus imp) {
        this.imp = imp;
    }
    
    public int getOffCenterRadius(){
        return off_center_radius;
    }
    
    public void setOffCenterRadius(int off_center_radius){
        this.off_center_radius = off_center_radius;
    }
    
    public void setOffCenterPoints(int[] xpoints, int[] ypoints){
        int npoints = xpoints.length;
        float[] px = new float[npoints];
        float[] py = new float[npoints];        
        for (int i=0; i<npoints; i++) {                
            int x = xpoints[i];
            int y = ypoints[i];
            px[i] = x;
            py[i] = y;
        }        
        points = new FloatPolygon(px, py);        
    }
    
    public boolean hasOffCenterPoints(){
        return (points!=null && points.npoints>0 && off_center_radius>0);
    }
    
    
    public void run(String arg) {        
        imp = IJ.getImage();
        ImageProcessor ip = imp.getProcessor();
        Object obj = imp.getProperty("FHT");
        FHT fht = (obj instanceof FHT)?(FHT)obj:null;
        boolean inverse;
        if (fht!=null) {
            inverse = true;            
        } else {            
            fht = newFHT(ip);
            inverse = false;
        }
        if (inverse){
            ImageProcessor ip2 = doInverseTransform(fht, imp.getProcessor());
            String title = imp.getTitle();
            if (title.startsWith("FFT of "))
                title = title.substring(7, title.length());
            ImagePlus imp2 = new ImagePlus("Inverse FFT of "+title, ip2);
            imp2.setCalibration(imp.getCalibration());
            imp2.show();    
        }
        else {            
            ImagePlus fft_ip = doForwardTransform(fht);   
            fft_ip.show();
        }    
        IJ.showProgress(1.0);
    }
    
    public void process(){
        imp = IJ.getImage();
        if(isInFrequencyDomain(imp)){
            IJ.showMessage("The spatial domain image is required");
            return;
        }
        
        if(imp.getStackSize()<2) process_all_frames = false;
        MaskFilterDialogResults rst = showProcessDialog(imp);
        if(!rst.dialog_return) return;         
        
        if(!process_all_frames){
            if(use_phase_correlation){
                ImageProcessor first_img = imp.getStack().getProcessor(1);
                first_frame_fht = newFHT(first_img);        
                first_frame_fht.setShowProgress(false);
                ImagePlus fft_ip = doForwardTransform(first_frame_fht);  
                int xc = first_img.getWidth()/2;
                int yc = first_img.getHeight()/2;
                first_frame_fht = phase_shift(first_frame_fht, -xc, -yc);
            }
            boolean has_off_center_mask = hasOffCenterPoints();//FFT_Overlay.size()>0;
            //if(has_off_center_mask) process_current_frame(off_center_radius, end_radius_off_center, step_radius_off_center);
            if(has_off_center_mask) process_current_frame(rst.radius_start, rst.radius_end, rst.radius_step);            
            else process_current_frame();
            return;
        }
        
        mask_filter_stack(use_phase_correlation);
    }
    
    public boolean isEmptyMask(){
        if(off_center_radius<1 && center_radius<1) return true;        
        return (!hasOffCenterPoints() && center_radius<1);
    }
    
    public static boolean isInFrequencyDomain(ImagePlus imp){
        //ImagePlus imp = IJ.getImage();        
        Object obj = imp.getProperty("FHT");
        FHT fht = (obj instanceof FHT)?(FHT)obj:null;
        return fht!=null;
    }
    
    private void process_current_frame(int start_radius, int end_radius, int step_radius){        
        if(step_radius<1 || start_radius<1 || end_radius<start_radius){
            process_current_frame();
            return;
        }
        
        ImageProcessor img = imp.getProcessor();     
        FHT fht = newFHT(img);        
        fht.setShowProgress(false);
        ImagePlus fft_ip = doForwardTransform(fht); 
       
        int w = imp.getWidth();
        int h = imp.getHeight();
        ImageStack stack = new ImageStack(w, h);
        ImageStack phase_stack = new ImageStack(w, h);
        ImageStack phase_diff_stack = new ImageStack(fft_ip.getWidth(), fft_ip.getHeight());
        
        //int num_img = end_radius-start_radius+1;        
        for(int i=start_radius; i<=end_radius; i += step_radius){            
            FHT fft_img = getMaskedFHT(fht, i);
            if(this.use_phase_correlation){
                FHT first_fft_img = getMaskedFHT(first_frame_fht, i);
                //ImageProcessor phase_cc_img = this.phase_correlation(first_fft_img, fft_img);
                FHT phase_cc_fht = fft_img.conjugateMultiply(first_fft_img);
                phase_cc_fht.originalWidth = w;
                phase_cc_fht.originalHeight = h;
                phase_cc_fht.setShowProgress(false);
                //FHT norm_phase_fht = phase_cc_fht.conjugateMultiply(phase_cc_fht);
                //phase_cc_fht = phase_cc_fht.divide(norm_phase_fht);
                if(show_phase_difference){
                    FHT phase_diff_fht = phase_cc_fht.getCopy();//fft_img.divide(first_fft_img);//phase_cc_fht.getCopy();
                    ImageStack re_im = phase_diff_fht.getComplexTransform();
                    float[] re = (float[])re_im.getProcessor(1).getPixels();
                    float[] im = (float[])re_im.getProcessor(2).getPixels();
                    int img_w = re_im.getWidth();
                    int img_h = re_im.getHeight();
                    int img_size = img_w*img_h;
                    float[] phase = new float[img_size];
                    for(int p=0; p<img_size; p++) phase[p] = (float)Math.atan2(im[p], re[p]);
                    phase_diff_stack.addSlice("mask radius="+i, phase);
                    //re_im.addSlice("Phase", phase);
                    //FloatProcessor phase_img = new FloatProcessor(img_w, img_h, phase);
                    //ImagePlus phase_diff_ip = new ImagePlus("phase difference-" + imp.getShortTitle(), phase_img);
                    //phase_diff_ip.show();
                }
                else{
                    ImageProcessor phase_cc_img = doInverseTransform(phase_cc_fht,false);               
                    phase_stack.addSlice("mask radius="+i, phase_cc_img);
                }
            }
            
            ImageProcessor ifft_img = doInverseTransform(fft_img, false);                    
            stack.addSlice("mask radius="+i, ifft_img);
            
            IJ.showProgress(i, end_radius);
        }        
        ImagePlus new_ip = new ImagePlus("FFT_Filter-" + imp.getShortTitle(), stack);
        new_ip.show();
        
        if(this.use_phase_correlation){
            if(show_phase_difference){
                ImagePlus phase_ip = new ImagePlus("PhaseDiff-" + imp.getShortTitle(), phase_diff_stack);
                phase_ip.show();
            }
            else{
                ImagePlus phase_ip = new ImagePlus("PhaseCC-" + imp.getShortTitle(), phase_stack);
                phase_ip.show();
            }
        }
        
        IJ.showProgress(1.0);
    }
    
    private void process_current_frame(){
        ImageProcessor img = imp.getProcessor();            
        FHT fht = newFHT(img);        
        fht.setShowProgress(false);
        ImagePlus fft_ip = doForwardTransform(fht); 
        //ImageProcessor mask = fft_ip.getProcessor();         
        FHT fft_img = getMaskedFHT(fht, off_center_radius);
        
        if(this.use_phase_correlation){
            FHT first_fft_img = getMaskedFHT(first_frame_fht, off_center_radius);
            
            FHT phase_cc_fht = fft_img.conjugateMultiply(first_fft_img);
            phase_cc_fht.originalWidth = img.getWidth();
            phase_cc_fht.originalHeight = img.getHeight();
            phase_cc_fht.setShowProgress(false);
            //FHT norm_phase_fht = phase_cc_fht.conjugateMultiply(phase_cc_fht);
            //phase_cc_fht = phase_cc_fht.divide(norm_phase_fht);
            if(show_phase_difference){
                FHT phase_diff_fht = phase_cc_fht.getCopy();//fft_img.divide(first_fft_img);//phase_cc_fht.getCopy();
                ImageStack re_im = phase_diff_fht.getComplexTransform();
                float[] re = (float[])re_im.getProcessor(1).getPixels();
                float[] im = (float[])re_im.getProcessor(2).getPixels();
                int img_size = re_im.getWidth()*re_im.getHeight();
                float[] phase = new float[img_size];
                for(int i=0; i<img_size; i++) phase[i] = (float)Math.atan2(im[i], re[i]);
                //re_im.addSlice("Phase", phase);
                FloatProcessor phase_img = new FloatProcessor(re_im.getWidth(), re_im.getHeight(), phase);
                ImagePlus phase_diff_ip = new ImagePlus("phase_diff-" + imp.getShortTitle(), phase_img);
                phase_diff_ip.show();
            }
            else{
                ImageProcessor phase_cc_img = doInverseTransform(phase_cc_fht,false);   
                //ImageProcessor phase_cc_img = this.phase_correlation(first_fft_img, fft_img);
                ImagePlus phase_ip = new ImagePlus("PhaseCC-" + imp.getShortTitle(), phase_cc_img);
                phase_ip.show();
            }
            
        }
        
        ImageProcessor ifft_img = doInverseTransform(fft_img, false);        
        ImagePlus new_ip = new ImagePlus("FFT_Filter-" + imp.getShortTitle(), ifft_img);
        new_ip.show();        
        
        IJ.showProgress(1.0);
    }
    
    public FDResult process_frame(ImageProcessor img, FHT reference_fht, int off_center_radius){        
        FHT fht = forward_transform(img); 
        FHT fht_mask = GetOffCenterMaskedFHT(fht, off_center_radius);    
        //FHT ref_mask = GetOffCenterMaskedFHT(reference_fht, off_center_radius);            
        //ImageProcessor phase_img = phase_correlation(ref_mask, fht_mask);       
        ImageProcessor ifft_mask = inverse_fft_transform(fht_mask);        
        return new FDResult(fht, fht_mask, ifft_mask);
    }
    
    public FDResult process_frame(FHT fht, int off_center_radius){        
        //FHT fht = forward_transform(img); 
        FHT fht_mask = GetOffCenterMaskedFHT(fht, off_center_radius);         
        ImageProcessor ifft_mask = inverse_fft_transform(fht_mask);        
        return new FDResult(fht, fht_mask, ifft_mask);
    }
    
        
    private ImageProcessor process(ImageProcessor img, int off_center_radius){
        FHT fht = newFHT(img);        
        fht.setShowProgress(false);
        ImagePlus fft_ip = doForwardTransform(fht); 
        ImageProcessor mask = fft_ip.getProcessor(); 
        
        FHT fft_img = getMaskedFHT(fht, mask, off_center_radius);
        ImageProcessor ifft_img = doInverseTransform(fft_img, false);       

        return ifft_img;
    }
    
    private FHT getMaskedFHT(FHT fht, int off_center_radius){
        FHT fht1 = fht.getCopy();
        ImageProcessor mask = fht1.getPowerSpectrum(); 
        return getMaskedFHT(fht1, mask, off_center_radius);
    }
    
    private FHT getMaskedFHT(FHT fht, ImageProcessor mask, int off_center_radius){
        
        int img_w = mask.getWidth();
        int img_h = mask.getHeight();
        
        int xc = img_w/2;
        int yc = img_h/2;
        
        Color c_center = white_mask ? Color.WHITE : Color.BLACK;  
        Color c_offcenter = c_center;  
        Overlay off_center_oval = get_mask_off_center(xc, yc, off_center_radius, c_offcenter);

        mask.setColor(c_offcenter);
        for(int i=0; i<off_center_oval.size(); i++) mask.fill(off_center_oval.get(i));  
        
        if(center_radius>0){
            int w = center_radius*2+1;        
            OvalRoi center_oval = new OvalRoi(xc-center_radius, yc-center_radius, w, w);         
            mask.setColor(c_center);  
            mask.fill(center_oval);  
        }
        
        FHT fft_img = fht.getCopy();
        fft_img.setShowProgress(false);                          
        doMasking(fft_img, mask);

        return fft_img;
    }
    
    private FHT GetOffCenterMaskedFHT(FHT fht, int off_center_radius){
        //FHT fht1 = fht.getCopy();
        ImageProcessor mask = fht.getPowerSpectrum(); 
        return GetOffCenterMaskedFHT(fht, mask, off_center_radius);
    }
    
    private FHT GetOffCenterMaskedFHT(FHT fht, ImageProcessor mask, int off_center_radius){
        int img_w = mask.getWidth();
        int img_h = mask.getHeight();
        
        int xc = img_w/2;
        int yc = img_h/2;
        
        Color c_offcenter = Color.WHITE;  
        Overlay off_center_oval = get_mask_off_center(xc, yc, off_center_radius, c_offcenter);
        
        mask.setColor(c_offcenter);
        for(int i=0; i<off_center_oval.size(); i++) mask.fill(off_center_oval.get(i));  
                
        FHT fft_img = fht.getCopy();
        fft_img.setShowProgress(false);                          
        doMasking(fft_img, mask);
        
        return fft_img;
    }
    
//    void phase_correlation_stack(){
//        ImageStack stack = imp.getImageStack();
//        int nslice = stack.getSize();
//        int img_w = imp.getWidth();
//        int img_h = imp.getHeight();
//        ImageStack is = new ImageStack(img_w, img_h);
//        //Color c_offcenter = white_mask_off_center ? Color.WHITE : Color.BLACK;
//        //Color c_center = white_mask_center ? Color.WHITE : Color.BLACK;  
//        Color c_center = white_mask ? Color.WHITE : Color.BLACK;  
//        Color c_offcenter = c_center;
//        
//        int xc = img_w/2;
//        int yc = img_h/2;
//        int w = low_radius*2+1;        
//        OvalRoi center_oval = new OvalRoi(xc-low_radius, yc-low_radius, w, w);
//        FHT first_fft_img = null;
//        //ImageProcessor last_mask = null;
//        for(int s=1; s<=nslice; s++){            
//            ImageProcessor img = stack.getProcessor(s);       
//            FHT fht = newFHT(img);
//            fht.setShowProgress(false);
//            ImagePlus fft_ip = doForwardTransform(fht); 
//            ImageProcessor mask = fft_ip.getProcessor(); 
//            mask.setColor(c_offcenter);            
//            for(int i=0; i<FFT_Overlay.size(); i++) mask.fill(FFT_Overlay.get(i));  
//            mask.setColor(c_center);  
//            mask.fill(center_oval);                             
//            
//            FHT fft_img = fht.getCopy();
//            fft_img.setShowProgress(false); 
//            doMasking(fft_img, mask);
//            if(s==1) first_fft_img = fft_img;            
//            ImageProcessor ifft_phase_cc = phase_correlation(first_fft_img, fft_img);
//            is.addSlice("phase_cc"+(s-1)+"-"+s, ifft_phase_cc);
//            
//            IJ.showStatus("mask filtering:" + s + "/" + nslice);
//            IJ.showProgress(s,nslice);
//        }        
//        ImagePlus new_ip = new ImagePlus("FFT_Filter-" + imp.getShortTitle(), is);
//        new_ip.show();
//    }
    
    private float sqr(float x) {
        return x*x;
    }
    
    void normalize(FHT img){
        float[] fht = (float[])img.getPixels();
        int maxN = img.getWidth();
        for (int row=0; row<maxN; row++) {
            //amplitude(row, maxN, fht, amp);
            int base = row*maxN;
            int l;
            for (int c=0; c<maxN; c++) {
                l = ((maxN-row)%maxN) * maxN + (maxN-c)%maxN;
                double amp = Math.sqrt((sqr(fht[base+c]) + sqr(fht[l]))/2);
                fht[base+c] /= amp;
                fht[l] /= amp;
            }
        }
    }
    
    // translate an image with phase shift in Fourier space
    public static FHT phase_shift(FHT img, double dx, double dy){
        int rowMod, cMod, colMod;
        int maxN = img.getWidth();
        double h2e, h2o;
        float[] h1 = (float[])img.getPixels();
        
        double px = dx*Math.PI*2/maxN;
        double py = dy*Math.PI*2/maxN;
//        double px = dx*Math.PI*2/img.originalWidth;
//        double py = dy*Math.PI*2/img.originalHeight;
        
        FHT shift_fht = img.getCopy();
        float[] tmp = (float[])shift_fht.getPixels();
        int center = maxN/2;
        for (int r =0; r<maxN; r++) {
            rowMod = (maxN - r) % maxN;
            double phasey = (r<=center)? r*py : (-rowMod)*py;
            for (int c=0; c<maxN; c++) {
                colMod = (maxN - c) % maxN;
                double phasex = (c<=center) ? c*px : (-colMod)*px;
                double phase = phasex+phasey;
                h2e = Math.cos(phase);//(h2[r * maxN + c] + h2[rowMod * maxN + colMod]) / 2;
                h2o = Math.sin(phase);//(h2[r * maxN + c] - h2[rowMod * maxN + colMod]) / 2;                
                tmp[r * maxN + c] = (float)(h1[r * maxN + c] * h2e + h1[rowMod * maxN + colMod] * h2o);
            }
        }
        
        shift_fht.originalWidth = img.originalWidth;
        shift_fht.originalHeight = img.originalHeight;
        shift_fht.setShowProgress(false);
        return shift_fht;
    }
    
    // phase correlation between two images
    public ImageProcessor phase_correlation(FHT fht1, FHT fht2) {        
        FHT fht = fht2.conjugateMultiply(fht1);
        fht.originalWidth = fht1.originalWidth;
        fht.originalHeight = fht1.originalHeight;
        fht.setShowProgress(false);
        return doInverseTransform(fht,false);
    }
    
    public ImageProcessor phase_correlation(FHT reference_fht, ImageProcessor img){
        FHT fht = forward_transform(img);
        FHT phase_cc_fht = fht.conjugateMultiply(reference_fht);                     
        phase_cc_fht.originalWidth = img.getWidth();
        phase_cc_fht.originalHeight = img.getHeight();
        phase_cc_fht.setShowProgress(false);
        ImageProcessor phase_cc_img = doInverseTransform(phase_cc_fht,false);                                   
        return phase_cc_img;
    }
    
    void mask_filter_stack(boolean phase_cc){
        ImageStack stack = imp.getImageStack();
        int nslice = stack.getSize();
        int img_w = imp.getWidth();
        int img_h = imp.getHeight();
        ImageStack ifft_stack = new ImageStack(img_w, img_h);
        ImageStack phase_stack = new ImageStack(img_w, img_h);   
        ImageStack phase_diff_stack = null;
        //Color c_offcenter = white_mask_off_center ? Color.WHITE : Color.BLACK;
        //Color c_center = white_mask_center ? Color.WHITE : Color.BLACK;  
        Color c_center = white_mask ? Color.WHITE : Color.BLACK;  
        Color c_offcenter = c_center;    
        Overlay off_center_oval = null;
        FHT first_fft_img = null;
        for(int s=1; s<=nslice; s++){            
            ImageProcessor img = stack.getProcessor(s);       
            FHT fht = newFHT(img);
            fht.setShowProgress(false);
            ImagePlus fft_ip = doForwardTransform(fht); 
            ImageProcessor mask = fft_ip.getProcessor();             
            int xc = mask.getWidth()/2;
            int yc = mask.getWidth()/2;
            
            if(off_center_oval==null) off_center_oval = get_mask_off_center(xc, yc, off_center_radius, c_offcenter);                        
            mask.setColor(c_offcenter);  
            for(int i=0; i<off_center_oval.size(); i++) mask.fill(off_center_oval.get(i));  
            
            if(center_radius>0){                
                int w = center_radius*2+1;        
                OvalRoi center_oval = new OvalRoi(xc-center_radius, yc-center_radius, w, w);
                mask.setColor(c_center);  
                mask.fill(center_oval);
            }
            
            FHT fft_img = fht.getCopy();             
            fft_img.setShowProgress(false); 
            doMasking(fft_img, mask);            

            if(phase_cc){
                if(s==1){
                    first_fft_img = phase_shift(fft_img, -img_w/2, -img_h/2);//fft_img.getCopy();
                    phase_diff_stack = new ImageStack(fft_ip.getWidth(), fft_ip.getHeight());
                }            
                
                FHT phase_cc_fht = fft_img.conjugateMultiply(first_fft_img);                     
                phase_cc_fht.originalWidth = img_w;
                phase_cc_fht.originalHeight = img_h;
                phase_cc_fht.setShowProgress(false);
                //FHT norm_phase_fht = phase_cc_fht.conjugateMultiply(phase_cc_fht);
                //phase_cc_fht = phase_cc_fht.divide(norm_phase_fht);
                //ImageProcessor ifft_phase_cc = phase_correlation(first_fft_img, fft_img);
                if(show_phase_difference){
                    FHT phase_diff_fht = phase_cc_fht.getCopy();//fft_img.divide(first_fft_img);//
                    ImageStack re_im = phase_diff_fht.getComplexTransform();                    
                    float[] re = (float[])re_im.getProcessor(1).getPixels();
                    float[] im = (float[])re_im.getProcessor(2).getPixels();
                    int img_size = re_im.getWidth()*re_im.getHeight();
                    float[] phase = new float[img_size];
                    for(int p=0; p<img_size; p++) phase[p] = (float)Math.atan2(im[p], re[p]);
                    phase_diff_stack.addSlice("phase_diff_1-"+s, phase);
                }
                else{
                    ImageProcessor phase_cc_img = doInverseTransform(phase_cc_fht,false);                                   
                    phase_stack.addSlice("phase_cc_1-"+s, phase_cc_img);
                }
                //phase_stack.addSlice("phase_cc"+(s-1)+"-"+s, ifft_phase_cc);
            }            
                        
            ImageProcessor ifft_img = doInverseTransform(fft_img, false);
            ifft_stack.addSlice(ifft_img); 
            
            IJ.showStatus("mask filtering:" + s + "/" + nslice);
            IJ.showProgress(s,nslice);
        }        
        ImagePlus new_ip = new ImagePlus("FFT_Filter-" + imp.getShortTitle(), ifft_stack);
        new_ip.show();
        
        if(phase_cc){
            if(show_phase_difference){
                ImagePlus phase_ip = new ImagePlus("Phase-" + imp.getShortTitle(), phase_diff_stack);
                phase_ip.show();
            }
            else{
                ImagePlus phase_ip = new ImagePlus("Phase-" + imp.getShortTitle(), phase_stack);
                phase_ip.show();
            }
        }
    }
    
    ImageProcessor doInverseTransform(FHT fht, ImageProcessor mask) {
        fht = fht.getCopy();
        doMasking(fht, mask);
        return doInverseTransform(fht, false);
    }
    
    ImageProcessor doInverseTransform(FHT fht) {
          return doInverseTransform(fht, true);
    }

public static ImageProcessor inverse_fft_transform(FHT fht) {        
        //showStatus("Inverse transform");
        fht.inverseTransform();
        if (fht.quadrantSwapNeeded)
            fht.swapQuadrants();
        //fht.resetMinAndMax();
        ImageProcessor ip2 = fht;
        if (fht.originalWidth>0) {
            fht.setRoi(0, 0, fht.originalWidth, fht.originalHeight);
            ip2 = fht.crop();
        }
        return ip2;
}

ImageProcessor doInverseTransform(FHT fht, boolean convert_bitDepth) {        
        //showStatus("Inverse transform");
        fht.inverseTransform();
        if (fht.quadrantSwapNeeded)
            fht.swapQuadrants();
        fht.resetMinAndMax();
        ImageProcessor ip2 = fht;
        if (fht.originalWidth>0) {
            fht.setRoi(0, 0, fht.originalWidth, fht.originalHeight);
            ip2 = fht.crop();
        }
        
        if(convert_bitDepth){
            int bitDepth = fht.originalBitDepth>0?fht.originalBitDepth:imp.getBitDepth();
            switch (bitDepth) {
                case 8: ip2 = ip2.convertToByte(false); break;
                case 16: ip2 = ip2.convertToShort(false); break;
                case 24:
                    //showStatus("Setting brightness");
                    if (fht.rgb==null || ip2==null) {
                        IJ.error("FFT", "Unable to set brightness");
                        return null;
                    }
                    ColorProcessor rgb = (ColorProcessor)fht.rgb.duplicate();
                    rgb.setBrightness((FloatProcessor)ip2);
                    ip2 = rgb; 
                    fht.rgb = null;
                    break;
                case 32: break;
            }
            if (bitDepth!=24 && fht.originalColorModel!=null)
                ip2.setColorModel(fht.originalColorModel);
        }
        return ip2;
    }    
    
    public FHT forward_transform(ImageProcessor img){
        FHT fht = newFHT(img);
        fht.setShowProgress(false);
        fht.transform();
        return fht;
    }
    
    ImagePlus doForwardTransform(FHT fht) {
        showStatus("Forward transform");
        fht.transform();           
        showStatus("Calculating power spectrum");
        ImageProcessor ps = fht.getPowerSpectrum();
        ImagePlus imp2 = new ImagePlus("FFT of "+imp.getTitle(), ps);
        //imp2.show();
        imp2.setProperty("FHT", fht);
        imp2.setCalibration(imp.getCalibration());
        String properties = "Fast Hartley Transform\n";
        properties += "width: "+fht.originalWidth + "\n";
        properties += "height: "+fht.originalHeight + "\n";
        properties += "bitdepth: "+fht.originalBitDepth + "\n";
        imp2.setProperty("Info", properties);
        return imp2;
    }
    
    FHT newFHT(ImageProcessor ip) {
        FHT fht;
        if (ip instanceof ColorProcessor) {
            showStatus("Extracting brightness");
            ImageProcessor ip2 = ((ColorProcessor)ip).getBrightness();
            fht = new FHT(pad(ip2));
            fht.rgb = (ColorProcessor)ip.duplicate(); // save so we can later update the brightness
        } else
            fht = new FHT(pad(ip));
        if (padded) {
            fht.originalWidth = originalWidth;
            fht.originalHeight = originalHeight;
        }
        int bitDepth = ip.getBitDepth();
        fht.originalBitDepth = bitDepth;
        if (bitDepth!=24)
        	fht.originalColorModel = ip.getColorModel();
        return fht;
    }
    
    ImageProcessor pad(ImageProcessor ip) {
        originalWidth = ip.getWidth();
        originalHeight = ip.getHeight();
        int maxN = Math.max(originalWidth, originalHeight);
        int i = 2;
        while(i<maxN) i *= 2;
        if (i==maxN && originalWidth==originalHeight) {
            padded = false;
            return ip;
        }
        maxN = i;
        showStatus("Padding to "+ maxN + "x" + maxN);
        ImageStatistics stats = ImageStatistics.getStatistics(ip, MEAN, null);
        ImageProcessor ip2 = ip.createProcessor(maxN, maxN);
        ip2.setValue(stats.mean);
        ip2.fill();
        ip2.insert(ip, 0, 0);
        padded = true;
        Undo.reset();
        //new ImagePlus("padded", ip2.duplicate()).show();
        return ip2;
    }
    
    void showStatus(String msg) {
//        if (stackSize>1)
//            IJ.showStatus("FFT: " + slice+"/"+stackSize);
//        else
//            IJ.showStatus(msg);
        //IJ.showStatus(msg);
    }
    
    void doMasking(FHT ip, ImageProcessor mask) {
//        if (stackSize>1)
//            return;
        float[] fht = (float[])ip.getPixels();
        //ImageProcessor mask = ip.getCopy();//imp.getProcessor();
        mask = mask.convertToByte(false);        
        ImageStatistics stats = ImageStatistics.getStatistics(mask, MIN_MAX, null);
        if (stats.histogram[0]==0 && stats.histogram[255]==0)
            return;
        boolean passMode = stats.histogram[255]!=0;
        //IJ.showStatus("Masking: "+(passMode?"pass":"filter"));
        mask = mask.duplicate();
        if (passMode)
            changeValuesAndSymmetrize(mask, (byte)255, (byte)0); //0-254 become 0
        else
            changeValuesAndSymmetrize(mask, (byte)0, (byte)255); //1-255 become 255
        //long t0=System.currentTimeMillis();
        for (int i=0; i<3; i++)
            smooth(mask);
        //IJ.log("smoothing time:"+(System.currentTimeMillis()-t0));
        if (IJ.debugMode || IJ.altKeyDown())
        	new ImagePlus("mask", mask.duplicate()).show();
        ip.swapQuadrants(mask);
        byte[] maskPixels = (byte[])mask.getPixels();
        for (int i=0; i<fht.length; i++) {
            fht[i] = (float)(fht[i]*(maskPixels[i]&255)/255.0);
        }
        //FloatProcessor fht2 = new FloatProcessor(mask.getWidth(),mask.getHeight(),fht,null);
        //new ImagePlus("fht", fht2.duplicate()).show();
    }

    // Change pixels not equal to v1 to the new value v2.
    // For pixels equal to v1, also the symmetry-equivalent pixel is set to v1
    // Requires a quadratic 8-bit image.
    void changeValuesAndSymmetrize(ImageProcessor ip, byte v1, byte v2) {
        byte[] pixels = (byte[])ip.getPixels();
        int n = ip.getWidth();
        for (int i=0; i<pixels.length; i++) {
            if (pixels[i] == v1) {  //pixel has been edited for pass or filter, set symmetry-equivalent
                if (i%n==0) {       //left edge
                    if (i>0) pixels[n*n-i] = v1;
                } else if (i<n)     //top edge
                    pixels[n-i] = v1;
                else                //no edge
                    pixels[n*(n+1)-i] = v1;
            } else
                pixels[i] = v2;     //reset all other pixel values
        }
    }

    // Smooth an 8-bit square image with periodic boundary conditions
    // by averaging over 3x3 pixels
    // Requires a quadratic 8-bit image.
    static void smooth(ImageProcessor ip) {
        byte[] pixels = (byte[])ip.getPixels();
        byte[] pixels2 = (byte[])pixels.clone();
        int n = ip.getWidth();
        int[] iMinus = new int[n];  //table of previous index modulo n
        int[] iPlus = new int[n];   //table of next index modulo n
        for (int i=0; i<n; i++) {   //creating the tables in advance is faster calculating each time
            iMinus[i] = (i-1+n)%n;
            iPlus[i] = (i+1)%n;
        }
        for (int y=0; y<n; y++) {
            int offset1 = n*iMinus[y];
            int offset2 = n*y;
            int offset3 = n*iPlus[y];
            for (int x=0; x<n; x++) {
                int sum = (pixels2[offset1+iMinus[x]]&255)
                        + (pixels2[offset1+x]&255)
                        + (pixels2[offset1+iPlus[x]]&255)
                        + (pixels2[offset2+iMinus[x]]&255)
                        + (pixels2[offset2+x]&255)
                        + (pixels2[offset2+iPlus[x]]&255)
                        + (pixels2[offset3+iMinus[x]]&255)
                        + (pixels2[offset3+x]&255)
                        + (pixels2[offset3+iPlus[x]]&255);
                pixels[offset2 + x] = (byte)((sum+4)/9);
            }
        }
    }

    void redisplayPowerSpectrum() {
        FHT fht = (FHT)imp.getProperty("FHT");
        if (fht==null)
            {IJ.error("FFT", "Frequency domain image required"); return;}
        ImageProcessor ps = fht.getPowerSpectrum();
        imp.setProcessor(null, ps);
    }
    
    public boolean showHighPassDialog(ImagePlus ip){
        imp = ip;        
        Roi roi = ip.getRoi();
        if(roi!=null){        
            if(PointRoi.class.isInstance(roi)) points = roi.getFloatPolygon();            
        }        
//        if(points == null || points.npoints<1){            
//            IJ.showMessage("There is no point selections");
//            return false;
//        }
        
        redisplayPowerSpectrum();             
        boolean do_ifft = false;        
        GenericDialog gd = new GenericDialog("Settings for Mask Filter");            
        int radius = off_center_radius;
        int center_radius = this.center_radius;
        //Color c_offcenter = white_mask_off_center ? Color.WHITE : Color.BLACK;
        //Color c_center = white_mask_center ? Color.WHITE : Color.BLACK;        
        Color c_center = white_mask ? Color.WHITE : Color.BLACK;
        Color c_offcenter = c_center;        
        draw_mask(radius, center_radius, c_offcenter, c_center);   
        
        gd.addNumericField("off-center mask radius:",  radius, 0, 10, "");
        gd.addNumericField("    center mask radius:", center_radius, 0, 10, "");
        gd.addCheckbox("white mask?", white_mask);
        //gd.addCheckbox("Fill White for off-center mask?", white_mask_off_center);
        //gd.addCheckbox("Fill White for center mask?", white_mask_center);
        gd.addCheckbox("do inverse FFT?", do_ifft);
        //gd.addPreviewCheckbox(null, "Preview point selection");
        gd.addDialogListener(this);
        previewing = true;
        gd.showDialog();
        if (gd.wasCanceled()) return false;        
        
        radius = (int)gd.getNextNumber();    
        center_radius = (int)gd.getNextNumber();   
        //boolean white_offcen = gd.getNextBoolean();
        //boolean white_cen = gd.getNextBoolean();
        boolean is_white = gd.getNextBoolean();
        if(radius<0){
            IJ.log("the off-center radius must be larger than zero!");
            return false;
        }
        
        previewing = false;
        //if (!dialogItemChanged(gd, null))   //read parameters
        //    return true;
        off_center_radius = radius;
        this.center_radius = center_radius;
        do_ifft = gd.getNextBoolean();
        white_mask = is_white;
        //white_mask_off_center = white_offcen;
        //white_mask_center = white_cen;
        //FFT_Overlay = new Overlay();
        int xc = ip.getWidth()/2;
        int yc = ip.getHeight()/2;
        ImageProcessor img = ip.getProcessor(); 
        //c_offcenter = white_mask_off_center ? Color.WHITE : Color.BLACK;
        c_offcenter = white_mask ? Color.WHITE : Color.BLACK;
        img.setColor(c_offcenter);
        int w = radius*2+1;
        
        if(off_center_radius>0 && points!=null){
            for (int i=0; i<points.npoints; i++) {
                int x = Math.round(points.xpoints[i]);
                int y = Math.round(points.ypoints[i]);
                int dx = xc -x;
                int dy = yc- y;
                if(dx*dx+dy*dy>9){
                    OvalRoi oval = new OvalRoi(x-radius, y-radius, w, w);
                    img.fill(oval);   
                    //FFT_Overlay.add(oval);
                }
            }
        }
        
        if(this.center_radius>0){
            //c_center = white_mask_center ? Color.WHITE : Color.BLACK;
            c_center = white_mask ? Color.WHITE : Color.BLACK;
            img.setColor(c_center);
            w = this.center_radius*2+1;
            OvalRoi oval = new OvalRoi(xc-this.center_radius, yc-this.center_radius, w, w);
            //FFT_Overlay.add(oval);
            img.fill(oval);   
        }
        
        end_radius_off_center = off_center_radius;
        
        ip.updateAndDraw();
        return do_ifft;
    }
    
    public List<IntPoint> getOffCenterPoints(){
        List<IntPoint> p = new ArrayList();
        if(points==null) return p;        
        for (int i=0; i<points.npoints; i++) {
            int x = Math.round(points.xpoints[i]);
            int y = Math.round(points.ypoints[i]);
            p.add(new IntPoint(x,y));            
        }
        return p;
    }
    
    public Overlay get_mask_off_center(int xc, int yc, int off_center_radius, Color color_off_center){
        return get_mask_off_center(points, xc, yc, off_center_radius, color_off_center);
    }
    
    public static Overlay get_mask_off_center(FloatPolygon points, int xc, int yc, int off_center_radius, Color color_off_center){
        Overlay fft_overlay = new Overlay();

        int w = off_center_radius*2+1;
        //FloatPolygon new_points = new FloatPolygon();
        if(off_center_radius>0 && points!=null){
            for (int i=0; i<points.npoints; i++) {
                int x = Math.round(points.xpoints[i]);
                int y = Math.round(points.ypoints[i]);
                int dx = xc - x;
                int dy = yc - y;
                if(dx*dx+dy*dy>9){
                    OvalRoi oval = new OvalRoi(x-off_center_radius, y-off_center_radius, w, w);
                    oval.setFillColor(color_off_center);
                    //img.fill(oval);   
                    //oval.setStrokeColor(color_off_center);
                    fft_overlay.add(oval);
                    //new_points.addPoint(x, y);
                }
            }
        }
        //points = new_points;
        return fft_overlay;
    }
    
    private void draw_mask(int off_center_radius, int center_radius, Color color_off_center, Color color_center){
        
        int xc = imp.getWidth()/2;
        int yc = imp.getHeight()/2;

        Overlay fft_overlay = get_mask_off_center(xc, yc, off_center_radius, color_off_center);
        if(center_radius>0){
            //c = white_mask_center ? Color.WHITE : Color.BLACK;
            int w = center_radius*2+1;
            OvalRoi oval = new OvalRoi(xc-center_radius, yc-center_radius, w, w);
            oval.setFillColor(color_center);
            fft_overlay.add(oval);
            //img.fill(oval);   
        }
        imp.setOverlay(fft_overlay);
        imp.updateAndDraw();
    }
    
    /** Read the parameters (during preview or after showing the dialog) */
    public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
        if(previewing){
            int radius = (int)gd.getNextNumber();            
            if(radius<0) return false;        
            int center_radius = (int)gd.getNextNumber();  
            //boolean white_offcenter = gd.getNextBoolean();
            //boolean white_center = gd.getNextBoolean();
            //Color c_center = white_center ? Color.WHITE : Color.BLACK;
            //Color c_offcenter = white_offcenter ? Color.WHITE : Color.BLACK;
            boolean is_white = gd.getNextBoolean();
            Color c_center = is_white ? Color.WHITE : Color.BLACK;
            Color c_offcenter = c_center;
            draw_mask(radius, center_radius, c_offcenter, c_center);
        }
        
       return (!gd.invalidNumber());
    } // public boolean DialogItemChanged
    
    public MaskFilterDialogResults showProcessDialog(ImagePlus ip){                
        boolean has_off_center_mask = hasOffCenterPoints();
        boolean has_one_frame = (ip.getStackSize()==1);
        //if(has_one_frame && !has_off_center_mask) return true;
        MaskFilterDialogResults result = new MaskFilterDialogResults();
        result.dialog_return = false;
        
        GenericDialog gd = new GenericDialog("Process with Fourier Mask Filter");        
        gd.addMessage("Current radius of off-center mask: " + off_center_radius);
        if(!has_one_frame) gd.addCheckbox("Process all frames?", false);
          
        //gd.addCheckbox("show phase correlation with the 1st frame?", false);
        if(has_off_center_mask){
            gd.addMessage("----Iterate the radius of off-center mask ONLY if processing current frame");
            gd.addNumericField("    start_radius:", off_center_radius, 0, 10, ">0");            
            gd.addNumericField("    end_radius:", end_radius_off_center, 0, 10, ">0");
            gd.addNumericField("    step_radius:", step_radius_off_center, 0, 10, ">0");
        }
        gd.showDialog();
        if (gd.wasCanceled()) return result;
        
        if(!has_one_frame) process_all_frames = gd.getNextBoolean();                    
        
        //use_phase_correlation = gd.getNextBoolean();
        if(process_all_frames){
            result.dialog_return = true;
            return result;
        }
        
        if(has_off_center_mask){            
            int start_radius=(int)gd.getNextNumber();
            int end_radius=(int)gd.getNextNumber();
            int step_raiuds = (int)gd.getNextNumber();
            if(start_radius<1) {
                IJ.showMessage("the start radius must be larger than 0");
                return result;
            }

            if(end_radius<start_radius) {
                IJ.showMessage("the end radius must be larger than start radius");
                return result;
            }
            
            if(step_raiuds<1) {
                IJ.showMessage("the step radius must be larger than 0");
                return result;
            }
            
            result.setRaiuds(start_radius, end_radius, step_raiuds);
            //step_radius_off_center = start_radius;
            //end_radius_off_center = end_radius;
        }
        
        result.dialog_return = true;
        return result;
    }
}
