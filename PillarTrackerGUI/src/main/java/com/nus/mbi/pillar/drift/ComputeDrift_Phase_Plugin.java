package com.nus.mbi.pillar.drift;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.Roi;
import ij.measure.CurveFitter;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import java.awt.AWTEvent;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author xiaochun
 */
public class ComputeDrift_Phase_Plugin implements PlugIn, DialogListener{
    private ImagePlus imp;
    
    private int window_radius = 10;
    private int threshold_max_drift = 5; //in pixel
    //private int MAX_N;    
    
    private boolean previewing;
    
    public double convert2phase(double d, int N){
        return d*Math.PI*2/(N*Math.sqrt(2));
    }
    
    public double convert2space(double phase, int N){
        return phase*Math.sqrt(2)*N/Math.PI/2;
    }

    @Override
    public void run(String arg) {
        imp = IJ.getImage();
        if(imp.getWidth()!=imp.getHeight()){
            IJ.showMessage("require phase difference image in Frequency Domain");
            return;
        }
        if(!showDialog(imp)) return;
        
        //double threshold_phase = convert2phase(threshold_max_drift, MAX_N);
        //double[] dxy = process(imp.getProcessor(), threshold_phase);
        //IJ.log("dx=" + convert2space(dxy[0], MAX_N) + " dy=" + convert2space(dxy[1], MAX_N));
        
        process_stack(imp);
        
    }
    
    public void process_stack(ImagePlus ip){
        ImageStack stack = ip.getImageStack();
        int nslice = stack.getSize();
        int img_w = ip.getWidth();
        int img_h = ip.getHeight();
        int xc = img_w/2;
        int yc = img_h/2;
        
        int lowx = xc - window_radius;
        int lowy = yc - window_radius;
        int highx = xc + window_radius;
        int highy = yc + window_radius;
        
        if(lowx<0) lowx=0;
        if(lowy<0) lowy=0;
        if(highx>=img_w) highx = img_w-1;
        if(highy>=img_h) highy = img_h-1;
        
        int w = highx-lowx+1;
        int h = highy-lowy+1;
        Rectangle rect = new Rectangle(lowx, lowy, w, h);
        
        double[] dx = new double[nslice];
        double[] dy = new double[nslice];
        double[] xaxis = new double[nslice];     
        //IJ.log("dx dy");
        for(int s=0; s<nslice; s++){            
            ImageProcessor img = stack.getProcessor(s+1);  
            double[] dxy = process(img, rect);            
            dx[s] = convert2space(dxy[0], img_w);
            dy[s] = convert2space(dxy[1], img_w);
            xaxis[s] = s;
            //IJ.log(""+dx[s]+" "+dy[s]);
        }
               
        Plot ploter_driftXY = DriftAnalysisPlotter.get_plot_driftXY(xaxis, dx, dy, "drift from phase difference");
        ploter_driftXY.setXYLabels("frame", "displacement(pixel) refer to 1st frame");
        ploter_driftXY.setLegend("x-drift\ny-drift\ndrift",Plot.AUTO_POSITION|Plot.LEGEND_TRANSPARENT);
        ploter_driftXY.show();	
    }
    
    private double[] plane_slover_mean(List<MyPoint3> data){
        int N = data.size();
        double sumX, sumY, sumZ;
        sumX = sumY = sumZ = 0;
        for(int i=0; i<N; i++){
            MyPoint3 p = data.get(i);
            sumX += p.x;
            sumY += p.y;
            sumZ += p.z;            
        }
        
        double avgX = sumX/N;
        double avgY = sumY/N;
        double avgZ = sumZ/N;
        
        double sumX2, sumXY, sumY2, sumXZ, sumYZ;
        sumX2 = sumXY = sumY2 = sumXZ = sumYZ = 0;
        
        for(int i=0; i<N; i++){
            MyPoint3 p = data.get(i);
            double x = p.x - avgX;
            double y = p.y - avgY;
            double z = p.z - avgZ;
            sumX2 += x*x;
            sumY2 += y*y;
            sumXY += x*y;
            sumXZ += x*z;
            sumYZ += y*z;
        }
        double det = sumX2*sumY2-sumXY*sumXY;
        double dx = (sumY2*sumXZ-sumXY*sumYZ)/det;
        double dy = (sumX2*sumYZ-sumXY*sumXZ)/det;
        return new double[]{dx, dy};
    }
    
//    //simplest way to solve the plane, not good as average version. plane_slover_mean
//    public double[] plane_slover(List<MyPoint3> data){
//        int N = data.size();
//        
//        double sumX2, sumXY, sumY2, sumXZ, sumYZ;
//        sumX2 = sumXY = sumY2 = sumXZ = sumYZ = 0;
//        
//        for(int i=0; i<N; i++){
//            MyPoint3 p = data.get(i);
//            double x = p.x;
//            double y = p.y;
//            double z = p.z;
//            sumX2 += x*x;
//            sumY2 += y*y;
//            sumXY += x*y;
//            sumXZ += x*z;
//            sumYZ += y*z;
//        }
//        double det = sumX2*sumY2-sumXY*sumXY;
//        double dx = (sumY2*sumXZ-sumXY*sumYZ)/det;
//        double dy = (sumX2*sumYZ-sumXY*sumXZ)/det;
//        return new double[]{dx, dy};
//    }
    
    public double[] process(ImageProcessor img, Rectangle rect){
        int lowx = rect.x;
        int lowy = rect.y;
        int w = rect.width;
        int h = rect.height;
        int highx = w+lowx-1;
        int highy = h+lowy-1;
        
        List<MyPoint3> data = new ArrayList();
        for(int x=lowx; x<=highx; x++){
            for(int y=lowy; y<=highy; y++){
                double z = img.getf(x, y);
                data.add(new MyPoint3(x,y,z));
                //if(Math.abs(z)<threshold_phase){
                //    data.add(new MyPoint3(x,y,z));
                //}
            }
        }
        
        return plane_slover_mean(data);
    }
    
//    public double[] process(ImageProcessor img, Rectangle rect, double threshold_phase){
//        int lowx = rect.x;
//        int lowy = rect.y;
//        int w = rect.width;
//        int h = rect.height;
//        int highx = w+lowx-1;
//        int highy = h+lowy-1;
//        
////        List<Double> dataX = new ArrayList();
////        List<Double> dataY = new ArrayList();
////        List<Double> dataZ = new ArrayList();
//        double[] projectX = new double[w];
//        double[] projectY = new double[h];
//        
//        int i=0;
//        for(int x=lowx; x<=highx; x++){
//            double sumY = 0;
//            for(int y=lowy; y<=highy; y++){
//                double z = img.getf(x, y);
//                if(Math.abs(z)<threshold_phase){
//                    sumY += z;
//                }
//            }
//            projectX[i]=sumY;
//            i++;
//        }       
//        
//        int j=0;
//        for(int y=lowy; y<=highy; y++){
//            double sumX = 0;
//            for(int x=lowx; x<=highx; x++){
//                double z = img.getf(x, y);
//                if(Math.abs(z)<threshold_phase){
//                    sumX += z;
//                }
//            }
//            projectY[j]=sumX;
//            j++;
//        }
//        
//        double dx = slope(projectX);
//        double dy = slope(projectY);
//        double[] dxy = new double[]{dx, dy};
//        return dxy;
//    }
    
    private double slope(double[] data){
        int n = data.length;
        double[] x = new double[n];
        for(int i=0; i<n; i++) x[i] = i;
        CurveFitter cf = new CurveFitter(x, data);
        cf.doFit(CurveFitter.STRAIGHT_LINE);
        double[] paras = cf.getParams();
        return paras[1];
    }
    
    public boolean showDialog(ImagePlus ip){
        //imp = ip;        
        //MAX_N = ip.getWidth();
        GenericDialog gd = new GenericDialog("Settings for Drift Computation");            
        int radius = window_radius;
        drawRoi(radius);   
        //int max_drift = threshold_max_drift;
        gd.addNumericField("    center mask radius:", radius, 0, 10, ""); 
        //gd.addNumericField("    maximum drift:", max_drift, 0, 10, "");
        //gd.addPreviewCheckbox(null, "Preview point selection");
        gd.addDialogListener(this);
        previewing = true;
        gd.showDialog();
        if (gd.wasCanceled()) return false;        
        
        radius = (int)gd.getNextNumber();    
        if(radius<0){
            IJ.log("the center radius must be larger than zero!");
            return false;
        }
        
        //max_drift = (int)gd.getNextNumber();    
        //if(max_drift<0){
        //    IJ.log("the maximum drift must be larger than zero!");
         //   return false;
        //}
        
        previewing = false;
        window_radius = radius;
        //threshold_max_drift = max_drift;
        
        ip.updateAndDraw();
        return true;
    }
    
    private void drawRoi(int radius){
        int xc = imp.getWidth()/2;
        int yc = imp.getHeight()/2;
        int w = radius*2+1;
        Roi oval = new Roi(xc-radius, yc-radius, w, w);
        imp.setRoi(oval);
    }
    
     /** Read the parameters (during preview or after showing the dialog) */
    public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
        if(previewing){
            int radius = (int)gd.getNextNumber();            
            if(radius<0) return false;        
            //int max_drift = (int)gd.getNextNumber();    
            //if(max_drift<0) return false;
            
            drawRoi(radius);            
        }
        
       return (!gd.invalidNumber());
    } // public boolean DialogItemChanged
}
