package com.nus.mbi.pillar.tracker;

import com.nus.mbi.pillar.detection.HighPassFFT;
import com.nus.mbi.pillar.drift.DriftAnalysisPlotter;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.process.FHT;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.util.ArrayList;
import com.nus.mbi.pillar.stat.GaussianFitter;

/**
 *
 * @author xiaochun
 */
public class Compute_Drift_Phase_Correlation_Plugin implements PlugIn{
    private ImagePlus experiment;
    private int kernel_radius = 4;
    //private int window_length = 9;
    private boolean show_shift_img = false;
    public boolean show_registeration = false;
    public double sigma = 1.5;
    //public boolean dark_pillar = false;
    
    public void setup(ImagePlus ip, int radius){
        experiment = ip;
        if(radius>0) kernel_radius = radius;
    }
    
    private int[] search_maxima(int xc, int yc, double[] pixels, int width, int height, int radius){
        int r = radius;
        double max = Double.MIN_VALUE;
        int xm = -1;
        int ym = -1;
        for(int i=-r; i<=r; i++){
            int y = yc + i;
            for(int j=-r; j<=r; j++){
                int x = xc + j;
                if(x>=0 && x<width && y>=0 && y<height){
                    double v = pixels[x+y*width];
                    if(v>max){
                        max = v;
                        xm = x;
                        ym = y;
                    }
                }
            }
        }
        if(xm<0) return null;
        return new int[]{xm, ym};
    }
    
    public static double[] get_profile(int xc, int yc, double[] pixels, int width, int height, int radius){
        int r = radius;
        int w = 2*r+1;
        if(xc<r || xc>=width-r) return null;
        
        double[] profile = new double[w];
        int yoff = yc*width;
        for(int j=-r; j<=r; j++){
                int x = xc + j;
                profile[j+r] = pixels[x+yoff];        
        }
        return profile;
    }
    
    public void process(){                
        ImageStack stack = experiment.getImageStack();
        Roi work_roi = experiment.getRoi();
        boolean use_roi = (work_roi!=null && work_roi.isArea());
        
        ImageProcessor ref_img = stack.getProcessor(1);
        if(use_roi){
            ref_img.setRoi(work_roi);
            ref_img = ref_img.crop();
        }
        
        int w = ref_img.getWidth();
        int h = ref_img.getHeight();        
        int xc = w/2;
        int yc = h/2;
        int XM = xc;
        int YM = yc;
        
        HighPassFFT fft = new HighPassFFT();
        FHT fht = fft.forward_transform(ref_img);
        FHT reference_fht = HighPassFFT.phase_shift(fht, -xc, -yc);
        
        ImageProcessor phase_img = fft.phase_correlation(reference_fht, fht);//auto correlation    
        ImagePlus phase_ip = new ImagePlus("phase_cc_" + experiment.getShortTitle(), phase_img);
        if(show_shift_img){
            phase_ip.show();
//            ImageProcessor shift_img = HighPassFFT.inverse_fft_transform(reference_fht); 
//            ImagePlus shift_ip = new ImagePlus("phase_cc_" + experiment.getShortTitle(), shift_img);
//            shift_ip.show();
//            return;
        }
        
        float[] pixels = (float[])phase_img.getPixels();        
        int img_size = w*h;
        double[] pixels_double = new double[img_size];
        for(int k=0; k<img_size; k++) pixels_double[k] = pixels[k];
        
        double[] dataY = get_profile(XM,YM,pixels_double,w, h, kernel_radius);
        double[] gaussian_para = GaussianFitter.do_fit(dataY);
        double sigma = gaussian_para[3];
        IJ.log("estmated sigma:" + sigma);
        
        double[] seed_cx, seed_cy;
        seed_cx = new double[]{xc};
        seed_cy = new double[]{yc};
        int kw = 2*kernel_radius + 1;
//        double[][] estimated_xy;
//        double[][] estimated_xy = LocalizationFunction.ML_localization(1, pixels_double, w, h, seed_cx, seed_cy, kw, sigma, sigma, 2, false);  
//        if(estimated_xy==null){
//            IJ.showMessage("can not find a point to start");
//            return;            
//        }
//        double[] mlxy = estimated_xy[0];
//        if(mlxy==null){
//            IJ.showMessage("can not find a point to start");
//            return;    
//        }
        
        int nframes = stack.getSize();
        ImageStack stack_reg = new ImageStack(stack.getWidth(),stack.getHeight());         
        //ImagePlus stablized_ip = new ImagePlus("stablized_" + experiment.getShortTitle(), ref_img);        
        if(show_registeration){
            //stablized_ip.show();
            FloatProcessor first_img = stack.getProcessor(1).convertToFloatProcessor();
            stack_reg.addSlice(first_img);
            //stack_reg.addSlice("f="+1, ref_img);
        }
        
//        double[] dx = new double[nframes];
//        double[] dy = new double[nframes];
//        double[] xaxis = new double[nframes];     
//        //IJ.log("dx dy");
//        for(int s=0; s<nframes; s++){            
//            dx[s] = 0;
//            dy[s] = 0;
//            xaxis[s] = s;
//            //IJ.log(""+dx[s]+" "+dy[s]);
//        }
//               
//        Plot ploter_driftXY = DriftAnalysisPlotter.get_plot_driftXY(xaxis, dx, dy, "drift from phase correlation");
//        ploter_driftXY.setXYLabels("frame", "displacement(pixel) refer to 1st frame");
//        ploter_driftXY.setLegend("x-drift\ny-drift",Plot.AUTO_POSITION|Plot.LEGEND_TRANSPARENT);
//        ploter_driftXY.show();	
                
        double[] mlx = new double[nframes];
        double[] mly = new double[nframes];
        mlx[0] = XM;//mlxy[0]; 
        mly[0] = YM;//mlxy[1];
        
        boolean suc = true;          
        for(int s=1; s<nframes; s++){            
            ImageProcessor ip_s = stack.getProcessor(s+1);   
            if(use_roi){
                ip_s.setRoi(work_roi);
                ip_s = ip_s.crop();
            }
            //phase_img = fft.phase_correlation(reference_fht, ip);
            fht = fft.forward_transform(ip_s);
            FHT phase_cc_fht = fht.conjugateMultiply(reference_fht);                     
            phase_cc_fht.originalWidth = w;
            phase_cc_fht.originalHeight = h;
            phase_cc_fht.setShowProgress(false);
            phase_img = HighPassFFT.inverse_fft_transform(phase_cc_fht);
            if(show_shift_img){
                phase_ip.setProcessor(phase_img);
                phase_ip.setRoi(xc-kernel_radius, yc-kernel_radius, kw, kw);
                //phase_ip.updateAndRepaintWindow();
            }       
            
            pixels = (float[])phase_img.getPixels();
            for(int k=0; k<img_size; k++) pixels_double[k] = pixels[k];
            
            int[] xy = search_maxima(xc, yc, pixels_double, w, h, kernel_radius);
            if(xy==null){
                suc = false;
                break;
            }
            seed_cx[0] = xc = xy[0]; seed_cy[0] = yc = xy[1];
            double[][] estimated_xy = LocalizationFunction.ML_localization(1, pixels_double, w, h, seed_cx, seed_cy, kw, sigma, sigma, 2, false);  
            if(estimated_xy==null){                
                suc=false;
                break;
            }
            
            double[] mlxy = estimated_xy[0];
            if(mlxy==null){
                suc=false;
                break;
            }
            mlx[s] = mlxy[0];
            mly[s] = mlxy[1];
            
//            dx[s] = mlx[0] - mlx[s];
//            dy[s] = mly[0] - mly[s];
//            
//            ploter_driftXY.updateImage();
            
            if(show_registeration){
                if(use_roi) fht = fft.forward_transform(stack.getProcessor(s+1));
                FHT phase_fht = HighPassFFT.phase_shift(fht,XM-mlx[s],YM-mly[s]);
                ImageProcessor shift_img = HighPassFFT.inverse_fft_transform(phase_fht); 
                stack_reg.addSlice(shift_img);
                //stack_reg.setProcessor(shift_img, s+1);
            }
            IJ.showProgress(s, nframes);
        }
        IJ.showProgress(1.0);
        IJ.showStatus("tracking done");  
        
        if(!suc){
            IJ.showMessage("tracking is not sucessful");
            return;
        }
        
        if(show_registeration){
            ImagePlus stablized_stack = new ImagePlus("stablized_" + experiment.getShortTitle(), stack_reg);
            stablized_stack.show();
        }
        
        double[] dx = new double[nframes];
        double[] dy = new double[nframes];
        double[] xaxis = new double[nframes];     
        //IJ.log("dx dy");
        for(int s=0; s<nframes; s++){            
            dx[s] = mlx[0] - mlx[s];
            dy[s] = mly[0] - mly[s];
            xaxis[s] = s;
            //IJ.log(""+dx[s]+" "+dy[s]);
        }
               
        Plot ploter_driftXY = DriftAnalysisPlotter.get_plot_driftXY(xaxis, dx, dy, "drift from phase correlation");
        ploter_driftXY.setXYLabels("frame", "displacement(pixel) refer to 1st frame");
        ploter_driftXY.setLegend("x-drift\ny-drift\ndrift",Plot.AUTO_POSITION|Plot.LEGEND_TRANSPARENT);
        ploter_driftXY.show();	
    }
    /**
    * This method gets called by ImageJ / Fiji.
    *
    * @param arg can be specified in plugins.config
    * @see ij.plugin.PlugIn#run(java.lang.String)
    */
    @Override
    public void run(String arg) {
        boolean check_once = arg.contains("once");
        boolean consistent;        
        do {
                consistent=true;
                if (!showDialog()) return;					
                if (kernel_radius<0) {IJ.showMessage("the radius of kernel must greater than zero"); consistent=false;}                
                if (sigma<0) {IJ.showMessage("the gaussian simga must greater than zero"); consistent=false;} 
        } while (!consistent && !check_once);	
        if(!consistent) return;
        
        process();
    }
        
    boolean showDialog() {
	
	int[] wList = WindowManager.getIDList();	
	if (wList==null || wList.length<1) {
		IJ.showMessage("Phase Correlation", "one image are required");	
		return false;
	}

	int wlistlen = wList.length;
	
	ArrayList<String> titles_list = new ArrayList<String>();
	for (int i=0; i<wlistlen; i++) {
		ImagePlus imp = WindowManager.getImage(wList[i]);
		if(imp!=null) titles_list.add(imp.getTitle());
	}

	if(titles_list.isEmpty()){
		IJ.showMessage("Phase Correlation", "one image are required");	
		return false;
	}
	
	String titles[] = titles_list.toArray(new String[0]);
	GenericDialog gd = new GenericDialog("Compute Drift");
	gd.addChoice("Image:", titles, titles[0]);	
	gd.addNumericField("kernel radius:", kernel_radius, 0, 10, "pixel");
        //gd.addMessage("--------settings for Gaussian Fit localization------------");
        //gd.addNumericField("    Gaussian's sigma:", sigma, 2, 10, "pixel");
        //gd.addCheckbox("dark pillar?", dark_pillar);
        gd.addCheckbox("show intermediate image?", false);
        gd.addCheckbox("show stabilized image?", false);
	gd.showDialog();
	if (gd.wasCanceled()) return false;
	String choice_name1 = gd.getNextChoice();
	experiment = WindowManager.getImage(choice_name1);

	kernel_radius = (int) gd.getNextNumber();	        
        //sigma = gd.getNextNumber();
        //dark_pillar = gd.getNextBoolean();
        show_shift_img = gd.getNextBoolean();
        show_registeration=gd.getNextBoolean();
	return true;
    }
}

