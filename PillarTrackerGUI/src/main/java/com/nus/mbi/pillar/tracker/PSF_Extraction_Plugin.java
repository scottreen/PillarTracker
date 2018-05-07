package com.nus.mbi.pillar.tracker;

import com.nus.mbi.pillar.detection.HighPassFFT;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.process.FHT;
import ij.process.FloatPolygon;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author xiaochun
 */
public class PSF_Extraction_Plugin implements PlugIn{
    private ImagePlus experiment;
    private int search_window_radius = 4;
    private int kernel_radius = 4;
    
    private boolean show_shift_img = false;
    public boolean sub_pixel = false;
    public double sigma = 1.5;
    public boolean dark_pillar = false;
    public boolean use_FFT_translate = false;
    
    public void setup(ImagePlus ip, int radius){
        experiment = ip;
        if(radius>0){
            search_window_radius = radius;
            kernel_radius = radius;
        }
    }
    
    public void process()
    {                
        Roi point_roi= experiment.getRoi();			
        if(point_roi==null || point_roi.isArea() || point_roi.isLine()){
            //IJ.showMessage("there is no point selections");
            return;
        }
        
        FloatPolygon polygon_roi = point_roi.getFloatPolygon();
        float[] roix = polygon_roi.xpoints;
        float[] roiy = polygon_roi.ypoints;
        int nrois = polygon_roi.npoints; 				                
        
        ImageProcessor ip = experiment.getProcessor();
        int w = ip.getWidth();
        int h = ip.getHeight();
        float[][] pixels = ip.getFloatArray();
        
        int kw = 2*search_window_radius + 1;
        double[][] psf = new double[kw][kw];        
        for(int i=0; i<kw; i++) for(int j=0; j<kw; j++) psf[i][j] = 0;
        int r = search_window_radius;
        int xmax = w-r-1;
        int ymax = h-r-1;     
        int num = 0;
        List<Integer> list = new ArrayList();
        for(int k=0; k<nrois; k++){
            int x = (int) Math.round(roix[k]);
            int y = (int) Math.round(roiy[k]);
            if(x>r && x<xmax && y>r && y<ymax){
                num++;
                list.add(k);
                for(int i=-r; i<=r; i++){                   
                    int yy = y + i;
                    int ii = i + r;
                    int dy2 = i*i;
                    for(int j=-r; j<=r; j++){
                        int xx = x + j;
                        int jj = j + r;
                        //if(dy2+j*j<=r*r)
                            psf[jj][ii] += pixels[xx][yy];
                    }
                }
            }
        }
        
        if(!sub_pixel){  
            FloatProcessor fp = new FloatProcessor(kw, kw);
            for(int i=0; i<kw; i++) for(int j=0; j<kw; j++) fp.setf(i, j, (float)(psf[i][j]/num));
            ImagePlus psf_ip = new ImagePlus("PSF_" + experiment.getShortTitle(), fp);
            psf_ip.show();
            return;
        }
        
        //localization
        nrois = list.size();        
        double[] seed_cx = new double[nrois];
        double[] seed_cy = new double[nrois];
        for(int k=0; k<nrois; k++){
            int i = list.get(k);
            seed_cx[k] = roix[i];
            seed_cy[k] = roiy[i];
        }
        
        int img_size = w*h;
        double[] pixels_double = new double[img_size];
        for(int k=0; k<img_size; k++){
            int x = k%w;
            int y = k/w;
            pixels_double[k] = pixels[x][y];
        }
        
        double[][] estimated_xy = LocalizationFunction.ML_localization(1, pixels_double, w, h, seed_cx, seed_cy, 2*kernel_radius + 1, sigma, sigma, 2, dark_pillar);  
        if(estimated_xy==null){
            IJ.showMessage("localization faild!");
            return;
        }
        
        num = 0;
        for(int k=0; k<nrois; k++){
            double[] ml_xy = estimated_xy[k];
            if(ml_xy!=null)num++;            
        }
        if(num<1){
            IJ.showMessage("can not find a point to start!");
            return;
        }
        
        String title = "shift_" + experiment.getShortTitle();        
        ImagePlus shift_imp = new ImagePlus(title, ip);
        if(show_shift_img) shift_imp.show();
        for(int i=0; i<kw; i++) for(int j=0; j<kw; j++) psf[i][j] = 0;
        num = 0;
        
        if(use_FFT_translate){
            IJ.showStatus("phase shift started");
            //FFT forward
            HighPassFFT fft = new HighPassFFT();
            FHT fht = fft.forward_transform(ip);
            //fht.setShowProgress(false);
            //fht.transform();

            //ImageProcessor shift_ip = HighPassFFT.inverse_fft_transform(fht.getCopy());
            //FHT phase_fht0 = HighPassFFT.phase_shift(fht,10.5,15.5);               
            //ImageProcessor shift_ip0 = HighPassFFT.inverse_fft_transform(phase_fht0);
            //ImagePlus shift_imp = new ImagePlus("shift_" + experiment.getShortTitle(), shift_ip0);
            for(int k=0; k<nrois; k++){
                double[] ml_xy = estimated_xy[k];
                if(ml_xy!=null) {    
                    double xc = ml_xy[0];
                    double yc = ml_xy[1];                
                    int x = (int) Math.round(xc);
                    int y = (int) Math.round(yc);
                    //FHT phase_fht = HighPassFFT.phase_shift(fht,xc-x,yc-y);    
                    FHT phase_fht = HighPassFFT.phase_shift(fht,x-xc,y-yc);    
                    ImageProcessor shift_ip = HighPassFFT.inverse_fft_transform(phase_fht);
                    if(show_shift_img) {
                        shift_imp.setProcessor(shift_ip);
                        shift_imp.setRoi(x-r, y-r, kw, kw);                    
                        //shift_imp.updateAndRepaintWindow();
                    }
                    float[][] shift_pixels = shift_ip.getFloatArray();
                    for(int i=-r; i<=r; i++){                   
                        int yy = y + i;
                        int ii = i + r;                    
                        for(int j=-r; j<=r; j++){
                            int xx = x + j;
                            int jj = j + r;                        
                            psf[jj][ii] += shift_pixels[xx][yy];
                        }
                    }                
                    num++;
                    IJ.showProgress(k, nrois);
                }
            }
        }
        else{   
            for(int k=0; k<nrois; k++){
                double[] ml_xy = estimated_xy[k];
                if(ml_xy!=null) {    
                    double xc = ml_xy[0];
                    double yc = ml_xy[1];                
                    int x = (int) Math.round(xc);
                    int y = (int) Math.round(yc);
                    //FHT phase_fht = HighPassFFT.phase_shift(fht,xc-x,yc-y);    
                    //FHT phase_fht = HighPassFFT.phase_shift(fht,x-xc,y-yc);    
                    ImageProcessor shift_ip = ip.duplicate();
                    shift_ip.setInterpolationMethod(ImageProcessor.BICUBIC);
                    //shift_ip.translate(x-xc,y-yc);
                    shift_ip.translate(xc-x,yc-y);
                    if(show_shift_img) {
                        shift_imp.setProcessor(shift_ip);
                        shift_imp.setRoi(x-r, y-r, kw, kw);                    
                        //shift_imp.updateAndRepaintWindow();
                    }
                    float[][] shift_pixels = shift_ip.getFloatArray();
                    for(int i=-r; i<=r; i++){                   
                        int yy = y + i;
                        int ii = i + r;                    
                        for(int j=-r; j<=r; j++){
                            int xx = x + j;
                            int jj = j + r;                        
                            psf[jj][ii] += shift_pixels[xx][yy];
                        }
                    }                
                    num++;
                    IJ.showProgress(k, nrois);
                }
            }
        }
        
        FloatProcessor fp = new FloatProcessor(kw, kw);
        for(int i=0; i<kw; i++) for(int j=0; j<kw; j++) fp.setf(i, j, (float)(psf[i][j]/num));
        ImagePlus psf_ip = new ImagePlus("PSF_" + experiment.getShortTitle(), fp);
        psf_ip.show();
        
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
                if (search_window_radius<0) {IJ.showMessage("the radius of search window must greater than zero"); consistent=false;}                
                if (kernel_radius<0) {IJ.showMessage("the radius of kernel must greater than zero"); consistent=false;}                
                if (sigma<0) {IJ.showMessage("the gaussian simga must greater than zero"); consistent=false;} 
        } while (!consistent && !check_once);	
        if(!consistent) return;
        
        Roi point_roi= experiment.getRoi();			
        if(point_roi==null || point_roi.isArea() || point_roi.isLine()){
            IJ.showMessage("there is no point selections");
            return;
        }
        
        process();
    }
        
    boolean showDialog() {
	
	int[] wList = WindowManager.getIDList();	
	if (wList==null || wList.length<1) {
		IJ.showMessage("Kernel Extraction", "one image are required");	
		return false;
	}

	int wlistlen = wList.length;
	
	ArrayList<String> titles_list = new ArrayList<String>();
	for (int i=0; i<wlistlen; i++) {
		ImagePlus imp = WindowManager.getImage(wList[i]);
		if(imp!=null) titles_list.add(imp.getTitle());
	}

	if(titles_list.isEmpty()){
		IJ.showMessage("Kernel Extraction", "one image are required");	
		return false;
	}
	
	String titles[] = titles_list.toArray(new String[0]);
	GenericDialog gd = new GenericDialog("Kernel Extraction");
	gd.addChoice("Image:", titles, titles[0]);
	
        gd.addNumericField("window radius:", search_window_radius, 0, 10, "pixel");	
        
        gd.addCheckbox("subpixel registration by Gaussian Fit?", sub_pixel);        
        gd.addMessage("--------settings for Gaussian Fit------------");
        gd.addNumericField("    Kernel radius:",kernel_radius , 0, 10, "pixel");
        gd.addNumericField("    Gaussian's sigma:", sigma, 2, 10, "pixel");
        gd.addCheckbox("dark pillar?", dark_pillar);
        gd.addCheckbox("use_fft translation?(otherwise,interpolation)", use_FFT_translate);
        gd.addCheckbox("show intermediate shifted image?", false);
        
	gd.showDialog();
	if (gd.wasCanceled()) return false;
	String choice_name1 = gd.getNextChoice();
	experiment = WindowManager.getImage(choice_name1);
        search_window_radius = (int) gd.getNextNumber();		
        sub_pixel = gd.getNextBoolean();
        kernel_radius = (int) gd.getNextNumber();
        sigma = gd.getNextNumber();
        dark_pillar = gd.getNextBoolean();
        use_FFT_translate=gd.getNextBoolean();
        show_shift_img = gd.getNextBoolean();
	return true;
    }

}

