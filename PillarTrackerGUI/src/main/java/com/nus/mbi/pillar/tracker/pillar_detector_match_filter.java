package com.nus.mbi.pillar.tracker;
// updates:
//2016-10-13: 
//			 fix the problem of adding integel point ROI, only works for double points.
//			 check oblique angle, combined with spacing check
import com.nus.mbi.pillar.detection.PillarDetector;
import ij.*;
import ij.plugin.*;
import ij.gui.GenericDialog;
import ij.process.ImageProcessor;

/**
 * @author Xiaochun Xu <mbixxc@nus.edu.sg>
 * @version 1.0
 * @since 2015-07-21
 */
public class pillar_detector_match_filter implements PlugIn {	

private ImagePlus experiment;

private int kernelw;		// the kernels and their dimensions
private double[] kernel; 	// used when a single kernel is applied
private int[] kernel_dx;
private int[] kernel_dy;

//private double diameter=16.;
private double spacing=24;
private double psf_sigma=7;
private double precentage = 5; //within (0~100)
private boolean dark_pillars = true;
private double oblique = 0;
private double oblique_tol = 20;
private double grid_angle = 90;
private int channels_orig; // # of used channels

private boolean plotting = false;
private boolean show_dialog = true;

private PillarDetector pillar_detector = new PillarDetector();
//public int num_local_maxima = 0;
//public int num_detected_maxima = 0;
/**
 * Set whether show the dialog, set it to false to run plugin without show dialog
 * @param show				show dialog?
 */
public void setShowDialog(boolean show){
    show_dialog = show;
}

/**
 * Set up the parameters used for pillar detection, then use {@link #process()}
 * 
 * @param ip 				time-series image stack 
 * @param spacing			pillar spacing in pixel
 * @param sigma				Gaussian sigma in pixel
 * @param oblique			Grid oblique in degree
 * @param grid_angle			Grid angle in degree
 * @param dark				dark pillar?
 * @return 					true if successfuly setup
 */
public boolean setup(ImagePlus ip, double spacing, double sigma, double oblique, double grid_angle, boolean dark)
{
	boolean suc = true;
	if(ip!=null) experiment = ip;
	else suc = false;

	if(sigma>0) psf_sigma = sigma;
	else suc = false;
	
	//if(diameter>0) this.diameter = diameter;
	//else suc = false;

	if(spacing>0) this.spacing = spacing;
	else suc = false;
        
        if(oblique<=45 && oblique>=-45) this.oblique = oblique;
	else suc = false;
	
        if(grid_angle>0 && grid_angle<=90) this.grid_angle = grid_angle;
	else suc = false;
        
	dark_pillars = dark;

	IJ.log("set up successful? " + suc);	
	IJ.log("sigma :" + psf_sigma);
	IJ.log("pillar spacing:" + this.spacing);	
        IJ.log("grid oblique:" + this.oblique);
	IJ.log("grid angle:" + this.grid_angle);	
	IJ.log("dark pillars:" + dark_pillars);	
	
	return suc;
}

public void setDarkPillar(boolean is_dark_pilalr){
    dark_pillars = is_dark_pilalr;
}

public int getNumLocalMaximas(){
    return pillar_detector.num_local_maxima;
}

public int getNumDetectedMaximas(){
    return pillar_detector.num_detected_maxima;
}

/**
 *  process pillars localization for the current time-series image stack,
 *  after setup {@link #setup(ImagePlus ip, double diameter, double spacing, double sigma, boolean dark)}
 */	
	public void process(){
		
                //PillarDetector pillar_detector = new PillarDetector();
		//double sigma = Math.min(sigmax_PSF, sigmay_PSF);	
                pillar_detector = new PillarDetector();
		pillar_detector.setup(experiment, oblique, spacing, psf_sigma, dark_pillars);
                pillar_detector.setGridAngle(grid_angle);
                pillar_detector.process();
                //num_local_maxima = pillar_detector.num_local_maxima;
                //num_detected_maxima = pillar_detector.num_detected_maxima;
                
		if(plotting){
                    pillar_detector.plot();                    
		}			
	}
        
        public void drawLocalMaxima(ImageProcessor ip, double threshold){
            pillar_detector.drawLocalMaxima(ip, threshold);
        }
        
        public int getNumLocalMaximas(ImageProcessor ip, double threshold){
            return pillar_detector.getNumLocalMaximas(ip, threshold);
        }
        
/**
 *  process pillars localization for the current time-series image stack,
 *  after setup {@link #setup(ImagePlus ip, double diameter, double spacing, double sigma, boolean dark)}
 *  @param ep                   ImageProcessor
 *  @param use_minimum_error    pre-processing use the minimum error solver?
 */	
	public void process(ImageProcessor ep, boolean use_minimum_error){
                //PillarDetector pillar_detector = new PillarDetector();
		//double sigma = Math.min(sigmax_PSF, sigmay_PSF);
                pillar_detector = new PillarDetector();
		pillar_detector.setup(experiment, oblique, spacing, psf_sigma, dark_pillars);
                pillar_detector.setGridAngle(grid_angle);
                pillar_detector.process(ep, use_minimum_error);
                //num_local_maxima = pillar_detector.num_local_maxima;
                //num_detected_maxima = pillar_detector.num_detected_maxima;
		if(plotting){
                    pillar_detector.plot();                    
		}			
	}


public int[][] convert_map(boolean[] maxl, int width, int height){
	int centC=0; //centroid counter per frame
	for(int k=0; k<width*height; k++)	if(maxl[k]) centC++;	
	int[][] centroidsXY=new int[2][centC];			
	centC=0; 
	for(int x=0; x<width; x++) {
		for(int y=0; y<height; y++) {
			if(maxl[x+y*width]){
			 	centroidsXY[0][centC]=x;
				centroidsXY[1][centC]=y;				 	
			 	centC++;
			}
		}
	}		

	return centroidsXY;
}

/**
 * run the plugin with dialog.
 * @param arg	not used 
 *
 */
public void run(String arg) { 
		IJ.freeMemory();
		if(show_dialog){
                    boolean consistent;			
                    do {
                            consistent=true;
                            if (!showDialog()) return;
                            channels_orig=experiment.getNChannels();
                            if (channels_orig>1) {IJ.log("multi channel image, only first channel will be processed");}			
                            //else if(kernelw<1) {IJ.showMessage("pillar diameter must be greater than zero"); consistent=false;}	
                            //else if(diameter<0) {IJ.showMessage("pillar diameter must be greater than zero"); consistent=false;}	
                            else if(spacing <0) {IJ.showMessage("pillar spacing must be must be greater than zero"); consistent=false;}	
                            else if(psf_sigma<0) {IJ.showMessage("Gaussian psf sigma must be greater than zero"); consistent=false;}		
                            else if(oblique<-45 || oblique>45) {IJ.showMessage("Grid oblique must be among the range of -45 to 45 degree "); consistent=false;}	
                            else if(grid_angle<=0 || grid_angle>90) {IJ.showMessage("Grid angle must be among the range of 0 to 90 degree "); consistent=false;}	
                    } while (!consistent);
                }
		process();
	}


boolean showDialog() {
	int[] wList = WindowManager.getIDList();
	if (wList==null || wList.length<1) {
		IJ.showMessage("PILLAR DETECTOR SHAPE", "a single channel image or time series required");
		return false;
	}
	String[] titles = new String[wList.length];
	for (int i=0; i<wList.length; i++) {
		ImagePlus imp = WindowManager.getImage(wList[i]);
		titles[i] = imp!=null?imp.getTitle():"";
	}

	GenericDialog gd = new GenericDialog("PILLAR DETECTOR SHAPE");	
	gd.addChoice("image or time series", titles, titles[0]);
	
	//gd.addMessage("pillar settings");
	//gd.addNumericField("pillar diameter in pixel:",diameter, 0);
	gd.addNumericField("pillar_spacing in pixel:", spacing, 2);
	gd.addNumericField("Grid_oblique( [-45~45] degree)", oblique, 2);
        gd.addNumericField("Grid_angle( (0~90] degree)", grid_angle, 2);
	gd.addNumericField("Gaussian_PSF_sigma", psf_sigma, 2);	
	gd.addCheckbox("dark pillars", dark_pillars);		
	gd.addCheckbox("plotting intermediant images?",plotting); 	
	
	gd.showDialog();
	if (gd.wasCanceled()) return false;

	int index1 = gd.getNextChoiceIndex();
	experiment = WindowManager.getImage(wList[index1]);

	//diameter=(double)gd.getNextNumber();
	spacing=(double)gd.getNextNumber();
	oblique=(double)gd.getNextNumber();
        grid_angle=(double)gd.getNextNumber();
	psf_sigma=(double)gd.getNextNumber();
	
	dark_pillars = gd.getNextBoolean();
	plotting = gd.getNextBoolean();
	return true;
}
}

