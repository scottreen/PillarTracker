package com.nus.mbi.pillar.tracker;

import ij.*;
import ij.io.*;
import ij.plugin.*;
import ij.process.*;
import ij.gui.*;
import fiji.util.gui.*;
import java.io.*;
import java.io.File;
import ij.util.ThreadUtil; //using the mulit threading

import com.nus.mbi.pillar.detection.PillarDetector;
import com.nus.mbi.pillar.detection.CrossCorreclation_Plugin;
import com.nus.mbi.pillar.drift.DriftAnalysisPlotter;
import com.nus.mbi.pillar.drift.DriftCorrection;
import com.nus.mbi.pillar.stat.IntPoint;

import ij.plugin.filter.RankFilters;
import java.awt.Rectangle;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Xiaochun Xu <mbixxc@nus.edu.sg>
 * @version 2.3
 * @since 2015-07-09
 */
public class pillar_tracking implements PlugIn {

//measurment image stack
ImagePlus experiment;
ImageStack experiment_s;
public double lattice = 1.0;
public double diameter = Double.NaN;

double spacing = 24;   //in pixel
double oblique = 1.5;
double oblique_tol = 20;
double grid_angle = 90;
double sigmax_PSF=7; //in pixel
double sigmay_PSF=7; //in pixel

double catch_radius = 12;
Roi [] rois;
int nrois;

// the kernels and their dimensions
int kernel_w = 21;
int box_constrian_R = kernel_w/2;
// used when a single kernel is applied
double[] kernel;//[]=new float[kernel_w*kernel_w];
int kernel_size;
double[] window;


int channels_orig; // how man channels are originally present

int vector_zoom = 15; // how much zoomed in should the vector plot appear?

boolean dark_pillars = false;
boolean debug = false;
boolean plotting = false;
//private boolean otsu_force_8bits = true;

int numThread = 0;
String output_fname = "pillar_tracks.bin";

boolean use_minimum_std = false;

public static int fileversion1 = 589450;

//write into text file
FileWriter out_text_File = null;
PrintWriter out_text = null;
SimpleDateFormat time_formatter = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
String timestamp_start;

private boolean gpu_on = false;

public boolean apply_mean_rank_filter = false;

public static int localization_algorithm_CG = 0;
public static int localization_algorithm_Levmar = 1;
int localization_algorithm = localization_algorithm_Levmar;
boolean use_metric_CG = false;
boolean show_dialog = true;

public boolean apply_drift_creation = true;

CrossCorreclation_Plugin image_enhancer_plugin = null;
/**
 * Set Enhancement Solver
 * @param enhancer				Enhancement solver?
 */
public void setEnhancementSolver(CrossCorreclation_Plugin enhancer){
    image_enhancer_plugin = enhancer;    
}

/**
 * Set whether show the dialog, set it to false to run plugin without show dialog
 * @param show				show dialog?
 */
public void setShowDialog(boolean show){
    show_dialog = show;
}

/**
 * Set up the parameters used for pillar tracking, then use {@link #process()}
 * 
 * @param ip 				time-series image stack
 * @param sigma 			Gaussian PSF Sigma used for Gaussian Fit
 * @param window_width		window width  
 * @param limit 			constrain box limit used in the solver
 * @param dark				dark pillar?
 * @param spacing			pillar spacing in pixel
 * @param oblique			Grid oblique in degree
 * @param num_threads		num of threads used for multi-threading
 * @param output_filename	output file path of binary file
 * @return 					true if successfuly setup
 */
public boolean setup(ImagePlus ip, double sigma, int window_width, int limit, boolean dark, double spacing,  double oblique, double catch_radius, int num_threads, String output_filename){
	boolean suc = true;
	if(ip!=null) experiment = ip;
	else suc = false;
	
	if(sigma>0) {sigmax_PSF = sigma; sigmay_PSF = sigma; }
	else suc = false;
	
	if(window_width>0) kernel_w = (window_width/2)*2+1;
	else suc = false;
	
	if(limit>0)	box_constrian_R = limit;
	else suc = false;

	if(spacing>0)	this.spacing = spacing;
	else suc = false;

//	if(oblique>0)	this.oblique = oblique;
//	else suc = false;
// not check the oblique angle
        this.oblique = oblique;

	if(catch_radius>0)	this.catch_radius = catch_radius;
	else suc = false;
	
	if(num_threads>0) numThread = num_threads;
	else suc = false;
	
	output_fname = output_filename;

	dark_pillars = dark;
	
	IJ.log("set up successful for pillar tracking? " + suc);	
	IJ.log("sigma X:" + sigmax_PSF + "	sigma Y:" + sigmay_PSF);	
	IJ.log("pillar spacing:" + this.spacing);			
	IJ.log("grid oblique angle:" + this.oblique);
	IJ.log("search window:" + kernel_w);		
	IJ.log("constrian limit:" + box_constrian_R);
	IJ.log("maximum drift between frames:" + catch_radius);
	IJ.log("num threads:" + numThread);
	IJ.log("dark pillars:" + dark_pillars);
	IJ.log("save to binary file:" + output_fname);			
	
	return suc;
}

/**
* Set up the grid oblique tolerance used for pillar detection, then use {@link #process()}
* 	 
* @param grid_angle		grid angle ( 0~90 degree)
* @return 					true if successfuly setup
*/
public boolean setGridAngle(double grid_angle)
{		
       boolean suc = true;
       if(grid_angle>0 && grid_angle<=90){
           this.grid_angle = grid_angle;           
       }
       else suc = false;		

       IJ.log("grid angle:" + this.grid_angle);	

       return suc;
}

/**
 * Set up the parameters used for pillar tracking, then use {@link #process()}
 *  @param use_minimum		use minimum std?
 */
public void setDriftCorrectionMethod(boolean use_minimum){
	use_minimum_std = use_minimum;	
	IJ.log("use minimum std? " + use_minimum_std);		
}

/**
 * Set GPU on/off for pillar tracking, then use {@link #process()}, GPU is on by default.
 *  @param use_gpu		use GPU computing?
 */
public void setGPUOn(boolean use_gpu){
	gpu_on = use_gpu;	
	IJ.log("use gpu? " + gpu_on);		
}

/**
 * Set GPU on/off for pillar tracking, then use {@link #process()}, GPU is on by default.
 *  @param use_gpu		use GPU computing?
 */
public void setAlgorithm(int solver){
    localization_algorithm = solver;
    use_metric_CG = (solver == localization_algorithm_CG);
    IJ.log("localization algorithm: " + localization_algorithm);
    IJ.log("use metric CG? " + use_metric_CG);
}

/**
 *  process pillars localization for the current time-series image stack,
 *  after setup {@link #setup(ImagePlus ip, double sigma, int window_width, int limit, boolean dark, int num_threads, String output_filename)}
 */
public void process(){
		if(experiment==null) return;
		timestamp_start = time_formatter.format(System.currentTimeMillis());
                IJ.log("------------Start Tracking at "+ timestamp_start+"------------");
                experiment.setOverlay(null);
		
		int width = experiment.getWidth();
		int height = experiment.getHeight();		
		experiment_s = experiment.getStack();
		
		IJ.log("channels="+experiment.getNChannels()+"\nframes="+experiment.getNFrames()+"\nslices="+experiment.getNSlices());
		if (experiment.isHyperStack()) IJ.log("is hyperstack\n");
		IJ.log("stack equivalent slices="+experiment_s.getSize());
		experiment.changes = false;

		//show the file path
		FileInfo fi = experiment.getOriginalFileInfo();
		String info = experiment.getTitle();
		if(fi != null){
			info = fi.directory + fi.fileName;
			IJ.log("file path=" + info);
		}
		
		kernel_size = kernel_w*kernel_w;
		int sliceN = Math.max(experiment.getNSlices(),experiment.getNFrames());
                
                boolean use_enhancer = image_enhancer_plugin!=null;
                boolean dark_object = use_enhancer? false : dark_pillars;
                //if(image_enhancer_plugin!=null) dark_pillars = false;
                
                PillarDetector pillar_detector = new PillarDetector();
		double sigma = Math.min(sigmax_PSF, sigmay_PSF);                
		pillar_detector.setup(experiment, oblique, spacing, sigma, dark_object);
		pillar_detector.setGridAngle(grid_angle);
		//if (roimanager_null)
		{
			boolean suc_load_roi = true;			
			
			Roi point_roi= experiment.getRoi();			
			if(point_roi==null || point_roi.isArea() || point_roi.isLine()){
				//no point selections, try to creat new point roi firstly				
                                ImageProcessor ep = experiment_s.getProcessor(1);
                                                                
                                if(apply_mean_rank_filter){
                                    RankFilters mean_filter = new RankFilters();
                                    ImageProcessor slice_filter = ep.duplicate();
                                    mean_filter.rank(slice_filter, 1, RankFilters.MEAN);
                                    ep = slice_filter;
                                }
                                
                                if(image_enhancer_plugin!=null){
                                    ep = image_enhancer_plugin.process(ep);
                                    pillar_detector.process(ep, false);
                                }
                                else pillar_detector.process(ep, true);
				point_roi= experiment.getRoi();
			}
			
			if(point_roi==null) suc_load_roi = false;			
			else if(!point_roi.isArea() && !point_roi.isLine()){
				FloatPolygon polygon_roi = point_roi.getFloatPolygon();
				float[] roix = polygon_roi.xpoints;
				float[] roiy = polygon_roi.ypoints;
                                int margin = (int)Math.ceil(spacing);
                                int margin_w = 2*margin+1;
                                Roi roi_margin = new Roi(margin,margin,width-margin_w,height-margin_w);
                                int num_points = polygon_roi.npoints; 				
                                int num_valid = 0;
                                for(int i=0; i<num_points; i++){
                                    int x= (int)Math.floor(roix[i]);
                                    int y= (int)Math.floor(roiy[i]);
                                    if(roi_margin.contains(x, y))num_valid++;
                                }
				//nrois = polygon_roi.npoints; 				
                                nrois = num_valid;
				if(nrois>0){
					rois = new Roi[nrois];
					int halfw_roi = kernel_w/2;//roi_half_width;//
					int roi_w = 2*halfw_roi + 1;
                                        num_valid = 0;
					for(int i=0; i<num_points; i++){
                                            int x = (int)Math.floor(roix[i]);
                                            int y = (int)Math.floor(roiy[i]);    
                                            if(roi_margin.contains(x, y)){
                                                rois[num_valid] = new Roi(x-halfw_roi,y-halfw_roi,roi_w,roi_w);
                                                num_valid++;
                                            }
					}					
                                }
				else suc_load_roi = false; 
			}
			else suc_load_roi = false;
			
			if(!suc_load_roi) {
				IJ.log("There is no any Point selection, the centriod of image will be procesed!");									
				nrois = 1;
				rois = new Roi[nrois];				
				Roi area_roi= experiment.getRoi();
				if(area_roi==null || false==area_roi.isArea()) rois[0] = new Roi(0,0,width,height);						
				else rois[0] = new Roi(area_roi.getBounds());				
			}
			
		}
		
		double[] distance = new double[sliceN];
		double[] xaxis = new double[sliceN];
		for(int s=0;s<sliceN;s++) xaxis[s] = s;	
		
		Rectangle[] rect_rois = new Rectangle[nrois];
		int[] rcx_rois = new int[nrois];
		int[] rcy_rois = new int[nrois];
		
		for(int i=0; i<nrois;i++) {
			Rectangle rect =rois[i].getBounds();	
			if(debug) IJ.log(" rect x="+rect.x+" y="+ rect.y+ " w="+rect.width+" h="+ rect.height);
			rcx_rois[i] = (rect.x + rect.width/2);//rect.getCenterX();
			rcy_rois[i] = (rect.y + rect.height/2);//rect.getCenterY();
			rect_rois[i] = rect;
		}

		int size = width*height;
		double[] pixels = new double[size];

                int numCPUs = ThreadUtil.getNbCpus();
                if(numThread<=0) numThread = numCPUs;
                else{
                        if(numThread>nrois) numThread = nrois;
                        int maxthreads = 20*numCPUs;
                        if(numThread>maxthreads) numThread=maxthreads;			
                }

                if(image_enhancer_plugin!=null) image_enhancer_plugin.setNumThreads(numThread);
                
		IJ.log("num of threads:" + numThread);
		if(!pillar_detector.isPrepared()) pillar_detector.prepare();
		// prepare centroid buffers
		//tracing settings
		boolean[][] active=new boolean[sliceN][nrois];
		double[][] trackX=new double[sliceN][nrois];
		double[][] trackY=new double[sliceN][nrois];
                
                int[] tx = new int[nrois];
                int[] ty = new int[nrois];                   	
                // seed first frame		
                for(int i=0; i<nrois; i++) {
                        active[0][i]=true;
                        trackX[0][i]=tx[i]=rcx_rois[i];//centroidsX[0][i];
                        trackY[0][i]=ty[i]=rcy_rois[i];//centroidsY[0][i];                
                }
                
		//optimization
		double[][] cx = new double[nrois][sliceN];
		double[][] cy = new double[nrois][sliceN];
		double[][] disX = new double[nrois][sliceN];
		double[][] disY = new double[nrois][sliceN];

		double[] seed_cx = new double[nrois];
		double[] seed_cy = new double[nrois];	
                
                LocalizationWindowQueue queue = new LocalizationWindowQueue();

		//detection
                //IJ.log("frame#  detected    matched     forward link12  forward linkt2"); 
                //experiment.deleteRoi();
                IJ.log("frame#  detected    matched"); 

		for(int s=0;s<sliceN;s++){//pillar detection for each frame
			//detection
                        ImageProcessor slice = experiment_s.getProcessor(s + 1);  
                        pixels = getSlicePixels(slice);
                        
			if(s>0){
                            //detection
                            int[][] centroidsXY;
                            if(apply_mean_rank_filter){
                                RankFilters mean_filter = new RankFilters();
                                ImageProcessor slice_filter = slice.duplicate();
                                mean_filter.rank(slice_filter, 1, RankFilters.MEAN);
                                double[] filter_pixels = getSlicePixels(slice_filter);                                 
                                centroidsXY = pillar_detector.detect_fast(filter_pixels); 
                            }
                            else centroidsXY = pillar_detector.detect_fast(pixels); //use faster detection version, without minimum error check.

                            int ncent = centroidsXY[0].length;
       
                            int[] cx2 = centroidsXY[0];//centroidsX[s];                
                            int[] cy2 = centroidsXY[1];//centroidsY[s];
                            
                            boolean[] active_flag = new boolean[nrois];
                            for(int i=0; i<nrois; i++) active_flag[i]=false;
                            List<IntPoint> pair = cross_check(tx, ty, cx2, cy2, width, height, catch_radius);
                            for(int k=0; k<pair.size(); k++) {
                                IntPoint p = pair.get(k);
                                int i  = p.x;
                                int i2 = p.y;
                                tx[i]=cx2[i2];
                                ty[i]=cy2[i2];			
                                active_flag[i] = true;                            
                            }
                            
                            int matches=pair.size(); 
                            for(int i=0; i<nrois; i++) {
                                active[s][i]=active_flag[i];
                                trackX[s][i]=tx[i];
                                trackY[s][i]=ty[i]; 
                            }

                            IJ.log((s+1)+"  "+ncent+"   "+matches);
                        }
                        else{
                            IJ.log((s+1)+"  "+nrois);
                        }
			
                        //optimization
                        //for(int p=0; p<size;p++) pixels[p] = slice.getf(p);
			for(int y=0; y<nrois; y++) {
                            seed_cx[y] = trackX[s][y];
                            seed_cy[y] = trackY[s][y];
			}				

                        double[][] estimated_xy = null;
                        if(localization_algorithm == localization_algorithm_CG){
                            estimated_xy = LocalizationFunction.metric_centroid_safe(numThread, pixels, width, height, seed_cx, seed_cy, kernel_w/2, 1, box_constrian_R);
                            for(int i=0; i<nrois; i++){
                                cx[i][s] = Double.NaN;
                                cy[i][s] = Double.NaN;    
                                if(estimated_xy!=null){
                                    double[] ml_xy = estimated_xy[i];
                                    if(ml_xy!=null){
                                        cx[i][s] = ml_xy[0];
                                        cy[i][s] = ml_xy[1];
                                    }                                
                                }
                            }                      
                        } //else mulit_threading(numThread, pixels,width,height,s,cx,cy,seed_cx,seed_cy,rect_rois); 
                        else{
                            for(int k=0; k<nrois; k++) {
                                double xx = seed_cx[k];
                                double yy = seed_cy[k];
                                if(Double.isNaN(xx) || Double.isNaN(yy)){
                                    cx[k][s] = Double.NaN;
                                    cy[k][s] = Double.NaN;    
                                }
                                else{
                                    int x = (int)Math.round(xx);
                                    int y = (int)Math.round(yy);                                  
                                    LocalizationWindow win = LocalizationWindow.create(x, y, pixels, width, height, kernel_w);                                    
                                    if(queue.isFull()){
                                        estimated_xy = LocalizationFunction.ML_localization(numThread, queue, sigma, sigma, box_constrian_R, dark_object);                                        
                                        queue.reset();
                                        
                                        if(estimated_xy!=null){
                                            int num = estimated_xy.length;
                                            int[] pid = queue.getPillarID();
                                            int[] fid = queue.getFrameID();
                                            for(int p=0; p<num; p++){
                                                double[] ml_xy = estimated_xy[p];
                                                int i = pid[p];
                                                int f = fid[p];
                                                if(ml_xy!=null){                                        
                                                    cx[i][f] = ml_xy[0];
                                                    cy[i][f] = ml_xy[1];
                                                }
                                                else{
                                                    cx[i][f] = Double.NaN;
                                                    cy[i][f] = Double.NaN;    
                                                }
                                            }
                                        }
                                    } 
                                    queue.add(win, k, s);
                                }
                            }
                            //estimated_xy = LocalizationFunction.ML_localization(numThread, pixels, width, height, seed_cx, seed_cy, kernel_w, sigma, sigma, box_constrian_R, dark_object);                            
                        }
                                                
			IJ.showStatus("tracking in slice="+s);
			IJ.showProgress(s+1,sliceN);
		}
                
                if(queue.getCount()>0){
                    double[][] estimated_xy = LocalizationFunction.ML_localization(numThread, queue, sigma, sigma, box_constrian_R, dark_object);
                    if(estimated_xy!=null){
                        int num = estimated_xy.length;
                        int[] pid = queue.getPillarID();
                        int[] fid = queue.getFrameID();
                        for(int k=0; k<num; k++){
                            double[] ml_xy = estimated_xy[k];
                            int i = pid[k];
                            int f = fid[k];
                            if(ml_xy!=null){                                        
                                cx[i][f] = ml_xy[0];
                                cy[i][f] = ml_xy[1];
                            }
                            else{
                                cx[i][f] = Double.NaN;
                                cy[i][f] = Double.NaN;    
                            }    
                        }
                    }
                }
                
                for(int s=0; s<sliceN; s++){
                    for(int i=0; i<nrois; i++){
                            disX[i][s] = cx[i][0] - cx[i][s];
                            disY[i][s] = cy[i][0] - cy[i][s];	
                    }
                } 

		//saving the localization data.	
		FileOutputStream fos = null;
                FileChannel fch = null;				
		try{
			File output_file = new File(output_fname);			
			fos = new FileOutputStream(output_file);	
                        fch = fos.getChannel();								
			IJ.log("write raw data into binary file:" + output_fname);
			IJ.showStatus("writing raw data into binary file");                        	
                        writeFileHeader_ver1(fch,sliceN);
                        writeCentroids(fch,rcx_rois, rcy_rois);			
                        writePoints(fch, cx, cy);
			// Write to text file
			out_text_File = new FileWriter(output_fname+".txt");
			out_text = new PrintWriter(out_text_File);
                        out_text.println("Start at " + timestamp_start);
			out_text.println("image from:" + info);
			out_text.println("----parameters used-----");
			out_text.println("sigma X:" + sigmax_PSF + "	sigma Y:" + sigmay_PSF);
			out_text.println("pillar spacing:" + spacing);	
			out_text.println("grid oblique angle:" + oblique);	
                        out_text.println("grid shape angle:" + grid_angle);	
			out_text.println("search window:" + kernel_w);			
			out_text.println("searching radius:" + box_constrian_R);
			out_text.println("maximum drift between frames:" + catch_radius);
			out_text.println("dark pillars:" + dark_pillars);
			out_text.println("use minimum std:" + use_minimum_std);
                        out_text.println("use mean filter to find maxima? " + apply_mean_rank_filter);
			out_text.println("number of pillars:" + nrois);
			out_text.println("number of frames:" + sliceN);
			out_text.println("save to binary file:" + output_fname);			
			//out_text.println("\r\n--------Important Notice:Binary file is in the big-endian order.------\r\nWithout specification as before, the file was saved as little-endian order!\r\n");
			IJ.log("write logs into text file:" + output_fname+".txt");			
			OpenDialog.setDefaultDirectory(output_file.getAbsoluteFile().getParent());			
			IJ.showStatus("done with saving raw data");
		}
		catch(Exception ex){
			IJ.log("write file error!");
		}		
		
                double[][] drift_correctX = null;
                double[][] drift_correctY = null;
                
		if(sliceN<=1) plotting = false;		
		boolean workingOnWindows = IJ.isWindows();
		boolean[] reference_flags = null;
		//using the reduced std algorithm to correct drift if have the multiple pillars.	
		if(apply_drift_creation && sliceN>1 && nrois>30){
			reference_flags = new boolean[nrois];//null;
			double[] varsum_array = new double[nrois];//null;
			double[] driftXY = null;
                        if(use_minimum_std){
                            driftXY = DriftCorrection.get_driftXY_minimum_std(cx, cy, reference_flags, varsum_array);
                        }
                        else{
                            driftXY = DriftCorrection.get_driftXY_reduced_std(cx, cy, reference_flags, varsum_array);
                        }
                        if(driftXY!=null){
                            //double[] drift_correctXY = new double[nrois*2*sliceN];			
                            drift_correctX = new double[nrois][sliceN];			
                            drift_correctY = new double[nrois][sliceN];			
                            double[] driftX = new double[sliceN];
                            double[] driftY = new double[sliceN];
                            for(int i=0; i<sliceN; i++){					
                                    double dx = driftXY[2*i];
                                    double dy = driftXY[2*i+1];				
                                    driftX[i] = dx; 
                                    driftY[i] = dy;

                                    for(int k=0; k<nrois; k++) {
                                            double ccx = cx[k][i] + dx;
                                            double ccy = cy[k][i] + dy;
                                            
                                            drift_correctX[k][i] = ccx;
                                            drift_correctY[k][i] = ccy;

                                            cx[k][i] = ccx;
                                            cy[k][i] = ccy;

                                            disX[k][i] = cx[k][0]-cx[k][i];
                                            disY[k][i] = cy[k][0]-cy[k][i];
                                    }			
                            }

                            try{
                                    //IJ.log("write corrected data into binary file:" + output_fname);
                                    IJ.showStatus("writing corrected data into binary file");
                                    //writeFloatImage(fch,drift_correctXY);	
                                    writePoints(fch, drift_correctX, drift_correctY);
                                    writeFloatImage(fch,driftXY);	
                                    writeBooleanArray(fch,reference_flags);

                                    out_text.println("----drift XY----");
                                    out_text.println("X	Y");
                                    writeFloatImage2Text(out_text, driftXY, 2);

                                    out_text.println("----reference pillars----");
                                    out_text.println("	pillar_ID	quiet_pillar?");
                                    writeBooleanArray2Text(out_text,reference_flags);

                                    IJ.showStatus("done with saving corrected data");
                            }
                            catch(Exception ex){
                                    IJ.log("write drift file error!");
                            }						

                            if(plotting && workingOnWindows)
                            {					
                                    //plot the drift xy
                                    Plot ploter_driftXY = DriftAnalysisPlotter.get_plot_driftXY(xaxis,driftX,driftY,"drift x&y");
                                    ploter_driftXY.setLegend("x-drift\ny-drift\ndrift",Plot.AUTO_POSITION|Plot.LEGEND_TRANSPARENT);
                                    ploter_driftXY.show();	

                                    double [] plot_x = new double[varsum_array.length];
                                    double [] plot_y = new double[varsum_array.length];
                                    double [] plot_std = new double[varsum_array.length];
                                    for (int i = 0; i < varsum_array.length; i++){ 
                                            plot_x[i] = i + 1; 			
                                            plot_y[i] = Math.sqrt(varsum_array[i]);
                                            plot_std[i] = Math.sqrt(var(cx[i])+var(cy[i]));
                                    }		

                                    if(use_minimum_std){
                                            Plot ploter = new Plot("reduced std xy", "num subset", "mean xy std(in pixel)",plot_x,plot_y);
                                            ploter.draw();
                                            ploter.show();
                                    }

                                    Plot ploter_std = new Plot("std xy", "pillar #", "std xy(in pixel)",plot_x,plot_std);
                                    ploter_std.draw();
                                    ploter_std.show();
                            }
                        }
		}

		//close the file handle
		try{
                        if(fch!=null) fch.close();
                        if(fos!=null) fos.close();
                        String timestamp_end = time_formatter.format(System.currentTimeMillis());
                        IJ.log("----------Tracking End at "+ timestamp_end+"------------");
                        out_text.println("End at " + timestamp_end);
			out_text.close();		
		}
		catch(Exception ex){
			IJ.log("close file error!");
		}
                
                saveCSV(cx, cy, drift_correctX, drift_correctY);
}

public void saveCSV(double[][] cx, double[][] cy, double[][] corrected_CX, double[][] corrected_CY){
        if(cx==null) return;
        int np = cx.length;
        if(np<1) return;        
        int nf = cx[0].length;
        if(nf<1) return;   
        try{    
            out_text_File = new FileWriter(output_fname+".csv");// Write to text file
            out_text = new PrintWriter(out_text_File);
            IJ.showStatus("writing raw data into csv file");
            boolean corrected = (corrected_CX!=null);
            if(corrected) out_text.println("frame,trajectory,x,y,x_raw,y_raw");
            else out_text.println("frame,trajectory,x,y");
            for(int f=0; f<nf; f++) {        
                for(int c=0; c<np; c++){
                    out_text.print(""+ (f+1) + ","+ (c+1)+",");
                    if(corrected) out_text.print((corrected_CX[c][f] + ","+ corrected_CY[c][f]+","));
                    out_text.print((cx[c][f] + ","+ cy[c][f]));
                    //out_text.print((disX[c][f] + ","+ disY[c][f]));
                    out_text.println();
                }
            }
            IJ.log("write data into CSV file:" + output_fname+".csv"); 
            IJ.showStatus("done with saving CSV data");
        }
        catch(Exception ex){
            IJ.log("write file error!");
        }        
        //close the file handle
        try{	
            out_text.close();
        }
        catch(Exception ex){
            IJ.log("close file error!");
        }
    }

public void save_current_settings(String filename){        
        int sliceN = 0;        
        nrois = 0;
        int[] rcx_rois = null;
        int[] rcy_rois = null;
        if(experiment!=null){
            sliceN = Math.max(experiment.getNSlices(),experiment.getNFrames());        
            Roi point_roi= experiment.getRoi();
            if(point_roi!=null && !point_roi.isArea() && !point_roi.isLine()){
                FloatPolygon polygon_roi = point_roi.getFloatPolygon();
                float[] roix = polygon_roi.xpoints;
                float[] roiy = polygon_roi.ypoints;            
                int num_points = polygon_roi.npoints;            

                rcx_rois = new int[num_points];
                rcy_rois = new int[num_points];
                for(int i=0; i<num_points; i++){
                    rcy_rois[i] = (int)Math.floor(roix[i]);
                    rcy_rois[i] = (int)Math.floor(roiy[i]);
                }            
                nrois = num_points;        
            }
        }
        
        FileOutputStream fos = null;
        FileChannel fch = null;        
        try{
            File output_file = new File(filename);
            fos = new FileOutputStream(output_file);
            fch = fos.getChannel();
            IJ.log("write current settings into binary file:" + filename);
            IJ.showStatus("writing current settings into binary file");
            writeFileHeader_ver1(fch,sliceN);
            if(nrois>0) writeCentroids(fch,rcx_rois,rcy_rois);
        }
        catch(Exception ex){
            IJ.log("write file error!");
        }
        
        //close the file handle
        try{	
            if(fch!=null) fch.close();
            if(fos!=null) fos.close();            
        }
        catch(Exception ex){
            IJ.log("close file error!");
        }
    }

double[] getSlicePixels(ImageProcessor slice){    
    int width = slice.getWidth();
    int height = slice.getHeight();
    int size = width*height;
    
    float[] pixels_float = new float[size];
    double[] pixels = new double[size];
    for(int p=0; p<size;p++) pixels_float[p] = slice.getf(p);
    if(image_enhancer_plugin!=null) pixels_float = image_enhancer_plugin.process(pixels_float, width, height);    
    for(int p=0; p<size;p++) pixels[p] = pixels_float[p];    
   
    return pixels;
}

static List<IntPoint> cross_check(int[] cx1, int[] cy1, int[] cx2, int[] cy2, int width, int height, double catch_radius){
    int[][] mask1 = PillarDetector.create_mask(cx1, cy1, width, height);	
    int[][] mask2 = PillarDetector.create_mask(cx2, cy2, width, height);	
    int[][] map   = PillarDetector.cross_check_frames(cx1, cy1, cx2, cy2, mask1, mask2, catch_radius);
    int[] forward = map[0];    
    int npillars = cx1.length;
    List<IntPoint> active_list = new ArrayList();
    for(int i=0; i<npillars; i++) {        
        int i2 = forward[i];
        if(i2>=0) active_list.add(new IntPoint(i, i2));
    }             
    return active_list;
}

static List<IntPoint> cross_check(int[] cx1, int[] cy1, int[] cx2, int[] cy2, int[][] mask1, int[][] mask2, double catch_radius){	
    int[][] map   = PillarDetector.cross_check_frames(cx1, cy1, cx2, cy2, mask1, mask2, catch_radius);
    int[] forward = map[0];    
    int npillars = cx1.length;
    List<IntPoint> active_list = new ArrayList();
    for(int i=0; i<npillars; i++) {        
        int i2 = forward[i];
        if(i2>=0) active_list.add(new IntPoint(i, i2));
    }             
    return active_list;
}

/**
 * run the plugin with dialog.
 * Marco command example: 
 * {@code run("Pillar Localization", "only=75fps_cell2_dish2_large_X27.tif sigma=7 search_window_width=23 constrian_radius=10 zoom-in=15 number_of_threads=20 dark plotting output=pillar-tracks-big.bin");}
 * 
 * @param arg	not used 
 *
 */
public void run(String arg) { 
		
		IJ.freeMemory();
		if(show_dialog){
                    boolean checkonce = arg.contains("checkonce");
                    boolean consistent;	
                    do {
                            consistent=true;
                            if (!showDialog()) return;
                            channels_orig=experiment.getNChannels();
                            if (channels_orig!=1) {IJ.showMessage("stack must have ONLY one channels"); consistent=false;}
                            if(spacing <0) {IJ.showMessage("pillar spacing must be must be greater than zero"); consistent=false;}				
                            if(oblique<-45 || oblique>45) {IJ.showMessage("grid oblique must be among the range of -45 to 45 degree "); consistent=false;}
                            if(grid_angle<=0 || grid_angle>90) {IJ.showMessage("grid angle must be among the range of 0 to 90 degree "); consistent=false;}
                            if (box_constrian_R<=0) {IJ.showMessage("searching radius must greater than zero"); consistent=false;}
                            if (box_constrian_R>kernel_w){ IJ.showMessage("The searching radius is not valid, which must be smaller than half of window width"); consistent = false;}
                            if (catch_radius<=0){ IJ.showMessage("maximum drift between frames must be greater than zero"); consistent = false;}
                            if (kernel_w<=0) {IJ.showMessage("window length must greater than zero"); consistent=false;}
                            if (sigmax_PSF<1.0) {IJ.showMessage("sigmaX of PSF must be larger than 1"); consistent=false;}	
                            if (sigmay_PSF<1.0) {IJ.showMessage("sigmaY of PSF must be larger than 1"); consistent=false;}	
                            	
                    } while (!consistent && !checkonce);
                    if(!consistent) return;
                }
		if (catch_radius*2>spacing) IJ.log("maximum drift between frames is greater than half of spacing, may fail tracing the pillars"); 
				
		process();
	}


	public boolean IsNotNaN(double a, double b, double c, double d){
		if(Double.isNaN(a)) return false;
		if(Double.isNaN(b)) return false;
		if(Double.isNaN(c)) return false;
		if(Double.isNaN(d)) return false;
		return true;
	}

//support nan double
private double avg(double[] values) {
		int len = values.length;
		int num = 0;
	    double avg = 0;//values[0];
	    for(int k=0; k<len; k++) {                	
    		double v = values[k];
    		if(!Double.isNaN(v)){
            	avg += v;
            	num++;
    		}
	    }
	    
	    avg = num>0? avg/num : Double.NaN;
	    return avg;
} 

private double var(double[] values)
{
		int len = values.length;
		double a = avg(values);
		int num = 0;
		double var = 0;
		for(int k=0; k<len; k++){
			double v = values[k];
			if(!Double.isNaN(v)){				
				double da = v-a;
				var += da*da;
				num++;
			}
		}
		
		var = num>0? var/num : Double.NaN;
		return var;	
}

private void writeFileHeader_ver1(FileChannel out, int nframes)  throws IOException {        
        int header_size = Integer.SIZE*5 + Double.SIZE*7 + Byte.SIZE*4;
        ByteBuffer header_buffer = ByteBuffer.allocate(header_size);
        header_buffer.putInt(fileversion1);
        header_buffer.putInt(nrois);
        header_buffer.putInt(nframes);
        header_buffer.putDouble(lattice);   
        header_buffer.putDouble(diameter);  
        header_buffer.putDouble(spacing);
        header_buffer.putDouble(oblique);
        header_buffer.putDouble(grid_angle);
        header_buffer.putDouble(sigmax_PSF);
        header_buffer.putDouble(catch_radius);
        header_buffer.putInt(kernel_w);
        header_buffer.putInt(box_constrian_R);
        //header_buffer.putInt(num_start_points);
        header_buffer.put(dark_pillars?(byte)1:0);
        //header_buffer.put(use_minimum_std?(byte)1:0);
        header_buffer.put(apply_mean_rank_filter?(byte)1:0);
        header_buffer.put(use_metric_CG?(byte)1:0);
        header_buffer.put(image_enhancer_plugin!=null?(byte)1:0);        
        header_buffer.flip();
        out.write(header_buffer);                    
        if(image_enhancer_plugin!=null) writePSFImage(out, image_enhancer_plugin);        
}

void writeFloatImage(FileChannel out, double[] pixels)  throws IOException {
        int len = pixels.length;
        ByteBuffer buff = ByteBuffer.allocate(Double.SIZE*len); 
        for(int i=0; i<len; i++) buff.putDouble(pixels[i]);        
        buff.flip();
        out.write(buff);
}

void writePSFImage(FileChannel out, CrossCorreclation_Plugin cc)  throws IOException {
        int header_size = Integer.SIZE*3 + Double.SIZE*1 + Byte.SIZE*3;
        ByteBuffer header_buffer = ByteBuffer.allocate(header_size);
        header_buffer.put(cc.getMode()?(byte)1:0);
        header_buffer.put(cc.isUsingGaussianPSF()?(byte)1:0);
        header_buffer.put(cc.isDarkObject()?(byte)1:0);
        header_buffer.putInt(cc.getGaussianRaidus());
        header_buffer.putDouble(cc.getGaussianSigma());        
        ImageProcessor ip = cc.getPSF().getProcessor();
        int w = ip.getWidth();
        int h = ip.getHeight();
        header_buffer.putInt(w);
        header_buffer.putInt(h);
        header_buffer.flip();
        out.write(header_buffer);    
        
        int len = w*h;
        ByteBuffer buff = ByteBuffer.allocate(Double.SIZE*len); 
        for(int i=0; i<len; i++) buff.putDouble(ip.getf(i));        
        buff.flip();
        out.write(buff);
}

void writeCentroids(FileChannel out, int[] x, int[] y)  throws IOException {
        int len = x.length;
        ByteBuffer buff = ByteBuffer.allocate(Integer.SIZE*len*2); 
        for(int i=0; i<len; i++){
            buff.putInt(x[i]);
            buff.putInt(y[i]);
        }        
        buff.flip();
        out.write(buff);
}

void writePoints(FileChannel out, double[][] x, double[][] y)  throws IOException {
    int cens = x.length;
    int frames = x[0].length;
    for(int f=0; f<frames; f++) {
        ByteBuffer buff = ByteBuffer.allocate(Double.SIZE*2 * cens); 
        for(int c=0; c<cens; c++) {	
            buff.putDouble(x[c][f]);
            buff.putDouble(y[c][f]);                                                        
        }        
        buff.flip();
        out.write(buff);        
        //IJ.showProgress(c+1, cens);    
    }        
}

void writeBooleanArray(DataOutputStream out, boolean[] array)  throws IOException {
        int len = array.length;
        for(int i=0; i<len; i++) out.writeBoolean(array[i]);
} 

void writeBooleanArray(FileChannel out, boolean[][] array)  throws IOException {
        int w = array.length;
        int h = array[0].length;
        for(int i=0; i<w; i++){
            ByteBuffer buff = ByteBuffer.allocate(Byte.SIZE*h); 
            for(int j=0; j<h; j++) buff.put(array[i][j]?(byte)1:0);            
            buff.flip();
            out.write(buff);
        }
}

void writeBooleanArray(FileChannel out, boolean[] array)  throws IOException {
        int len = array.length;
        ByteBuffer buff = ByteBuffer.allocate(Byte.SIZE*len); 
        for(int i=0; i<len; i++){
            byte a = array[i] ? (byte)1 : 0;
            buff.put(a);
        }
        buff.flip();
        out.write(buff);
} 

void writeFloatImage2Text(PrintWriter out, double[] pixels, int ncols)  throws IOException {    
    int len = pixels.length;
    for (int i=0; i<len; i++) {
            if(i>0 && i%ncols == 0) out.println();   				       				
            out.print((pixels[i] + "	"));
    }
    out.println();                
}

void writeFloatPoints2Text(PrintWriter out, double[][] x, double[][] y)  throws IOException {    
    int cens = x.length;
    int frames = x[0].length;
    for(int f=0; f<frames; f++) {        
        for(int c=0; c<cens; c++) out.print((x[c][f] + "	"+ y[c][f]+"	"));        
        out.println();
    }        
}

void writeBooleanArray2Text(PrintWriter out, boolean[] array)  throws IOException {   		
    int len = array.length;
    for (int i=0; i<len; i++) {
            out.println("	" + (i+1) + "	" + array[i]);
    }                    
}

int indexof(String[] items,int len,  String item)
{
    int index = -1;
    for(int i=0; i<len; i++){
            if(item.equals(items[i])){
                    index = i;
                    break;
            }
    }
    return index;
}

boolean showDialog() {	

	int[] wList = WindowManager.getIDList();	
	if (wList==null || wList.length<1) {
		IJ.showMessage("Pillar Tracking", "1 channel time series are required");
		return false;
	}

	int wlistlen = wList.length;
	
	ArrayList<String> titles_list = new ArrayList<String>();
	for (int i=0; i<wlistlen; i++) {
		ImagePlus imp = WindowManager.getImage(wList[i]);
		//if(imp!=null && imp.getNDimensions()==3) titles_list.add(imp.getTitle());
                //if(imp!=null) titles_list.add(imp.getTitle());
		if(imp!=null && imp.getNDimensions()<=3) titles_list.add(imp.getTitle());
	}

	if(titles_list.isEmpty()){
		IJ.showMessage("Pillar Tracking", "1 channel time series are required");	
		return false;
	}
	
	String titles[] = titles_list.toArray(new String[0]);
	
	GenericDialogPlus gd = new GenericDialogPlus("Pillar Tracking");
	gd.addChoice("ONLY one channel time series", titles, titles[0]);

	gd.addMessage("----PILLAR AND GRID PROFILE----");
	gd.addNumericField("pilalr_spacing in pixel:", spacing, 1);
	gd.addNumericField("grid_oblique angle( [-45~45] degree)", oblique, 2);	
        gd.addNumericField("grid_shape angle( (0~90] degree)", grid_angle, 2);	
	gd.addNumericField("sigma of Gaussian PSF in pixel:", sigmax_PSF, 1);
	
        gd.addMessage("----TRACKING AND OPTIMIZATION ----");	
	gd.addNumericField("search_window_width in pixel:", kernel_w, 0);
	gd.addNumericField("constrian_radius in pixel:", box_constrian_R, 0);
	gd.addNumericField("maximum drift between frames in pixel:", catch_radius, 2);
	gd.addMessage("----Multi-Threading for Acceleration");
        
	numThread = ThreadUtil.getNbCpus();
	gd.addNumericField("number_of_threads:", numThread, 0);		

	gd.addMessage("----Plot Control for Super Resolved Image, Only Effective When Plotting");	
	gd.addNumericField("Zoom-in factor in super resolved image:", vector_zoom, 0);	
	gd.addCheckbox("dark pillars?", dark_pillars);
	gd.addCheckbox("use_center_of_mass?", use_metric_CG);		
        gd.addCheckbox("use_mean_filter to find maxima?", apply_mean_rank_filter);
	gd.addCheckbox("plotting movement? (require large memory)",plotting); 
	gd.addCheckbox("output debug information?",debug); 
		
	String lastname = OpenDialog.getDefaultDirectory();//OpenDialog.getLastDirectory();
	if(lastname.isEmpty()) lastname = OpenDialog.getLastDirectory();
	output_fname = lastname + "pillar-tracks.bin";
	int str_len = output_fname.length();
	if(str_len<50) str_len = 50;
	gd.addFileField("output file name:", output_fname, str_len);
	
	gd.showDialog();
	if (gd.wasCanceled()) return false;
	String choice_name1 = gd.getNextChoice();
	experiment = WindowManager.getImage(choice_name1);

	
        spacing=(double)gd.getNextNumber();
	oblique=(double)gd.getNextNumber();
        grid_angle=(double)gd.getNextNumber();
	sigmax_PSF=(double)gd.getNextNumber();	
	sigmay_PSF=sigmax_PSF;//use the same sigma with sigmax
	
	//kernel_w = (int)((Math.round(spacing)/2)*2 + 1);
	kernel_w=(int)gd.getNextNumber();
	if(kernel_w%2==0){
		kernel_w++;
		//IJ.showMessage("Warning", "The window width is not valid, which must be odd number");
		//return false;
	}
	box_constrian_R=(int)gd.getNextNumber();
	
	//IJ.log("kernel_w=" + kernel_w + "	searching radius=" + box_constrian_R);
	catch_radius = (double)gd.getNextNumber();
	
	numThread = (int)gd.getNextNumber();
	vector_zoom=(int)gd.getNextNumber();
	dark_pillars = gd.getNextBoolean();
	use_metric_CG = gd.getNextBoolean();        
        if(use_metric_CG) localization_algorithm = localization_algorithm_CG;
        else localization_algorithm = localization_algorithm_Levmar;
        
        apply_mean_rank_filter = gd.getNextBoolean();	
	plotting = gd.getNextBoolean();	
	debug = gd.getNextBoolean();	
	output_fname = gd.getNextString();
	
	return true;
}

}
