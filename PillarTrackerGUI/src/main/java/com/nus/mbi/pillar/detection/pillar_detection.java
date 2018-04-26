package com.nus.mbi.pillar.detection;

//Xu Xiaochun @ Mechnobiology Institute, Singapore
//contact:mbixxc@nus.edu.sg
//update history:
//2015-07-09: add the javadoc
//2015-06-20: Todo: add 2D label for the pillars. 
//2015-06-12: add the try-all for auto threshold and auto local threshold, when showing in the original image, don't show table and threshold image
//2015-06-09: make a point selection use different thresholding metholds,
//			  such as global thrsholding, auto thrsholding, and auto local thresholding

import ij.*;
import ij.io.*;
import ij.plugin.*;
import ij.plugin.filter.*;
import ij.measure.*;
import ij.process.*;
import ij.gui.*;
import ij.gui.GenericDialog;
import java.util.ArrayList;

/**
 * @author Xiaochun Xu <mbixxc@nus.edu.sg>
 * @version 1.4
 * @since 2015-07-09
 */
public class pillar_detection implements PlugIn {

//measurment image stack
private ImagePlus experiment;
private ImageStack experiment_s;

//ROI 
//RoiManager rmanager;
Roi [] rois;
int nrois;

// the kernels and their dimensions
private double auto_local_radius = 21;
private double auto_local_para1 = 0;
private double auto_local_para2 = 0;

//calculation method
//0=>global thresholding, 
//1=>auto thresholding
//2=>auto local thresholding
private int calc_method = 0;

private String method_auto_threshold = "Default";
private String method_auto_local_threshold = "";
private double global_threshold = 20000;

private double min_area = 30;
private double max_area = 500;
private double min_circ = 0.5;
private double max_circ = 1.0;

private boolean dark_pillars = true;
private boolean otsu_force_8bits = false;
private boolean debug = false;
private boolean showInOrigin = false;

private String str_tryall = "[Try all]";

/**
 * Set up the parameters used for pillar detection, then use {@link #process()}
 * 
 * @param ip 				time-series image stack 
 * @param dark				dark pillar?
 * @param creat				creat the selections in original image?
 * @return 					true if successfuly setup
 */
public boolean setup(ImagePlus ip, boolean dark, boolean creat)
{
	boolean suc = true;
	if(ip!=null)experiment = ip;
	else suc = false;
	dark_pillars = dark;
	showInOrigin = creat;
	return suc;
}

/**
 * Set to Global Threshold
 * 
 * @param threshold			global threshold
 */
public void setGlobalThreshold(double threshold)
{
	calc_method = 0;
	global_threshold = threshold;
}

/**
 * Set to Auto Threshold.
 * please see: http://fiji.sc/Auto_Threshold.
 * The method can be one of {@code "[Try all]","Default", "Huang", "Intermodes", "IsoData", "Li", "MaxEntropy","Mean", "MinError(I)", "Minimum","Moments","Otsu", "Percentile", "RenyiEntropy", "Shanbhag", "Triangle","Yen"};
 * @param method			auto threshold methold
 */
public void setAutoThreshold(String method, boolean force_8bits)
{
	calc_method = 1;
	method_auto_threshold = method;
	otsu_force_8bits = force_8bits;
}

/**
 * Set to Auto Local Threshold.
 * please see: http://fiji.sc/Auto_Local_Threshold.
 * The method can be one of {@code "[Try all]","Bernsen", "Contrast", "Mean", "Median", "MidGrey", "Niblack","Otsu", "Phansalkar", "Sauvola"};
 * @param method			auto local threshold threshold
 * @param radius			radius used for the current method
 * @param para1				para1 used for the current method
 * @param para2				para2 used for the current method
 * 
 */
public void setAutoLocalThreshold(String method, double radius, double para1, double para2)
{
	calc_method = 2;
	method_auto_local_threshold = method;
	auto_local_radius = radius;
	auto_local_para1 = para1;
	auto_local_para2 = para2;
}

/**
 * Set Partical Filter.
 * please see: http://fiji.sc/Particle_Analysis
 * @param min_area			minimum area
 * @param max_area			maximum area
 * @param min_circ			minimum circulrity
 * @param max_area			maximum circulrity
 */
public void setParticalFilter(double min_area, double max_area, double min_circ, double max_circ){
	this.min_area = min_area;
	this.max_area = max_area;
	this.min_circ = min_circ;
	this.max_circ = max_circ;
}

/**
 *  process pillars detection for the current time-series image stack,
 *  after setup {@link #setup(ImagePlus ip, boolean dark, boolean creat)}
 */
public void process(){
		if(experiment==null) return;
		int width = experiment.getWidth();
		int height = experiment.getHeight();		
		experiment_s = experiment.getStack();
				
		IJ.log("channels="+experiment.getNChannels()+"\nframes="+experiment.getNFrames()+"\nslices="+experiment.getNSlices());
		if (experiment.isHyperStack()) IJ.log("is hyperstack\n");
		IJ.log("stack equivalent slices="+experiment_s.getSize());
		experiment.changes = false;	
		
		FileInfo fi = experiment.getOriginalFileInfo();
		if(fi != null){
			String info = fi.directory + fi.fileName;
			IJ.log("path=" + info);
		}
		
		int sliceN = Math.max(experiment.getNSlices(),experiment.getNFrames());
		
		//if (roimanager_null)
		{
			boolean suc_load_roi = false;			
			//experiment.deleteRoi();
			Roi point_roi= experiment.getRoi();
			if(point_roi==null) suc_load_roi = false; //else if(!point_roi.isArea() && !point_roi.isLine()){			
			else if(!point_roi.isArea()){
				//IJ.log("points has been selected, clear the selctions!");					
				//FloatPolygon polygon_roi = point_roi.getFloatPolygon();				
				//nrois = polygon_roi.npoints; 				
				//if(nrois<1) suc_load_roi = false; 
				experiment.deleteRoi();				
				point_roi=null;
				//experiment.updateAndDraw();
			}
			else{
				suc_load_roi = false;
			}
			
			if(!suc_load_roi) {
				//creat the point-like selection	
				ImageProcessor slice = experiment.getProcessor();	
				ImagePlus slice_imp = duplicateImage(slice);
				byte[] binaryImg = new byte[width*height];
				boolean dowhite = !dark_pillars;
				boolean tryall = false;
				if(calc_method == 1) //auto thresholding
				{	
					ImagePlus thrshold_imgp = otsu_force_8bits? duplicateImage(slice.convertToByteProcessor(true)) : duplicateImage(slice);													
					String options = "method=" + method_auto_threshold;
					if(dowhite) options  = options + " white";
					options = options + " show";
					IJ.log("Auto Threshold->" + options);
					IJ.run(thrshold_imgp,"Auto Threshold",options);
					while (IJ.macroRunning()) {
						// Nothing needed here. Just wait for macro to end.
					}	
							
					if(method_auto_threshold == str_tryall) tryall = true;
					else{
						ImageProcessor thresholdImage = thrshold_imgp.getProcessor();
						for(int x=0; x<width; x++) for(int y=0; y<height; y++) binaryImg[y*width + x] = (byte)(thresholdImage.getf(x,y)>0 ? 255 : 0);															
					}
				}
				else if(calc_method == 2) //auto local thresholding
				{					
					ImagePlus thrsholdImage8bits = duplicateImage(slice.convertToByteProcessor(true));
					//thrsholdImage8bits.show();
					String options = "method="+method_auto_local_threshold+" radius="+auto_local_radius+" parameter_1="+auto_local_para1+" parameter_2="+auto_local_para2;					
					if(dowhite) options  = options + " white";
					IJ.run(thrsholdImage8bits,"Auto Local Threshold",options);
					IJ.log("Auto Local Threshold->"+options);
					
					while (IJ.macroRunning()) {
						// Nothing needed here. Just wait for macro to end.
					}
					
					if(method_auto_local_threshold == str_tryall) tryall = true;
					else{
						ImageProcessor thresholdImage8bits = thrsholdImage8bits.getProcessor();
						for(int x=0; x<width; x++) for(int y=0; y<height; y++) binaryImg[y*width + x] = (byte) (thresholdImage8bits.getf(x,y)>0 ? 255 : 0);	
					}		
				}
				else //global thrsholding
				{					
					ImageProcessor thrshold_img = otsu_force_8bits? slice.convertToByteProcessor(true) : slice;			
					double threshold = global_threshold;//dowhite ? slice.getMinThreshold() : slice.getMaxThreshold();										
					for(int x=0; x<width; x++) for(int y=0; y<height; y++){
						if(dowhite) binaryImg[y*width + x] = (byte)(thrshold_img.getf(x,y)>threshold ? 255 : 0);
						else binaryImg[y*width + x] = (byte)(thrshold_img.getf(x,y)>threshold ? 0 : 255);
					}
					
					IJ.log("Global Threshold->" + threshold);
				}
				//thrshold_imgp.show();
				if(tryall) return;
				
				ByteProcessor thresholdImage = new ByteProcessor(width, height, binaryImg);
				ImagePlus thrshold_imgp = new ImagePlus("thresholding for the first frame", thresholdImage);
				//thrshold_imgp.show();								
				
				thresholdImage.setThreshold(128,255,0);
				ResultsTable rt = new ResultsTable();
				ParticleAnalyzer.setFontColor("yellow");
				ParticleAnalyzer pa = new ParticleAnalyzer(ParticleAnalyzer.EXCLUDE_EDGE_PARTICLES|ParticleAnalyzer.SHOW_OVERLAY_MASKS,Measurements.CENTROID,rt,min_area,max_area,min_circ,max_circ);				
				//pa.analyze(thrshold_imgp);
				pa.analyze(slice_imp,thresholdImage);
				//pa.run(thresholdImage);
				//rt.show("centroids for the first frame");

				int ncol = rt. getLastColumn();
				int rtsize = rt.size();
				IJ.log("rtsize=" + rtsize + "	 ncol = " + ncol);				
				
				if(rtsize>0){					
					if(!showInOrigin){
						thrshold_imgp.show();
						rt.show("centroids for the first frame");			
					}
					double[] cetroidsx = rt.getColumnAsDoubles(ncol-1);
					double[] cetroidsy = rt.getColumnAsDoubles(ncol);
				
					PointRoi mulit_proi = null;						
					nrois = 0;
					for(int i=0; i<rtsize; i++){
						double cx = cetroidsx[i];
						double cy = cetroidsy[i];	
						int roisx = (int)Math.round(cx);
						int roisy = (int)Math.round(cy);									
						if(point_roi==null || point_roi.contains(roisx,roisy)){																
							if(mulit_proi == null) mulit_proi = new PointRoi(cx, cy);
							else mulit_proi = mulit_proi.addPoint(cx, cy);
							//IJ.log("x=" + cx +"	 y=" + cy);
							nrois++;
						}
					}
					
					if(mulit_proi != null){
						Overlay overlay_draw = new Overlay();					
						//mulit_proi.setFillColor(Color.CYAN);
						overlay_draw.add(mulit_proi);											
						thrshold_imgp.setOverlay(overlay_draw);
						thrshold_imgp.updateAndDraw();	
						IJ.log("npoints=" + mulit_proi.getNCoordinates());
						if(showInOrigin){
							//experiment.setOverlay(overlay_draw);
							experiment.deleteRoi();
							experiment.setRoi(mulit_proi);
							//experiment.updateAndDraw();			
						}
						
					}		
				}
				/*
				if(point_roi.isArea()){ //area selection					
					IJ.log("process area selection");
					
				}*/
				else{ //no selection
					//IJ.showMessage("There is no any Point selection"); 
					IJ.log("There is no any Point selection!");					
				}
				
			}
			
		}	
						
		if(debug){
			IJ.log("This software is implemented by Xu Xiaochun in MechnoBilogy Institute,Sinapore(MBI)");
			IJ.log("Bug report to: xuxiaochun27@gmail.com");					
		}		
}

/**
 * run the plugin with dialog.
 * Marco command example: 
 * {@code run("Creat selections", "only=75fps_cell2_dish2_large_X27.tif choose=[Global Threshold] global=24000 auto=[[Try all]] auto=[[Try all]] radius=21 para_1=0 para_2=0 min_area=30 max_area=500 min_circularity=0.50 max_circularity=1 dark creat");}
 * 
 * @param arg	not used 
 *
 */
public void run(String arg) { 
		boolean checkonce = arg.contains("checkonce");
		boolean consistent;	
		do {
			consistent=true;
			if (!showDialog()) return;
			int channels_orig=experiment.getNChannels();
			if (channels_orig!=1) {IJ.showMessage("stack must have ONLY one channels"); consistent=false;}			
			if (calc_method == 2 && auto_local_radius<0) {IJ.showMessage("the radius for auto local thresholding must greater than zero"); consistent=false;}
		} while (!consistent && !checkonce);
		
		if(!consistent) return;

		process();
	}	

	private ImagePlus duplicateImage(ByteProcessor iProcessor){
		int w=iProcessor.getWidth();
		int h=iProcessor.getHeight();
		ImagePlus iPlus=NewImage.createByteImage("Image", w, h, 1, NewImage.FILL_BLACK);
		ImageProcessor imageProcessor=iPlus.getProcessor();
		imageProcessor.copyBits(iProcessor, 0,0, Blitter.COPY);
		return iPlus;
	}

	 private ImagePlus duplicateImage(ShortProcessor iProcessor){
		int w=iProcessor.getWidth();
		int h=iProcessor.getHeight();
		ImagePlus iPlus=NewImage.createShortImage("Image", w, h, 1, NewImage.FILL_BLACK);
		ImageProcessor imageProcessor=iPlus.getProcessor();
		imageProcessor.copyBits(iProcessor, 0,0, Blitter.COPY);
		return iPlus;
	}

	private ImagePlus duplicateImage(ImageProcessor iProcessor){
		int w=iProcessor.getWidth();
		int h=iProcessor.getHeight();
		ImagePlus iPlus=NewImage.createShortImage("Image", w, h, 1, NewImage.FILL_BLACK);
		ImageProcessor imageProcessor=iPlus.getProcessor();
		imageProcessor.copyBits(iProcessor, 0,0, Blitter.COPY);
		return iPlus;
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
		IJ.showMessage("Pillar Detection", "1 channel time series are required");
		return false;
	}

	int wlistlen = wList.length;
	
	ArrayList<String> titles_list = new ArrayList<String>();
	for (int i=0; i<wlistlen; i++) {
		ImagePlus imp = WindowManager.getImage(wList[i]);
		if(imp!=null) titles_list.add(imp.getTitle());
	}

	if(titles_list.isEmpty()){
		IJ.showMessage("Pillar Detection", "1 channel time series are required");	
		return false;
	}
	
	String titles[] = titles_list.toArray(new String[0]);
	GenericDialog gd = new GenericDialog("Pillar Detection");
	gd.addChoice("ONLY one channel time series", titles, titles[0]);
	
	gd.addMessage("--Segmentation Process Control");
	String[] items = {"Global Threshold", "Auto Threshold", "Auto Local Threshold"}; 
	//gd.addRadioButtonGroup("Choose an algorithm:",items, 3, 1, items[calc_method]); 	
	gd.addChoice("Choose an algorithm:",items, items[calc_method]); 

	gd.addMessage("(1)Parameters @ Global Threshold");
	gd.addNumericField("   global threshold:", global_threshold, 1);

	gd.addMessage("(2)Parameters @ Auto Threshold");
	//String[] methods_auto = AutoThresholder.getMethods();
	String[] methods_auto = {str_tryall,"Default", "Huang", "Intermodes", "IsoData", "Li", "MaxEntropy","Mean", "MinError(I)", "Minimum","Moments","Otsu", "Percentile", "RenyiEntropy", "Shanbhag", "Triangle","Yen"};
	gd.addChoice("Auto Threshold Algorithm:",methods_auto, methods_auto[0]); 

	gd.addMessage("(3)Parameters @ Auto Local Threshold");			
	String[] methods_auto_local = {str_tryall,"Bernsen", "Contrast", "Mean", "Median", "MidGrey", "Niblack","Otsu", "Phansalkar", "Sauvola"};
	gd.addChoice("Auto Local Threshold Algorithm:",methods_auto_local, methods_auto_local[0]); 
	gd.addNumericField("radius:", auto_local_radius, 0);
	gd.addNumericField("para_1:", auto_local_para1, 0);
	gd.addNumericField("para_2:", auto_local_para2, 0);

	gd.addMessage("--Partical Shape Control");	
	gd.addNumericField("min_area:", min_area, 0);
	gd.addNumericField("max_area:", max_area, 0);
	gd.addNumericField("min_circularity:", min_circ, 2);
	gd.addNumericField("max_circularity:", max_circ, 2);
	
	gd.addCheckbox("dark pillars", dark_pillars);	
	gd.addCheckbox("convert to 8 bit when segmentaion?", otsu_force_8bits);
	gd.addCheckbox("output debug information?",debug); 	
	gd.addCheckbox("creat the selections in original image?",showInOrigin); 	

	gd.showDialog();
	if (gd.wasCanceled()) return false;
	String choice_name1 = gd.getNextChoice();
	experiment = WindowManager.getImage(choice_name1);

	//String method = gd.getNextRadioButton();		
	//calc_method = indexof(items,3,method);
	calc_method = gd.getNextChoiceIndex();	
	if(calc_method<0 ||calc_method>2){
		IJ.showMessage("Pillar Detection", "The choosen mothod is not valid:" + calc_method);
		return false;
	}
	
	global_threshold = gd.getNextNumber();
	method_auto_threshold = gd.getNextChoice();	
	method_auto_local_threshold = gd.getNextChoice();	
	
	auto_local_radius = gd.getNextNumber();
	auto_local_para1 = gd.getNextNumber();
	auto_local_para2 = gd.getNextNumber();

	min_area = gd.getNextNumber();
	max_area = gd.getNextNumber();
	min_circ = gd.getNextNumber();
	max_circ = gd.getNextNumber();	
	
	dark_pillars = gd.getNextBoolean();
	otsu_force_8bits = gd.getNextBoolean();
	debug = gd.getNextBoolean();
	showInOrigin = gd.getNextBoolean();
	
	return true;
}

}

