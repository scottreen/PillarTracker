package com.nus.mbi.pillar.detection;

import com.nus.mbi.pillar.tracker.NearestNeigbour;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.gui.ShapeRoi;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.util.ArrayList;
import java.util.List;
import com.nus.mbi.pillar.solver.MinimumErrorSolver;

/**
 *
 * @author xiaochun
 */
public class PillarDetector {
        private ImagePlus experiment;		
	private int channels_orig; // # of used channels
	
	private int kernelw;		// the kernels and their dimensions
	private double[] kernel; 	// used when a single kernel is applied
	private int[] kernel_dx;
	private int[] kernel_dy;
	
	//private double diameter=16.;
	private double spacing=23.;
	private double psf_sigma=7;
	public double precentage = 5; //within (0~100)
	private boolean dark_pillars = true;
	
	private double oblique = 0;
	private double oblique_tol = 20;        
       
        private double grid_angle = 90;
        private boolean is_square_grid = true;
        
	int width = 0;
	int height = 0;	
	int NN=0;	
	int kernel_psf_size = 0;
	ShapeRoi expand_roi = null; 
	com.nus.mbi.pillar.grid.SearchAreaCreation area_creator = new com.nus.mbi.pillar.grid.SearchAreaCreation(); 
        
        double[] quality   = null;
        double[] error     = null;
        double[] amplitude = null;
        double[] offset    = null;
	boolean[] maxl     = null; // is an amplitude pixel a maximum?
        int[][] centroisXY = null;
        
        boolean prepared = false;
        
	//output
	public boolean[] valid_map = null;
        public int num_local_maxima = 0;
        public int num_detected_maxima = 0;

	public boolean isPrepared()
	{
		return prepared;
	}
        
        public double[] getKernel(){
            return kernel;
        }
        
        public int[] getKernelDX(){
            return kernel_dx;
        }
        
        public int[] getKernelDY(){
            return kernel_dy;
        }
        
        public int getKernelWidth(){
            return kernelw;
        }
        
        public int getKernelSize(){
            return kernel_psf_size;
        }
        
        public ByteProcessor getValidMap(){                        
            //if(!prepared) return null;
            ByteProcessor bp = new ByteProcessor(width, height);
            
            if(valid_map==null){;
                byte[] pixels = (byte[])bp.getPixels();
                for(int i=0; i<NN; i++) pixels[i] = valid_map[i] ? (byte) 255 : 0;            
            }
            
            return bp;
        }
	
	/**
	 * Set up the parameters used for pillar detection, then use {@link #process()}
	 * 
	 * @param ip 				time-series image stack 
	 * @param oblique_angle		grid oblique angle( -45~45 degree)	 
	 * @param spacing			pillar spacing in pixel
	 * @param sigma				Gaussian sigma in pixel
	 * @param dark				dark pillar?
	 
	 * @return 					true if successfuly setup
	 */
	public boolean setup(ImagePlus ip, double oblique_angle, double spacing, double sigma, boolean dark)
	{
		boolean suc = true;
		if(ip!=null) experiment = ip;
		else suc = false;
	
		if(sigma>0) psf_sigma = sigma;
		else suc = false;		
	
		if(spacing>0) this.spacing = spacing;
		else suc = false;

		if(oblique_angle>-45 && oblique_angle<45) this.oblique = oblique_angle;
		else suc = false;		
		
		dark_pillars = dark;
	
//		IJ.log("set up successful? " + suc);	
//		IJ.log("sigma :" + psf_sigma);
//		IJ.log("pillar spacing:" + this.spacing);	
//		IJ.log("dark pillars:" + dark_pillars);	
//		IJ.log("oblique angle:" + oblique);	
		
		return suc;
	}

	/**
	 * Set up the grid oblique tolerance used for pillar detection, then use {@link #process()}
	 * 	 
	 * @param oblique_tol		grid oblique angle toleracne( 0~45 degree)
	 * @return 					true if successfuly setup
	 */
	public boolean setObliqueTolerance(double oblique_tol)
	{		
		boolean suc = true;
		if(oblique_tol>0 && oblique_tol<45) this.oblique_tol = oblique_tol;
		else suc = false;		
		
		IJ.log("oblique angle tol:" + this.oblique_tol);	
		
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
                    is_square_grid = (grid_angle==90);
                }
		else suc = false;		
		
//		IJ.log("grid angle:" + this.grid_angle);	
		
		return suc;
	}

	public void prepare(){
		if(experiment==null) return;
		width = experiment.getWidth();
		height = experiment.getHeight();	
		NN=width*height;	
//		ImageStack experiment_s = experiment.getStack();
				
//		IJ.log("channels="+experiment.getNChannels()+"\nframes="+experiment.getNFrames()+"\nslices="+experiment.getNSlices());
//		if (experiment.isHyperStack()) IJ.log("is hyperstack\n");
//		IJ.log("stack equivalent slices="+experiment_s.getSize());
		experiment.changes = false;	
		
//		FileInfo fi = experiment.getOriginalFileInfo();
//		if(fi != null){
//			String info = fi.directory + fi.fileName;
//			IJ.log("path=" + info);
//		}
		
		int sliceN = Math.max(experiment.getNSlices(),experiment.getNFrames());
		kernelw = (int)Math.round(spacing);
		if (kernelw%2 == 0) kernelw++; // round up to oddnumber
		int kernel_size  = kernelw*kernelw;
		int kernal_halfw = kernelw/2;

		kernel=new double[kernel_size];
		kernel_dx=new int[kernel_size];
		kernel_dy=new int[kernel_size];
		
		//int kernel_psf_size = make_Gaussian_PSF(psf_sigma, diameter, kernel, kernel_dx, kernel_dy, kernal_halfw);
		kernel_psf_size = make_Gaussian_PSF(psf_sigma, kernel, kernel_dx, kernel_dy, kernal_halfw);
//		IJ.log("kernel_psf_size = " + kernel_psf_size);

		Roi work_roi=experiment.getRoi();
		Roi limit_roi = new Roi(kernal_halfw,kernal_halfw,width-kernelw,height-kernelw);
		if(work_roi == null || false==work_roi.isArea()) work_roi = limit_roi;			
		
		if(work_roi != null && work_roi.isArea()){
			java.awt.Rectangle rect=work_roi.getBounds();			
			int lowx  = rect.x - kernal_halfw; 
			int width_new = rect.width + 2*kernal_halfw;
			int lowy  = rect.y - kernal_halfw; 
			int height_new = rect.height + 2*kernal_halfw;

			if(lowx<0) lowx=0;
			if(lowy<0) lowy=0;
			if (width_new>width) width_new = width;
			if (height_new>height) height_new = height;			
                        
                        Roi limitRoi = new Roi(lowx,lowy,width_new,height_new);
                        
			ShapeRoi sr = new ShapeRoi(limitRoi);
                        ShapeRoi wr = new ShapeRoi(work_roi);
                        expand_roi = sr.and(wr);
                                
		}
		else{
                    Roi limitRoi = new Roi(0,0,width,height);
                    expand_roi = new ShapeRoi(limitRoi);
		}	
                
                area_creator.create_search_area(spacing, oblique, grid_angle);                
                
		prepared = true;
	}
	
	/**
	 *  process pillars localization for the current time-series image stack,
	 *  after setup {@link #setup(ImagePlus ip, double diameter, double spacing, double sigma, boolean dark)}         
	 */	
	public void process(){		
		ImageProcessor ep = experiment.getProcessor();	
                process(ep, true);
	}        
        
        /**
	 *  process pillars localization for the current time-series image stack,
	 *  after setup {@link #setup(ImagePlus ip, double diameter, double spacing, double sigma, boolean dark)}
         *  @param ep		ImageProcessor
         *  @param use_minimum_error    pre-processing use the minimum error solver?
	 */	
	public void process(ImageProcessor ep, boolean use_minimum_error){		
		prepare();
		
		//valid_max     = detect(ep);
		centroisXY = detect_slow(ep, use_minimum_error);
                draw_maxima(centroisXY, experiment);
//		
//		//valid_max     = new boolean[NN]; // is an amplitude pixel a maximum?
//		PointRoi mulit_proi = null;						
//		int nrois = 0;
//		if(centroisXY!=null && centroisXY[0]!=null){
//			int numcent = centroisXY[0].length;
//			for(int c0=0; c0<numcent; c0++) {
//				int x = (int)Math.round(centroisXY[0][c0]);
//				int y = (int)Math.round(centroisXY[1][c0]);
//			//	valid_max[x+y*width]=true;
//				if(mulit_proi == null) mulit_proi = new PointRoi(x, y);
//				else mulit_proi = mulit_proi.addPoint((double)x,(double)y);
//				//IJ.log("x=" + cx +"	 y=" + cy);
//				nrois++;
//			}
//		}
//						
//		if(mulit_proi != null){
//			Overlay overlay_draw = new Overlay();					
//			overlay_draw.add(mulit_proi);	
//                        //num_detected_maxima = mulit_proi.getNCoordinates();
//			IJ.log("npoints=" + mulit_proi.getNCoordinates());
//			experiment.deleteRoi();
//			experiment.setRoi(mulit_proi);		
//		}				
	}     
        
        private static void draw_maxima(int[][] centroisXY, ImagePlus ip){
            //valid_max     = new boolean[NN]; // is an amplitude pixel a maximum?
		PointRoi mulit_proi = null;						
		int nrois = 0;
		if(centroisXY!=null && centroisXY[0]!=null){
			int numcent = centroisXY[0].length;
			for(int c0=0; c0<numcent; c0++) {
				int x = (int)Math.round(centroisXY[0][c0]);
				int y = (int)Math.round(centroisXY[1][c0]);
			//	valid_max[x+y*width]=true;
				if(mulit_proi == null) mulit_proi = new PointRoi(x, y);
				else mulit_proi = mulit_proi.addPoint((double)x,(double)y);
				//IJ.log("x=" + cx +"	 y=" + cy);
				nrois++;
			}
		}
						
		if(mulit_proi != null){
			Overlay overlay_draw = new Overlay();					
			overlay_draw.add(mulit_proi);	
                        //num_detected_maxima = mulit_proi.getNCoordinates();
			IJ.log("npoints=" + mulit_proi.getNCoordinates());
			ip.deleteRoi();
			ip.setRoi(mulit_proi);		
		}				
        }
        
        public void drawLocalMaxima(ImageProcessor ip, double threshold){
                int[][] centroids = find_local_maximas(ip, threshold);
		draw_maxima(centroids, experiment);
        }
        
        public int getNumLocalMaximas(ImageProcessor ip, double threshold){
            find_local_maximas(ip, threshold);
            return num_local_maxima;
        }
        
        public void plot(){
            FloatProcessor errfp = new FloatProcessor(width, height);
            FloatProcessor ampfp = new FloatProcessor(width, height);
            FloatProcessor offsetfp = new FloatProcessor(width, height);
            FloatProcessor qualityfp = new FloatProcessor(width, height);
            ByteProcessor maxminp = new ByteProcessor(width, height);
            int num_maxmin = 0;
            for(int k=0; k<NN; k++) {
                    errfp.setf(k,(float)error[k]);
                    ampfp.setf(k,(float)amplitude[k]);				
                    offsetfp.setf(k,(float)offset[k]);	
                    if(error[k]>0) qualityfp.setf(k,(float)(amplitude[k]/Math.sqrt(error[k])));		
                    if(maxl[k]){
                            maxminp.set(k, 255);
                            num_maxmin++;
                    }
            }
            ImagePlus amplitudei = new ImagePlus("amplitude ", ampfp);
            amplitudei.show();		

            ImagePlus offseti = new ImagePlus("offset ", offsetfp);
            offseti.show();

            ImagePlus errori = new ImagePlus("error ", errfp);
            errori.show();

            ImagePlus qualityi = new ImagePlus("quality ", qualityfp);
            qualityi.show();

            String title = dark_pillars ? "local minimum " : "local maximum ";
            ImagePlus maxmini = new ImagePlus(title + num_maxmin, maxminp);
            maxmini.show();

            ByteProcessor valid_maxminp = new ByteProcessor(width, height);
            int numcent = centroisXY[0].length;
            if(numcent>0){
                    for(int c0=0; c0<numcent; c0++) {
                            int x = (int)Math.round(centroisXY[0][c0]);
                            int y = (int)Math.round(centroisXY[1][c0]);
                            valid_maxminp.set(x,y,255);
                    }
            }

            title = dark_pillars ? "valid local minimum " : "valid local maximum ";
            ImagePlus valid_maxmini = new ImagePlus(title + numcent, valid_maxminp);
            valid_maxmini.show();
        }
        
	public int[][] detect_slow(ImageProcessor ep, boolean use_minimum_error){
		if(!prepared) return null;
		
		quality   = new double[NN];
		error     = new double[NN];
		amplitude = new double[NN];
		offset    = new double[NN];
		maxl     = new boolean[NN]; // is an amplitude pixel a maximum?
		
		for(int i=0; i<NN; i++) {
			quality[i]=0.; error[i]=0.;
			amplitude[i]=0.; maxl[i]=false;
		}

		//ImageProcessor ep = experiment.getProcessor();
		double[] slice = new double[NN];		
		for(int k=0; k<NN; k++) slice[k] = ep.getf(k);
                int kernal_halfw = kernelw/2;
                if(use_minimum_error){               
                    MinimumErrorSolver.MinimumError(slice,amplitude,offset,error, width, height, kernel, kernel_dx, kernel_dy, kernel_psf_size, kernal_halfw);                                
//                    for(int k=0; k<NN; k++){
//                            if(error[k]>0) quality[k] = (float)(amplitude[k]/Math.sqrt(error[k]));			
//                    }
                }
                else{
                    for(int k=0; k<NN; k++) amplitude[k] = slice[k];			
                }

		num_local_maxima = local_max_min(maxl, amplitude, amplitude, width, height, expand_roi, 0, dark_pillars);	
		//IJ.log("num of local maxmin=" + num);
		double[] score_map = getScoreMap(maxl, width, height);
                double orphan_score = getOrphanScore();
                deleteOrphan(maxl, score_map, orphan_score);
                //int[][] centroidsXY = getPoints(maxl, width, height);
                //int[][] centroidsXY = check_local_maxmin(maxl, score_map, width, height, spacing*1.5, kernal_halfw);	
		check_local_maxmin(maxl, amplitude,score_map, width, height, spacing*1.5, kernal_halfw, dark_pillars);	
                
                score_map = getScoreMap(maxl, width, height);                
                deleteOrphan(maxl, score_map, orphan_score);
		check_local_maxmin(maxl, score_map, width, height, spacing*1.5, kernal_halfw);                     
                //if(is_square_grid) check_angle(maxl, score_map, width, height, spacing*2.0, oblique, oblique_tol);	
                check_local_maxmin(maxl, amplitude, width, height, spacing*1.5, kernal_halfw, dark_pillars);	
                score_map = getScoreMap(maxl, width, height);        
                deleteOrphan(maxl, score_map, orphan_score);
                int[][] centroidsXY = getPoints(maxl, width, height);
                num_detected_maxima = centroidsXY[0].length;
		//seting the output max map;
		valid_map     = new boolean[NN]; 
		for(int i=0; i<NN; i++) valid_map[i] = maxl[i];
		
		return centroidsXY;
	}
        
        public int[][] detect(ImageProcessor ep){
		if(!prepared) return null;
		
//		quality   = new double[NN];
		error     = new double[NN];
		amplitude = new double[NN];
		offset    = new double[NN];
		maxl     = new boolean[NN]; // is an amplitude pixel a maximum?
		
		for(int i=0; i<NN; i++) {
//			quality[i]=0.;
                        error[i]=0.;
			amplitude[i]=0.; maxl[i]=false;
		}

		//ImageProcessor ep = experiment.getProcessor();
		double[] slice = new double[NN];		
		for(int k=0; k<NN; k++) slice[k] = ep.getf(k);
		int kernal_halfw = kernelw/2;
		MinimumErrorSolver.MinimumError(slice,amplitude,offset,error, width, height, kernel, kernel_dx, kernel_dy, kernel_psf_size, kernal_halfw);
//		for(int k=0; k<NN; k++){
//			if(error[k]>0) quality[k] = (float)(amplitude[k]/Math.sqrt(error[k]));			
//		}

		num_local_maxima = local_max_min(maxl, amplitude, amplitude, width, height, expand_roi, 0, dark_pillars);	
		//IJ.log("num of local maxmin=" + num);
		double[] score_map = getScoreMap(maxl, width, height);
                
                //deleteOrphan(maxl, score_map, getOrphanScore());
                //int[][] centroidsXY = getPoints(maxl, width, height);
                //int[][] centroidsXY = check_local_maxmin(maxl, score_map, width, height, spacing*1.5, kernal_halfw);	
		int[][] centroidsXY = check_local_maxmin(maxl, amplitude,score_map, width, height, spacing*1.5, kernal_halfw, dark_pillars);	
                
                //score_map = getScoreMap(maxl, width, height);                
                //deleteOrphan(maxl, score_map, getOrphanScore());
		//centroidsXY = check_local_maxmin(maxl, score_map, width, height, spacing*1.5, kernal_halfw);
                num_detected_maxima = centroidsXY[0].length;
		//seting the output max map;
		valid_map     = new boolean[NN]; 
		for(int i=0; i<NN; i++) valid_map[i] = maxl[i];
		
		return centroidsXY;
	}
        
        public int[][] find_local_maximas(ImageProcessor ep, double threshold){
		if(!prepared) return null;
		boolean[] maxl     = new boolean[NN]; // is an amplitude pixel a maximum?		
		double[] slice = new double[NN];		
		for(int k=0; k<NN; k++) slice[k] = ep.getf(k);
                
		int kernal_halfw = kernelw/2;
		num_local_maxima = local_max_min(maxl, slice, slice, width, height, expand_roi, threshold, dark_pillars);
                //int num = local_max_min(maxl, slice, width, height, expand_roi, dark_pillars);	
		//IJ.log("num of local maxmin=" + num);
                //int[][] centroidsXY = getPoints(maxl, width, height);
		double[] score_map = getScoreMap(maxl, width, height);
                
		//int[][] centroidsXY = check_local_maxmin(maxl, score_map, width, height, spacing*1.5, kernal_halfw);
                int[][] centroidsXY = check_local_maxmin(maxl, slice, score_map, width, height, spacing*1.5, kernal_halfw, dark_pillars);			
		
		return centroidsXY;
	}
        
        
	public int[][] detect_fast(ImageProcessor ep){
		if(!prepared) return null;
		double[] slice = new double[NN];		
		for(int k=0; k<NN; k++) slice[k] = ep.getf(k);
                return detect_fast(slice);
	}
        
        public int[][] detect_fast(double[] slice){
		if(!prepared) return null;
		boolean[] maxl     = new boolean[NN]; // is an amplitude pixel a maximum?	                
		int kernal_halfw = kernelw/2;
		//num_local_maxima = local_max_min(maxl, slice, slice, width, height, expand_roi, 0, dark_pillars);	
                num_local_maxima = local_max_min(maxl, slice, width, height, expand_roi, dark_pillars);	
                //int[][] centroidsXY = getPoints(maxl, width, height);                
		double[] score_map = getScoreMap(maxl, width, height);
                int[][] centroidsXY = check_local_maxmin(maxl, slice, score_map, width, height, spacing*1.5, kernal_halfw, dark_pillars);
                num_detected_maxima = centroidsXY[0].length;
		//seting the output max map;
//		valid_map     = new boolean[NN]; 
//		for(int i=0; i<NN; i++) valid_map[i] = maxl[i];
		valid_map = maxl;
		return centroidsXY;
	}        
        
        public int[][] detect_fastest(ImageProcessor ep){
		if(!prepared) return null;
		double[] slice = new double[NN];		
		for(int k=0; k<NN; k++) slice[k] = ep.getf(k);
                return detect_fastest(slice);
	}
        
        public int[][] detect_fastest(double[] slice){
		if(!prepared) return null;
		boolean[] maxl   = new boolean[NN]; // is an amplitude pixel a maximum?	                		
                num_local_maxima = local_max_min(maxl, slice, width, height, expand_roi, dark_pillars);	
                int[][] centroidsXY = getPoints(maxl, width, height);                
		num_detected_maxima = centroidsXY[0].length;
		//seting the output max map;
		//valid_map     = new boolean[NN]; 
                //System.arraycopy(maxl, 0, valid_map, 0, NN);		
                valid_map = maxl;
		return centroidsXY;
	}       
                
        public double[][] create_map(int[] x, int[] y, int map_width, int map_height){
            double[][] map = new double[map_height][map_width];
            for(int i=0; i<map_width; i++){
                for(int j=0; j<map_height; j++){
                    map[j][i] = 0;
                }
            }
            int npillars = x.length;
            
            for(int i=0; i<npillars; i++){
                int rxx = x[i];
                int ryy = y[i];
                if(rxx>=0 && ryy>=0 && rxx<map_width && ryy<map_height)
                    map[ryy][rxx] = i+1;
            }
            
            return map;
        }
        
        public double[] create_score_map(int[] x, int[] y, int map_width, int map_height, double[] scores){
            int img_size = map_height*map_width;
            double[] map = new double[img_size];
            for(int i=0; i<img_size; i++) map[i] = Double.MAX_VALUE;
            
            int npillars = x.length;
            
            for(int i=0; i<npillars; i++){
                int rxx = x[i];
                int ryy = y[i];
                if(rxx>=0 && ryy>=0 && rxx<map_width && ryy<map_height)
                    map[rxx+ryy*map_width] = scores[i];
            }
            
            return map;
        }
        
        private int[][] getPoints(boolean[] max, int width, int height){            
            int img_size = width*height;
            int n = 0;
            for(int i=0; i<img_size; i++) { if(max[i]) n++;}
            if(n<1) return null;
            
            int[][] centroidsXY=new int[2][n];		
            int cen=0;
            for(int x=0; x<width; x++) {
                for(int y=0; y<height; y++) {
                    if(max[x+y*width]) {
                        centroidsXY[0][cen]=x;
                        centroidsXY[1][cen]=y;
                        cen++;
                    }
                }
            }
            
            return centroidsXY;
        }
        
        private void deleteOrphan(boolean[] max, double[] score, double score_orphan){
            int img_size = max.length;
            for(int i=0; i<img_size; i++){
                if(max[i]){
                    if(score[i]<score_orphan){ }
                    else{
                        score[i] = Double.MAX_VALUE;
                        max[i] = false;
                    }
                }
            }
        }
        
        private double getScoreWall(){
            double r = area_creator.getOuterRadius();
            double score_wall = r*r;
            return score_wall;
        }
        
        private double getOrphanScore(){
            double score_wall = getScoreWall();
            return 4*score_wall;
        }
        
        private double[] getScoreMap(boolean[] max, int width, int height){
            double[] nx = area_creator.getNX();
            double[] ny = area_creator.getNY();
            double score_wall = getScoreWall();
            int[][] xy = getPoints(max, width, height);
            if(xy==null) return null;
            double[][] map = create_map(xy[0], xy[1], width, height);
            double[][] np = getScores(map, area_creator.getArea(), xy[0], xy[1], nx,ny, score_wall);
            return create_score_map(xy[0], xy[1], width, height, np[4]);
        }
        
        
        private double[][] getScores(double[][] map_image, double[][] search_area, int[] rx, int[] ry, double[] nx, double[] ny, double score_wall){
            int narea = 4;
            int npillars = rx.length;
            //double[] match_error = new double[npillars];
            double[][] np = new double[narea+1][npillars];
            double r2 = score_wall;//radius*radius;
            int map_w = map_image[0].length;
            int map_h = map_image.length;
            int area_size = search_area.length;
            for(int k=0; k<npillars; k++){           
                double x = rx[k];
                double y = ry[k];  
                double score = 0;
                for(int i=0; i<narea; i++){            
                    //% search the neighbours along x forward
                    double[] min_p_dis = search_neighbor(x,y,map_image,map_w,map_h,search_area,area_size,nx[i],ny[i]);
                    double min_p = min_p_dis[0];
                    double min_dis = min_p_dis[1];
                    np[i][k] = min_p;
                    if (min_p>0)
                        score = score + min_dis;
                    else
                        score = score + r2;
                }            
                //match_error[k] = score;
                np[narea][k] = score;
            }
            return np;
        }
    
    public double[] search_neighbor(double x, double y, double[][] map_image, int map_w, int map_h, double[][] area, int area_size, double nx, double ny){               
        double min_dis = -1;
        double min_p = -1;

        int xc = (int)Math.round(x + nx);
        int yc = (int)Math.round(y + ny);
        
        for(int k=0; k<area_size; k++){
           double d2 = area[k][0];       
           int dx = (int)area[k][1];
           int dy = (int)area[k][2];
           int xx = xc + dx;
           int yy = yc + dy;           
           if (xx>=0 && yy>=0 && xx<map_w && yy<map_h){
                double p = map_image[yy][xx];
                if(p>0){
                    double dis = d2;
                    if(min_dis<0 || min_dis>dis){
                        min_dis = dis;
                        min_p = p;
                    }
                }
           }
        }   
        
        double[] ret = {min_p, min_dis};
        return ret;
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

public static int[][] create_mask(int[] xc, int[] yc, int width, int height){
	int num=xc.length; //centroid counter per frame
        int[][] mask = new int[width][height];        
        for(int x=0; x<width; x++) for(int y=0; y<height; y++) mask[x][y] = -1;                
        for(int k=0; k<num; k++){
            int x = xc[k];
            int y = yc[k];
//            if(x>=0 && x<width && y>=0 && y<height)
            {
                mask[x][y] = k;
            }
        }
        return mask;
}

public static int search_neighbor(int x, int y, double radius, int[][] mask){
    int r = (int)Math.round(radius);
    double r2 = radius*radius;
    int min_d2 = Integer.MAX_VALUE;
    int min_nk = -1;
    for(int i=-r; i<=r; i++){
        int ny = y + i;
        int i2 = i*i;
        for(int j=-r; j<=r; j++){                
            int d2 = i2 + j*j;
            if(d2<=r2){ 
                int nx = x + j;
                int nk = mask[nx][ny];
                if(nk>=0){
                    if(d2<min_d2){
                        min_d2 = d2;
                        min_nk = nk;
                    }
                }
            }
        }            
    }
    
    return min_nk;
}

public static NearestNeigbour search_nearest_neighbors(int x, int y, double radius, int[][] mask){
    int r = (int)Math.ceil(radius);
    List<Integer> neigbors = new ArrayList();
    double r2 = radius*radius;
    int min_d2 = Integer.MAX_VALUE;
    //int min_nk = -1;
    for(int i=-r; i<=r; i++){
        int ny = y + i;
        int i2 = i*i;
        for(int j=-r; j<=r; j++){                
            int d2 = i2 + j*j;
            if(d2<=r2){ 
                int nx = x + j;
                int nk = mask[nx][ny];
                if(nk>=0){
                    if(d2<min_d2){
                        neigbors.clear();
                        neigbors.add(nk);
                        min_d2 = d2;
                        //min_nk = nk;
                    }
                    else if(d2==min_d2){
                        neigbors.add(nk);
                    }
                }
            }
        }            
    }
    
    return new NearestNeigbour(neigbors, min_d2);
}

public static int[][] cross_check_frames(int[] xc1, int[] yc1, int[] xc2, int[] yc2, int[][] I_mask_1, int[][] I_mask_2, double catch_radius)
{
    double r = catch_radius;        
    int num1 = xc1.length;
    int[] map12 = new int[num1];
    for(int k=0; k<num1; k++) map12[k] = -1;
    
    int num2 = xc2.length;
    int[] map21 = new int[num2];
    for(int k=0; k<num2; k++) map21[k] = -1; 
    
    boolean found_new = true;
    while(found_new){
        found_new = false;
        for(int k=0; k<num1; k++){
            if(map12[k]<0){
                NearestNeigbour nn1 = search_nearest_neighbors(xc1[k], yc1[k], r, I_mask_2);
                List<Integer> nb = nn1.neigbours;                
                if(nb.size()>0){
                    double shrink_radius = Math.sqrt(nn1.distance);
                    for(int i=0; i<nb.size(); i++){
                        int kk = nb.get(i);                        
                        NearestNeigbour nn2 = search_nearest_neighbors(xc2[kk], yc2[kk], shrink_radius, I_mask_1);
                        if(nn2.distance>=nn1.distance){                            
                            map12[k] = kk;
                            map21[kk] = k;
                            int x1 = xc1[k];
                            int y1 = yc1[k];
                            int x2 = xc2[kk];
                            int y2 = yc2[kk];
                            I_mask_1[x1][y1] = -1;
                            I_mask_2[x2][y2] = -1;
                            found_new = true;
                            break;
                        }
                        //else map12[k] = -1; 
                    }
                }
            }
        }
    }
    
    int[][] map = new int[2][];
    map[0] = map12;
    map[1] = map21;
    return map;
}

public static int[][] cross_check_frames_dead_lock(int[] xc1, int[] yc1, int[] xc2, int[] yc2, int[][] I_mask_1, int[][] I_mask_2, double catch_radius)
{
    double r = catch_radius;        
    int num1 = xc1.length;
    int[] map12 = new int[num1];
    for(int k=0; k<num1; k++) map12[k] = -1;
    //for(int k=0; k<num1; k++) map12[k] = search_neighbor(xc1[k], yc1[k], r, I_mask_2);
    
    int num2 = xc2.length;
    int[] map21 = new int[num2];
    for(int k=0; k<num2; k++) map21[k] = -1; 
//    for(int k=0; k<num1; k++){
//        int kk = map12[k];
//        if(kk>=0){
//            int min_nk = search_neighbor(xc2[kk], yc2[kk], r, I_mask_1);
//            if(min_nk == k){
//                map21[kk] = k;
//                int x1 = xc1[k];
//                int y1 = yc1[k];
//                int x2 = xc2[kk];
//                int y2 = yc2[kk];
//                I_mask_1[x1][y1] = -1;
//                I_mask_2[x2][y2] = -1;
//            }
//            else map12[k] = -1;            
//        }
//    }
    
    boolean found_new = true;
    while(found_new){
        found_new = false;
        for(int k=0; k<num1; k++){
            if(map12[k]<0){
                int kk = map12[k] = search_neighbor(xc1[k], yc1[k], r, I_mask_2);
                if(kk>=0){
                    int min_nk = search_neighbor(xc2[kk], yc2[kk], r, I_mask_1);
                    if(min_nk == k){
                        map21[kk] = k;
                        int x1 = xc1[k];
                        int y1 = yc1[k];
                        int x2 = xc2[kk];
                        int y2 = yc2[kk];
                        I_mask_1[x1][y1] = -1;
                        I_mask_2[x2][y2] = -1;
                        found_new = true;
                    }
                    else map12[k] = -1;       
                }
            }
        }
    }
    
    int[][] map = new int[2][];
    map[0] = map12;
    map[1] = map21;
    return map;
}

public static int[][] cross_check_frames_competitive(int[] xc1, int[] yc1, int[] xc2, int[] yc2, int[][] I_mask_1, int[][] I_mask_2, double catch_radius)
{
    double r = catch_radius;        
    int num1 = xc1.length;
    int[] map12 = new int[num1];
    for(int k=0; k<num1; k++) map12[k] = search_neighbor(xc1[k], yc1[k], r, I_mask_2);
    
    int num2 = xc2.length;
    int[] map21 = new int[num2];
    for(int k=0; k<num2; k++) map21[k] = -1; 
    for(int k=0; k<num1; k++){
        int kk = map12[k];
        if(kk>=0){
            int min_nk = search_neighbor(xc2[kk], yc2[kk], r, I_mask_1);
            if(min_nk == k) map21[kk] = k;
            else map12[k] = -1;            
        }
    }
    
    int[][] map = new int[2][];
    map[0] = map12;
    map[1] = map21;
    return map;
//    r = catch_radius;
//    nrange = -r:r;
//    
//    num1 = numel(xc1);
//    map12 = zeros(num1, 1);        
//    for i=1:num1
//        x = xc1(i);
//        y = yc1(i);
//        I_neigbor = I_mask_2(y+nrange,x+nrange);
//        [ny, nx, ni] = find(I_neigbor);
//        if(numel(ni)>0)
//            dis = (nx-r-1).^2 + (ny-r-1).^2;
//            [~, k] = min(dis);
//            map12(i) = ni(k);
//        end
//    end
//    
//    num2 = numel(xc2);
//    map21 = zeros(num2, 1);        
//    for j=1:num1
//        i = map12(j);
//        if(i>0)           
//            x = xc2(i);
//            y = yc2(i);
//            I_neigbor = I_mask_1(y+nrange,x+nrange);
//            [ny, nx, ni] = find(I_neigbor);
//            if(numel(ni)>0)
//                dis = (nx-r-1).^2 + (ny-r-1).^2;
//                [~, k] = min(dis);
//                jj = ni(k);                
//                if(jj == j)
//                    map21(i) = j;
//                else
//                    map12(j) = 0;
//                end
//            else
//                map12(j) = 0;
//            end
//        end
//    end
}


//int[][] check_angle(boolean[]maxl, double[] amp, int width, int height, double diameter, double angle, double tol, boolean dark_objects)
//{	
//		double a = -angle*Math.PI/180;
//		double sina = Math.sin(a);
//		double cosa = Math.cos(a);
//		
//		int circler=(int)Math.round(diameter/2);
//		int circlew=2*circler+1;
//		int[] circlex=new int[circlew*circlew];
//		int[] circley=new int[circlew*circlew];
//		int cic=0; // circle element counter
//		for(int y=-circler; y<=circler; y++) {
//			for(int x=-circler; x<=circler; x++) {
//				if((x*x+y*y<=circler*circler) && ((x!=0) || (y!=0))) {
//					double x0 = x*cosa - y*sina;
//					double y0 = x*sina + y*cosa;
//					double dx = Math.abs(x0);
//					double dy = Math.abs(y0);
//					if(dy>dx){
//						double t = dx;
//						dx = dy;
//						dy = t;
//					}
//					double aa = Math.atan(dy/dx)*180/Math.PI;
//					
//					if(aa>tol){
//						circlex[cic]=x;
//						circley[cic]=y;						
//						cic++;
//					}
//				}
//			}
//		}
//
//		int img_size = width*height;
//		boolean[] maxl_c = new boolean[img_size];
//		for(int i=0; i<img_size; i++) maxl_c[i] = maxl[i];
//                
//                for(int x=0; x<width; x++) {
//			for(int y=0; y<height; y++) {
//                                int p = x+y*width;
//				if(maxl_c[p]) {
//					for(int i=0; i<cic; i++) {
//						int x1=x+circlex[i];
//						int y1=y+circley[i];
//						if((x1>=0) && (x1<width) && (y1>=0) && (y1<height)) {
//                                                    int np = x1+y1*width;
//                                                    if(maxl_c[np]) {								
//                                                            if(dark_objects && amp[np]>amp[p]){
//                                                                maxl[np] = false;
//                                                            }
//                                                            else if(!dark_objects && amp[np]<amp[p]){
//                                                                maxl[np] = false;
//                                                            }
//							}
//						}
//					}
//				}			
//			}
//		}
//              
//		int centC=0; //centroid counter per frame
//		for(int k=0; k<width*height; k++)	if(maxl[k]) centC++;	
//		int[][] centroidsXY=new int[2][centC];			
//		centC=0; 
//		for(int x=0; x<width; x++) {
//			for(int y=0; y<height; y++) {
//				if(maxl[x+y*width]){
//				 	centroidsXY[0][centC]=x;
//	  				centroidsXY[1][centC]=y;				 	
//				 	centC++;
//				}
//			}
//		}		
//
//		return centroidsXY;
//}
//
//int[][] check_angle(boolean[]maxl, double[] score, int width, int height, double diameter, double angle, double tol)
//{	
//		double a = -angle*Math.PI/180;
//		double sina = Math.sin(a);
//		double cosa = Math.cos(a);
//		
//		int circler=(int)Math.round(diameter/2);
//		int circlew=2*circler+1;
//		int[] circlex=new int[circlew*circlew];
//		int[] circley=new int[circlew*circlew];
//		int cic=0; // circle element counter
//		for(int y=-circler; y<=circler; y++) {
//			for(int x=-circler; x<=circler; x++) {
//				if((x*x+y*y<=circler*circler) && ((x!=0) || (y!=0))) {
//					double x0 = x*cosa - y*sina;
//					double y0 = x*sina + y*cosa;
//					double dx = Math.abs(x0);
//					double dy = Math.abs(y0);
//					if(dy>dx){
//						double t = dx;
//						dx = dy;
//						dy = t;
//					}
//					double aa = Math.atan(dy/dx)*180/Math.PI;
//					
//					if(aa>tol){
//						circlex[cic]=x;
//						circley[cic]=y;						
//						cic++;
//					}
//				}
//			}
//		}
//
//		int img_size = width*height;
//		boolean[] maxl_c = new boolean[img_size];
//		for(int i=0; i<img_size; i++) maxl_c[i] = maxl[i];
//                
//                for(int x=0; x<width; x++) {
//			for(int y=0; y<height; y++) {
//                                int p = x+y*width;
//				if(maxl_c[p]) {
//					for(int i=0; i<cic; i++) {
//						int x1=x+circlex[i];
//						int y1=y+circley[i];
//						if((x1>=0) && (x1<width) && (y1>=0) && (y1<height)) {
//                                                    int np = x1+y1*width;
//                                                    if(maxl_c[np] && score[np]>score[p]) maxl[np] = false;
//						}
//					}
//				}			
//			}
//		}
//                
//		int centC=0; //centroid counter per frame
//		for(int k=0; k<width*height; k++)	if(maxl[k]) centC++;	
//		int[][] centroidsXY=new int[2][centC];			
//		centC=0; 
//		for(int x=0; x<width; x++) {
//			for(int y=0; y<height; y++) {
//				if(maxl[x+y*width]){
//				 	centroidsXY[0][centC]=x;
//	  				centroidsXY[1][centC]=y;				 	
//				 	centC++;
//				}
//			}
//		}		
//
//		return centroidsXY;
//}


int[][] check_local_maxmin(boolean[]maxl, double[] score, int width, int height, double diameter, int margin)
{
		// contract them to point-like maxima
		int marginsr = margin;
		
		int circler=(int)Math.round(diameter*0.5);
		int circlew=2*circler+1;
		int[] circlex=new int[circlew*circlew];
		int[] circley=new int[circlew*circlew];
		int cic=0; // circle element counter
		for(int y=-circler; y<=circler; y++) {
			for(int x=-circler; x<=circler; x++) {
				if((x*x+y*y<=circler*circler) && ((x!=0) || (y!=0))) {
					circlex[cic]=x;
					circley[cic]=y;
					cic++;
				}
			}
		}

		int img_size = width*height;
		boolean[] maxl_c = new boolean[img_size];
		for(int i=0; i<img_size; i++) maxl_c[i] = maxl[i];
		// delete lesser maxima within the diameter of a pillar
		
                for(int x=0; x<width; x++) {
			for(int y=0; y<height; y++) {
                                int p = x+y*width;
				if(maxl_c[p]) {
					for(int i=0; i<cic; i++) {
						int x1=x+circlex[i];
						int y1=y+circley[i];
						if((x1>=0) && (x1<width) && (y1>=0) && (y1<height)) {
                                                    int np = x1+y1*width;
                                                    if(maxl_c[np]) {
								if(score[np]>score[p]){
                                                                    maxl[np] = false;
                                                                }   
							}
						}
					}
				}			
			}
		}
                
		int centC=0; //centroid counter per frame
		for(int k=0; k<width*height; k++)	if(maxl[k]) centC++;		
		int[] cx=new int[centC];
		int[] cy=new int[centC];
		centC=0; 
		for(int x=marginsr; x<width-marginsr; x++) {
			for(int y=marginsr; y<height-marginsr; y++) {
				if(maxl[x+y*width]){
				 	cx[centC] = x;
				 	cy[centC] = y;
				 	centC++;
				}
			}
		}
                
                int[][] centroidsXY=new int[2][];		
                centroidsXY[0] = cx;
                centroidsXY[1] = cy;
		
//		// delete the colliding centroids
//		boolean[] active=new boolean[centC];
//		for(int c=0; c<centC; c++) active[c]=true;
//		
//		int cent2nd=0;
//		for(int c0=0; c0<centC; c0++) {if(active[c0]) {cent2nd++;}}
//		int[][] centroidsXY=new int[2][cent2nd];		
//		for(int i=0; i<img_size; i++) maxl[i] = false;
//		int cen=0;
//		for(int c0=0; c0<centC; c0++) {
//			if(active[c0]) {
//				int x = cx[c0];
//				int y = cy[c0];	  			
//				maxl[x+y*width]=true;
//				centroidsXY[0][cen]=x;
//	  			centroidsXY[1][cen]=y;
//	  			cen++;
//			}
//		}
		
		return centroidsXY;
}

int[][] check_local_maxmin(boolean[]maxl, double[]amp, int width, int height, double diameter, int margin, boolean dark)
{
		// contract them to point-like maxima
		int marginsr = margin;
		
		int circler=(int)Math.round(diameter*0.5);
		int circlew=2*circler+1;
		int[] circlex=new int[circlew*circlew];
		int[] circley=new int[circlew*circlew];
		int cic=0; // circle element counter
		for(int y=-circler; y<=circler; y++) {
			for(int x=-circler; x<=circler; x++) {
				if((x*x+y*y<=circler*circler) && ((x!=0) || (y!=0))) {
					circlex[cic]=x;
					circley[cic]=y;
					cic++;
				}
			}
		}

		int img_size = width*height;
		boolean[] maxl_c = new boolean[img_size];
		for(int i=0; i<img_size; i++) maxl_c[i] = maxl[i];
		// delete lesser maxima within the diameter of a pillar
		for(int x=0; x<width; x++) {
			for(int y=0; y<height; y++) {
                                int p = x+y*width;
				if(maxl_c[p]) {
					for(int i=0; i<cic; i++) {
						int x1=x+circlex[i];
						int y1=y+circley[i];
						if((x1>=0) && (x1<width) && (y1>=0) && (y1<height)) {
                                                    int np = x1+y1*width;
                                                    if(maxl_c[np]) {								
                                                            if(dark && amp[np]>amp[p]){
                                                                maxl[np] = false;
                                                            }
                                                            else if(!dark && amp[np]<amp[p]){
                                                                maxl[np] = false;
                                                            }
							}
						}
					}
				}			
			}
		}
		
		int centC=0; //centroid counter per frame
		for(int k=0; k<width*height; k++)	if(maxl[k]) centC++;		
		int[] cx=new int[centC];
		int[] cy=new int[centC];
		centC=0; 
		for(int x=marginsr; x<width-marginsr; x++) {
			for(int y=marginsr; y<height-marginsr; y++) {
				if(maxl[x+y*width]){
				 	cx[centC] = x;
				 	cy[centC] = y;
				 	centC++;
				}
			}
		}
                
                int[][] centroidsXY=new int[2][];		
                centroidsXY[0] = cx;
                centroidsXY[1] = cy;
		
//		// delete the colliding centroids
//		boolean[] active=new boolean[centC];
//		for(int c=0; c<centC; c++) active[c]=true;
//		
//		int cent2nd=0;
//		for(int c0=0; c0<centC; c0++) {if(active[c0]) {cent2nd++;}}
//		int[][] centroidsXY=new int[2][cent2nd];		
//		for(int i=0; i<img_size; i++) maxl[i] = false;
//		int cen=0;
//		for(int c0=0; c0<centC; c0++) {
//			if(active[c0]) {
//				int x = cx[c0];
//				int y = cy[c0];	  			
//				maxl[x+y*width]=true;
//				centroidsXY[0][cen]=x;
//	  			centroidsXY[1][cen]=y;
//	  			cen++;
//			}
//		}
		
		return centroidsXY;
}

int[][] check_local_maxmin(boolean[]maxl, double[]amp, double[] score, int width, int height, double diameter, int margin, boolean dark)
{
		// contract them to point-like maxima
		int marginsr = margin;
		
		int circler=(int)Math.round(diameter*0.5);
		int circlew=2*circler+1;
		int[] circlex=new int[circlew*circlew];
		int[] circley=new int[circlew*circlew];
		int cic=0; // circle element counter
		for(int y=-circler; y<=circler; y++) {
			for(int x=-circler; x<=circler; x++) {
				if((x*x+y*y<=circler*circler) && ((x!=0) || (y!=0))) {
					circlex[cic]=x;
					circley[cic]=y;
					cic++;
				}
			}
		}

		int img_size = width*height;
		boolean[] maxl_c = new boolean[img_size];
		for(int i=0; i<img_size; i++) maxl_c[i] = maxl[i];
		// delete lesser maxima within the diameter of a pillar
		for(int x=0; x<width; x++) {
			for(int y=0; y<height; y++) {
                                int p = x+y*width;
				if(maxl_c[p]) {
					for(int i=0; i<cic; i++) {
						int x1=x+circlex[i];
						int y1=y+circley[i];
						if((x1>=0) && (x1<width) && (y1>=0) && (y1<height)) {
                                                    int np = x1+y1*width;
                                                    if(maxl_c[np]) {
								if(score[np]>score[p]){
                                                                    if(dark && amp[np]>=amp[p]){
                                                                        maxl[np] = false;
                                                                    }
                                                                    else if(!dark && amp[np]<=amp[p]){
                                                                        maxl[np] = false;
                                                                    }
                                                                }   
							}
						}
					}
				}			
			}
		}
		
		int centC=0; //centroid counter per frame
		for(int k=0; k<width*height; k++)	if(maxl[k]) centC++;		
		int[] cx=new int[centC];
		int[] cy=new int[centC];
		centC=0; 
		for(int x=marginsr; x<width-marginsr; x++) {
			for(int y=marginsr; y<height-marginsr; y++) {
				if(maxl[x+y*width]){
				 	cx[centC] = x;
				 	cy[centC] = y;
				 	centC++;
				}
			}
		}
		
                int[][] centroidsXY=new int[2][];		
                centroidsXY[0] = cx;
                centroidsXY[1] = cy;
                
//		// delete the colliding centroids
//		boolean[] active=new boolean[centC];
//		for(int c=0; c<centC; c++) active[c]=true;
//		
//		int cent2nd=0;
//		for(int c0=0; c0<centC; c0++) {if(active[c0]) {cent2nd++;}}
//		int[][] centroidsXY=new int[2][cent2nd];		
//		for(int i=0; i<img_size; i++) maxl[i] = false;
//		int cen=0;
//		for(int c0=0; c0<centC; c0++) {
//			if(active[c0]) {
//				int x = cx[c0];
//				int y = cy[c0];	  			
//				maxl[x+y*width]=true;
//				centroidsXY[0][cen]=x;
//	  			centroidsXY[1][cen]=y;
//	  			cen++;
//			}
//		}
		
		return centroidsXY;
}
	
	int local_maximum(boolean[] max, double[] amplitute_map, double[] quality_map, int width, int height, Roi work_roi, double threshold){
		int yy1=1; 
		int yyh=height-1; 
		int xx1=1;
		int xxw=width-1;  
		
		if(work_roi != null && work_roi.isArea()) {		
			java.awt.Rectangle rect=work_roi.getBounds();
			int lowx=rect.x; 
			int highx=lowx+rect.width;
			int lowy=rect.y; 
			int highy=lowy+rect.height;
			
			if(lowy>1) yy1=lowy;
			if (highy<yyh) yyh=highy;
			if(lowx>1) xx1=lowx;
			if(highx<xxw) xxw=highx;
		}
		else{
			work_roi = new Roi(0,0,width,height);	
		}
	
		for(int ii=0; ii<width*height; ii++) max[ii]=false;
		int[] p_neighbors = {-width-1, -width, -width+1, 1, width+1, width, width-1, -1};
		int neighbors = p_neighbors.length;
		int num=0;
		for(int y=yy1; y<yyh;y++) {
			for(int x=xx1; x<xxw;x++) {
				if(work_roi.contains(x,y)){
					 int p = y*width + x;			 
					 double t=quality_map[p];
					 double e=amplitute_map[p];
					 if(e>threshold){ //only one channel is used
					 	boolean is_local_maxima = true;
                                                int num_nbs = 0;
					 	for(int i=0; i<neighbors; i++){
					 		int np = p + p_neighbors[i];
                                                        double quality_nb = quality_map[np];
					 		if(t<quality_nb){
					 			is_local_maxima = false;
					 			break;
					 		}else if(t>quality_nb) num_nbs++;		 		
					 	}
					 	
					 	if(is_local_maxima && num_nbs>4){
					 		num++;
                                                        max[p]=true;
					 	}
					 	else max[p]=false;
					 }		
				}	 
			}
		}
		return num;
	}
	
	int local_minimum(boolean[] max, double[] amplitute_map, double[] quality_map, int width, int height, Roi work_roi, double threshold) {
		int yy1=1; 
		int yyh=height-1; 
		int xx1=1;
		int xxw=width-1;  
		
		if(work_roi != null && work_roi.isArea()) {		
			java.awt.Rectangle rect=work_roi.getBounds();
			int lowx=rect.x; 
			int highx=lowx+rect.width;
			int lowy=rect.y; 
			int highy=lowy+rect.height;
			
			if(lowy>1) yy1=lowy;
			if (highy<yyh) yyh=highy;
			if(lowx>1) xx1=lowx;
			if(highx<xxw) xxw=highx;
		}
		else{
			work_roi = new Roi(0,0,width,height);	
		}
		
		//int yy1=1; if(lowy>1) yy1=lowy;
		//int yyh=height-1; if (highy<yyh) yyh=highy;
		//int xx1=1; if(lowx>1) xx1=lowx;
		//int xxw=width-1; if(highx<xxw) xxw=highx;
	
		for(int ii=0; ii<width*height; ii++) max[ii]=false;
		int[] p_neighbors = {-width-1, -width, -width+1, 1, width+1, width, width-1, -1};
		int neighbors = p_neighbors.length;
		int num=0;
		for(int y=yy1; y<yyh;y++) {
			for(int x=xx1; x<xxw;x++) {
				 if(work_roi.contains(x,y)){
				 	 int p = y*width + x;			 
					 double t=quality_map[p];
					 double e=amplitute_map[p];
					 if(e<threshold){ //only one channel is used
					 	boolean is_local_minimum = true;
                                                int num_nbs = 0;
					 	for(int i=0; i<neighbors; i++){
					 		int np = p + p_neighbors[i];
                                                        double quality_nb = quality_map[np];
					 		if(t>quality_nb){
					 			is_local_minimum = false;
					 			break;
                                                        }else if(t<quality_nb) num_nbs++;
					 	}
					 	
					 	if(is_local_minimum && num_nbs>4) {
                                                    num++;
                                                    max[p]=true;
					 	}
					 	else max[p]=false;
					 }		
				 }	 
			}
		}
		return num;
	}
	
	public int local_max_min(boolean[] max, double[] amplitute_map, double[] quality_map, int width, int height, Roi work_roi, double threshold, boolean dark_pillars)
	{	
		return dark_pillars ? local_minimum(max,amplitute_map, quality_map, width, height, work_roi, threshold) : local_maximum(max,amplitute_map, quality_map, width, height, work_roi, threshold);	
	}


	static int local_maximum(boolean[] max, double[] amplitute_map, int width, int height, Roi work_roi){
		int yy1=1; 
		int yyh=height-1; 
		int xx1=1;
		int xxw=width-1;  
		
		if(work_roi != null && work_roi.isArea()) {		
			java.awt.Rectangle rect=work_roi.getBounds();
			int lowx=rect.x; 
			int highx=lowx+rect.width;
			int lowy=rect.y; 
			int highy=lowy+rect.height;
			
			if(lowy>1) yy1=lowy;
			if (highy<yyh) yyh=highy;
			if(lowx>1) xx1=lowx;
			if(highx<xxw) xxw=highx;
		}
		else{
			work_roi = new Roi(0,0,width,height);	
		}
	
		for(int ii=0; ii<width*height; ii++) max[ii]=false;
		int[] p_neighbors = {-width-1, -width, -width+1, 1, width+1, width, width-1, -1};
		int neighbors = p_neighbors.length;
		int num=0;
		for(int y=yy1; y<yyh;y++) {
			for(int x=xx1; x<xxw;x++) {
				if(work_roi.contains(x,y)){
					int p = y*width + x;			 
					double t=amplitute_map[p];
					boolean is_local_maxima = true;
                                        int num_nbs = 0;
					for(int i=0; i<neighbors; i++){
						int np = p + p_neighbors[i];
                                                double amp_np = amplitute_map[np];
						if(t<amp_np){
							is_local_maxima = false;
							break;
						}
                                                else if(t>amp_np) num_nbs++;
					}
					
					if(is_local_maxima && num_nbs>4) {
						num++;
						max[p]=true;
					}
					else max[p]=false;						
				}	 
			}
		}
		return num;
	}
	
	static int local_minimum(boolean[] max, double[] amplitute_map, int width, int height, Roi work_roi) {
		int yy1=1; 
		int yyh=height-1; 
		int xx1=1;
		int xxw=width-1;  
		
		if(work_roi != null && work_roi.isArea()) {		
			java.awt.Rectangle rect=work_roi.getBounds();
			int lowx=rect.x; 
			int highx=lowx+rect.width;
			int lowy=rect.y; 
			int highy=lowy+rect.height;
			
			if(lowy>1) yy1=lowy;
			if (highy<yyh) yyh=highy;
			if(lowx>1) xx1=lowx;
			if(highx<xxw) xxw=highx;
		}
		else{
			work_roi = new Roi(0,0,width,height);	
		}		
			
		for(int ii=0; ii<width*height; ii++) max[ii]=false;
		int[] p_neighbors = {-width-1, -width, -width+1, 1, width+1, width, width-1, -1};
		int neighbors = p_neighbors.length;
		int num=0;
		for(int y=yy1; y<yyh;y++) {
			for(int x=xx1; x<xxw;x++) {
				 if(work_roi.contains(x,y)){
					int p = y*width + x;			 
					double t=amplitute_map[p];
					boolean is_local_minimum = true;
                                        int num_nbs = 0;
					for(int i=0; i<neighbors; i++){
						int np = p + p_neighbors[i];
                                                double amp_np = amplitute_map[np];
						if(t>amp_np){
							is_local_minimum = false;
							break;
						}
                                                else if(t<amp_np) num_nbs++;
					}
					
					if(is_local_minimum && num_nbs>4) {
						num++;
						max[p]=true;
					}
					else max[p]=false;					 		
				 }	 
			}
		}
		return num;
	}
	
	public static int local_max_min(boolean[] max, double[] amplitute_map,int width, int height, Roi work_roi, boolean dark_pillars)
	{	
		return dark_pillars ? local_minimum(max,amplitute_map, width, height, work_roi) : local_maximum(max,amplitute_map, width, height, work_roi);	
	}
	
//	// local maxima centroid counters
//	private int cmx;
//	private int cmy;
//	private int cc;
//	
//	private void paint_max(boolean[] maxl, int width, int height, int x, int y) {
//		cmx=0; cmy=0; cc=0;
//		paint_mc(maxl,width, height,x,y);
//		if(cc>0) {
//			int x1=cmx/cc;
//			int y1=cmy/cc;
//			maxl[x1+y1*width]=true;
//		}
//	}
//	
//	private void paint_mc(boolean[] maxl,int width, int height, int x, int y) {
//		if((x>=0) && (x<width) && (y>=0) && (y<height)) {
//			int ind=x+y*width;
//			if(maxl[ind]) {
//				maxl[ind]=false;
//				cc++;
//				cmx+=x;
//				cmy+=y;
//				paint_mc(maxl,width,height,x-1,y  );
//				paint_mc(maxl,width,height,x+1,y  );
//				paint_mc(maxl,width,height,x-1,y-1);
//				paint_mc(maxl,width,height,x  ,y-1);
//				paint_mc(maxl,width,height,x+1,y-1);
//				paint_mc(maxl,width,height,x-1,y+1);
//				paint_mc(maxl,width,height,x  ,y+1);
//				paint_mc(maxl,width,height,x+1,y+1);
//			}
//		}
//	}
	/*
	private double[] calc_sigma_psf(double[] psf, int size)
	{
		double[] array = new double[3];
		
		double sigmaY = 0.0;
		double sigmaY2 = 0.0;
	
		for(int i=0; i<size; i++) 
		{
			sigmaY += psf[i];
			sigmaY2 += psf[i] * psf[i];
		}
		
		double a_scale = size/(sigmaY2*size-sigmaY*sigmaY);
	
		array[0] = sigmaY;
		array[1] = sigmaY2;
		array[2] = a_scale;
	
		return array;
	}
	
	private void MinimumError(double[] imagef, double[] amp_map, double[] offset_map, double[] err_map, int width, int height, double[] psf, int[] psf_dx, int[] psf_dy, int psf_size, int psf_halfw)
	{
		double[] array = calc_sigma_psf(psf, psf_size);
		double sigmaY = array[0];
		double sigmaY2 = array[1];
		double a_scale = array[2];
		
		//IJ.log("sigmaY=" + sigmaY);
		//IJ.log("sigmaY2=" + sigmaY2);
		//IJ.log("a_scale=" + a_scale);
		
		for(int x=psf_halfw; x<width-psf_halfw; x++)
		{
			for(int y=psf_halfw; y<height-psf_halfw; y++)
			{
				int p= x + y*width;
				
				double sigmaXY = 0.0;
				double sigmaX = 0.0;
				double sigmaX2 = 0.0;
				
				for(int i=0; i<psf_size; i++)
				{
					int x1 = x + psf_dx[i];
					int y1 = y + psf_dy[i];
		
					double mi = imagef[x1 + y1*width];
					double dy = psf[i];
		
					sigmaX  += mi;//sigmaX + mi;
					sigmaX2 += mi*mi;//sigmaX2 + mi*mi;
					sigmaXY += dy*mi;//sigmaXY + //psf[i]*mi;
				}
	
				double lxy = sigmaXY-sigmaX*sigmaY/psf_size;		
				double amp = lxy*a_scale;//sigmaXY/sigmaY2;//sigmaXY;//
				double o = (sigmaX - amp*sigmaY)/psf_size;
				double err = 0;
				
				double absdet = Math.abs(a_scale);
				if(absdet<1e-20){			
					amp = 0;
					o= 0;
					err = sigmaX2;			
				}				
				else{
					 err = amp*amp*sigmaY2 + psf_size*o*o + sigmaX2 + 2*amp*o*sigmaY - 2*amp*sigmaXY - 2*o*sigmaX;
				}
							
				amp_map[p] = amp;
				err_map[p] = err;
				offset_map[p] = o;			
			}
		}
	
		// do the margins only
		int margin=psf_halfw;
		int kernelw = psf_halfw*2 + 1;
		int kernelw2=kernelw*kernelw;
		double[] mask=new double[kernelw*kernelw];
		for(int x=0;x<width;x++) {
			for(int y=0;y<height;y++) {
				if((x<margin) || (x>=width-margin) || (y<margin) | (y>=height-margin)) {
					int i=0;
	  				for(int y1=y-margin; y1<=y+margin; y1++) {
			  			for(int x1=x-margin; x1<=x+margin; x1++) {
		  					//if((x1>=0) && (x1<width) && (y1>=0) && (y1<height)) mask[i]=slice[x1+width*y1]; else mask[i]=0;
		  					int x2=x1;
		  					if(x1<0) x2=0; else if (x1>=width) x2=width-1;
		  					int y2=y1;
		  					if(y1<0) y2=0; else if (y1>=height) y2=height-1;
		  					mask[i]=imagef[x2+width*y2];
		  					i++;
		  				}
		  			}
		  			double sigmaX=0.;
		  			for(i=0; i<kernelw2; i++) sigmaX=sigmaX+mask[i];
		  			double sigmaX2=0.;
		  			for(i=0; i<kernelw2; i++) sigmaX2=sigmaX2+mask[i]*mask[i];
		  			// calculate the variable solver parameters
				  	double sigmaXY=0.;
				  	for(i=0; i<kernelw2; i++) sigmaXY=sigmaXY+kernel[i]*mask[i];
				  	double amp=(sigmaXY-(1./(double)kernelw2)*sigmaX*sigmaY)*a_scale;
				  	double offset  =(sigmaX-amp*sigmaY)/(double)kernelw2;
				  	double err=amp*amp*sigmaY2+(double)kernelw2*offset*offset+sigmaX2+2.*amp*offset*sigmaY-2.*amp*sigmaXY-2.*offset*sigmaX;
				  	double absdet = Math.abs(a_scale);
				  	if(absdet<1e-20){			
						amp = 0;
						offset= 0;
						err = sigmaX2;			
					}				
				  	int sindex=x+y*width;			  	
		  			amp_map[sindex]=amp;
		  			offset_map[sindex] = offset;			  	
				  	err_map[sindex]=err;
		  		}	  		
			}		
		}
	
		array = null;
	}
	*/
	private int make_Gaussian_PSF(double sigma, double radius, double[] psf, int[] dx, int[] dy, int halfw)
	{	
		int count = 0;
		double radius2 = radius*radius;
		double sigma2  = 2*sigma*sigma;
		for(int x=-halfw; x<=halfw; x++)
		{
			for(int y=-halfw; y<=halfw; y++)
			{
	  			double r2 = x*x + y*y;
	  			if(r2<=radius2){
					psf[count]= Math.exp(-r2/sigma2);
					dx[count] = x;
					dy[count] = y;
					count++;		
	  			}	
			}
	 	}	
	 	
	 	return count;
	}
	
	public static int make_Gaussian_PSF(double sigma, double[] psf, int[] dx, int[] dy, int halfw)
	{	
		int count = 0;	
		double sigma2  = 2*sigma*sigma;
		for(int x=-halfw; x<=halfw; x++)
		{
			for(int y=-halfw; y<=halfw; y++)
			{
	  			double r2 = x*x + y*y;  			
	  			psf[count]= Math.exp(-r2/sigma2);
				dx[count] = x;
				dy[count] = y;
				count++;						
			}
	 	}	
	 	
	 	return count;
	}
	
	private double[] make_Gaussian_PSF(double sigma, int halfw)
	{	
		int kernelw = 2*halfw+1;
		double[] psf =new double[kernelw*kernelw];	
		double sigma2  = 2*sigma*sigma;
		int count=0;
		for(int x=-halfw; x<=halfw; x++)
		{
			for(int y=-halfw; y<=halfw; y++)
			{
	  			double r2 = x*x + y*y;  			
	  			psf[count]= Math.exp(-r2/sigma2);
				count++;		  			
			}
	 	}	
	 	
	 	return psf;
	}
	
	
	public double precentile(double[] values, double precentage){
		if (values ==null || values.length==0) return Double.NaN;
		
	    int len = values.length;    
		double[] a = new double[len];
		a[0] = values[0];
		double min = a[0];
		double max = a[0];
		for(int i=1; i<len; i++){
			double v = values[i];
			a[i] = v;
			if(min>v) min=v;
			if(max<v) max=v;
		}	
		double x = precentage*len/100;
		if(x<=0) return min;
		else if(x>=len) return max;
			
		quicksort(a);
		
		int x1 = (int)Math.floor(x);		
		double y1 = a[x1];
		double y = y1;
		
		int x2 = x1+1;
		if(x2<len){
			double y2 =  a[x2];
			y = y1+(x-x1)*(y2-y1);
		}
	
		a = null;
		return y;
	}
	
	public void quicksort(double[] values) {
	    // Check for empty or null array
	    if (values ==null || values.length==0){
	      return;
	    }
	    qsort(values, 0, values.length-1);
	}
	
	private void qsort(double[] numbers, int low, int high) {
	    int i = low, j = high;
	    // Get the pivot element from the middle of the list
	    double pivot = numbers[low + (high-low)/2];
	
	    // Divide into two lists
	    while (i <= j) {
	      // If the current value from the left list is smaller then the pivot
	      // element then get the next element from the left list
	      while (numbers[i] < pivot) {
	        i++;
	      }
	      // If the current value from the right list is larger then the pivot
	      // element then get the next element from the right list
	      while (numbers[j] > pivot) {
	        j--;
	      }
	
	      // If we have found a values in the left list which is larger then
	      // the pivot element and if we have found a value in the right list
	      // which is smaller then the pivot element then we exchange the
	      // values.
	      // As we are done we can increase i and j
	      if (i <= j) {
	      	double temp=numbers[i]; numbers[i]=numbers[j]; numbers[j]=temp; //exchange(i, j)
	        i++;
	        j--;
	      }
	    }
	    // Recursion
	    if (low < j)  qsort(numbers, low, j);
	    if (i < high) qsort(numbers, i, high);
	}
}
