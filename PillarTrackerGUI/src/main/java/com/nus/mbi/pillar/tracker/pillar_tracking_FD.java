package com.nus.mbi.pillar.tracker;

//Xu Xiaochun @ Mechnobiology Institute, Singapore
//contact:mbixxc@nus.edu.sg
//update history:
//2017-02-08 (V1.1.3):
//  apply the pillar reconstrucation in frequency domain(FD)

import ij.*;
import ij.io.*;
import ij.process.*;
import ij.gui.*;
import java.io.*;
import java.io.File;
import ij.util.ThreadUtil; //using the mulit threading

import com.nus.mbi.pillar.detection.PillarDetector;
import com.nus.mbi.pillar.detection.FDResult;
import com.nus.mbi.pillar.detection.HighPassFFT;
import com.nus.mbi.pillar.drift.DriftCorrection;
import com.nus.mbi.pillar.drift.DriftDataLoader;
import com.nus.mbi.pillar.drift.FileHeaderDrift;
import fiji.util.gui.GenericDialogPlus;
import com.nus.mbi.pillar.grid.Grid_Creation_Plugin;
import ij.plugin.filter.RankFilters;
import java.awt.Color;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import com.nus.mbi.pillar.stat.IntPoint;

/**
 * @author Xiaochun Xu <mbixxc@nus.edu.sg>
 * @version 1.1.3
 * @since 2017-02-08
 */
public class pillar_tracking_FD extends pillar_tracking{
//public static int fileversion1 = 589450;
public static int fileversion2 = 589451;
public HighPassFFT fft = null;
public boolean show_reconstruction = true;
public boolean show_FFT = false;
private String input_settings_fname;
public boolean apply_mask_filter = false;

public pillar_tracking_FD(){
    super();
}

    @Override
    public void process() {
        if(experiment==null) return;
        if(fft==null || !fft.hasOffCenterPoints()) return;
        timestamp_start = time_formatter.format(System.currentTimeMillis());
        IJ.log("------------Start Tracking at "+ timestamp_start+"------------");
        
        experiment.setOverlay(null);
        fft.setup(experiment);
        
        int width = experiment.getWidth();
        int height = experiment.getHeight();
        int size = width*height;
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
        //if(image_enhancer_plugin!=null) dark_pillars = false;
        boolean use_enhancer = image_enhancer_plugin!=null;
        //boolean fft_valid = (fft!=null && fft.hasOffCenterPoints());
        
        PillarDetector pillar_detector = new PillarDetector();
        PillarDetector ref_pillar_detector = new PillarDetector();
        double sigma = Math.min(sigmax_PSF, sigmay_PSF);        
        
        if(use_enhancer) pillar_detector.setup(experiment, oblique, spacing, sigma, false);
        else pillar_detector.setup(experiment, oblique, spacing, sigma, dark_pillars);        
        pillar_detector.setGridAngle(grid_angle);
        
        ref_pillar_detector.setup(experiment, oblique, spacing, sigma, dark_pillars);
        ref_pillar_detector.setGridAngle(grid_angle);        
        
        List<IntPoint> list_roi = new ArrayList();
        
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
                else{
                    pillar_detector.process(ep, true);
                }
//                pillar_detector.process(ep, true);
                point_roi= experiment.getRoi();
            }
            
            if(point_roi==null) suc_load_roi = false;
            else if(!point_roi.isArea() && !point_roi.isLine()){
                FloatPolygon polygon_roi = point_roi.getFloatPolygon();
                float[] roix = polygon_roi.xpoints;
                float[] roiy = polygon_roi.ypoints;
                int margin = (int)Math.ceil(1.5*spacing);
                int margin_w = 2*margin+1;
                Roi roi_margin = new Roi(margin,margin,width-margin_w,height-margin_w);
                int num_points = polygon_roi.npoints;
                int num_valid = 0;
                for(int i=0; i<num_points; i++){
                    int x= (int)Math.floor(roix[i]);
                    int y= (int)Math.floor(roiy[i]);
                    if(roi_margin.contains(x, y)){
                        list_roi.add(new IntPoint(x,y));
                        num_valid++;
                    }
                }
                if(num_valid<1) suc_load_roi = false;
            }
            else suc_load_roi = false;
            
            if(!suc_load_roi) {
                IJ.log("There is no any Point selection, the centriod of image will be procesed!");
                return;
            }
            
        }
        int num_roi = list_roi.size();
        int numCPUs = ThreadUtil.getNbCpus();
        if(numThread<=0) numThread = numCPUs;
        else{
            if(numThread>num_roi) numThread = num_roi;
            int maxthreads = 20*numCPUs;
            if(numThread>maxthreads) numThread=maxthreads;
        }
        if(image_enhancer_plugin!=null) image_enhancer_plugin.setNumThreads(numThread);
        IJ.log("num of threads:" + numThread);
        
        if(!pillar_detector.isPrepared()) pillar_detector.prepare();
        if(!ref_pillar_detector.isPrepared()) ref_pillar_detector.prepare();
        
        int[] rcx = new int[num_roi];
        int[] rcy = new int[num_roi];
        for(int i=0; i<num_roi; i++) {
            IntPoint p = list_roi.get(i);
            rcx[i]=p.x;
            rcy[i]=p.y;
        }
        //check the first frame;
        ImageProcessor first_slice = experiment_s.getProcessor(1);        
//        SlicePixels slice_pixels = getSlicePixels(first_slice);
        FDResult ref_fft_data = fft.process_frame(first_slice, null, fft.off_center_radius);
        
        if(show_reconstruction){
            ImageProcessor img_ref = ref_fft_data.ifft_img;
            ImageStack stack_ref = new ImageStack(img_ref.getWidth(), img_ref.getHeight());
            stack_ref.addSlice("r="+fft.off_center_radius, img_ref);              
            ImagePlus ip = new ImagePlus("Reconstructed_"+experiment.getTitle(), stack_ref);
            ip.show();
        }
        
        if(show_FFT){
            ImageProcessor fft_ref = ref_fft_data.fft_img.getPowerSpectrum();
            SetMaskFilterRadius_Plugin.draw_mask(fft_ref, fft.points, fft.off_center_radius, Color.WHITE);
            ImageStack stack_ref = new ImageStack(fft_ref.getWidth(), fft_ref.getHeight());
            stack_ref.addSlice("r="+fft.off_center_radius, fft_ref);                        
            ImagePlus fft_ip = new ImagePlus("FFT_"+experiment.getTitle(), stack_ref);            
            fft_ip.show();                 
        }        
        
        double[] ref_pixels = getReferencePixels(ref_fft_data.ifft_img);
        int[][] ref_centroidsXY = ref_pillar_detector.detect_fast(ref_pixels);        
        int[] ref_cx = ref_centroidsXY[0];
        int[] ref_cy = ref_centroidsXY[1];        
        //int[][] txy = tracking_relaxiation(ref_cx, ref_cy, width, height, slice_pixels.ref_fft_data.fft_img, ref_pillar_detector);
        //List<IntPoint> pair = cross_check(txy[0],txy[1], rcx, rcy, width, height, catch_radius);
        List<IntPoint> pair = cross_check(ref_cx, ref_cy, rcx, rcy, width, height, catch_radius);
        
        nrois = pair.size();
        if(nrois<1) {
            IJ.showMessage("can not find a point");
            return;
        }

        // seed first frame
        int[] rcx_rois = new int[nrois];
        int[] rcy_rois = new int[nrois];
        int[] gx = new int[nrois];
        int[] gy = new int[nrois];
        int[] tx = new int[nrois];
        int[] ty = new int[nrois];
        for(int i=0; i<nrois; i++) {
            IntPoint p = pair.get(i);
            gx[i]=ref_cx[p.x];
            gy[i]=ref_cy[p.x];
            tx[i]=rcx_rois[i]=rcx[p.y];
            ty[i]=rcy_rois[i]=rcy[p.y];
        }
        
        IJ.log("npillars="+nrois);
        
        double[] pixels = getSlicePixels(first_slice);         
        if(!apply_mean_rank_filter){                                    
            int[][] centroidsXY = pillar_detector.detect_fastest(pixels);
            int num_detected = centroidsXY[0].length;
            if(num_detected>nrois*4){
                IJ.log("The local maximas are overwhelming, will denoise the image");
                IJ.log("detected in first frame:" + num_detected);
                //using_rank_filter = true;
                apply_mean_rank_filter = true;
                RankFilters mean_filter = new RankFilters();
                ImageProcessor slice_filter = first_slice.duplicate();
                mean_filter.rank(slice_filter, 1, RankFilters.MEAN);
                double[] filter_pixels = getSlicePixels(slice_filter); 
                centroidsXY = pillar_detector.detect_fastest(filter_pixels);
                num_detected = centroidsXY[0].length;
                IJ.log("check again after denoising:" + num_detected);
            }
        }        
        
        //optimization
        double[][] cx = new double[nrois][sliceN];
        double[][] cy = new double[nrois][sliceN];
        double[][] grid_x = new double[nrois][sliceN];
        double[][] grid_y = new double[nrois][sliceN];
        
        double[][] disX = new double[nrois][sliceN];
        double[][] disY = new double[nrois][sliceN];
        LocalizationWindowQueue queue_ref = new LocalizationWindowQueue();
        LocalizationWindowQueue queue = new LocalizationWindowQueue();
        
        int num_threads_FFT = ThreadUtil.getNbCpus();
        if(sliceN<num_threads_FFT) num_threads_FFT = sliceN;
        if(numThread<num_threads_FFT) num_threads_FFT = numThread;
        
        FDResult[] FFT_DATA = new FDResult[num_threads_FFT];
        
        IJ.log("frame#  detected_ref    matched_ref  detected    matched");
        for(int s=0;s<sliceN;s++){//pillar detection for each frame            
            if(s>0){
                if((s-1)%num_threads_FFT==0){
                    //FFT
                    ExecutorService exec=Executors.newFixedThreadPool(num_threads_FFT);  
                    List<Callable<Integer>> callList=new ArrayList<>();  
                    List<Integer> list = new ArrayList<>();
                    for(int f=s; f<(s+num_threads_FFT); f++) if(f<sliceN) list.add(f);
                    
                    for(Integer f : list) {
                        ImageProcessor slice = experiment_s.getProcessor(f + 1);                        
                        int k = f-s;
                        callList.add((Callable<Integer>) () -> {                                                    
                            FFT_DATA[k] = fft.process_frame(slice, null, fft.off_center_radius);
                            return 1;
                        });
                    }

                    try{
                        exec.invokeAll(callList);  
                    }
                    catch(Exception ex){ }
                }
                //detection
                ImageProcessor slice = experiment_s.getProcessor(s + 1);
                //slice_pixels = getSlicePixels(slice);
                pixels = getSlicePixels(slice);
                ref_fft_data = FFT_DATA[(s-1)%num_threads_FFT];
                ref_pixels = getReferencePixels(ref_fft_data.ifft_img);
                
                //ip.setProcessor(slice_pixels.fft_data.ifft_img);
                ref_centroidsXY = ref_pillar_detector.detect_fastest(ref_pixels);
                ref_cx = ref_centroidsXY[0];
                ref_cy = ref_centroidsXY[1];
                pair = cross_check(gx, gy, ref_cx, ref_cy, width, height, catch_radius);
                int grid_matches = pair.size();
                int[][] new_gxy = tracking(gx, gy, ref_cx, ref_cy, pair);
                gx=new_gxy[0];
                gy=new_gxy[1];
                
                int[][] txy = tracking_relaxiation(gx, gy, width, height, ref_fft_data.fft_img, ref_pillar_detector);
                
                int[][] centroidsXY;
                if(apply_mean_rank_filter){
                    RankFilters mean_filter = new RankFilters();
                    ImageProcessor slice_filter = slice.duplicate();
                    mean_filter.rank(slice_filter, 1, RankFilters.MEAN);
                    double[] filter_pixels = getSlicePixels(slice_filter); 
                    centroidsXY = pillar_detector.detect_fastest(filter_pixels);
                }else centroidsXY = pillar_detector.detect_fastest(pixels);
                
                int[] org_cx = centroidsXY[0];
                int[] org_cy = centroidsXY[1];
                pair = cross_check(txy[0], txy[1], org_cx, org_cy, width, height, catch_radius);
                for(int k=0; k<pair.size(); k++) {
                    IntPoint p = pair.get(k);
                    int i  = p.x;
                    int i2 = p.y;
                    tx[i]=org_cx[i2];
                    ty[i]=org_cy[i2];			
                }
                int nrefs = ref_cx.length;
                int ncent = org_cx.length;
                int matches = pair.size();
                IJ.log((s+1)+"  "+nrefs+"  "+grid_matches+"  "+ncent+"   "+matches);                
            }
            else{
                int nrefs = ref_cx.length;
                int ncent = num_roi;
                int matches = pair.size();
                IJ.log((s+1)+"  "+nrefs+"  "+ncent+"   "+matches);    
            }
            
            //optimization
            boolean dark_object = dark_pillars;
            if(use_enhancer) dark_object = false;
            if(localization_algorithm == localization_algorithm_CG)
            {
                double[][] estimated_xy = optimize(tx, ty, pixels, width, height, sigma, dark_object);
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
                //optimization reference, pillar reconstruction
                estimated_xy = optimize(gx, gy, ref_pixels, width, height, sigma, dark_pillars);
                for(int i=0; i<nrois; i++){
                    grid_x[i][s] = Double.NaN;
                    grid_y[i][s] = Double.NaN;
                    if(estimated_xy!=null){
                        double[] ml_xy = estimated_xy[i];
                        if(ml_xy!=null){
                            grid_x[i][s] = ml_xy[0];
                            grid_y[i][s] = ml_xy[1];
                        }
                    }
                }
            }
            else{                
                optimize_queue(cx, cy, queue, s, tx, ty, pixels, width, height, sigma, dark_object);
                optimize_queue(grid_x, grid_y, queue_ref, s, gx, gy, ref_pixels, width, height, sigma, dark_pillars);                				
            }

            IJ.showStatus("tracking in slice="+s);
            IJ.showProgress(s+1,sliceN);
        }
        
        if(queue.getCount()>0){
            boolean dark_object = dark_pillars;
            if(use_enhancer) dark_object = false;
            optimize_queue(cx, cy, queue, sigma, dark_object);              
        }
        
        if(queue_ref.getCount()>0){
            optimize_queue(grid_x, grid_y, queue_ref, sigma, dark_pillars); 
        }
        
        for(int s=0; s<sliceN; s++){
            for(int i=0; i<nrois; i++){
                disX[i][s] = cx[i][s] - grid_x[i][s];
                disY[i][s] = cy[i][s] - grid_y[i][s];
            }
        }           
        
        //saving the localization data.
        FileOutputStream fos = null;
        FileChannel fch = null;
        //String text_header = " ";
        //for(int i=0; i<nrois; i++) text_header += ("X" + i + "	Y" + i + "	");
        try{
            File output_file = new File(output_fname);
            fos = new FileOutputStream(output_file);
            fch = fos.getChannel();
            IJ.showStatus("writing raw data into binary file");
            writeFileHeader_ver2(fch,sliceN);
            writeCentroids(fch,rcx_rois,rcy_rois);
            writePoints(fch,cx,cy);
            writePoints(fch,disX,disY);
            IJ.log("write raw data into binary file:" + output_fname);
            
            out_text_File = new FileWriter(output_fname+".txt");// Write to text file
            out_text = new PrintWriter(out_text_File);
            write_header_text(info, sliceN);
            
            IJ.log("write logs into text file:" + output_fname+".txt");            
            //OpenDialog.setDefaultDirectory(output_file.getAbsoluteFile().getParent());
            IJ.showStatus("done with saving raw data");
        }
        catch(Exception ex){
            IJ.log("write file error!");
        }
        
        double[][] corrected_CX = null;
        double[][] corrected_CY = null;
        
        if(sliceN>1){
            double[][] mxy = Grid_Creation_Plugin.average_tracks(grid_x, grid_y);            
            double[][] drift = Grid_Creation_Plugin.compute_drift(mxy[0], mxy[1]);
            
            boolean[] reference_flags = new boolean[nrois];            
            for(int i=0; i<nrois; i++) reference_flags[i] = false;    
            
            double[] driftXY = new double[sliceN*2];
            corrected_CX = new double[nrois][sliceN];
            corrected_CY = new double[nrois][sliceN];
            for(int i=0; i<sliceN; i++){
                double dx = driftXY[2*i]   = drift[0][i];
                double dy = driftXY[2*i+1] = drift[1][i];
                for(int k=0; k<nrois; k++) {
                    corrected_CX[k][i] = cx[k][i] + dx;
                    corrected_CY[k][i] = cy[k][i] + dy;
                }
            }
            
            DriftCorrection.search_slient_pillars(corrected_CX, corrected_CY, reference_flags);
            
            try{
                //IJ.log("write corrected data into binary file:" + output_fname);
                IJ.showStatus("writing corrected data into binary file");                
                writePoints(fch,corrected_CX,corrected_CY);
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
        }
        
        //close the file handle
        try{	//dos.close();
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
        
        saveCSV(cx, cy, disX, disY, corrected_CX, corrected_CY);
        
        IJ.showStatus("tracking done");
    }
    
    public void saveCSV(double[][] cx, double[][] cy, double[][] disX, double[][] disY, double[][] corrected_CX, double[][] corrected_CY){
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
            if(corrected) out_text.println("frame,trajectory,x,y,x_raw,y_raw,dx,dy");
            else out_text.println("frame,trajectory,x,y,dx,dy");
            for(int f=0; f<nf; f++) {        
                for(int c=0; c<np; c++){
                    out_text.print(""+ (f+1) + ","+ (c+1)+",");
                    if(corrected) out_text.print((corrected_CX[c][f] + ","+ corrected_CY[c][f]+","));
                    out_text.print((cx[c][f] + ","+ cy[c][f]+","));
                    out_text.print((disX[c][f] + ","+ disY[c][f]));
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
    
    private void optimize_queue(double[][] cx, double[][] cy, LocalizationWindowQueue queue, double sigma,boolean dark_object){
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
    
    private void optimize_queue(double[][] cx, double[][] cy, LocalizationWindowQueue queue, int frame,int[] tx, int[] ty, double[] pixels, int width, int height, double sigma, boolean dark_object){
        for(int k=0; k<nrois; k++) {                    
            int x = tx[k];
            int y = ty[k];                                  
            LocalizationWindow win = LocalizationWindow.create(x, y, pixels, width, height, kernel_w);                                    
            if(queue.isFull()){
                optimize_queue(cx, cy, queue, sigma, dark_object);
                queue.reset();
            } 
            queue.add(win, k, frame);
        }
    }
    
    @Override
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
            writeFileHeader_ver2(fch,sliceN);
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
    
public static double[][] getXYOnFrame(double[][] X, double[][] Y, int iframe){                            
    int npillars = X.length;        
    double[][] xy = new double[2][npillars];
    for(int i=0; i<npillars; i++){
        xy[0][i] = X[i][iframe];
        xy[1][i] = Y[i][iframe];
    }
    return xy;
}
    
private double[][] optimize(int[] x, int[] y, double[] pixels, int width, int height, double sigma, boolean dark_object){
    int npillars = x.length;
    double[] seed_cx = new double[npillars];
    double[] seed_cy = new double[npillars];
    for(int i=0; i<npillars; i++) {
        seed_cx[i] = x[i];
        seed_cy[i] = y[i];
    }
    double[][] estimated_xy=null;
    if(localization_algorithm == localization_algorithm_CG){
        if(dark_object){    //invert the image, CG need the dark background       
            int size = width*height;
            double[] invert_pixels = new double[size];
            for(int i=0; i<size; i++) invert_pixels[i] = -pixels[i];
            estimated_xy = LocalizationFunction.metric_centroid_safe(numThread, invert_pixels, width, height, seed_cx, seed_cy, kernel_w/2, 1, box_constrian_R);                        
        }
        else
         estimated_xy= LocalizationFunction.metric_centroid_safe(numThread, pixels, width, height, seed_cx, seed_cy, kernel_w/2, 1, box_constrian_R);                        
    }
    else estimated_xy= LocalizationFunction.ML_localization(     numThread, pixels, width, height, seed_cx, seed_cy, kernel_w, sigma, sigma, box_constrian_R, dark_object);
    return estimated_xy;
}



private static int[][] tracking(int[] cx1, int[] cy1, int[] cx2, int[] cy2, List<IntPoint> pair){        
    int npillars = cx1.length;
    boolean[] flag = new boolean[npillars];    
    int[][] new_p = new int[2][npillars];
    for(int i=0; i<npillars; i++){
        flag[i]=false;
        new_p[0][i] = cx1[i];
        new_p[1][i] = cy1[i];
    }
    int num = pair.size();
    if(num<1) return new_p;
    
    if(num==npillars){
        for(int k=0; k<num; k++) {                                
            IntPoint p = pair.get(k);
            int i  = p.x;
            int i2 = p.y;
            int newx = cx2[i2];
            int newy = cy2[i2];
            new_p[0][i] = newx;
            new_p[1][i] = newy;                                                    
        }
        return new_p;
    }
    
    double dx = 0;
    double dy = 0;    
    for(int k=0; k<num; k++) {                                
        IntPoint p = pair.get(k);
        int i  = p.x;
        int i2 = p.y;
        int newx = cx2[i2];
        int newy = cy2[i2];
        dx += (newx-cx1[i]);
        dy += (newy-cy1[i]);
        new_p[0][i] = newx;
        new_p[1][i] = newy;                                        
        flag[i] = true;
    }           
    dx /= num;
    dy /= num;    
    for(int i=0; i<npillars; i++){
        if(flag[i]==false){
            new_p[0][i] = (int)Math.round(cx1[i] + dx);
            new_p[1][i] = (int)Math.round(cy1[i] + dy);
        }
    }    
    return new_p;
}

private double[] getReferencePixels(ImageProcessor ifft_img){
    float[] ref_pixels_float = (float[])ifft_img.getPixels();
    int size = ref_pixels_float.length;    
    double[] ref_pixels = new double[size];
    for(int p=0; p<size;p++) ref_pixels[p] = ref_pixels_float[p];
    return ref_pixels;
}

private int[][] tracking_relaxiation(int[] gx, int[] gy, int width, int height, FHT fht, PillarDetector detector){
    int npillars = gx.length;
    int[] tx = new int[npillars];
    int[] ty = new int[npillars];
    System.arraycopy(gx, 0, tx, 0, npillars);
    System.arraycopy(gy, 0, ty, 0, npillars);
    
    int start_r=fft.off_center_radius+fft.step_radius_off_center;    
    int step_r = fft.step_radius_off_center;
    List<Integer> list = new ArrayList<>();
    for(int r=start_r; r<=fft.end_radius_off_center; r += step_r) list.add(r);
    int list_size = list.size();
    if(list_size>1){
        int[][] cx = new int[list_size][];
        int[][] cy = new int[list_size][];
        
        int num_threads = ThreadUtil.getNbCpus();
        if(list_size<num_threads) num_threads = list_size;
        if(numThread<num_threads) num_threads = numThread;
        
        ExecutorService exec=Executors.newFixedThreadPool(num_threads);  
        List<Callable<Integer>> callList=new ArrayList<>();          
        int len = list_size/num_threads;
        if(len==0){
            num_threads = list_size;
            len = 1;
        }

        for(int i=0; i<num_threads; i++){
            final List<Integer> sub_list;
            if(i==num_threads-1) sub_list = list.subList(i*len, list_size);
            else sub_list=list.subList(i*len, len*(i+1)>list.size()?list.size():len*(i+1));

            // 
            callList.add((Callable<Integer>) () -> {
                for(Integer r : sub_list) {
                    double[] pixels = getMaskedSlicePixels(fht, r);
                    int[][] centroidsXY = detector.detect_fastest(pixels);
                    int k = (r-start_r)/step_r;
                    cx[k] = centroidsXY[0];
                    cy[k] = centroidsXY[1];
                }
                return 1;
            });
        }

        try{
            exec.invokeAll(callList);  
        }
        catch(Exception ex){ }
        
        for(int j=0; j<list_size; j++){
            int[] org_cx = cx[j];
            int[] org_cy = cy[j];
            List<IntPoint> pair = cross_check(tx, ty, org_cx, org_cy, width, height, catch_radius);
            for(int k=0; k<pair.size(); k++) {
                IntPoint p = pair.get(k);
                int i  = p.x;
                int i2 = p.y;
                tx[i]=org_cx[i2];
                ty[i]=org_cy[i2];
            }
        }
    }
    else{    
        for(int r=start_r; r<=fft.end_radius_off_center; r += step_r){
            double[] pixels = getMaskedSlicePixels(fht, r);
            int[][] centroidsXY = detector.detect_fastest(pixels);
            int[] org_cx = centroidsXY[0];
            int[] org_cy = centroidsXY[1];
            List<IntPoint> pair = cross_check(tx, ty, org_cx, org_cy, width, height, catch_radius);
            for(int k=0; k<pair.size(); k++) {
                IntPoint p = pair.get(k);
                int i  = p.x;
                int i2 = p.y;
                tx[i]=org_cx[i2];
                ty[i]=org_cy[i2];
            }
        }
    }
    
    int[][] txy = new int[2][];
    txy[0] = tx;
    txy[1] = ty;
    return txy;
}

private double[] getMaskedSlicePixels(FHT fht, int off_center_radius){    
    FDResult fft_data = fft.process_frame(fht, off_center_radius);
    ImageProcessor ifft = fft_data.ifft_img;
    float[] pixels_float = (float[]) ifft.getPixels();
    int width = ifft.getWidth();
    int height = ifft.getHeight();
    int size = width*height;
    double[] pixels = new double[size];    
    for(int p=0; p<size;p++) pixels[p] = pixels_float[p];  
    return pixels;
}

private void write_header_text(String info, int sliceN){
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
    out_text.println("use mean filter to find maxima:" + apply_mean_rank_filter);
    out_text.println("number of pillars:" + nrois);
    out_text.println("number of frames:" + sliceN);
    out_text.println("save to binary file:" + output_fname);			
    //out_text.println("----raw coordinates of tracks----");
    //out_text.println(text_header);
}

private void writeFileHeader_ver2(FileChannel out, int nframes)  throws IOException {        
        int header_size = Integer.SIZE*5 + Double.SIZE*7 + Byte.SIZE*5;
        ByteBuffer header_buffer = ByteBuffer.allocate(header_size);
        header_buffer.putInt(fileversion2);
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
        header_buffer.put(dark_pillars?(byte)1:0);
        //header_buffer.put(use_minimum_std?(byte)1:0);
        header_buffer.put(apply_mean_rank_filter?(byte)1:0);
        header_buffer.put(use_metric_CG?(byte)1:0);
        header_buffer.put(image_enhancer_plugin!=null?(byte)1:0);
        //header_buffer.put(false?(byte)1:0); //dont use enhancer
        boolean has_fft = (fft!=null && fft.hasOffCenterPoints());
        header_buffer.put(has_fft?(byte)1:0);        
        header_buffer.flip();
        out.write(header_buffer);                    
        if(image_enhancer_plugin!=null) writePSFImage(out, image_enhancer_plugin);        
        if(has_fft) writeFFT(out, fft);
}

private void writeFFT(FileChannel out, HighPassFFT myfft)  throws IOException {
    List<IntPoint> list = myfft.getOffCenterPoints();    
    int npoints = list.size();
    int header_size = Integer.SIZE*(5+npoints*2);
    
    ByteBuffer header_buffer = ByteBuffer.allocate(header_size);    
    header_buffer.putInt(myfft.getOffCenterRadius());
    header_buffer.putInt(myfft.center_radius);
    header_buffer.putInt(myfft.step_radius_off_center);
    header_buffer.putInt(myfft.end_radius_off_center);
    header_buffer.putInt(list.size());
    for(int i=0; i<npoints; i++){
        IntPoint p = list.get(i);
        header_buffer.putInt(p.x);
        header_buffer.putInt(p.y);
    }    
    header_buffer.flip();
    out.write(header_buffer);    
}

private int[][] load_settings(){
        DriftDataLoader drift_loader = new DriftDataLoader();
        int[][] start_points = null;
        if(drift_loader.load_settings(input_settings_fname)){
            FileHeaderDrift header = drift_loader.file_header;                                
            lattice = header.lattice;
            spacing = header.spacing;
            oblique = header.oblique;
            grid_angle = header.grid_angle;
            sigmax_PSF = sigmay_PSF = header.sigma_PSF;
            diameter = header.diameter;                                
            dark_pillars = header.dark_pillars;
            localization_algorithm = header.use_metric_CG ? pillar_tracking.localization_algorithm_CG : pillar_tracking.localization_algorithm_Levmar;                                
            kernel_w = header.kernel_w;
            box_constrian_R = header.box_constrian_R;
            catch_radius = header.catch_radius;                                
            if(header.use_enhancer)
                image_enhancer_plugin = drift_loader.enhancer;
            if(header.use_fft){
                if(fft==null) fft = new HighPassFFT();
                fft.setOffCenterRadius(drift_loader.mask_radius);
                fft.step_radius_off_center = drift_loader.start_radius;
                fft.end_radius_off_center = drift_loader.end_radius;                
                int[][] points = drift_loader.fft_points;
                fft.setOffCenterPoints(points[0], points[1]);
            }  
            start_points = drift_loader.startXY;
            if(start_points!=null) experiment.setRoi(new PointRoi(start_points[0], start_points[1], start_points[0].length));
        }
        return start_points;
}

@Override
public void run(String arg){
                IJ.freeMemory();
		if(show_dialog){
                    boolean checkonce = arg.contains("checkonce");
                    boolean consistent;	
                    do {
                            consistent=true;
                            if (!showDialog()) return;                            
                            load_settings();
                            if(fft==null || !fft.hasOffCenterPoints()){IJ.showMessage("the mask filter for pillar reconstrucation are not avaliable"); consistent=false;} 
                            
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

@Override
boolean showDialog() {
	int[] wList = WindowManager.getIDList();	
	if (wList==null || wList.length<1) {
		IJ.showMessage("Pillar Tracking With Grid Reconstruction", "1 channel time series are required");
		return false;
	}

	int wlistlen = wList.length;
	
	ArrayList<String> titles_list = new ArrayList<String>();
	for (int i=0; i<wlistlen; i++) {
		ImagePlus imp = WindowManager.getImage(wList[i]);		
		if(imp!=null && imp.getNDimensions()<=3) titles_list.add(imp.getTitle());
	}

	if(titles_list.isEmpty()){
		IJ.showMessage("Pillar Tracking With Grid Reconstruction", "1 channel time series are required");	
		return false;
	}
	
	String titles[] = titles_list.toArray(new String[0]);
	
	GenericDialogPlus gd = new GenericDialogPlus("Pillar Tracking With Grid Reconstruction");
        gd.addImageChoice("ONLY one channel time series:", titles[0]);
        output_fname = "pillar_tracks.fd";
        if(input_settings_fname!=null) output_fname = input_settings_fname + ".fd";
	int str_len = output_fname.length();
	if(str_len<50) str_len = 50;
	gd.addFileField("input settings file:", input_settings_fname, str_len);
        gd.addFileField("output file name:", output_fname, str_len);
	
	gd.showDialog();
	if (gd.wasCanceled()) return false;
	experiment = gd.getNextImage();
	input_settings_fname = gd.getNextString();
	output_fname = gd.getNextString();
	
	return true;
}
}
