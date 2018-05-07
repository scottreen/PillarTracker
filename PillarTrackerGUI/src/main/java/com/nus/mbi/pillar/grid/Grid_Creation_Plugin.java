package com.nus.mbi.pillar.grid;

import com.nus.mbi.pillar.drift.DriftDataLoader;
import ij.plugin.PlugIn;
import fiji.util.gui.*;
//import static com.nus.mbi.pillar.grid.GridDataLoader.*;
import ij.IJ;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Files;
import java.nio.file.Path;
import com.nus.mbi.pillar.stat.BasicStatisitic;

/**
 *
 * @author xiaochun
 */
public class Grid_Creation_Plugin implements PlugIn{
    private String centroid_name = "H:\\pillar_tracking_slover\\grid_analyse\\korean_data\\pillar-tracks.bin";
    private double lattice = 208;
    private double spacing = 8.37;
    private double oblique = 3.75;
    private double grid_angle = 90;
    private double diameter = 4;
        
    private String output_name = "H:\\pillar_tracking_slover\\grid_analyse\\korean_data\\pillar-tracks.bin.dxy";
    private boolean use_restored_pillars = false;
    private String restored_fname = "H:\\pillar\\korea\\test_fft\\FFT_Filter-MDA_pillar_10sec_15min__R3D_1.bin";
    
    double[][] raw_tracksX;
    double[][] raw_tracksY;
    
    double[][] corrected_tracksX;
    double[][] corrected_tracksY;    
    
    private double[][] DX;
    private double[][] DY;
    
    private int[][] IX;
    private int[][] IY;
    
    private boolean need_show_dialog = true;
    private boolean dialog_canceld = false;
    private boolean use_1st_frame_pillars = false;
    
    public void setShowDialog(boolean show){
        need_show_dialog = show;
    }
    
    public boolean IsDialogCanceled(){
        return dialog_canceld;
    }
    
    private boolean check_writable_file(String fname){
        boolean can_write = false;
        try{
            Path file = new File(fname).toPath();
            if(Files.exists(file)){
                boolean iw = Files.isWritable(file);
                if(iw) can_write = true;
            }
            else if(Files.notExists(file)){
                Files.createFile(file);
                can_write = true;
            }            
        }
        catch(Exception e){
            IJ.log("file is not writable->" + fname);
        }
        
        return can_write;
    }
    
    private boolean check_exist_file(String fname){
        boolean exist = false;
        try{
            Path file = new File(fname).toPath();
            if(Files.exists(file)){
                exist = true;
            }            
        }
        catch(Exception e){
            IJ.log("file is not accessable->" + fname);
        }
        
        return exist;
    }
    
    public void run(String arg) {
        if(need_show_dialog){
            dialog_canceld = false; 
            boolean consistent;			
            do {
                    consistent=true;
                    if (!showDialogFD()){
                        dialog_canceld = true; 
                        return;
                    } 
                    else if(!check_writable_file(output_name)) {IJ.showMessage("output file is not writable"); consistent=false;}
                    else if(!check_exist_file(centroid_name)) {IJ.showMessage("input file doesn't exist"); consistent=false;}	
                    else if(spacing <0) {IJ.showMessage("pillar spacing must be must be greater than zero"); consistent=false;}	
                    else if(lattice<0) {IJ.showMessage("lattice be greater than zero"); consistent=false;}	
                    else if(diameter<0) {IJ.showMessage("lattice be greater than zero"); consistent=false;}	
                    else if(oblique<-45 || oblique>45) {IJ.showMessage("Grid oblique must be among the range of -45 to 45 degree "); consistent=false;}	
                    else if(grid_angle<=0 || grid_angle>90) {IJ.showMessage("Grid angle must be among the range of 0 to 90 degree "); consistent=false;}
            } while (!consistent);        
        }
        if(use_restored_pillars){
            if(!check_exist_file(restored_fname)){
                use_restored_pillars = false;
                IJ.log("the reference file is not existed");
            }
        }
        print_current_settings();        
        process();
    }
    
    private void print_current_settings(){
        IJ.log("input fname:" + this.centroid_name);
        if(use_restored_pillars) IJ.log("input reference fname:" + this.restored_fname);
	IJ.log("lattice:" + this.lattice);	
	IJ.log("pillar spacing:" + this.spacing);	
        IJ.log("grid oblique:" + this.oblique);
        IJ.log("grid angle:" + this.grid_angle);
	IJ.log("pillar diameter:" + this.diameter);	
    }
    
    public boolean setup(String input, double lattice, double spacing, double oblique, double grid_angle, double diameter){
        boolean suc = true;
        if(check_exist_file(input)) this.centroid_name = input;
        else suc = false;
        
        if(lattice>0) this.lattice = lattice;
        else suc = false;
        
        if(spacing>0) this.spacing = spacing;
        else suc = false;
        
        if(oblique>-45 && oblique<45) this.oblique = oblique;
        else suc = false;        
        
        if(grid_angle>0 && grid_angle<=90) this.grid_angle = grid_angle;
        else suc = false;
        
        if(diameter>0) this.diameter = diameter;
        else suc = false;
        
        IJ.log("set up successful for grid creation? " + suc);	
        print_current_settings();
        //IJ.log("input fnmae:" + this.centroid_name);
	//IJ.log("lattice:" + this.lattice);	
	//IJ.log("pillar spacing:" + this.spacing);	
        //IJ.log("grid oblique:" + this.oblique);
	//IJ.log("pillar diameter:" + this.diameter);			
        
        this.output_name = input + ".dxy";
        
	return suc;
    }
    
    public boolean setupReferenceFilename(String ouput){
        boolean suc = true;        
        if(check_exist_file(ouput)){
            use_restored_pillars = true;
            this.restored_fname = ouput;
        }
        else suc = false;
        
        return suc;
    }
    
    public boolean setupOutputFilename(String ouput){
        boolean suc = true;
        if(check_writable_file(ouput)) this.output_name = ouput;
        else suc = false;
        
        return suc;
    }
    
    public String getOutputFilename(){
        return this.output_name;
    }
    
    public String getInputFilename(){
        return this.centroid_name;
    }
    
    public void process(){
        DriftDataLoader drift_data = new DriftDataLoader();
        drift_data.load_tracks(centroid_name, lattice);   
        DriftDataLoader reference_data = new DriftDataLoader();
        if(use_restored_pillars) reference_data.load_tracks(restored_fname, lattice);           
        if(drift_data.load_data_suc){
            IJ.showStatus("Grid creation");            
            int npillars = drift_data.npillars;
            int nframes = drift_data.nframes;
            DX = new double[npillars][nframes];
            DY = new double[npillars][nframes];
            if(reference_data.load_data_suc){
                raw_tracksX = drift_data.raw_tracksX;
                raw_tracksY = drift_data.raw_tracksY;  
                corrected_tracksX = drift_data.correct_tracksX;
                corrected_tracksY = drift_data.correct_tracksY; 
                double[][] GX = reference_data.raw_tracksX;
                double[][] GY = reference_data.raw_tracksY;
                //update the correct data;         
                double[][] mxy = average_tracks(GX, GY);
                double[][] drift = compute_drift(mxy[0], mxy[1]);
                if(drift!=null){
                    if(corrected_tracksX==null){
                        corrected_tracksX = new double[npillars][nframes];
                        corrected_tracksY = new double[npillars][nframes];
                    }
                    for(int f=0; f<nframes; f++){                    
                        for(int c=0; c<npillars; c++){
                            corrected_tracksX[c][f] = raw_tracksX[c][f] + drift[0][f];
                            corrected_tracksY[c][f] = raw_tracksY[c][f] + drift[1][f];
                        }
                    }
                }
                
                double catch_radius = spacing/2.0;
                double max_x = reference_data.X_MAX/lattice;//max(xpoints);
                double max_y = reference_data.Y_MAX/lattice;//max(ypoints);
                int map_width = (int)Math.ceil(max_x);
                int map_height = (int)Math.ceil(max_y);
                for(int f=0; f<nframes; f++){
                    double[][] xyc = drift_data.getXYOnFrame(f);                    
//                    double[][] gxy = correct_drift(GX, GY, mxy[0],mxy[1],f);
//                    for(int c=0; c<npillars; c++){
//                            gxy[0][c] /= lattice;
//                            gxy[1][c] /= lattice;
//                    }
                    double[][] gxy = reference_data.getXYOnFrame(f); 
                    if(xyc!=null && gxy!=null){
                        double[][] dxy = GridCreation.computeDeflections(xyc[0], xyc[1], gxy[0], gxy[1], map_width, map_height, catch_radius);
                        for(int c=0; c<npillars; c++){
                            DX[c][f] = dxy[0][c]*lattice;
                            DY[c][f] = dxy[1][c]*lattice;
                        }
                    }
                    else{
                        for(int c=0; c<npillars; c++){
                            DX[c][f] = Double.NaN;
                            DY[c][f] = Double.NaN;
                        }
                    }
                    IJ.showProgress(f+1, nframes);   
                }               
                
                IJ.showStatus("saving the deflection data");
                writeDXY_version3(output_name);
            }
            else if(use_1st_frame_pillars){
                raw_tracksX = drift_data.raw_tracksX;
                raw_tracksY = drift_data.raw_tracksY;  
                corrected_tracksX = drift_data.correct_tracksX;
                corrected_tracksY = drift_data.correct_tracksY; 
                
                double catch_radius = spacing/2.0;
                double max_x = drift_data.X_MAX/lattice;//max(xpoints);
                double max_y = drift_data.Y_MAX/lattice;//max(ypoints);
                int map_width = (int)Math.ceil(max_x);
                int map_height = (int)Math.ceil(max_y);
                double[][] gxy = drift_data.getXYOnFrame(0); 
                for(int f=0; f<nframes; f++){
                    double[][] xyc = drift_data.getXYOnFrame(f);                    
                    if(xyc!=null && gxy!=null){
                        double[][] dxy = GridCreation.computeDeflections(xyc[0], xyc[1], gxy[0], gxy[1], map_width, map_height, catch_radius);
                        for(int c=0; c<npillars; c++){
                            DX[c][f] = dxy[0][c]*lattice;
                            DY[c][f] = dxy[1][c]*lattice;
                        }
                    }
                    else{
                        for(int c=0; c<npillars; c++){
                            DX[c][f] = Double.NaN;
                            DY[c][f] = Double.NaN;
                        }
                    }
                    IJ.showProgress(f+1, nframes);   
                }               
                
                IJ.showStatus("saving the deflection data");
                writeDXY_version3(output_name);
            }
            else{
                IX = new int[npillars][nframes];
                IY = new int[npillars][nframes];            
                for(int f=0; f<nframes; f++){
                    double[][] xy = drift_data.getXYOnFrame(f);
                    GridCreation grid_creator = new GridCreation(xy[0], xy[1]);
                    grid_creator.label_grid(spacing, oblique, grid_angle);
                    grid_creator.create_grid_points(spacing/2.0);
                    double[] dx = grid_creator.getDX();
                    double[] dy = grid_creator.getDY();
                    int[] ix = grid_creator.getIX();
                    int[] iy = grid_creator.getIY();

                    for(int c=0; c<npillars; c++){
                        DX[c][f] = dx[c]*lattice;
                        DY[c][f] = dy[c]*lattice;
                        IX[c][f] = ix[c];
                        IY[c][f] = iy[c];
                    }
                    IJ.showProgress(f+1, nframes);                
                }

                raw_tracksX = drift_data.raw_tracksX;
                raw_tracksY = drift_data.raw_tracksY;                   
                corrected_tracksX = drift_data.correct_tracksX;
                corrected_tracksY = drift_data.correct_tracksY; 

                IJ.showStatus("saving the deflection data");
                writeDXY_version2(output_name);
            }
            //writeDXY(output_name, GridDataLoader.fileformat_1, lattice, diameter, spacing, tracksX, tracksY, DX, DY);
            IJ.showStatus("Grid Creation Done");
            IJ.log("deflection file saved to->" + output_name);
        }
    }
    
    private double[][] correct_drift(double[][] cx, double[][] cy, double[] mx, double[] my, int frame){
        int nrois = cx.length;
        int sliceN = cx[0].length;                
        double[][] drift = compute_drift(mx, my, frame);
        double[] ccx = new double[sliceN];
        double[] ccy = new double[sliceN];
        double[][] mxy = new double[2][nrois];
        for(int c=0; c<nrois; c++){            
            for(int f=0; f<sliceN; f++){    
                ccx[f] = cx[c][f] + drift[0][f];
                ccy[f] = cy[c][f] + drift[1][f];
            }
            mxy[0][c] = BasicStatisitic.avg(ccx);
            mxy[1][c] = BasicStatisitic.avg(ccy);
        }
        return mxy;
    }       
    
    public static double[][] compute_drift(double[] mx, double[] my){        
        return compute_drift(mx, my, 0);
    }
    
    public static double[][] compute_drift(double[] mx, double[] my, int reference_frame){        
        int sliceN = mx.length;
        double[][] driftXY = new double[2][sliceN];
        for(int i=0; i<sliceN; i++){	
            double dx = mx[reference_frame] - mx[i];
            double dy = my[reference_frame] - my[i];				
            driftXY[0][i] = dx;
            driftXY[1][i] = dy;            
        }
        return driftXY;
    }
    
    public static int check_continuous_tracks(double[][] cx, double[][] cy, boolean[] flags){
        int nrois = cx.length;
        int sliceN = cx[0].length;                
        //boolean[] flags = new boolean[nrois];
        int nrois_checked = 0; 		
        for(int k=0; k<nrois; k++){			
                flags[k] = true;			
                for(int i=0; i<sliceN; i++){	
                    if(flags[k]){
                       if(Double.isNaN(cx[k][i]) || Double.isNaN(cy[k][i])) flags[k] = false;
                    }
                    else break;
                }			
                //best_flags[k] = flags[k];
                if(flags[k]) nrois_checked++;
        }
        return nrois_checked;
    }
    
    public static double[][] average_tracks(double[][] cx, double[][] cy, boolean[] flags){
        int nrois = cx.length;
        int sliceN = cx[0].length;                
        double[] cxt = new double[nrois];
        double[] cyt = new double[nrois];
        double[][] mxy = new double[2][sliceN];
        for(int i=0; i<sliceN; i++){	
            for(int k=0; k<nrois; k++){
                cxt[k] = cx[k][i];
                cyt[k] = cy[k][i];
            }			
            mxy[0][i] = BasicStatisitic.avg(cxt,flags);
            mxy[1][i] = BasicStatisitic.avg(cyt,flags);				
        }

        return mxy;
    }
    
    public static double[][] average_tracks(double[][] cx, double[][] cy){
        int nrois = cx.length;        
        boolean[] flags = new boolean[nrois];
        int nrois_checked = check_continuous_tracks(cx, cy, flags);
        if(nrois_checked<1) return null;   
        return average_tracks(cx, cy, flags);
    }
    
    public static double[] var_tracks(double[][] cx, double[][] cy){
        int nrois = cx.length;        
        double[] var_xy = new double[nrois];
        for(int k=0; k<nrois; k++){
            double varx = BasicStatisitic.mean_sqr(cx[k]);
            double vary = BasicStatisitic.mean_sqr(cy[k]);				
            var_xy[k] = varx + vary;
        }
        return var_xy;
    }
    
    // extort some information from the user
    boolean showDialog() {
            GenericDialogPlus gd = new GenericDialogPlus("Grid Generation");
            gd.addFileField("centroid_file_name:", centroid_name, 30);
            
            //gd.addFileField("centroid_file_name:", centroid_name, 30);
            gd.addNumericField("pixels_size in nm:", lattice, 2);		
            gd.addNumericField("pillar_spacing in pixel:", spacing, 2);
            gd.addNumericField("Grid_oblique ([-45~45] degree)", oblique, 2);	
            gd.addNumericField("Grid_angle ((0~90] degree)", grid_angle, 2);	
            gd.addNumericField("pillar_diameter:", diameter, 2);	
            //output_name = centroid_name + ".dxy";
            gd.addFileField("output_file_name:", output_name, 30);
            gd.showDialog();

            if (gd.wasCanceled()) return false;

            centroid_name=gd.getNextString();
            lattice=gd.getNextNumber();		
            spacing=gd.getNextNumber();		
            oblique=gd.getNextNumber();	
            grid_angle=gd.getNextNumber();	
            diameter=gd.getNextNumber();	
            output_name=gd.getNextString();

            return true;
    }    
    
    // extort some information from the user
    boolean showDialogFD() {
            GenericDialogPlus gd = new GenericDialogPlus("Grid Generation");
            gd.addFileField("input_centroid_file_name:", centroid_name, 30);
            String[] items = {"use reconstructed pillars as reference",
                              "use the pllars in 1st frame as reference",
                              "use intersections of 2nd-order polynomial fit"};
            gd.addChoice("reference choice:", items, items[0]);
            //gd.addCheckbox("use reconstructed pillars as reference?", use_restored_pillars);            
            gd.addMessage("-----settings for reconstructed pillars---------");
            gd.addFileField("input_reconstructed_pillars_file_name:", restored_fname, 30); 
            gd.addMessage("-----settings for grid creation by polynomial fit---------");
            gd.addNumericField("pixels_size in nm:", lattice, 2);		
            gd.addNumericField("pillar_spacing in pixel:", spacing, 2);
            gd.addNumericField("Grid_oblique ([-45~45] degree)", oblique, 2);	
            gd.addNumericField("Grid_angle ((0~90] degree)", grid_angle, 2);	
            gd.addNumericField("pillar_diameter:", diameter, 2);	
            //output_name = centroid_name + ".dxy";
            gd.addMessage("-----output binary file---------");
            gd.addFileField("output_file_name:", output_name, 30);
            gd.showDialog();

            if (gd.wasCanceled()) return false;
            
            centroid_name=gd.getNextString();
            int ref_choice = gd.getNextChoiceIndex();
            use_restored_pillars=(ref_choice==0);
            use_1st_frame_pillars=(ref_choice==1);
            restored_fname=gd.getNextString();            
            //use_restored_pillars = gd.getNextBoolean();
            lattice=gd.getNextNumber();		
            spacing=gd.getNextNumber();		
            oblique=gd.getNextNumber();	
            grid_angle=gd.getNextNumber();	
            diameter=gd.getNextNumber();	
            output_name=gd.getNextString();

            return true;
    }    
    
    private void writeDXY_version2(String fname){
        GridFileHeader header = new GridFileHeader(lattice,diameter, spacing, oblique, grid_angle, Double.NaN);
        writeDXY_version2(fname, header, 
             raw_tracksX, raw_tracksY, 
             corrected_tracksX, corrected_tracksY,
             DX, DY,
             IX, IY);
    }
    
    private void writeDXY_version3(String fname){
        GridFileHeader header = new GridFileHeader(lattice,diameter, spacing, oblique, grid_angle, Double.NaN);
        writeDXY_version3(fname, header, 
             raw_tracksX, raw_tracksY, 
             corrected_tracksX, corrected_tracksY,
             DX, DY);
    }
    
    private void writeDXY(String fname, int version){
        int frames = raw_tracksX[0].length;
        int cens = raw_tracksX.length;
        boolean[][] active = new boolean[cens][frames];
        byte[] active_lin = new byte[frames*cens];
        int j=0;
        for(int c=0; c<cens; c++) {
            for(int f=0; f<frames; f++) {
                   if(Double.isNaN(raw_tracksX[c][f]) || Double.isNaN(raw_tracksY[c][f]) || Double.isNaN(DX[c][f]) || Double.isNaN(DY[c][f])){
                        active_lin[j] = 0;
                        active[c][f] = false;
                   }
                   else{
                       active_lin[j] = 1;
                       active[c][f] = true;
                   }
                   j++;
            }
        }
        //write out result        
        // output file
        try {
                FileOutputStream pout = new FileOutputStream(fname);
                DataOutputStream printout = new DataOutputStream(pout);
                printout.writeInt(version);
                printout.writeDouble(lattice);
                printout.writeDouble(diameter);
                printout.writeDouble(spacing);
                printout.writeDouble(oblique);
                printout.writeDouble(grid_angle);                
                printout.writeInt(frames);
                printout.writeInt(cens);
                printout.write(active_lin, 0, cens*frames);
                for(int c=0; c<cens; c++) {
                        for(int f=0; f<frames; f++) {	
                                if (active[c][f]) {
                                        printout.writeDouble(raw_tracksX[c][f]);
                                        printout.writeDouble(raw_tracksY[c][f]);                                        
                                        printout.writeDouble(corrected_tracksX[c][f]);
                                        printout.writeDouble(corrected_tracksY[c][f]);                                                                                
                                        printout.writeDouble(DX[c][f]);
                                        printout.writeDouble(DY[c][f]);
                                }
                        }
                        IJ.showProgress(c+1, cens);    
                }
                
                for(int c=0; c<cens; c++) {
                        for(int f=0; f<frames; f++) {	
                                if (active[c][f]) {
                                        printout.writeInt(IX[c][f]);
                                        printout.writeInt(IY[c][f]);
                                }
                        }
                }
                printout.close();
        }
        catch(FileNotFoundException fe)
        {
                IJ.log("FileNotFoundException - writing failed : " + fe);
        }
        catch(IOException ioe)
        {
                IJ.log("IOException - writing failed  : " + ioe);
        }				
    }
    
    private static void writeDXY(String fname, int version, double lattice, double diameter, double spacing, double[][] tracksX, double[][] tracksY, double[][] DX, double[][] DY){
        int frames = tracksX[0].length;
        int cens = tracksX.length;
        boolean[][] active = new boolean[cens][frames];
        byte[] active_lin = new byte[frames*cens];
        int j=0;
        for(int c=0; c<cens; c++) {
            for(int f=0; f<frames; f++) {
                   if(Double.isNaN(tracksX[c][f]) || Double.isNaN(tracksY[c][f]) || Double.isNaN(DX[c][f]) || Double.isNaN(DY[c][f])){
                        active_lin[j] = 0;
                        active[c][f] = false;
                   }
                   else{
                       active_lin[j] = 1;
                       active[c][f] = true;
                   }
                   j++;
            }
        }
        //write out result        
        // output file
        try {
                FileOutputStream pout = new FileOutputStream(fname);
                DataOutputStream printout = new DataOutputStream(pout);
                printout.writeInt(version);
                printout.writeDouble(lattice);
                printout.writeDouble(diameter);
                printout.writeDouble(spacing);
                printout.writeInt(frames);
                printout.writeInt(cens);
                printout.write(active_lin, 0, cens*frames);
                for(int c=0; c<cens; c++) {
                        for(int f=0; f<frames; f++) {	
                                if (active[c][f]) {
                                        printout.writeDouble(tracksX[c][f]);
                                        printout.writeDouble(tracksY[c][f]);                                        
                                        printout.writeDouble(DX[c][f]);
                                        printout.writeDouble(DY[c][f]);
                                }
                        }
                }
                printout.close();
        }
        catch(FileNotFoundException fe)
        {
                IJ.log("FileNotFoundException - writing failed : " + fe);
        }
        catch(IOException ioe)
        {
                IJ.log("IOException - writing failed  : " + ioe);
        }				
    }
    
    public static void writeDXY_version2(String fname, GridFileHeader header, 
             double[][] raw_tracksX, double[][] raw_tracksY, 
             double[][] corrected_tracksX, double[][] corrected_tracksY,
             double[][] DX, double[][] DY,
             int[][] IX, int[][] IY)
    {
        int frames = raw_tracksX[0].length;
        int cens = raw_tracksX.length;
        boolean[][] active = new boolean[cens][frames];
        byte[] active_lin = new byte[frames*cens];
        int[] num_active = new int[cens];
        int j=0;
        for(int c=0; c<cens; c++) {
            num_active[c] = 0;
            for(int f=0; f<frames; f++) {
                   if(Double.isNaN(raw_tracksX[c][f]) || Double.isNaN(raw_tracksY[c][f]) || Double.isNaN(DX[c][f]) || Double.isNaN(DY[c][f])){
                        active_lin[j] = 0;
                        active[c][f] = false;
                   }
                   else{
                       active_lin[j] = 1;
                       active[c][f] = true;
                       num_active[c]++;
                   }
                   j++;
            }
        }
        //write out result        
        // output file
        try {
                FileOutputStream pout = new FileOutputStream(fname);
                FileChannel file_ch = pout.getChannel();
                
                int header_size = Double.SIZE*5+Integer.SIZE*3;
                ByteBuffer header_buffer = ByteBuffer.allocate(header_size);
                header_buffer.putInt(GridDataLoader.fileformat_2);
                header_buffer.putDouble(header.lattice);
                header_buffer.putDouble(header.diameter);
                header_buffer.putDouble(header.spacing);
                header_buffer.putDouble(header.grid_oblique);
                header_buffer.putDouble(header.grid_angle);      
                //header_buffer.putDouble(header.deflection_threshold);  
                header_buffer.putInt(frames);
                header_buffer.putInt(cens);
                //header_buffer.rewind();
                header_buffer.flip();
                file_ch.write(header_buffer);
                ByteBuffer active_buffer = ByteBuffer.wrap(active_lin);
                file_ch.write(active_buffer);
                
                for(int c=0; c<cens; c++) {
                    int num = num_active[c];
                    if(num>0){
                        ByteBuffer buff = ByteBuffer.allocate(Double.SIZE*6 * num); 
                        for(int f=0; f<frames; f++) {	
                                    if (active[c][f]) {
                                            buff.putDouble(raw_tracksX[c][f]);
                                            buff.putDouble(raw_tracksY[c][f]);
                                            buff.putDouble(corrected_tracksX[c][f]);
                                            buff.putDouble(corrected_tracksY[c][f]);
                                            buff.putDouble(DX[c][f]);
                                            buff.putDouble(DY[c][f]);                                            
                                    }
                        }
                        //buff.rewind();
                        buff.flip();
                        file_ch.write(buff);
                    }
                    //IJ.showProgress(c+1, cens);    
                }
                
                for(int c=0; c<cens; c++) {
                    int num = num_active[c];
                    if(num>0){
                        ByteBuffer buff = ByteBuffer.allocate(Integer.SIZE*2 * num); 
                        for(int f=0; f<frames; f++) {	
                                    if (active[c][f]) {
                                            buff.putInt(IX[c][f]);
                                            buff.putInt(IY[c][f]);                                             
                                    }
                        }
                        //buff.rewind();
                        buff.flip();
                        file_ch.write(buff);
                    }
                    IJ.showProgress(c+1, cens);    
                }
                file_ch.close();
//                DataOutputStream printout = new DataOutputStream(pout);
//                printout.writeInt(fileformat_2);
//                printout.writeDouble(header.lattice);
//                printout.writeDouble(header.diameter);
//                printout.writeDouble(header.spacing);
//                printout.writeDouble(header.grid_oblique);
//                printout.writeDouble(header.grid_angle);                
//                //printout.writeDouble(header.deflection_threshold);                
//                printout.writeInt(frames);
//                printout.writeInt(cens);
//                printout.write(active_lin, 0, cens*frames);
//                for(int c=0; c<cens; c++) {
//                        for(int f=0; f<frames; f++) {	
//                                if (active[c][f]) {
//                                        printout.writeDouble(raw_tracksX[c][f]);
//                                        printout.writeDouble(raw_tracksY[c][f]);                                        
//                                        printout.writeDouble(corrected_tracksX[c][f]);
//                                        printout.writeDouble(corrected_tracksY[c][f]);                                                                                
//                                        printout.writeDouble(DX[c][f]);
//                                        printout.writeDouble(DY[c][f]);
//                                }
//                        }
//                        IJ.showProgress(c+1, cens);    
//                }
//                
//                for(int c=0; c<cens; c++) {
//                        for(int f=0; f<frames; f++) {	
//                                if (active[c][f]) {
//                                        printout.writeInt(IX[c][f]);
//                                        printout.writeInt(IY[c][f]);
//                                }
//                        }
//                }
//                printout.close();
        }
        catch(FileNotFoundException fe)
        {
                IJ.log("FileNotFoundException - writing failed : " + fe);
        }
        catch(IOException ioe)
        {
                IJ.log("IOException - writing failed  : " + ioe);
        }				
    }
    
    public static void writeDXY_version3(String fname, GridFileHeader header, 
             double[][] raw_tracksX, double[][] raw_tracksY, 
             double[][] corrected_tracksX, double[][] corrected_tracksY,
             double[][] DX, double[][] DY
             )
    {
        int frames = raw_tracksX[0].length;
        int cens = raw_tracksX.length;
        boolean[][] active = new boolean[cens][frames];
        byte[] active_lin = new byte[frames*cens];
        int[] num_active = new int[cens];
        int j=0;
        for(int c=0; c<cens; c++) {
            num_active[c] = 0;
            for(int f=0; f<frames; f++) {
                   if(Double.isNaN(raw_tracksX[c][f]) || Double.isNaN(raw_tracksY[c][f]) || Double.isNaN(DX[c][f]) || Double.isNaN(DY[c][f])){
                        active_lin[j] = 0;
                        active[c][f] = false;
                   }
                   else{
                       active_lin[j] = 1;
                       active[c][f] = true;
                       num_active[c]++;
                   }
                   j++;
            }
        }
        //write out result        
        // output file
        try {
                FileOutputStream pout = new FileOutputStream(fname);
                FileChannel file_ch = pout.getChannel();
                
                int header_size = Double.SIZE*5+Integer.SIZE*3;
                ByteBuffer header_buffer = ByteBuffer.allocate(header_size);
                header_buffer.putInt(GridDataLoader.fileformat_3);
                header_buffer.putDouble(header.lattice);
                header_buffer.putDouble(header.diameter);
                header_buffer.putDouble(header.spacing);
                header_buffer.putDouble(header.grid_oblique);
                header_buffer.putDouble(header.grid_angle);      
                //header_buffer.putDouble(header.deflection_threshold);  
                header_buffer.putInt(frames);
                header_buffer.putInt(cens);
                //header_buffer.rewind();
                header_buffer.flip();
                file_ch.write(header_buffer);
                ByteBuffer active_buffer = ByteBuffer.wrap(active_lin);
                file_ch.write(active_buffer);
                
                for(int c=0; c<cens; c++) {
                    int num = num_active[c];
                    if(num>0){
                        ByteBuffer buff = ByteBuffer.allocate(Double.SIZE*6 * num); 
                        for(int f=0; f<frames; f++) {	
                                    if (active[c][f]) {
                                            buff.putDouble(raw_tracksX[c][f]);
                                            buff.putDouble(raw_tracksY[c][f]);
                                            buff.putDouble(corrected_tracksX[c][f]);
                                            buff.putDouble(corrected_tracksY[c][f]);
                                            buff.putDouble(DX[c][f]);
                                            buff.putDouble(DY[c][f]);                                            
                                    }
                        }
                        //buff.rewind();
                        buff.flip();
                        file_ch.write(buff);
                    }
                    //IJ.showProgress(c+1, cens);    
                }                
                file_ch.close();
        }
        catch(FileNotFoundException fe)
        {
                IJ.log("FileNotFoundException - writing failed : " + fe);
        }
        catch(IOException ioe)
        {
                IJ.log("IOException - writing failed  : " + ioe);
        }				
    }
    
    public static void writeDXY_version4(String fname, GridFileHeader header,              
             boolean[][] active, double[][] raw_CX, double[][] raw_CY, 
             double[][] CX, double[][] CY,
             double[][] DX, double[][] DY,
             int[][] IX, int[][] IY)
     {
         writeDXY_version4(fname, header, active, raw_CX, raw_CY, CX, CY, active, DX, DY, active, IX, IY);
     }
         
     public static void writeDXY_version4(String fname, GridFileHeader header,              
             boolean[][] active, double[][] raw_CX, double[][] raw_CY,
             double[][] CX, double[][] CY,
             boolean[][] DXY_active, double[][] DX, double[][] DY, 
             boolean[][] IXY_active, int[][] IX, int[][] IY)
     {
        int frames = raw_CX[0].length;
        int cens = raw_CX.length;
        //write out result        
        // output file
        try {
                FileOutputStream pout = new FileOutputStream(fname);
                //FileChannel file_ch = pout.getChannel();
                int header_size = Double.SIZE*6+Integer.SIZE*3;
                ByteBuffer header_buffer = ByteBuffer.allocate(header_size);
                header_buffer.putInt(GridDataLoader.fileformat_3);
                header_buffer.putDouble(header.lattice);
                header_buffer.putDouble(header.diameter);
                header_buffer.putDouble(header.spacing);
                header_buffer.putDouble(header.grid_oblique);
                header_buffer.putDouble(header.grid_angle);      
                header_buffer.putDouble(header.deflection_threshold);  
                header_buffer.putInt(frames);
                header_buffer.putInt(cens);   
                //byte[] headByte = new byte[header_buffer.limit()];      
                //header_buffer.get(headByte);                
                //pout.write(headByte);                
                pout.write(header_buffer.array());
                //header_buffer.rewind();
                //file_ch.write(header_buffer);
//                DataOutputStream printout = new DataOutputStream(pout);                
//                printout.writeInt(fileformat_3);
//                printout.writeDouble(header.lattice);
//                printout.writeDouble(header.diameter);
//                printout.writeDouble(header.spacing);
//                printout.writeDouble(header.grid_oblique);
//                printout.writeDouble(header.grid_angle);      
//                printout.writeDouble(header.deflection_threshold);  
//                printout.writeInt(frames);
//                printout.writeInt(cens);                
//                printout.close();
//                printout.write(active_lin, 0, cens*frames);
                
                //int byte_size = Double.SIZE*6 + Integer.SIZE*2 + Byte.SIZE*3;                     
                int byte_size = Double.SIZE*6 + Integer.SIZE*5;                     
                for(int c=0; c<cens; c++) {
                    ByteBuffer buff = ByteBuffer.allocate(byte_size * frames); 
                    for(int f=0; f<frames; f++) {	
//                        byte ta = active[c][f] ? (byte)1 : 0;
//                        byte da = DXY_active[c][f] ? (byte)1 : 0;
//                        byte ia = IXY_active[c][f] ? (byte)1 : 0;
                        int ta = active[c][f] ? 1 : 0;
                        int da = DXY_active[c][f] ? 1 : 0;
                        int ia = IXY_active[c][f] ? 1 : 0;
                       
                        buff.putInt(ta);                        
                        buff.putInt(da);                        
                        buff.putInt(ia);                        
                        buff.putDouble(raw_CX[c][f]);
                        buff.putDouble(raw_CY[c][f]);
                        buff.putDouble(CX[c][f]);
                        buff.putDouble(CY[c][f]);
                        buff.putDouble(DX[c][f]);
                        buff.putDouble(DY[c][f]);
                        buff.putInt(IX[c][f]);
                        buff.putInt(IY[c][f]);                          
                    }      
                    //buff.get(dataByte);    
                    //pout.write(dataByte);
                    //buff.rewind();
                    //file_ch.write(buff);
                    pout.write(buff.array());
                    IJ.showProgress(c+1, cens);  
                }
                //file_ch.close();
                pout.close();
        }
        catch(FileNotFoundException fe)
        {
                IJ.log("FileNotFoundException - writing failed : " + fe);
        }
        catch(IOException ioe)
        {
                IJ.log("IOException - writing failed  : " + ioe);
        }				
        
    }
}
