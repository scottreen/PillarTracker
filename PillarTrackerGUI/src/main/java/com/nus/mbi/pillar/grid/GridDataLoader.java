package com.nus.mbi.pillar.grid;

import com.nus.mbi.pillar.tracker.pillar_tracking_FD;
import com.nus.mbi.pillar.drift.DriftDataLoader;
import com.nus.mbi.pillar.drift.FileHeaderDrift;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.gui.TextRoi;
import ij.process.ColorProcessor;
import ij.process.FloatPolygon;
import java.awt.Color;
import java.awt.Font;
import java.awt.Polygon;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.ByteBuffer;
import static com.nus.mbi.pillar.stat.BasicStatisitic.max;

/**
 *
 * @author xiaochun
 */
public class GridDataLoader {
    /******** File Structure ****************
    int outchecksum
    double lattice
    double diameter
    double spacing
    int frames
    int cens
    byte[cens*frames] active_lin

    repeat for each active[c][f]:
            double trackX[c][f]
            double trackY[c][f]
            double dx[c][f]
            double dy[c][f]          
            
    ************************************************/
    boolean suc_load = false;
    public static final int fileformat_1=200700; // version 1 for deflection file
    public static final int fileformat_2=200701; // version 2 for deflection file
    public static final int fileformat_3=200702; // version 3 for deflection file
    
    public final double empty=1.0e20; // version 1
    
    // image and pillar params
    public int   frames;
    public int   cens; // # og centroids
    public double lattice;
    public double diameter;
    public double spacing;
    public int pillars_selected;

    public boolean[][] active;
    public double[][] trackX;
    public double[][] trackY;
    public double[][] CX;
    public double[][] CY;
    public boolean[][] DXY_active;
    public double[][] DX;
    public double[][] DY;   
    public boolean[][] IXY_active;
    public int[][] IX = null;
    public int[][] IY = null;    
    public int num_active = 0;

    public int label_size = 8;
    public Color label_color = Color.yellow;
    
    public double grid_oblique = Double.NaN;
    public double grid_angle = 90;
    public double deflection_threshold = Double.NaN;
    private int version;
    private String file_name;
    public String getFileName(){
        return file_name;
    }
    public int getFileVersion(){
        return version;
    }
    
    public boolean data_loader(String intext){        
        suc_load = false;
        grid_oblique = Double.NaN;
        grid_angle = 90;
        CX = null;
        CY = null;
        active = null;        
        DX = null;
        DY = null;
        DXY_active = null;        
        IX = null;
        IY = null;
        IXY_active = null;        
        FileInputStream pin;
        //DataInputStream readbin;
        try {
                pin = new FileInputStream(intext);                
                version=readInt(pin);
                if (version==fileformat_1) {
                        lattice=readDouble(pin);
                        diameter=readDouble(pin);
                        spacing=readDouble(pin);
                        frames=readInt(pin);
                        cens=readInt(pin);
                        
                        if(spacing>lattice) spacing = spacing/lattice;
                        if(diameter>lattice) diameter = diameter/lattice;
                        IJ.log("lattice=" + lattice + "	dia=" + diameter + " spacing=" + spacing + " frames=" + frames + " cents= " + cens);
                        
                        byte[] active_lin = readBytes(pin,cens*frames);                        
                        
                        active=new boolean[cens][frames];
                        trackX=new double[cens][frames];
                        trackY=new double[cens][frames];
                        DX=new double[cens][frames];
                        DY=new double[cens][frames];
                        //arrow_valid=new boolean[cens][frames];       
                        
                        num_active=0;
                        int j=0;                          
                        for(int c=0; c<cens; c++) {
                            int n_active = 0;
                            for(int f=0; f<frames; f++) {
                                        active[c][f]=(active_lin[j]==1); j++;                                        
                                        if(active[c][f]) n_active++;
                            }
                            num_active += n_active;
                            double[][] matrix = readfile2matrix(pin, 4, n_active);
                            n_active = 0;
                            for(int f=0; f<frames; f++) {                                        
                                    if(active[c][f]) {
                                            trackX[c][f]=matrix[0][n_active];
                                            trackY[c][f]=matrix[1][n_active];
                                            DX[c][f]=matrix[2][n_active];
                                            DY[c][f]=matrix[3][n_active];
                                            n_active++;
                                            //num_active++;
                                    }else{
                                        trackX[c][f]=Double.NaN;
                                        trackY[c][f]=Double.NaN;
                                        DX[c][f]=Double.NaN;
                                        DY[c][f]=Double.NaN;
                                    }
                            }
                            IJ.showProgress(c, cens);
                            IJ.showStatus("loading grid data... " + Math.round(c*100/cens) + "%");
                        }               
                        IJ.showStatus("grid data loaded");
                        suc_load = true;                                                
                }
                else if (version==fileformat_2) {
                        lattice=readDouble(pin);
                        diameter=readDouble(pin);
                        spacing=readDouble(pin);
                        grid_oblique=readDouble(pin);
                        grid_angle=readDouble(pin);
                        frames=readInt(pin);
                        cens=readInt(pin);
                        
                        if(spacing>lattice) spacing = spacing/lattice;
                        if(diameter>lattice) diameter = diameter/lattice;
                        IJ.log("lattice=" + lattice + "	dia=" + diameter + " spacing=" + spacing + " frames=" + frames + " cents= " + cens);
                        
                        byte[] active_lin = readBytes(pin,cens*frames);                        
                        
                        active=new boolean[cens][frames];
                        trackX=new double[cens][frames];
                        trackY=new double[cens][frames];
                        CX=new double[cens][frames];
                        CY=new double[cens][frames];
                        DX=new double[cens][frames];
                        DY=new double[cens][frames];
                        IX=new int[cens][frames];
                        IY=new int[cens][frames];
                        //arrow_valid=new boolean[cens][frames];       
                        
                        num_active=0;
                        int j=0;                          
                        for(int c=0; c<cens; c++) {
                            int n_active = 0;
                            for(int f=0; f<frames; f++) {
                                        active[c][f]=(active_lin[j]==1); j++;                                        
                                        if(active[c][f]) n_active++;
                            }
                            num_active += n_active;
                            double[][] matrix = readfile2matrix(pin, 6, n_active);
                            n_active = 0;
                            for(int f=0; f<frames; f++) {                                        
                                    if(active[c][f]) {
                                            trackX[c][f]=matrix[0][n_active];
                                            trackY[c][f]=matrix[1][n_active]; 
                                            CX[c][f]=matrix[2][n_active];
                                            CY[c][f]=matrix[3][n_active];
                                            DX[c][f]=matrix[4][n_active];
                                            DY[c][f]=matrix[5][n_active];
                                            n_active++;
                                            //num_active++;
                                    }else{
                                        trackX[c][f]=Double.NaN;
                                        trackY[c][f]=Double.NaN;
                                        CX[c][f]=Double.NaN;
                                        CY[c][f]=Double.NaN;
                                        DX[c][f]=Double.NaN;
                                        DY[c][f]=Double.NaN;
                                    }
                            }
                            IJ.showProgress(c, cens);
                            IJ.showStatus("loading grid data... " + Math.round(c*100/cens) + "%");
                        }
                                                
                        for(int c=0; c<cens; c++) {
                            int n_active = 0;
                            for(int f=0; f<frames; f++) {                                                                             
                                if(active[c][f]) n_active++;
                            }
                            int[][] matrix = readfile2Int_matrix(pin, 2, n_active);
                            n_active = 0;
                            for(int f=0; f<frames; f++) {	
                                if (active[c][f]) {
                                        IX[c][f]=matrix[0][n_active];
                                        IY[c][f]=matrix[1][n_active];
                                        n_active++;
                                }else{
                                        IX[c][f]=Integer.MAX_VALUE;
                                        IY[c][f]=Integer.MAX_VALUE;
                                    }
                            }
                        }
                        IJ.showStatus("grid data loaded");
                        suc_load = true;                                                
                }
                else if (version==fileformat_3) {
                        lattice=readDouble(pin);
                        diameter=readDouble(pin);
                        spacing=readDouble(pin);
                        grid_oblique=readDouble(pin);
                        grid_angle=readDouble(pin);
                        frames=readInt(pin);
                        cens=readInt(pin);
                        
                        if(spacing>lattice) spacing = spacing/lattice;
                        if(diameter>lattice) diameter = diameter/lattice;
                        IJ.log("lattice=" + lattice + "	dia=" + diameter + " spacing=" + spacing + " frames=" + frames + " cents= " + cens);
                        
                        byte[] active_lin = readBytes(pin,cens*frames);                        
                        
                        active=new boolean[cens][frames];
                        trackX=new double[cens][frames];
                        trackY=new double[cens][frames];
                        CX=new double[cens][frames];
                        CY=new double[cens][frames];
                        DX=new double[cens][frames];
                        DY=new double[cens][frames];                        
                        //arrow_valid=new boolean[cens][frames];       
                        
                        num_active=0;
                        int j=0;                          
                        for(int c=0; c<cens; c++) {
                            int n_active = 0;
                            for(int f=0; f<frames; f++) {
                                        active[c][f]=(active_lin[j]==1); j++;                                        
                                        if(active[c][f]) n_active++;
                            }
                            num_active += n_active;
                            double[][] matrix = readfile2matrix(pin, 6, n_active);
                            n_active = 0;
                            for(int f=0; f<frames; f++) {                                        
                                    if(active[c][f]) {
                                            trackX[c][f]=matrix[0][n_active];
                                            trackY[c][f]=matrix[1][n_active]; 
                                            CX[c][f]=matrix[2][n_active];
                                            CY[c][f]=matrix[3][n_active];
                                            DX[c][f]=matrix[4][n_active];
                                            DY[c][f]=matrix[5][n_active];
                                            n_active++;
                                            //num_active++;
                                    }else{
                                        trackX[c][f]=Double.NaN;
                                        trackY[c][f]=Double.NaN;
                                        CX[c][f]=Double.NaN;
                                        CY[c][f]=Double.NaN;
                                        DX[c][f]=Double.NaN;
                                        DY[c][f]=Double.NaN;
                                    }
                            }
                            IJ.showProgress(c, cens);
                            IJ.showStatus("loading grid data... " + Math.round(c*100/cens) + "%");
                        }                        
                        IJ.showStatus("grid data loaded");
                        suc_load = true;                                                             
                }
                else if (version==pillar_tracking_FD.fileversion2) {
                    FileHeaderDrift header = DriftDataLoader.read_header_ver2(pin);
                    lattice=header.lattice;
                    diameter=header.diameter;
                    spacing=header.spacing;
                    grid_oblique=header.oblique;
                    grid_angle=header.grid_angle;
                    frames=header.nframes;
                    cens=header.npillars;
                        
                    if(spacing>lattice) spacing = spacing/lattice;
                    if(diameter>lattice) diameter = diameter/lattice;
                    IJ.log("lattice=" + lattice + "	dia=" + diameter + " spacing=" + spacing + " frames=" + frames + " cents= " + cens);

                    active=new boolean[cens][frames];
                    trackX=new double[cens][frames];
                    trackY=new double[cens][frames];                        
                    DX=new double[cens][frames];
                    DY=new double[cens][frames];                        
                    
                    double[][] matrix_raw = readfile2matrix(pin, cens*2, frames);
                    double[][] matrix_DXY = readfile2matrix(pin, cens*2, frames);                    
                    for(int c=0; c<cens; c++) {
                        int n_active = 0;
                        for(int f=0; f<frames; f++) {                                        
                            trackX[c][f]=lattice*matrix_raw[2*c][f];
                            trackY[c][f]=lattice*matrix_raw[2*c+1][f];                                 
                            DX[c][f]=lattice*matrix_DXY[2*c][f];
                            DY[c][f]=lattice*matrix_DXY[2*c+1][f];    

                            active[c][f] = false;
                            if(!Double.isNaN(trackX[c][f]) && !Double.isNaN(DX[c][f])) {
                                active[c][f] = true;
                                n_active++;                                            
                            }
                        }
                        num_active += n_active;        
                        IJ.showProgress(c, cens);
                        IJ.showStatus("loading grid data... " + Math.round(c*100/cens) + "%");
                    }

                    double[][] matrix_corrected = readfile2matrix(pin, cens*2, frames);
                    //double[][] matrix_driftXY = readfile2matrix(pin, cens*2, frames);
                    if(matrix_corrected!=null){
                        CX=new double[cens][frames];
                        CY=new double[cens][frames];                                                                    
                        for(int c=0; c<cens; c++) {                        
                            for(int f=0; f<frames; f++) {                                        
                                CX[c][f]=lattice*matrix_corrected[2*c][f];
                                CY[c][f]=lattice*matrix_corrected[2*c+1][f];                                
                            }
                        }
                    }

                    IJ.showStatus("grid data loaded");
                    suc_load = true;                                                             
                }
                else{
                    IJ.log("file version error");
                }
                pin.close();   
        }                
        catch(FileNotFoundException fe) {
                //consistent=false;
                IJ.log("file path not found->" + intext);
        }
        catch(IOException ioe) {
                //consistent=false;
                IJ.log("file IO error");
        }
        catch(Exception ioe) {
                IJ.log("file reading error");
        }
        finally{
            
        }
        
        if(suc_load) file_name = intext;
        
        return suc_load;
    }
    
    public boolean load_settings(String intext){        
        boolean suc_load_settings = false;
        grid_oblique = Double.NaN;
        grid_angle = 90; 
        FileInputStream pin;       
        try {
                pin = new FileInputStream(intext);                
                version=readInt(pin);
                if(version==fileformat_1){
                        lattice=readDouble(pin);
                        diameter=readDouble(pin);
                        spacing=readDouble(pin);
                        frames=readInt(pin);
                        cens=readInt(pin);
                        
                        if(spacing>lattice) spacing = spacing/lattice;
                        if(diameter>lattice) diameter = diameter/lattice;
                        IJ.log("lattice=" + lattice + "	dia=" + diameter + " spacing=" + spacing + " frames=" + frames + " cents= " + cens);
                        
                        IJ.showStatus("grid settings loaded");
                        suc_load_settings = true;                                                
                }
                else if(version==fileformat_2){
                        lattice=readDouble(pin);
                        diameter=readDouble(pin);
                        spacing=readDouble(pin);
                        grid_oblique=readDouble(pin);
                        grid_angle=readDouble(pin);
                        frames=readInt(pin);
                        cens=readInt(pin);
                        
                        if(spacing>lattice) spacing = spacing/lattice;
                        if(diameter>lattice) diameter = diameter/lattice;
                        IJ.log("lattice=" + lattice + "	dia=" + diameter + " spacing=" + spacing + " frames=" + frames + " cents= " + cens);
                        IJ.showStatus("grid settings loaded");
                        suc_load_settings = true;                                                
                }
                else if(version==fileformat_3){
                        lattice=readDouble(pin);
                        diameter=readDouble(pin);
                        spacing=readDouble(pin);
                        grid_oblique=readDouble(pin);
                        grid_angle=readDouble(pin);
                        deflection_threshold = readDouble(pin);
                        frames=readInt(pin);
                        cens=readInt(pin);
                        
                        IJ.log("lattice=" + lattice + "	dia=" + diameter + " spacing=" + spacing + " frames=" + frames + " cents= " + cens);
                        IJ.showStatus("grid settings loaded");
                        suc_load_settings = true;                                                
                }
                else{
                    IJ.log("file version error");
                }
                pin.close();
        }                
        catch(FileNotFoundException fe) {
                //consistent=false;
                IJ.log("file path not found->" + intext);
        }
        catch(IOException ioe) {
                //consistent=false;
                IJ.log("file IO error");
        }
        catch(Exception ioe) {
                IJ.log("file reading error");
        }
        finally{
            
        }
        
        IJ.log("try loading settings from grid file sucessful? " + suc_load_settings);  
        return suc_load_settings;
    }
//    
//    public boolean data_loader_slow(String intext){        
//        suc_load = false;
//        grid_oblique = Double.NaN;
//        grid_angle = 90;
//        FileInputStream pin;
//        DataInputStream readbin;
//        try {
//                pin = new FileInputStream(intext);
//                readbin = new DataInputStream(pin);
//                int checksum=readbin.readInt();
//                if (checksum==fileformat_1) {
//                        lattice=readbin.readDouble();
//                        diameter=readbin.readDouble();
//                        spacing=readbin.readDouble();
//                        frames=readbin.readInt();
//                        cens=readbin.readInt();
//                        IJ.log("lattice=" + lattice + "	dia=" + diameter + " spacing=" + spacing + " frames=" + frames + " cents= " + cens);
//                        byte[] active_lin = new byte[cens*frames];
//                        readbin.read(active_lin,0,frames*cens);
//                        
//                        active=new boolean[cens][frames];
//                        trackX=new double[cens][frames];
//                        trackY=new double[cens][frames];
//                        DX=new double[cens][frames];
//                        DY=new double[cens][frames];
//                        //arrow_valid=new boolean[cens][frames];
//                        
//                        int j=0; 
//                        num_active=0;
//                        for(int c=0; c<cens; c++) {
//                                for(int f=0; f<frames; f++) {
//                                        active[c][f]=(active_lin[j]==1); j++;
//                                        //arrow_valid[c][f]=false;
//                                        if(active[c][f]) {
//                                                trackX[c][f]=readbin.readDouble();
//                                                trackY[c][f]=readbin.readDouble();
//                                                DX[c][f]=readbin.readDouble();
//                                                DY[c][f]=readbin.readDouble();
//                                                num_active++;
//                                        }else{
//                                            trackX[c][f]=Double.NaN;
//                                            trackY[c][f]=Double.NaN;
//                                            DX[c][f]=Double.NaN;
//                                            DY[c][f]=Double.NaN;
//                                        }
//                                }
//                                IJ.showProgress(c, cens);
//                                IJ.showStatus("loading grid data... " + Math.round(c*100/cens) + "%");
//                        }               
//                        IJ.showStatus("grid data loaded");
//                        suc_load = true;
//                }
//                readbin.close();
//                pin.close();
//        }                
//        catch(FileNotFoundException fe) {
//                //consistent=false;
//                IJ.log("file path not found->" + intext);
//        }
//        catch(IOException ioe) {
//                //consistent=false;
//                IJ.log("file IO error");
//        }
//        finally{
//            
//        }
//        
//        return suc_load;
//    }
//    
    public static int readInt(FileInputStream fis) throws Exception
    {	            
            int size_int = Integer.BYTES;
            byte[] header = new byte[size_int];
            int bytes = fis.read(header);
            if(bytes<=0) return Integer.MIN_VALUE;
            //  create a byte buffer and wrap the array
            ByteBuffer bb = ByteBuffer.wrap(header); 
            return bb.getInt();
    }
    
    public static byte[] readBytes(FileInputStream fis, int size) throws Exception
    {	    
            byte[] header = new byte[size];
            int bytes = fis.read(header);
            IJ.log("reading bytes="+bytes);
            if(bytes<=0) return null;           
            return header;
    }
    
    public static double readDouble(FileInputStream fis) throws Exception
    {	            
            int size_double = Double.BYTES;
            byte[] header = new byte[size_double];
            int bytes = fis.read(header);
            if(bytes<=0) return Double.NaN;
            //  create a byte buffer and wrap the array
            ByteBuffer bb = ByteBuffer.wrap(header); 
            return bb.getDouble();
    }
    
    public static double[][] readfile2matrix(FileInputStream fis, int ncols, int nrows) throws Exception
    {
            int size_double = Double.BYTES;		
            byte[] data = new byte[ncols*nrows*size_double];		
            int bytes = fis.read(data);	
            //IJ.log("reading bytes="+bytes);
            if(bytes<=0) return null;
            ByteBuffer bb = ByteBuffer.wrap(data);           

            double[][] matrix = new double[ncols][nrows];
            for(int i=0; i<nrows; i++){
                    for(int j=0; j<ncols; j++){							
                        matrix[j][i] = bb.getDouble();                        
                        //IJ.log("	 " + matrix[j][i]);			                        
                    }
            }

            //fis.close();
            //IJ.log("loading matrix file done");
            return matrix;
    }
    
    public static int[][] readfile2Int_matrix(FileInputStream fis, int ncols, int nrows) throws Exception
    {
            int size_int = Integer.BYTES;		
            byte[] data = new byte[ncols*nrows*size_int];		
            int bytes = fis.read(data);	
            //IJ.log("reading bytes="+bytes);
            if(bytes<=0) return null;
            ByteBuffer bb = ByteBuffer.wrap(data);           

            int[][] matrix = new int[ncols][nrows];
            for(int i=0; i<nrows; i++){
                    for(int j=0; j<ncols; j++){							
                        matrix[j][i] = bb.getInt();                        
                        //IJ.log("	 " + matrix[j][i]);			                        
                    }
            }

            //fis.close();
            //IJ.log("loading matrix file done");
            return matrix;
    }
    
    public int[] getSelections(Roi pr0){        
        int[] candidates=new int[cens];
        int ccount=0;
        if(PointRoi.class.isInstance(pr0)) {// handle points differently
                Polygon po=pr0.getPolygon();
                double sep=0.5*spacing;
                int RR2=(int)(sep*sep+0.5);
                for(int c=0; c<cens; c++) {
                        int f=0;
                        while ((!active[c][f]) && (f<frames-1)) f++;
                        //if(!active[c][f]) {IJ.showMessage("deflection file erroneous!"); return;}
                        if(!active[c][f]) continue;
                        int X0=(int)(trackX[c][f]/lattice);
                        int Y0=(int)(trackY[c][f]/lattice);
                        int pp=0;
                        boolean found=false;
                        while((pp<po.npoints) && (!found)) {
                                int dx=X0-po.xpoints[pp];
                                int dy=Y0-po.ypoints[pp];
                                int d2=dx*dx+dy*dy;
                                found=(d2<RR2);
                                pp++;
                        }
                        if(found) {//found a pillar
                                candidates[ccount]=c;
                                ccount++;
                        }
                }
        } else {
                for(int c=0; c<cens; c++) {
                        int f=0;
                        while ((!active[c][f]) && (f<frames-1)) f++;
                        //if(!active[c][f]) {IJ.showMessage("deflection file erroneous!"); return;}	
                        if(!active[c][f]) continue;
                        if(pr0.contains((int)(trackX[c][f]/lattice),(int)(trackY[c][f]/lattice))) {//found a pillar
                                candidates[ccount]=c;
                                ccount++;
                        }
                }
        }
        //IJ.log("PILLAR DEFLECTION PLOT f="+frames+" pillars="+cens+" lattice="+lattice+" diameter="+diameter+" spacing="+spacing+" total pillars in all frames="+num_active+" candidates="+ccount);
        if(ccount<1) return null;
        
        int [] subset = new int[ccount];
        for(int i=0; i<ccount; i++) subset[i] = candidates[i];
        return subset;
    }
    
    public int[] getSelections(Roi pr0, double zoom){        
        int[] candidates=new int[cens];
        int ccount=0;
        if(PointRoi.class.isInstance(pr0)) {// handle points differently
                FloatPolygon po=pr0.getFloatPolygon();                
                double sep=0.5*spacing;
                double RR2=sep*sep;
                for(int c=0; c<cens; c++) {
                        int f=0;
                        while ((!active[c][f]) && (f<frames-1)) f++;
                        //if(!active[c][f]) {IJ.showMessage("deflection file erroneous!"); return;}
                        if(!active[c][f]) continue;
                        double X0=trackX[c][f]/lattice;
                        double Y0=trackY[c][f]/lattice;
                        int pp=0;
                        boolean found=false;
                        while((pp<po.npoints) && (!found)) {
                                double dx=X0-po.xpoints[pp]/zoom;
                                double dy=Y0-po.ypoints[pp]/zoom;
                                double d2=dx*dx+dy*dy;
                                found=(d2<RR2);
                                pp++;
                        }
                        if(found) {//found a pillar
                                candidates[ccount]=c;
                                ccount++;
                        }
                }
        } else {
                for(int c=0; c<cens; c++) {
                        int f=0;
                        while ((!active[c][f]) && (f<frames-1)) f++;
                        //if(!active[c][f]) {IJ.showMessage("deflection file erroneous!"); return;}	
                        if(!active[c][f]) continue;
                        if(pr0.contains((int)(trackX[c][f]*zoom/lattice),(int)(trackY[c][f]*zoom/lattice))) {//found a pillar
                                candidates[ccount]=c;
                                ccount++;
                        }
                }
        }
        //IJ.log("PILLAR DEFLECTION PLOT f="+frames+" pillars="+cens+" lattice="+lattice+" diameter="+diameter+" spacing="+spacing+" total pillars in all frames="+num_active+" candidates="+ccount);
        if(ccount<1) return null;
        
        int [] subset = new int[ccount];
        for(int i=0; i<ccount; i++) subset[i] = candidates[i];
        return subset;
    }
    
    public int[] remove_large_deflection(int[] candidates, int frame, double threshold){        
        if(candidates==null || frame<0 || frame>=frames) return null;
        int n = candidates.length;
        int f = frame;
        int na = 0;
        double thrshold2 = threshold*threshold;
        
        boolean[] flag = new boolean[n];
        for(int jj=0; jj<n; jj++) {
            flag[jj] = false;    
            int c=candidates[jj];               
            if(active[c][f]){
                double x = DX[c][f];
                double y = DY[c][f];
                double d2 = x*x + y*y;
                if(d2<thrshold2){
                    flag[jj] = true;
                    na++;
                }
            }                
        }
        
        if(na<1) return null;
        
        int [] subset = new int[na];  
        int ccount=0;
        for(int i=0; i<n; i++){
            if(flag[i]){
                subset[ccount] = candidates[i];
                ccount++;
            }
        }
        
        return subset;
    }
    
    public int[] remove_large_deflection(int[] candidates, double[] dx, double[] dy, boolean[] active_d, double threshold){        
        if(candidates==null) return null;
        int n = candidates.length;        
        int na = 0;
        double thrshold2 = threshold*threshold;        
        boolean[] flag = new boolean[n];
        for(int jj=0; jj<n; jj++) {
            flag[jj] = false;    
            int c=candidates[jj];               
            if(active_d[c]){
                double x = dx[c];
                double y = dy[c];
                double d2 = x*x + y*y;
                if(d2<thrshold2){
                    flag[jj] = true;
                    na++;
                }
            }                
        }
        
        if(na<1) return null;
        
        int [] subset = new int[na];  
        int ccount=0;
        for(int i=0; i<n; i++){
            if(flag[i]){
                subset[ccount] = candidates[i];
                ccount++;
            }
        }
        
        return subset;
    }
    
    public boolean[] get_CU_flag(int[] candidates, boolean[] cu_flag, int frame){
        if(candidates==null || frame<0 || frame>=frames) return null;
        int n = candidates.length;
        int f = frame;
        int na = 0;        
        for(int jj=0; jj<n; jj++) {
                int c=candidates[jj];
                if(active[c][f]) na++;                
        }
        if(na<1) return null;
        
        boolean[] new_cu_flag = new boolean[na];
        na = 0;
        for(int jj=0; jj<n; jj++) {
                int c=candidates[jj];
                if(active[c][f]) {
                        new_cu_flag[na] = cu_flag[jj];
                        na++;
                }
        }
        
        return new_cu_flag;
    }
    
    public double[][] getGrids(int[] candidates, int frame){
        if(candidates==null || frame<0 || frame>=frames) return null;
        int n = candidates.length;
        int f = frame;
        int na = 0;        
        for(int jj=0; jj<n; jj++) {
                int c=candidates[jj];
                if(active[c][f]) na++;                
        }
        if(na<1) return null;
        
        double[][] gridXY = new double[2][na];
        na = 0;
        for(int jj=0; jj<n; jj++) {
                int c=candidates[jj];
                if(active[c][f]) {
                        double x = (trackX[c][f]-DX[c][f])/lattice;
                        double y = (trackY[c][f]-DY[c][f])/lattice;                        
                        gridXY[0][na] = x;
                        gridXY[1][na] = y;
                        na++;
                }
        }
        
        return gridXY;
    }
    
    public double[][] getGrids(int frame){
        if(frame<0 || frame>=frames) return null;
        int n = cens;
        int f = frame;
        int na = 0;        
        for(int c=0; c<n; c++) {                
                if(active[c][f]) na++;                
        }
        if(na<1) return null;
        
        double[][] gridXY = new double[2][na];
        na = 0;
        for(int c=0; c<n; c++) {                
                if(active[c][f]) {
                        double x = (trackX[c][f]-DX[c][f])/lattice;
                        double y = (trackY[c][f]-DY[c][f])/lattice;                        
                        gridXY[0][na] = x;
                        gridXY[1][na] = y;
                        na++;
                }
        }
        
        return gridXY;
    }
    
    public double[][] getGridXY(int frame){
        if(frame<0 || frame>=frames) return null;
        int n = cens;
        int f = frame;
        
        double[][] gridXY = new double[2][n];        
        for(int c=0; c<n; c++) {                
            if(active[c][f]) {
                    double x = (trackX[c][f]-DX[c][f])/lattice;
                    double y = (trackY[c][f]-DY[c][f])/lattice;                        
                    gridXY[0][c] = x;
                    gridXY[1][c] = y;
            }
            else{
                gridXY[0][c] = Double.NaN;
                gridXY[1][c] = Double.NaN;
            }
        }
        
        return gridXY;
    }
    
    public static double[][] getGrids(double[] xc, double[] yc, double[] ddx, double[] ddy){
        int n = xc.length;
        
        double[][] gridXY = new double[2][n];        
        for(int c=0; c<n; c++) {                
            double x = (xc[c]-ddx[c]);
            double y = (yc[c]-ddy[c]);                        
            gridXY[0][c] = x;
            gridXY[1][c] = y;
        }
        
        return gridXY;
    }
    
    public int[][] getGridsIXY(int frame){
        if(frame<0 || frame>=frames) return null;
        int n = cens;
        int f = frame;
        int na = 0;        
        for(int c=0; c<n; c++) {                
                if(active[c][f]) na++;                
        }
        if(na<1) return null;
        
        int[][] gridIXY = new int[2][na];
        na = 0;
        for(int c=0; c<n; c++) {                
                if(active[c][f]) {
                        gridIXY[0][na] = IX[c][f];
                        gridIXY[1][na] = IY[c][f];
                        na++;
                }
        }
        
        return gridIXY;
    }
    
    public double[][] getFrameXY(int frame){        
        int f = frame;
        double[][] XY = new double[2][cens];        
        for(int c=0; c<cens; c++) {                                    
            XY[0][c] = trackX[c][f]/lattice;
            XY[1][c] = trackY[c][f]/lattice;
        }
        
        return XY;
    }
    
    public double[][] getFrameDXY(int frame){        
        int f = frame;
        double[][] XY = new double[2][cens];        
        for(int c=0; c<cens; c++) {                                    
            XY[0][c] = DX[c][f]/lattice;
            XY[1][c] = DY[c][f]/lattice;
        }
        
        return XY;
    }
    
    public void updateFrameDXY(int frame, double[] dx, double[] dy, boolean[] active_d){        
        int f = frame;
        for(int c=0; c<cens; c++) {                                    
            DX[c][f] = dx[c];
            DY[c][f] = dy[c];
            active[c][f] = active_d[c];
        }    
    }
    
    public int[][] getFrameIXY(int frame){        
        int f = frame;
        int[][] XY = new int[2][cens];        
        for(int c=0; c<cens; c++) {                                    
            XY[0][c] = IX[c][f];
            XY[1][c] = IY[c][f];
        }
        
        return XY;
    }
    
    public boolean[] getFrameActive(int frame){        
        int f = frame;
        boolean[] is_labeled = new boolean[cens];        
        for(int c=0; c<cens; c++) is_labeled[c] = active[c][f];
        return is_labeled;
    }
    
    public double[][] getPillarXY(int pillar){        
        int c = pillar;
        double[][] XY = new double[2][frames];        
        for(int f=0; f<frames; f++) {                                    
            XY[0][f] = trackX[c][f];
            XY[1][f] = trackY[c][f];
        }
        
        return XY;
    }
    
    public double[][] getPillarDXY(int pillar){        
        int c = pillar;
        double[][] XY = new double[3][frames];        
        for(int f=0; f<frames; f++) {                                    
            double dx = DX[c][f];
            double dy = DY[c][f];
            XY[0][f] = dx;
            XY[1][f] = dy;
            XY[2][f] = Math.sqrt(dx*dx + dy*dy);
        }
        
        return XY;
    }
    
    public double[] getFrameSeries(){
        double[] xaxis = new double[frames];
        for(int s=0;s<frames;s++) xaxis[s] = s+1;	
        return xaxis;
    } 
    
    public double[] getXY(int ipillar, int frame){            
        if(ipillar<0 || ipillar>=cens) return null;
        if(frame<0 || frame>=frames) return null;
        double[] xy = new double[2];
        xy[0] = trackX[ipillar][frame]/lattice;
        xy[1] = trackY[ipillar][frame]/lattice;
        
        return xy;
    }
    
    public double[] getXY(int ipillar){ 
        return getXY(ipillar, 0);
    }
    
    public double[] getGridXY(int ipillar, int frame){            
        if(ipillar<0 || ipillar>=cens) return null;
        if(frame<0 || frame>=frames) return null;
        double[] xy = new double[2];
        xy[0] = (trackX[ipillar][frame]-DX[ipillar][frame])/lattice;
        xy[1] = (trackY[ipillar][frame]-DY[ipillar][frame])/lattice;
        
        return xy;
    }
    
    public double[][] getXY(int[] candidates, int frame){
        if(candidates==null || frame<0 || frame>=frames) return null;
        int n = candidates.length;
        int f = frame;
        int na = 0;        
        for(int jj=0; jj<n; jj++) {
                int c=candidates[jj];
                if(active[c][f]) na++;                
        }
        if(na<1) return null;
        
        double[][] XY = new double[2][na];
        na = 0;
        for(int jj=0; jj<n; jj++) {
                int c=candidates[jj];
                if(active[c][f]) {
                        double x = trackX[c][f]/lattice;
                        double y = trackY[c][f]/lattice;                        
                        XY[0][na] = x;
                        XY[1][na] = y;
                        na++;
                }
        }
        
        return XY;
    }
    
    public double[] getXYMAX(){
        double xmax = 0;
        double ymax = 0;
        for(int i=0; i<cens; i++){
                double x = max(trackX[i]);
                double y = max(trackY[i]);
                if(!Double.isNaN(x) && !Double.isNaN(y)){
                        if(x>xmax) xmax = x;
                        if(y>ymax) ymax = y;
                }
        }
        double[] xy = {xmax, ymax};
        return xy;
    }
    
    public double getXMAX(){
        double xmax = 0;       
        for(int i=0; i<cens; i++){
                double x = max(trackX[i]);                
                if(!Double.isNaN(x)){
                    if(x>xmax) xmax = x;
                }                        
                
        }        
        return xmax;
    }
    
    public double getYMAX(){
        double ymax = 0;       
        for(int i=0; i<cens; i++){
                double y = max(trackY[i]);                
                if(!Double.isNaN(y)){
                    if(y>ymax) ymax = y;
                }                        
                
        }        
        return ymax;
    }
    
    public ImagePlus create_empty_image(int width, int height, String title){
		ColorProcessor tracks_overlays_image = new ColorProcessor(width, height); 
		ImagePlus tracks_overlays_plus = new ImagePlus(title, tracks_overlays_image);				
		return tracks_overlays_plus;
    }
    
    public ImagePlus create_empty_map(int margin, String title){
            return create_empty_map(1, margin, title);
    }
    
    public ImagePlus create_empty_map(double zoomin_image, int margin, String title){
            double[] xymax = getXYMAX();
            double xmax = xymax[0];
            double ymax = xymax[1];
            int map_img_width = (int)Math.ceil(zoomin_image*xmax/lattice + margin);
            int map_img_height = (int)Math.ceil(zoomin_image*ymax/lattice + margin);
            IJ.log("xmax="+ xmax + " ymax=" + ymax + " map_img_width="+ map_img_width + " map_img_height=" + map_img_height);
            return create_empty_image(map_img_width, map_img_height, title);
    }
    
    public Overlay getLabels(int frame, double image_zoom, boolean showdxy){
            int nrois = cens;
            int font_size = label_size;//(int)Math.ceil(5);
            Font font = new Font("Arial",Font.PLAIN, font_size);
            Overlay overlay = new Overlay();            
            for(int i=0; i<nrois; i++){
                //double rcx = trackX[i][0];//avg(cx[i])/pixel_size;			
                //double rcy = trackY[i][0];//avg(cy[i])/pixel_size;
                double rcx = trackX[i][frame];			
                double rcy = trackY[i][frame];
                if(Double.isNaN(rcx) || Double.isNaN(rcy)) continue;
                
                double zoom_rcx = image_zoom*rcx/lattice;
                double zoom_rcy = image_zoom*rcy/lattice;
                double dx = DX[i][frame];
                double dy = DY[i][frame];
                double dis = Math.sqrt(dx*dx+dy*dy);
                String text = showdxy ? String.format("(%.0f)", dis) : String.format("%d", (i+1));
                TextRoi label = new TextRoi(zoom_rcx, zoom_rcy, text, font);
                label.setStrokeColor(label_color);
                overlay.add(label);                
            }
            return overlay;
    }
    
    public Overlay getLabels(int[] points, int frame, double image_zoom, boolean showdxy){
        if(points==null) return getLabels(frame, image_zoom, showdxy);
        
        int nrois = points.length;
        int font_size = label_size;//(int)Math.ceil(5);
        Font font = new Font("Arial",Font.PLAIN, font_size);
        Overlay overlay = new Overlay();            
        for(int k=0; k<nrois; k++){
            int i= points[k];
            //double rcx = trackX[i][0];//avg(cx[i])/pixel_size;			
            //double rcy = trackY[i][0];//avg(cy[i])/pixel_size;
            double rcx = trackX[i][frame];			
            double rcy = trackY[i][frame];            
            if(Double.isNaN(rcx) || Double.isNaN(rcy)) continue;

            double zoom_rcx = image_zoom*rcx/lattice;
            double zoom_rcy = image_zoom*rcy/lattice;
            double dx = DX[i][frame];
            double dy = DY[i][frame];
            double dis = Math.sqrt(dx*dx+dy*dy);
            String text = showdxy ? String.format("(%.0f)", dis) : String.format("%d", (i+1));
            TextRoi label = new TextRoi(zoom_rcx, zoom_rcy, text, font);
            label.setStrokeColor(label_color);
            overlay.add(label);                
        }
        return overlay;
    }
    
    public Overlay getLabels(int[] points, double[] xc, double[] yc, double[] ddx, double[] ddy, double image_zoom, boolean showdxy){
        int nrois = xc.length;
        int font_size = label_size;//(int)Math.ceil(5);
        Font font = new Font("Arial",Font.PLAIN, font_size);
        Overlay overlay = new Overlay();            
        for(int k=0; k<nrois; k++){
            int i= points[k];
            double rcx = xc[i];			
            double rcy = yc[i];            
            if(Double.isNaN(rcx) || Double.isNaN(rcy)) continue;

            double zoom_rcx = image_zoom*rcx;
            double zoom_rcy = image_zoom*rcy;
            double dx = ddx[i];
            double dy = ddy[i];
            double dis = Math.sqrt(dx*dx+dy*dy);
            String text = showdxy ? String.format("(%.0f)", dis) : String.format("%d", (i+1));
            TextRoi label = new TextRoi(zoom_rcx, zoom_rcy, text, font);
            label.setStrokeColor(label_color);
            overlay.add(label);                
        }
        return overlay;
    }
    
    public Overlay getLabels(int frame){
        return getLabels(frame, 1, false);
    }

     public Overlay getLabels(int frame, boolean showdxy){
        return getLabels(frame, 1, showdxy);
    }
    

}
