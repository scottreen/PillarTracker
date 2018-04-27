
package com.nus.mbi.pillar.drift;

//import static ;
import com.nus.mbi.pillar.detection.CrossCorreclation_Plugin;
import ij.*;
import ij.gui.Roi;
import ij.process.FloatProcessor;
import java.io.*;
import java.nio.*;
import com.nus.mbi.pillar.stat.BasicStatisitic;
import static com.nus.mbi.pillar.stat.BasicStatisitic.*;
import com.nus.mbi.pillar.stat.StatisitcCenterDistance;
/**
 *
 * @author xiaochun
 */
public class DriftDataLoader {

// image and pillar params
public double pixel_size = 1.0; // pixel size in nm
public int npillars = 0;
public int nframes = 0;
public double[][] raw_tracksX = null;
public double[][] raw_tracksY = null;
public double[][] correct_tracksX = null;
public double[][] correct_tracksY = null;
public double[][] DX = null;
public double[][] DY = null;
public double[]   driftX = null;
public double[]   driftY = null;
public boolean correted_data = false;
public boolean load_data_suc = false;
private boolean debug = true;
private static boolean little_endian = false; //Only used for the old binary format

public double[][] raw_tracksX_T = null;
public double[][] raw_tracksY_T = null;
public double[][] correct_tracksX_T = null;
public double[][] correct_tracksY_T = null;
public boolean[] flags_still_pillars = null;

public FileHeaderDrift file_header = null;
public CrossCorreclation_Plugin  enhancer = null;
public int[][] startXY;

public double X_MAX = -1;
public double Y_MAX = -1;

public boolean use_fft_mask = false;
public int mask_radius = 0;
public int start_radius = 0;
public int end_radius = 0;
public int[][] fft_points = null; 
public int center_radius;

private String file_name;
private int version;
public int getFileVersion(){
    return version;
}

public String getFileName() {
        return file_name;
}
    

public boolean load_tracks(String filename, double pixel_size){		
        load_data_suc = false;
        //int [] info=new int[2];
        double[][] matrix_raw = null;
        double[][] matrix_corrected = null;
        double[][] matrix_DXY = null;
        double[][] matrix_driftXY = null;
        boolean[] flags = null;
        file_header = null;
        enhancer = null;
        startXY=null;     
        use_fft_mask = false;
        fft_points = null;
        int num_pillars = 0;
        int num_frames = 0;
        try{
                FileInputStream fis = new FileInputStream(filename);
                version = readInteger(fis,little_endian);                                
                if(version==com.nus.mbi.pillar.tracker.pillar_tracking.fileversion1){                    
                    file_header = readfileinfo_ver1(fis,little_endian);
                    if(file_header!=null){
                        boolean enhancer_valid = true;
                        if(file_header.use_enhancer){
                            enhancer = readPSFImage(fis,little_endian);
                            if(enhancer==null) enhancer_valid=false;
                        }
                        if(enhancer_valid){
                            startXY = readCentroids(fis,file_header.npillars,little_endian);
                            if(startXY!=null){
                                num_pillars = file_header.npillars;
                                num_frames = file_header.nframes;                                                
                            }
                        }
                    }                    
                }
                else if(version==com.nus.mbi.pillar.tracker.pillar_tracking_FD.fileversion2){
                    read_file_header_ver2(fis);
                    if(startXY!=null){
                        num_pillars = file_header.npillars;
                        num_frames = file_header.nframes;                                                
                    }
                }
                else{
                    if(version%2==0){
                        num_pillars = version/2;
                        num_frames = readInteger(fis,little_endian);       
//                        if(num_pillars>0 && num_frames>0){                              
//                            matrix_raw = readfile2matrix(fis, info[0], info[1], little_endian);
//                            matrix_corrected = readfile2matrix(fis, info[0], info[1], little_endian);
//                            matrix_driftXY = readfile2matrix(fis, 2, info[1], little_endian);
//                            flags = readfile2BooleanArray(fis, info[0]/2);
//                        }
                    }
                }
                if(num_pillars>0 && num_frames>0){                    
                    matrix_raw = readfile2matrix(fis, num_pillars*2, num_frames, little_endian);
                    if(version==com.nus.mbi.pillar.tracker.pillar_tracking_FD.fileversion2)
                        matrix_DXY = readfile2matrix(fis, num_pillars*2, num_frames, little_endian);                    
                    matrix_corrected = readfile2matrix(fis, num_pillars*2, num_frames, little_endian);
                    matrix_driftXY = readfile2matrix(fis, 2, num_frames, little_endian);
                    flags = readfile2BooleanArray(fis, num_pillars, little_endian);
                }
                fis.close();		
        }
        catch(Exception e){
                matrix_corrected = null;
                IJ.log("read file failed->" + e.getMessage());
        }

        //if(matrix_raw==null || info[0]%2 != 0) return false;
        if(matrix_raw==null) return false;
        correted_data = (matrix_corrected!=null);
        
        if(file_header!=null && file_header.lattice>1) pixel_size = file_header.lattice;
        //tracks_file = filename;
        this.pixel_size = pixel_size;       
        npillars = num_pillars;
        nframes = num_frames;
        X_MAX = -1;
        Y_MAX = -1;
        raw_tracksX = new double[npillars][nframes];
        raw_tracksY = new double[npillars][nframes];		
        correct_tracksX = new double[npillars][nframes];
        correct_tracksY = new double[npillars][nframes];
        driftX = new double[nframes];
        driftY = new double[nframes];
        flags_still_pillars = new boolean[npillars];
        for(int j=0; j<npillars; j++) flags_still_pillars[j] = false;
        IJ.showStatus("loading the tracks");
        for(int j=0; j<npillars; j++){
            for(int i=0; i<nframes; i++){                
                    double x = raw_tracksX[j][i] = pixel_size*matrix_raw[2*j][i];
                    double y = raw_tracksY[j][i] = pixel_size*matrix_raw[2*j+1][i];
                    if(x>X_MAX) X_MAX = x;
                    if(y>Y_MAX) Y_MAX = y;
                    if(correted_data){
                            correct_tracksX[j][i] = pixel_size*matrix_corrected[2*j][i];
                            correct_tracksY[j][i] = pixel_size*matrix_corrected[2*j+1][i];			
                    }
            }
            IJ.showProgress(j, npillars);
        }
        
        if(correted_data){
            for(int i=0; i<nframes; i++){                
                driftX[i] = pixel_size*matrix_driftXY[0][i];
                driftY[i] = pixel_size*matrix_driftXY[1][i];
            }                
        }
        
        if(matrix_DXY!=null){
            DX = new double[npillars][nframes];
            DY = new double[npillars][nframes];
            IJ.showStatus("loading the deflections");
            for(int j=0; j<npillars; j++){
                for(int i=0; i<nframes; i++){                
                    DX[j][i] = pixel_size*matrix_DXY[2*j][i];
                    DY[j][i] = pixel_size*matrix_DXY[2*j+1][i];                        
                }
                IJ.showProgress(j, npillars);
            }            
        }
        
        if(correted_data){
            for(int j=0; j<npillars; j++) flags_still_pillars[j] = flags[j];
        }
        
        boolean[] discontinuous_flags = getDiscontinuousPillars();
        flags_still_pillars=cross_check_or_flags(flags_still_pillars,discontinuous_flags,npillars);        
        //getLargeMotionPillars();
        
        IJ.log("check the last values");
        IJ.log("	" + raw_tracksX[npillars-1][nframes-1] + "	 " + raw_tracksY[npillars-1][nframes-1]);
        if(correted_data){
                IJ.log("	" + (correct_tracksX[npillars-1][nframes-1] - raw_tracksX[npillars-1][nframes-1]) + "	 " + (correct_tracksY[npillars-1][nframes-1]-raw_tracksY[npillars-1][nframes-1]));
                IJ.log("	" + driftX[nframes-1] + "	 " + driftY[nframes-1]);
        }	
        
        load_data_suc = true;
        file_name = filename;
        return true;
}

public void read_file_header_ver2(FileInputStream fis) throws Exception{
    file_header = readfileinfo_ver2(fis,little_endian);
    if(file_header!=null){
        boolean enhancer_valid = true;
        if(file_header.use_enhancer){
            enhancer = readPSFImage(fis,little_endian);
            if(enhancer==null) enhancer_valid=false;
        }                        
        if(file_header.use_fft){
            use_fft_mask = true;
            mask_radius = readInteger(fis,little_endian); 
            center_radius = readInteger(fis,little_endian); 
            start_radius = readInteger(fis,little_endian); 
            end_radius = readInteger(fis,little_endian); 
            int npoints = readInteger(fis,little_endian); 
            fft_points = readCentroids(fis, npoints, little_endian);
        }

        if(enhancer_valid){
            startXY = readCentroids(fis,file_header.npillars,little_endian);            
        }
    }                     
}

public static FileHeaderDrift read_header_ver2(FileInputStream fis) throws Exception{
    FileHeaderDrift file_header = readfileinfo_ver2(fis,little_endian);
    if(file_header!=null){        
        if(file_header.use_enhancer){
            CrossCorreclation_Plugin enhancer = readPSFImage(fis,little_endian);            
        }                        
        if(file_header.use_fft){
            //use_fft_mask = true;
            int mask_radius = readInteger(fis,little_endian); 
            int center_radius = readInteger(fis,little_endian); 
            int start_radius = readInteger(fis,little_endian); 
            int end_radius = readInteger(fis,little_endian); 
            int npoints = readInteger(fis,little_endian); 
            readCentroids(fis, npoints, little_endian);
        }
        readCentroids(fis,file_header.npillars,little_endian);
    }
    return file_header;
}

public boolean[] getLargeMotionPillars(){    
    if(!correted_data) return null;
    if(file_header==null) return null;
    double s = file_header.spacing*pixel_size;
    double s2 = s*s;
    boolean[] flags = new boolean[npillars];
    for(int j=0; j<npillars; j++){
        flags[j] = false;
        if(!flags_still_pillars[j]){
            boolean found = false;
            for(int f=0; f<nframes-1; f++){
                for(int f1=f+1; f1<nframes; f1++){
                    double dx = correct_tracksX[j][f1] - correct_tracksX[j][f];  
                    double dy = correct_tracksY[j][f1] - correct_tracksY[j][f];  
                    if(dx*dx+dy*dy>s2){
                        found = true;
                        break;
                    }
                }
                if(found) break;
            }
            if(found){
                flags[j] = true;
                flags_still_pillars[j] = true;
            }
        }
    }
    return flags;
}

public boolean[] getDiscontinuousPillars(){        
    boolean[] flags = new boolean[npillars];
    for(int j=0; j<npillars; j++){        
        boolean is_continuous = true;
        for(int f=0; f<nframes; f++){            
            double v = raw_tracksX[j][f];
            if(Double.isNaN(v)||Double.isInfinite(v)){
                is_continuous = false;
                break;
            }
        }
        flags[j] = !is_continuous;            
    }    
    return flags;
}

private boolean[] cross_check_or_flags(boolean[] flag1, boolean[] flag2, int len){
    boolean[] flag = new boolean[len];
    for(int j=0; j<len; j++) flag[j] = (flag1[j]||flag2[j]);        
    return flag;
}

private boolean[] cross_check_and_flags(boolean[] flag1, boolean[] flag2, int len){
    boolean[] flag = new boolean[len];
    for(int j=0; j<len; j++) flag[j] = (flag1[j]&&flag2[j]);        
    return flag;
}

public boolean load_settings(String filename){		
        boolean load_setting_suc = false;        
        file_header = null;
        enhancer = null;
        startXY=null;        
        try{
              
            FileInputStream fis = new FileInputStream(filename);
                version = readInteger(fis,little_endian);                                
                if(version==com.nus.mbi.pillar.tracker.pillar_tracking.fileversion1){                    
                    file_header = readfileinfo_ver1(fis,little_endian);
                    if(file_header!=null){
                        boolean enhancer_valid = true;
                        if(file_header.use_enhancer){
                            enhancer = readPSFImage(fis,little_endian);
                            if(enhancer==null) enhancer_valid=false;
                        }
                        if(enhancer_valid){
                            startXY = readCentroids(fis,file_header.npillars,little_endian);
                            load_setting_suc = true;   
                            //if(startXY!=null) load_setting_suc = true;                            
                        }
                    }
                }
                else if(version==com.nus.mbi.pillar.tracker.pillar_tracking_FD.fileversion2){
                    read_file_header_ver2(fis);
                    load_setting_suc = true;   
                    //if(startXY!=null) load_setting_suc = true;                                            
                }
                fis.close();		
        }
        catch(Exception e){
                IJ.log("read file failed->" + e.getMessage());
        }
        IJ.showStatus("drfit settings loaded");
        IJ.log("try loading settings from drfit file sucessful? " + load_setting_suc);  
        return load_setting_suc;
}
/*
public boolean load_tracks()
{
        return load_tracks(tracks_file, pixel_size);		
}
*/

//public boolean read_settings(FileInputStream fis) throws Exception
//{	
//    DataInputStream din = new DataInputStream(fis);
//    int version = din.readInt();
//    int file_ver = PillarTracker.pillar_tracking.fileversion1;
//    if(version == file_ver){
//        double sigmax_PSF = din.readDouble();
//        double spacing = din.readDouble();
//        double oblique = din.readDouble();
//        double grid_angle = din.readDouble();
//        int kernel_w = din.readInt();
//        int box_constrian_R = din.readInt();
//        double catch_radius = din.readDouble();
//        boolean dark_pillars = din.readBoolean();
//        boolean use_minimum_std = din.readBoolean();        
//    }
//    din.close();
//    return (version == file_ver);
//}

public int[] readfileinfo(FileInputStream fis, boolean use_little_endian) throws Exception
{	
        //double [][] matrix = null;
        //FileInputStream fis = new FileInputStream(tracks_file);
        int[]info = new int[2];
        byte[] header = new byte[Integer.BYTES*2];
        int bytes = fis.read(header);
        if(bytes<=0) return null;
        //  create a byte buffer and wrap the array
        ByteBuffer bb = ByteBuffer.wrap(header);
        //  if the file uses little endian as apposed to network
        //  (big endian, Java's native) format,
        //  then set the byte order of the ByteBuffer
        if(use_little_endian) bb.order(ByteOrder.LITTLE_ENDIAN);

        for (int r = 0; r < info.length; r++){
                info[r] = bb.getInt();
                IJ.log("	" + info[r]);
        }

        return info;
}

public static int readInteger(FileInputStream fis, boolean use_little_endian) throws Exception
{	
        byte[] header = new byte[Integer.BYTES];
        int bytes = fis.read(header);
        if(bytes<=0) return 0;
        ByteBuffer bb = ByteBuffer.wrap(header);
        //  if the file uses little endian as apposed to network
        //  (big endian, Java's native) format,
        //  then set the byte order of the ByteBuffer
        if(use_little_endian) bb.order(ByteOrder.LITTLE_ENDIAN);
        int version = bb.getInt();
        return version;
}

public FileHeaderDrift readfileinfo_ver1(FileInputStream fis, boolean use_little_endian) throws Exception
{	
        int header_size = Integer.BYTES*4 + Double.BYTES*7 + Byte.BYTES*4;
//    header_buffer.putInt(fileversion1);
//    header_buffer.putInt(nrois);
//    header_buffer.putInt(nframes);
//    header_buffer.putDouble(lattice);   
//    header_buffer.putDouble(diameter);  
//    header_buffer.putDouble(spacing);
//    header_buffer.putDouble(oblique);
//    header_buffer.putDouble(grid_angle);
//    header_buffer.putDouble(sigmax_PSF);
//    header_buffer.putDouble(catch_radius);
//    header_buffer.putInt(kernel_w);
//    header_buffer.putInt(box_constrian_R);
//    header_buffer.put(dark_pillars?(byte)1:0);
//    header_buffer.put(use_minimum_std?(byte)1:0);
//    header_buffer.put(use_metric_CG?(byte)1:0);
//    header_buffer.put(image_enhancer_plugin!=null?(byte)1:0);            
        byte[] header = new byte[header_size];
        int bytes = fis.read(header);
        if(bytes<=0) return null;
        //  create a byte buffer and wrap the array
        ByteBuffer bb = ByteBuffer.wrap(header);
        //  if the file uses little endian as apposed to network
        //  (big endian, Java's native) format,
        //  then set the byte order of the ByteBuffer
        if(use_little_endian) bb.order(ByteOrder.LITTLE_ENDIAN);
        
        FileHeaderDrift file_header = new FileHeaderDrift();
        //file_header.file_version = bb.getInt();
        file_header.npillars = bb.getInt();
        file_header.nframes = bb.getInt();
        file_header.lattice = bb.getDouble();
        file_header.diameter = bb.getDouble();
        file_header.spacing = bb.getDouble();
        file_header.oblique = bb.getDouble();
        file_header.grid_angle = bb.getDouble();
        file_header.sigma_PSF = bb.getDouble();
        file_header.catch_radius = bb.getDouble();
        file_header.kernel_w = bb.getInt();
        file_header.box_constrian_R = bb.getInt();
        //file_header.num_start_points= bb.getInt();
        file_header.dark_pillars = bb.get()>0;
        file_header.apply_mean_rank_filter = bb.get()>0;
        file_header.use_metric_CG = bb.get()>0;
        file_header.use_enhancer = bb.get()>0;
 
        return file_header;
}

public static FileHeaderDrift readfileinfo_ver2(FileInputStream fis, boolean use_little_endian) throws Exception
{	
        int header_size = Integer.BYTES*4 + Double.BYTES*7 + Byte.BYTES*5;
        byte[] header = new byte[header_size];
        int bytes = fis.read(header);
        if(bytes<=0) return null;
        ByteBuffer bb = ByteBuffer.wrap(header);
        if(use_little_endian) bb.order(ByteOrder.LITTLE_ENDIAN);
        
        FileHeaderDrift file_header = new FileHeaderDrift();
        file_header.npillars = bb.getInt();
        file_header.nframes = bb.getInt();
        file_header.lattice = bb.getDouble();
        file_header.diameter = bb.getDouble();
        file_header.spacing = bb.getDouble();
        file_header.oblique = bb.getDouble();
        file_header.grid_angle = bb.getDouble();
        file_header.sigma_PSF = bb.getDouble();
        file_header.catch_radius = bb.getDouble();
        file_header.kernel_w = bb.getInt();
        file_header.box_constrian_R = bb.getInt();
        file_header.dark_pillars = bb.get()>0;
        file_header.apply_mean_rank_filter = bb.get()>0;
        file_header.use_metric_CG = bb.get()>0;
        file_header.use_enhancer = bb.get()>0;
        file_header.use_fft = bb.get()>0;
        return file_header;
}


private static CrossCorreclation_Plugin readPSFImage(FileInputStream fis, boolean use_little_endian)  throws IOException {
        int header_size = Integer.BYTES*3 + Double.BYTES*1 + Byte.BYTES*3;
        byte[] header = new byte[header_size];
        int bytes = fis.read(header);
        if(bytes<=0) return null;
        ByteBuffer bb = ByteBuffer.wrap(header);
        if(use_little_endian) bb.order(ByteOrder.LITTLE_ENDIAN);
        CrossCorreclation_Plugin cc = new CrossCorreclation_Plugin();
        cc.setMode(bb.get()>0);
        boolean use_gauss = bb.get()>0;
        boolean is_dark = bb.get()>0;
        int radius = bb.getInt();
        double sigma = bb.getDouble();        
        int w = bb.getInt();
        int h = bb.getInt();
        if(w<1 || h<1) return null;
        
        if(use_gauss) cc.setGaussianPSF(sigma, radius, is_dark);
        int len = w*h;
        byte[] psf_data = new byte[Double.BYTES*len];
        bytes = fis.read(psf_data);
        if(bytes<=0) return null;
        ByteBuffer buff = ByteBuffer.wrap(psf_data);        
        FloatProcessor fp = new FloatProcessor(w,h);
        for(int i=0; i<len; i++) fp.setf(i, (float)buff.getDouble());        
        ImagePlus psf = new ImagePlus("PSF", fp);        
        cc.setCustomPSF(psf);
        return cc;
}

public double[][] readfile2matrix(FileInputStream fis, int ncols, int nrows, boolean use_little_endian) throws Exception
{
        //int ncols = info[0];
        //int nrows = info[1];		
        byte[] data = new byte[ncols*nrows*Double.BYTES];		
        int bytes = fis.read(data);	
        IJ.log("reading bytes="+bytes);
        if(bytes<=0) return null;
        ByteBuffer bb = ByteBuffer.wrap(data);	
        if(use_little_endian) bb.order(ByteOrder.LITTLE_ENDIAN);

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

public static int[][] readCentroids(FileInputStream fis, int num, boolean use_little_endian) throws Exception
{
        byte[] data = new byte[num*2*Integer.BYTES];		
        int bytes = fis.read(data);	
        IJ.log("reading bytes="+bytes);
        if(bytes<=0) return null;
        ByteBuffer bb = ByteBuffer.wrap(data);	
        if(use_little_endian) bb.order(ByteOrder.LITTLE_ENDIAN);

        int[][] matrix = new int[2][num];
        for(int i=0; i<num; i++){
            matrix[0][i] = bb.getInt();
            matrix[1][i] = bb.getInt();                
        }
        return matrix;
}

public boolean[] readfile2BooleanArray(FileInputStream fis, int size, boolean use_little_endian) throws Exception
{
        byte[] data = new byte[size*Byte.BYTES];		
        int bytes = fis.read(data);	
        IJ.log("reading bytes="+bytes);
        if(bytes<=0) return null;
        ByteBuffer bb = ByteBuffer.wrap(data);	
        if(use_little_endian) bb.order(ByteOrder.LITTLE_ENDIAN);

        boolean[] array = new boolean[size];
        for(int i=0; i<size; i++) array[i] = bb.get()>0;        
        return array;
}

//public boolean[] readfile2BooleanArray(FileInputStream fis, int size) throws Exception
//{   
//        DataInputStream dis = new DataInputStream(fis);
//        boolean[] array = new boolean[size];
//        for(int i=0; i<size; i++) array[i] = dis.readBoolean();        
//        return array;
//}

    public double[] getFrameSeries(){
        double[] xaxis = new double[nframes];
        for(int s=0;s<nframes;s++) xaxis[s] = s+1;	
        return xaxis;
    }
    
    public double[][] getSubset(int ipillar, double[][] data){            
        if(ipillar<0 || ipillar>=npillars) return null;
        
        double[][] subset = new double[1][nframes];
        for(int s=0;s<nframes;s++) subset[0][s] = data[ipillar][s];                   
        return subset;
    }
    
    public double[] getXY(int ipillar){            
        if(ipillar<0 || ipillar>=npillars) return null;
        
        double[] xy = new double[2];
        xy[0] = raw_tracksX[ipillar][0]/pixel_size;
        xy[1] = raw_tracksY[ipillar][0]/pixel_size;
        
        return xy;
    }
    
    public double[][] getXYOnFrame(int iframe){            
        if(iframe<0 || iframe>=nframes) return null;
        
        double[][] xy = new double[2][npillars];
        for(int i=0; i<npillars; i++){
            xy[0][i] = raw_tracksX[i][iframe]/pixel_size;
            xy[1][i] = raw_tracksY[i][iframe]/pixel_size;
        }
        return xy;
    }
    
    public double[][] getSubset(int[] indexs, double[][] data)
    {        
        if(indexs==null) return null;
        
        int n = indexs.length;   
        int count = 0;
        for(int i=0;i<n;i++){
            if(indexs[i]>=0 && indexs[i]<npillars){
                count++;
            }
        }
        double[][] subset = new double[count][nframes];
        count = 0;
        for(int i=0;i<n;i++){
            int ipillar = indexs[i];
            if(ipillar>=0 && ipillar<npillars){                
                for(int s=0;s<nframes;s++) subset[count][s] = data[ipillar][s];
                count++;
            }
        }                   
        return subset;
    }
    
    public int[] getSelectedPillars(Roi roi){
        if(roi==null || !roi.isArea()) return null;
        
        int count = 0;
        for(int i=0;i<npillars;i++){
            int x = (int)Math.round(raw_tracksX[i][0]/pixel_size);
            int y = (int)Math.round(raw_tracksY[i][0]/pixel_size);
            if(roi.contains(x, y)){
                count++;
            }
        }

        int[] indexs = new int[count];        
        count = 0;
        for(int i=0;i<npillars;i++){
            int x = (int)Math.round(raw_tracksX[i][0]/pixel_size);
            int y = (int)Math.round(raw_tracksY[i][0]/pixel_size);
            if(roi.contains(x, y)){
                indexs[count] = i;
                count++;
            }
        }

        return indexs;
    }
    
    public int[] getSelectedPillars(Roi roi, int[] index){
        if(roi==null || !roi.isArea()) return null;
        if(index==null || index.length<1) return null;
        int num = index.length;
        
        int count = 0;
        for(int k=0;k<num;k++){
            int i = index[k];
            if(i>=0 && i<npillars){
                int x = (int)Math.round(raw_tracksX[i][0]/pixel_size);
                int y = (int)Math.round(raw_tracksY[i][0]/pixel_size);
                if(roi.contains(x, y)){
                    count++;
                }
            }
        }

        int[] new_indexs = new int[count];        
        count = 0;
        for(int k=0;k<num;k++){
            int i = index[k];
            if(i>=0 && i<npillars){
                int x = (int)Math.round(raw_tracksX[i][0]/pixel_size);
                int y = (int)Math.round(raw_tracksY[i][0]/pixel_size);
                if(roi.contains(x, y)){
                    new_indexs[count] = i;
                    count++;
                }
            }
        }

        return new_indexs;
    }
    
    public double[][] getDeflectionsRAW(int ipillar, int ref_frame)
    {        
        if(ref_frame<0) ref_frame = 0;
        else if(ref_frame>=nframes) ref_frame = nframes - 1;
     
        double[][] deflectionXYR = new double[3][nframes];
        for(int s=0;s<nframes;s++){
            double dx = raw_tracksX[ipillar][s] - raw_tracksX[ipillar][ref_frame];
            double dy = raw_tracksY[ipillar][s] - raw_tracksY[ipillar][ref_frame];			
            deflectionXYR[0][s] = dx;//raw_tracksX[ipillar][0] - raw_tracksX[ipillar][s];
            deflectionXYR[1][s] = dy;//raw_tracksY[ipillar][0] - raw_tracksY[ipillar][s];			
            deflectionXYR[2][s] = Math.sqrt(dx*dx + dy*dy);				
        }		
        
        return deflectionXYR;
    }   
    
    public double[][] getDeflections(int ipillar, int ref_frame)
    {
        if(ref_frame<0) ref_frame = 0;
        else if(ref_frame>=nframes) ref_frame = nframes - 1;
     
        double[][] deflectionXYR = new double[3][nframes];
        
        for(int s=0;s<nframes;s++){
            double dx = correct_tracksX[ipillar][s] - correct_tracksX[ipillar][ref_frame];
            double dy = correct_tracksY[ipillar][s] - correct_tracksY[ipillar][ref_frame];			
            deflectionXYR[0][s] = dx;//raw_tracksX[ipillar][0] - raw_tracksX[ipillar][s];
            deflectionXYR[1][s] = dy;//raw_tracksY[ipillar][0] - raw_tracksY[ipillar][s];			
            deflectionXYR[2][s] = Math.sqrt(dx*dx + dy*dy);				
        }		
        
        return deflectionXYR;
    }
    
    public double[] getDeflections(int ipillar, int frame, int ref_frame)
    {
        if(ref_frame<0) ref_frame = 0;
        else if(ref_frame>=nframes) ref_frame = nframes - 1;
     
        double[] deflectionXYR = new double[3];
        double dx = correct_tracksX[ipillar][frame] - correct_tracksX[ipillar][ref_frame];
        double dy = correct_tracksY[ipillar][frame] - correct_tracksY[ipillar][ref_frame];			
        deflectionXYR[0] = dx;//raw_tracksX[ipillar][0] - raw_tracksX[ipillar][s];
        deflectionXYR[1] = dy;//raw_tracksY[ipillar][0] - raw_tracksY[ipillar][s];			
        deflectionXYR[2] = Math.sqrt(dx*dx + dy*dy);				
        		
        
        return deflectionXYR;
    }
    
    
    public double[][] getDeflectionsRAW(int ipillar){
        return getDeflectionsRAW(ipillar, 0);
    }
    
    public double[][] getDeflections(int ipillar){
        return getDeflections(ipillar, 0);
    }    
    
    public double[][] getAbsDeflections(int ipillar){
        double[][] deflectionXYR = new double[3][nframes];
        
        for(int s=0;s<nframes;s++){
            double dx = DX[ipillar][s];
            double dy = DY[ipillar][s];			
            deflectionXYR[0][s] = dx;//raw_tracksX[ipillar][0] - raw_tracksX[ipillar][s];
            deflectionXYR[1][s] = dy;//raw_tracksY[ipillar][0] - raw_tracksY[ipillar][s];			
            deflectionXYR[2][s] = Math.sqrt(dx*dx + dy*dy);				
        }		
        
        return deflectionXYR;
    }  
    
     public float[][] getAbsDeflectionsFloat(int ipillar){
        float[][] deflectionXYR = new float[3][nframes];
        
        for(int s=0;s<nframes;s++){
            double dx = DX[ipillar][s];
            double dy = DY[ipillar][s];			
            deflectionXYR[0][s] = (float)dx;//raw_tracksX[ipillar][0] - raw_tracksX[ipillar][s];
            deflectionXYR[1][s] = (float)dy;//raw_tracksY[ipillar][0] - raw_tracksY[ipillar][s];			
            deflectionXYR[2][s] = (float)Math.sqrt(dx*dx + dy*dy);				
        }		
        
        return deflectionXYR;
    } 
     
    public StatisitcCenterDistance getDistanceCenterRAW(int ipillar){
        //double[] result = new double[4];
        double avg_raw_x = BasicStatisitic.avg(raw_tracksX[ipillar]);
        double avg_raw_y = BasicStatisitic.avg(raw_tracksY[ipillar]);
        double[] dis = new double[nframes];
        for(int s=0;s<nframes;s++){
            double dx = raw_tracksX[ipillar][s] - avg_raw_x;
            double dy = raw_tracksY[ipillar][s] - avg_raw_y;			            			
            dis[s] = Math.sqrt(dx*dx + dy*dy);				
        }
        
        double avg_dis = BasicStatisitic.avg(dis);
        double max_dis = BasicStatisitic.max(dis);
        double min_dis = BasicStatisitic.min(dis);
        double var_dis = BasicStatisitic.var(dis);
        double std_dis = Math.sqrt(var_dis);
        
        StatisitcCenterDistance result = new StatisitcCenterDistance(avg_raw_x, avg_raw_y, avg_dis, max_dis, min_dis, std_dis);
        
        return result;
    } 
    
    public StatisitcCenterDistance getDistanceCenter(int ipillar){
        //double[] result = new double[4];
        double avg_raw_x = BasicStatisitic.avg(correct_tracksX[ipillar]);
        double avg_raw_y = BasicStatisitic.avg(correct_tracksY[ipillar]);
        double[] dis = new double[nframes];
        for(int s=0;s<nframes;s++){
            double dx = correct_tracksX[ipillar][s] - avg_raw_x;
            double dy = correct_tracksY[ipillar][s] - avg_raw_y;			            			
            dis[s] = Math.sqrt(dx*dx + dy*dy);				
        }
        
        double avg_dis = BasicStatisitic.avg(dis);
        double max_dis = BasicStatisitic.max(dis);
        double min_dis = BasicStatisitic.min(dis);
        double var_dis = BasicStatisitic.var(dis);
        double std_dis = Math.sqrt(var_dis);
        
        StatisitcCenterDistance result = new StatisitcCenterDistance(avg_raw_x, avg_raw_y, avg_dis, max_dis, min_dis, std_dis);
        
        return result;
    } 
    
    public double[] getPillarSeries(){
        double [] plot_x = new double[npillars];        
        for (int i = 0; i < npillars; i++) plot_x[i] = i + 1; 			                    
        return plot_x;
    }
    
    public double[] getStdXY(){        
        double [] std = new double[npillars];
        for (int i = 0; i < npillars; i++){                 			
            std[i] = Math.sqrt(var(correct_tracksX[i])+var(correct_tracksY[i]));
        }
        return std;
    }
    
    public double[] getJumpLength(int ipillar){        
        double [] jump = new double[nframes-1];
        for (int f = 0; f < nframes-1; f++){                 			
            double dx = correct_tracksX[ipillar][f+1] - correct_tracksX[ipillar][f];
            double dy = correct_tracksY[ipillar][f+1] - correct_tracksY[ipillar][f];
            jump[f] = Math.sqrt(dx*dx+dy*dy);
        }
        return jump;
    }
    
    public double[] getPairDistance(int p1, int p2){        
        if(p1<0 || p1>=npillars || p2<0 || p2>=npillars)
            return null;
        
        double [] dis = new double[nframes];        
        for (int i = 0; i < nframes; i++){                 			
            double x1 = raw_tracksX[p1][i];
            double y1 = raw_tracksY[p1][i];
            double x2 = raw_tracksX[p2][i];
            double y2 = raw_tracksY[p2][i];
            
            double dx = x1-x2;
            double dy = y1-y2;
            dis[i] = Math.sqrt(dx*dx + dy*dy);
        }
        return dis;
    }
}
