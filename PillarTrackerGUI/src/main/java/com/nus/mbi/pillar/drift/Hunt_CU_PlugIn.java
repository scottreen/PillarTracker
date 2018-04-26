package com.nus.mbi.pillar.drift;

import com.nus.mbi.pillar.detection.ContractionUnit;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.plugin.PlugIn;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.List;
import com.nus.mbi.pillar.solver.MovingCrossCorrelationNaN;
import static com.nus.mbi.pillar.stat.BasicStatisitic.max;
import com.nus.mbi.pillar.stat.MyPoint;

/**
 *
 * @author xiaochun
 */
public class Hunt_CU_PlugIn implements PlugIn{
    public double catch_radius = 0;
    public double threshold_rs = 1;
    public double threshold_rp = -0.8;
    public double threshold_deflection = 20; //nm;
    public boolean merging_overlaps = true;
    public int frame_window;
    
    public boolean show_table = true;
    public static boolean using_deflections = true;
    public String output_fname = "";
    
    private List[] neigbors = null;
    private List<CU_Projection>[][] proj_list = null;
    private List<ContractionUnit> cu_list = null;
    
    public static final int fileversion1 = 332571;
    public static final int CHECK_NOT = 0;
    public static final int CHECK_CONTRACTING = 1;
    public static final int CHECK_STRECHING = 2;   
    public int check_contracting = CHECK_CONTRACTING;
    
    public Hunt_CU_PlugIn() {
    }

    public Hunt_CU_PlugIn(double catch_radius) {
        this.catch_radius = catch_radius;
    }
    
    public Hunt_CU_PlugIn(double catch_radius, double threshold_rs, double threshold_rp, int frame_window) {
        this.catch_radius = catch_radius;
        this.threshold_rs = threshold_rs;
        this.threshold_rp = threshold_rp;
        this.frame_window = frame_window;
    }
    
    public void copy_parameters(Hunt_CU_PlugIn hunter){
        this.catch_radius = hunter.catch_radius;
        this.threshold_rs = hunter.threshold_rs;
        this.threshold_rp = hunter.threshold_rp;
        this.threshold_deflection = hunter.threshold_deflection;
        this.frame_window = hunter.frame_window;
        this.merging_overlaps = hunter.merging_overlaps;
    }
    
    public boolean showDialog(int nframes){
        GenericDialogPlus gd = new GenericDialogPlus("Hunt Contractile Unit(CU)");
        String[] items = {"None", "Contraction", "Stretching"};
        gd.addChoice("Type:", items, items[check_contracting]);
        
        gd.addNumericField("catch_radius:", catch_radius, 2, 10, "pixel");
        gd.addNumericField("Rp_Threshold (perpendicular):", threshold_rp, 2, 10, "(-1~0)");
        gd.addNumericField("Rs_Threshold (orthogonal):", threshold_rs, 2, 10, "(0~1)");
        gd.addNumericField("Deflection_Threshold:", threshold_deflection, 2, 10, "nm");
        gd.addNumericField("frame_window", frame_window, 0, 10, "(10~"+nframes+")");
        gd.addCheckbox("merging overlaps?", merging_overlaps);        
        gd.addCheckbox("show CU table?", true);        
        gd.addFileField("save to:", output_fname, 50);
        
        gd.showDialog();
        if (gd.wasCanceled()) return false;
        
        check_contracting = gd.getNextChoiceIndex();
        catch_radius = (double)gd.getNextNumber();            
        threshold_rp = (double)gd.getNextNumber();
        threshold_rs = (double)gd.getNextNumber();
        threshold_deflection = (double)gd.getNextNumber();
        frame_window = (int)gd.getNextNumber();
        merging_overlaps = gd.getNextBoolean();
        show_table = gd.getNextBoolean();        
        
        output_fname = gd.getNextString();
        //hunt_single_frame = (frame_window==1);
        if (catch_radius<=0) {IJ.showMessage("the catch radius must be greater than zero");  return false;}
        if (threshold_rp<-1 || threshold_rp>0) {IJ.showMessage("threshold_Rp(perpendicular) must be (-1~0)"); return false;}        
        if (threshold_rs<0 || threshold_rs>1) {IJ.showMessage("threshold_Rs(orthogonal) must be (0~1)"); return false;}
        if (threshold_deflection<0) {IJ.showMessage("deflection threshold_ must be greater than zero)"); return false;}
        return true;
    }
    
    @Override
    public void run(String arg) {
        
    }
    
    public void save(String fname, List<ContractionUnit> cu){
        File output_file = new File(fname);
        //if(!output_file.isFile()|| output_file.isDirectory()) return;
        IJ.log("saving the contraction units->" + fname);
        FileOutputStream fos = null;
        FileChannel fch = null;
        try{            
            fos = new FileOutputStream(output_file);
            fch = fos.getChannel();
            IJ.showStatus("writing raw data into binary file");
            writeFileHeader_ver1(fch);
            write_CU(fch, cu);
            //IJ.showStatus("done with saving raw data");
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
        
        IJ.showStatus("saving done");
    }
    
    public boolean readfileinfo_ver1(FileInputStream fis, boolean use_little_endian) throws Exception
    {	
            int header_size = Integer.BYTES*1 + Double.BYTES*4 + Byte.BYTES*1;

            byte[] header = new byte[header_size];
            int bytes = fis.read(header);
            if(bytes<=0) return false;
            //  create a byte buffer and wrap the array
            ByteBuffer bb = ByteBuffer.wrap(header);        
            if(use_little_endian) bb.order(ByteOrder.LITTLE_ENDIAN);

            catch_radius = bb.getDouble();            
            threshold_rp = bb.getDouble();
            threshold_rs = bb.getDouble();
            threshold_deflection = bb.getDouble();
            frame_window = bb.getInt();
            merging_overlaps = bb.get()>0;

            return true;
    }
    
    public static boolean skip_fileinfo_ver1(FileInputStream fis, boolean use_little_endian) throws Exception
    {	
            int header_size = Integer.BYTES*1 + Double.BYTES*4 + Byte.BYTES*1;

            byte[] header = new byte[header_size];
            int bytes = fis.read(header);
            if(bytes<=0) return false;
            
            return true;
    }
    
    public static List<ContractionUnit> load(String fname){
        List<ContractionUnit> cu = null;
        
        FileInputStream pin;
        //DataInputStream readbin;
        try {
                pin = new FileInputStream(fname);                
                int version=DriftDataLoader.readInteger(pin, false);
                if (version==fileversion1) {
                    IJ.log("loading the contraction units->" + fname);    
                    if(skip_fileinfo_ver1(pin, false)){
                            int n = DriftDataLoader.readInteger(pin, false);
                            if(n>0) cu = readfile_CU(pin, n);                            
                    }
                }
                pin.close();   
        }
        catch(FileNotFoundException fe) {
                //consistent=false;
                IJ.log("file path not found->" + fname);
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
        
        IJ.showStatus("loading done");
        
        return cu;
    }
    
    private void writeFileHeader_ver1(FileChannel out)  throws IOException {        
        int header_size = Integer.SIZE*2 + Double.SIZE*4 + Byte.SIZE*1;
        ByteBuffer header_buffer = ByteBuffer.allocate(header_size);
        header_buffer.putInt(fileversion1);        
        header_buffer.putDouble(catch_radius);
        header_buffer.putDouble(threshold_rp);
        header_buffer.putDouble(threshold_rs);        
        header_buffer.putDouble(threshold_deflection);
        header_buffer.putInt(frame_window);        
        header_buffer.put(merging_overlaps?(byte)1:0);           
        header_buffer.flip();
        out.write(header_buffer);                    
    }
    
    private static void write_CU(FileChannel out, List<ContractionUnit> cu)  throws IOException {        
        int n = 0;
        if(cu!=null) n = cu.size();                
        int header_size = Integer.SIZE;
        ByteBuffer header_buffer = ByteBuffer.allocate(header_size);
        header_buffer.putInt(n);        
        header_buffer.flip();
        out.write(header_buffer);               
        if(n>0){
            ByteBuffer buff = ByteBuffer.allocate(ContractionUnit.SIZE*n); 
            for(int i=0; i<n; i++){
                ContractionUnit  c = cu.get(i);
                buff.putInt(c.i1);
                buff.putInt(c.i2);
                buff.putInt(c.f1);
                buff.putInt(c.f2);
                buff.putDouble(c.rp);
                buff.putDouble(c.rs);
            }
            buff.flip();
            out.write(buff);
        }
    }
    
    public static List<ContractionUnit> readfile_CU(FileInputStream fis, int n) throws Exception
    {
            List<ContractionUnit> cu = new ArrayList();
            int size_cu = ContractionUnit.BYTES;		
            byte[] data = new byte[n*size_cu];		
            int bytes = fis.read(data);	
            //IJ.log("reading bytes="+bytes);
            if(bytes<=0) return null;
            ByteBuffer bb = ByteBuffer.wrap(data);           

            for(int i=0; i<n; i++){
                ContractionUnit c = new ContractionUnit(0,0,0,0,0,0);
                c.i1 = bb.getInt();
                c.i2 = bb.getInt();
                c.f1 = bb.getInt();
                c.f2 = bb.getInt();
                c.rp = bb.getDouble();
                c.rs = bb.getDouble();
                cu.add(c);
            }

            return cu;
    }
    
    public static int[][] create_map(double[] x, double[] y, int width, int height){
        int n = x.length;
        int[][] map = new int[width][height];
        for(int i=0; i<width; i++) for(int j=0; j<height; j++) map[i][j] = -1;        
        for(int k=0; k<n; k++){
            double xc = x[k];
            double yc = y[k];
            if(!Double.isNaN(xc) && !Double.isNaN(yc)){
                int ix = (int)Math.round(xc);
                int iy = (int)Math.round(yc);
                if(iy>=0 && iy<height && ix>=0 && ix<width) map[ix][iy] = k;                
            }
        }
        return map;
    }
    
    public static int[][] create_map(int[] x, int[] y, int width, int height){
        int n = x.length;
        int[][] map = new int[width][height];
        for(int i=0; i<width; i++) for(int j=0; j<height; j++) map[i][j] = -1;        
        for(int k=0; k<n; k++){
            int ix = x[k];
            int iy = y[k];
            if(iy>=0 && iy<height && ix>=0 && ix<width) map[ix][iy] = k;                            
        }
        return map;
    }
    
     public static int search_nearest(double x, double y, double catch_radius, int[][] map){         
        int nearest_id = -1;
        if(!Double.isNaN(x) && !Double.isNaN(y)){
            int w = map.length;
            int h = map[0].length;
            int r=(int)Math.ceil(catch_radius);            
            double r2 = catch_radius*catch_radius;
            double maxr2 = r2*2;
            int xx = (int)Math.round(x);
            int yy = (int)Math.round(y);
            map[xx][yy] = -1;
            for(int i=-r; i<=r; i++){                
                int iy = yy + i;
                for(int j=-r; j<=r; j++){
                    int jx = xx + j;
                    if(iy>=0 && iy<h && jx>=0 && jx<w){
                        int d2 = i*i+j*j;
                        if(d2<=r2){
                            int index = map[jx][iy];
                            if(index>=0 && d2<maxr2){
                                maxr2 = d2;
                                nearest_id=index;
                            }
                        }
                    }
                }
            }
        }
        return nearest_id;
    }
    
    public static int[] search_nearest(double[] x, double[] y, double catch_radius, int map_w, int map_h){
        int n = x.length;
        int[] list = new int[n]; 
        int[][] map = create_map(x, y, map_w, map_h);
        for(int k=0; k<n; k++) list[k] = search_nearest(x[k], y[k], catch_radius, map);        
        return list;
    }    
    
    public static List<Integer> search_neighbors(double x, double y, double catch_radius, int[][] map){
        List<Integer> list = new ArrayList();        
        if(!Double.isNaN(x) && !Double.isNaN(y)){
            int w = map.length;
            int h = map[0].length;
            int r=(int)Math.ceil(catch_radius);            
            double r2 = catch_radius*catch_radius;
            int xx = (int)Math.round(x);
            int yy = (int)Math.round(y);
            map[xx][yy] = -1;
            for(int i=-r; i<=r; i++){                
                int iy = yy + i;
                for(int j=-r; j<=r; j++){
                    int jx = xx + j;
                    if(iy>=0 && iy<h && jx>=0 && jx<w){
                        int ij2 = i*i+j*j;
                        if(ij2>0 && ij2<=r2){
                            int index = map[jx][iy];
                            if(index>=0) list.add(index);
                        }
                    }
                }
            }
        }
        return list;
    }
    
    public static List[] search_neighbors(double[] x, double[] y, double catch_radius){
        int n = x.length;
        List<Integer>[] lists = new List[n];
        int r=(int)Math.ceil(catch_radius);     
        double maxx = max(x);
        double maxy = max(y);
        int map_w = (int)Math.ceil(maxx) + 2*r + 1;
        int map_h = (int)Math.ceil(maxy) + 2*r + 1;
        int[][] map = create_map(x, y, map_w, map_h);
        for(int k=0; k<n; k++) lists[k] = search_neighbors(x[k], y[k], catch_radius, map);        
        return lists;
    }    
    
    public static List[] get_distance(List<Integer>[] neigbors, MyPoint[] points){
        int n = neigbors.length;
        List<Double>[] list = new List[n];        
        for(int k=0; k<n; k++){
            List<Integer> nb = neigbors[k];
            List<Double> dist = new ArrayList();
            MyPoint g1 = points[k];            
            for(int i=0; i<nb.size(); i++){
                int j = nb.get(i);
                MyPoint g2 = points[j];
                double dx = g1.x-g2.x;
                double dy = g1.y-g2.y;
                double d = Math.sqrt(dx*dx+dy*dy);
                dist.add(d);
            }
            list[k] = dist;
        }        
        return list;
    }
    
//    public static List[] project(List<Integer>[] neigbors, List<Double>[] distance, MyPoint[] anchors, MyPoint[] centers){        
//        int n = neigbors.length;
//        List<CU_Projection>[] list = new List[n];        
//        for(int k=0; k<n; k++){
//            List<Integer> list_nb = neigbors[k];
//            List<Double> list_dis = distance[k];
//            List<CU_Projection> list_proj = new ArrayList();
//            MyPoint g1 = anchors[k];
//            MyPoint c1 = centers[k];            
//            for(int i=0; i<list_nb.size(); i++){
//                int j = list_nb.get(i);
//                double dis = list_dis.get(i);
//                MyPoint g2 = anchors[j];
//                MyPoint c2 = centers[j];
//                CU_Projection proj = new CU_Projection(0,0,0,0,0,0);
//                if(!Double.isNaN(c1.x) && !Double.isNaN(c2.x))
//                    proj = CU_Projection.create(g1, g2, dis, c1, c2);
//                list_proj.add(proj);                            
//            }
//            list[k] = list_proj;
//        }
//        
//        return list;
//    }
    
    public static List[] project(List<Integer>[] neigbors, List<Double>[] distance, MyPoint[] anchors, MyPoint[] centers){        
        int n = neigbors.length;
        List<CU_Projection>[] list = new List[n];        
        for(int k=0; k<n; k++){
            List<Integer> list_nb = neigbors[k];
            List<Double> list_dis = distance[k];
            List<CU_Projection> list_proj = new ArrayList();
            MyPoint g1 = anchors[k];
            MyPoint c1 = centers[k];            
            for(int i=0; i<list_nb.size(); i++){
                int j = list_nb.get(i);
                double dis = list_dis.get(i);
                MyPoint g2 = anchors[j];
                MyPoint c2 = centers[j];
                CU_Projection proj = new CU_Projection();//new CU_Projection(0,0,0,0,0,0);                        
                if(!Double.isNaN(c1.x) && !Double.isNaN(c2.x)){
                    if(using_deflections){
                        MyPoint gg = new MyPoint(g2.x-g1.x, g2.y-g1.y);
                        proj = CU_Projection.create(gg, dis, c1, c2);
                    }
                    else proj = CU_Projection.create(g1, g2, dis, c1, c2);
                }
                list_proj.add(proj);
            }
            list[k] = list_proj;
        }
        
        return list;
    }
    
//    public static double cross_correlation(double[] a, double[]b){
//        int n = a.length;
//        double sumI = 0;
//        double sumK = 0;
//        double sumI2 = 0;
//        double sumK2 = 0;
//        double sumIK = 0;
//        int counter = n;
//
//        for(int i=0; i<n; i++){
//            double mi = a[i];
//            double ki = b[i];
//            sumI += mi;
//            sumK += ki;
//            sumIK += mi*ki;
//            sumI2 += mi*mi;
//            sumK2 += ki*ki;
//        }
//
//        double std_img = sumI2 - sumI*sumI/counter;
//        double std_psf = sumK2 - sumK*sumK/counter;
//        double sum_img_psf = sumIK-sumI*sumK/counter;						
//        double cc = sum_img_psf*(counter-1)/(counter*Math.sqrt(std_img*std_psf));
//        return cc;
//    }
    
    public static MyPoint cross_correlation(CU_Projection[] projections, int start_f, int end_f){
        double rs = Double.NaN;
        double rp = Double.NaN;
        int frames = projections.length;
        if(start_f>=0 && start_f<frames && end_f>start_f && end_f<frames){
            int num = end_f-start_f+1;
            double[] a = new double[num];
            double[] b = new double[num];
            boolean[] nan_flag = new boolean[num];
            int n = 0;
            for(int f=start_f; f<=end_f; f++){
                a[n] = projections[f].proj_p1;
                b[n] = projections[f].proj_p2;
                nan_flag[n] = projections[f].isNaN;
                n++;
            }            
            rp = MovingCrossCorrelationNaN.cross_correlation(a,b,nan_flag);
            
            n = 0;
            for(int f=start_f; f<=end_f; f++){
                a[n] = projections[f].proj_s1;
                b[n] = projections[f].proj_s2;
                n++;
            }
            rs = MovingCrossCorrelationNaN.cross_correlation(a,b,nan_flag);
        }
        return new MyPoint(rs,rp);
    }
    
    public static MovingCrossCorrelationNaN[] cross_correlation(CU_Projection[] projections, int window){
        int frames = projections.length;
        double[] a = new double[frames];
        double[] b = new double[frames];   
        boolean[] nan_flag = new boolean[frames];
        for(int f=0; f<frames; f++){
            a[f] = projections[f].proj_p1;
            b[f] = projections[f].proj_p2;     
            nan_flag[f] = projections[f].isNaN;
        }
        MovingCrossCorrelationNaN rp = new MovingCrossCorrelationNaN(nan_flag, a, b, window);
        rp.moving_cross_correlation();
        
        for(int f=0; f<frames; f++){
            a[f] = projections[f].proj_s1;
            b[f] = projections[f].proj_s2;            
        }        
        MovingCrossCorrelationNaN rs = new MovingCrossCorrelationNaN(nan_flag, a, b, window);        
        rs.moving_cross_correlation();
        
        return new MovingCrossCorrelationNaN[]{rp,rs};
    }
    
//    public static MyPoint get_maximum(CU_Projection[] projections, int start_f, int end_f){
//        double max_1 = -1;
//        double max_2 = -1;
//        int frames = projections.length;
//        if(start_f>=0 && start_f<frames && end_f>start_f && end_f<frames){                        
//            for(int f=start_f; f<=end_f; f++){
//                double p1 = Math.abs(projections[f].proj_p1);
//                double p2 = Math.abs(projections[f].proj_p2);
//                if(p1>max_1) max_1 = p1;
//                if(p2>max_2) max_2 = p2;
//            }              
//        }
//        return new MyPoint(max_1,max_2);
//    }
    
    
    public static MyPoint get_length_perpendicular(CU_Projection[] projections, int start_f, int end_f){
        double length1 = 0;
        double length2 = 0;
        int frames = projections.length;
        if(start_f>=0 && start_f<frames && end_f>start_f && end_f<frames){                                    
            double p1 = projections[start_f].proj_p1;
            double p2 = projections[start_f].proj_p2;
            
            double min_1 = p1;
            double min_2 = p2;
            double max_1 = p1;
            double max_2 = p2;
            for(int f=start_f+1; f<=end_f; f++){
                p1 = projections[f].proj_p1;
                p2 = projections[f].proj_p2;
                if(p1>max_1) max_1 = p1;
                if(p2>max_2) max_2 = p2;
                if(p1<min_1) min_1 = p1;
                if(p2<min_2) min_2 = p2;
            }
            length1 = max_1-min_1;
            length2 = max_2-min_2;
        }
        return new MyPoint(length1,length2);
    }
    
    public static MyPoint get_length_orthogonal(CU_Projection[] projections, int start_f, int end_f){
        double length1 = 0;
        double length2 = 0;
        int frames = projections.length;
        if(start_f>=0 && start_f<frames && end_f>start_f && end_f<frames){                                    
            double p1 = projections[start_f].proj_s1;
            double p2 = projections[start_f].proj_s2;
            
            double min_1 = p1;
            double min_2 = p2;
            double max_1 = p1;
            double max_2 = p2;
            for(int f=start_f+1; f<=end_f; f++){
                p1 = projections[f].proj_s1;
                p2 = projections[f].proj_s2;
                if(p1>max_1) max_1 = p1;
                if(p2>max_2) max_2 = p2;
                if(p1<min_1) min_1 = p1;
                if(p2<min_2) min_2 = p2;
            }
            length1 = max_1-min_1;
            length2 = max_2-min_2;
        }
        return new MyPoint(length1,length2);
    }
    
    public static double percentage_contracting(CU_Projection[] projections, int start_f, int end_f){
        double percentage = 0;
        int frames = projections.length;
        if(start_f>=0 && start_f<frames && end_f>start_f && end_f<frames){
            int num = end_f-start_f+1;            
            int n = 0;
            for(int f=start_f; f<=end_f; f++){                
                double d = projections[f].proj_p1-projections[f].proj_p2;
                if(d>0) n++;
            }
            percentage = n*100/num;
        }
        return percentage;
    }
    
//    public static double line_least_square(CU_Projection[] projections, int start_f, int end_f){
//        double slop = 0;
//        int frames = projections.length;
//        if(start_f>=0 && start_f<frames && end_f>start_f && end_f<frames){
//            int num = end_f-start_f+1;
//            double[] a = new double[num];
//            double[] b = new double[num];
//            int n = 0;
//            for(int f=start_f; f<=end_f; f++){
//                a[n] = n;
//                b[n] = projections[f].proj_p1-projections[f].proj_p2;
//                n++;
//            }  
//            
//            if(n>2){
//                CurveFitter cf = new CurveFitter(a, b);
//                cf.doFit(CurveFitter.STRAIGHT_LINE);
//                double[] paras = cf.getParams();
//                slop = paras[1];
//            }
//        }
//        return slop;
//    }
    
//    public static List<ContractionUnit> hunt(double[]gx, double[]gy, double catch_radius, double[][] cx, double[][] cy, double pixel_size, double threshold_rs, double threshold_rp, double threshold_deflection){        
//        List[] neigbors = search_neighbors(gx, gy, catch_radius);
//        int n = neigbors.length;
//        MyPoint[] anchors = new MyPoint[n];
//        for(int k=0; k<n; k++) anchors[k] = new MyPoint(gx[k], gy[k]);
//        List[] distance = get_distance(neigbors, anchors);                        
//        int frames = cx[0].length;
//        IJ.showStatus("projecting");
//        List<CU_Projection>[][] proj_list = new List[frames][];     
//        for(int f=0; f<frames; f++){
//            MyPoint[] centers = new MyPoint[n];
//            for(int k=0; k<n; k++) centers[k] = new MyPoint(cx[k][f]/pixel_size, cy[k][f]/pixel_size);                
//            proj_list[f] = project(neigbors, distance, anchors, centers);
//            IJ.showProgress(f, frames);
//        }
//        
//        IJ.showStatus("hunting contractile unit");
//        List<ContractionUnit> cu_list = new ArrayList();    
//        for(int k=0; k<n; k++){
//            List<Integer> list_nb = neigbors[k];            
//            for(int i=0; i<list_nb.size(); i++){
//                int j = list_nb.get(i);
//                CU_Projection[] proj_f = new CU_Projection[frames];
//                for(int f=0; f<frames; f++) proj_f[f] = proj_list[f][k].get(i);                
//                MyPoint R = cross_correlation(proj_f, 0, frames-1);
//                double rs = R.x;
//                double rp = R.y;
//                if(Math.abs(rs)<threshold_rs && rp<threshold_rp){
//                    //if(percentage_contracting(proj_f, 0, frames-1)>80)
//                    {                    
//                        double deflection_threshold = threshold_deflection/pixel_size; //20/pixel_size;
//                        MyPoint length = get_length(proj_f, 0, frames-1);
//                        if(length.x>deflection_threshold && length.y>deflection_threshold){
//                            ContractionUnit cu = new ContractionUnit(k,j,0,frames-1,rs,rp);
//                            cu_list.add(cu);
//                        }
//                    }
//                }
//            }
//            IJ.showProgress(k, n);
//        }
//        IJ.showStatus("hunting done");
//        return cu_list;
//    }
    
    public void setToReProjecting(){
        proj_list=null;
    }
    
    public List<ContractionUnit> hunt(double[]gx, double[]gy, double[][] dx, double[][] dy, double pixel_size){    
        if(proj_list==null) projecting(gx, gy, catch_radius, dx, dy, pixel_size);
        //projecting(gx, gy, catch_radius, dx, dy, pixel_size);
        int npillars = neigbors.length;
        int nframes = dx[0].length;
        return hunt_cu(npillars, nframes, pixel_size);
        //return hunt(gx, gy, catch_radius, dx, dy, pixel_size,threshold_rs,threshold_rp, threshold_deflection,frame_window, merging_overlaps);
    }
    
    public List<CU_Projection>[][] projecting(double[]gx, double[]gy, double catch_radius, double[][] dx, double[][] dy, double pixel_size){        
        neigbors = search_neighbors(gx, gy, catch_radius);
        int n = neigbors.length;
        MyPoint[] anchors = new MyPoint[n];
        for(int k=0; k<n; k++) anchors[k] = new MyPoint(gx[k], gy[k]);
        List[] distance = get_distance(neigbors, anchors);                        
        int frames = dx[0].length;
        IJ.showStatus("projecting");
        proj_list = new List[frames][];     
        for(int f=0; f<frames; f++){
            MyPoint[] centers = new MyPoint[n];
            for(int k=0; k<n; k++) centers[k] = new MyPoint(dx[k][f]/pixel_size, dy[k][f]/pixel_size);                
            proj_list[f] = project(neigbors, distance, anchors, centers);
            IJ.showProgress(f, frames);
        }
        return proj_list;
    }
    
    public List<ContractionUnit> hunt_cu(int npillars, int frames, double pixel_size){ 
        int n = npillars;//neigbors.length;
        IJ.showStatus("hunting contractile unit");
        cu_list = new ArrayList();          
        if(frame_window<10 || frame_window>frames) frame_window = frames;  
        int num_Rsp = frames-frame_window+1;
        double deflection_threshold = threshold_deflection/pixel_size;//20/pixel_size;
        double threshold_d2 = deflection_threshold*deflection_threshold;                  
        //int N1 = frame_window-1;
        for(int k=0; k<n; k++){
            List<Integer> list_nb = neigbors[k];            
            for(int i=0; i<list_nb.size(); i++){
                int j = list_nb.get(i);
                CU_Projection[] proj_f = new CU_Projection[frames];
                for(int f=0; f<frames; f++) proj_f[f] = proj_list[f][k].get(i);
                MovingCrossCorrelationNaN[] window_R = cross_correlation(proj_f, frame_window);
                MovingCrossCorrelationNaN window_Rp = window_R[0];
                MovingCrossCorrelationNaN window_Rs = window_R[1];
                
                //search minimum rp
                List<Integer> minimas = new ArrayList();
                for(int f=0; f<num_Rsp; f++){
                    double rs = window_Rs.R[f];
                    double rp = window_Rp.R[f];
                    boolean is_accept = false;
                    if(Math.abs(rs)<threshold_rs && rp<threshold_rp){                        
                        int N = window_Rp.countN[f];
                        int N1 = N - 1;
                        double std_p1 = window_Rp.varA[f]/N1;
                        double std_p2 = window_Rp.varB[f]/N1;
                        double std_s1 = window_Rs.varA[f]/N1;
                        double std_s2 = window_Rs.varB[f]/N1;                    
                        if(Math.min(std_p1, std_p2)>Math.max(std_s1, std_s2))
                        //if(std_p1+std_p2>std_s1+std_s2)
                        {
                            if(check_contracting==CHECK_CONTRACTING || check_contracting==CHECK_STRECHING){
                                double avg_p1 = window_Rp.sumA[f]/N;
                                double avg_p2 = window_Rp.sumB[f]/N;
                                if(check_contracting==CHECK_CONTRACTING){
                                    if(avg_p1>deflection_threshold && avg_p2<-deflection_threshold){
                                        is_accept = true;
                                    }
                                }
                                else{
                                    if(avg_p1<-deflection_threshold && avg_p2>deflection_threshold){
                                        is_accept = true;
                                    }
                                }
                            }
                            else{
                                if(std_p1>threshold_d2 && std_p2>threshold_d2){
                                    is_accept = true;
                                }
                            }
                        }    
                    }
                    
                    if(is_accept) minimas.add(f);                    
                }
                
                int num_minimas = minimas.size();
                for(int m=0; m<num_minimas; m++)
                {                        
                    int f_start = minimas.get(m);          
                    int f_end = f_start+frame_window-1;                    
                    MyPoint R = new MyPoint(window_Rs.R[f_start], window_Rp.R[f_start]);//Rsp[f_start];
                    
                    if(merging_overlaps){                        
//                        CrossCorrelationNaN ccRs = window_Rs.get(f_start);
//                        CrossCorrelationNaN ccRp = window_Rp.get(f_start);   
//                        int ff_start = f_start;
//                        while((m+1)<num_minimas && minimas.get(m+1)<f_end){                        
//                            int f = minimas.get(m+1);                            
//                            int shift = f - ff_start;
//                            ccRs = ccRs.merge_overlaps(window_Rs.get(f), shift);
//                            ccRp = ccRp.merge_overlaps(window_Rp.get(f), shift);                            
//                            MyPoint RR = new MyPoint(ccRs.R,ccRp.R);
//                            double rs = RR.x;
//                            double rp = RR.y;
//                            if(Math.abs(rs)<threshold_rs && rp<threshold_rp){
//                                R = RR;   
//                                ff_start = f;
//                                f_end = f+frame_window-1;
//                                m++;                            
//                            }
//                            else break;
//                        }
                        while((m+1)<num_minimas && minimas.get(m+1)<f_end){                        
                            int f = minimas.get(m+1);
                            int ff_end = f+frame_window-1;

                            MyPoint RR = cross_correlation(proj_f, f_start, ff_end);
                            double rs = RR.x;
                            double rp = RR.y;
                            if(Math.abs(rs)<threshold_rs && rp<threshold_rp){
                                R = RR;                            
                                f_end = ff_end;
                                m++;                            
                            }
                            else break;
                        }
                    }
                    
                    double rs = R.x;
                    double rp = R.y;
                    ContractionUnit cu = new ContractionUnit(k,j,f_start,f_end,rs,rp);
                    cu_list.add(cu);
                }
            }
            IJ.showProgress(k, n);
        }
        IJ.showStatus("hunting done");
        return cu_list;
    }
    
    public List<ContractionUnit> hunt_cu_maxima(int npillars, int frames, double pixel_size){ 
        int n = npillars;//neigbors.length;
        IJ.showStatus("hunting contractile unit");
        cu_list = new ArrayList();          
        if(frame_window<10 || frame_window>frames) frame_window = frames;  
        int num_Rsp = frames-frame_window+1;
        double deflection_threshold = threshold_deflection/pixel_size;//20/pixel_size;
        double threshold_d2 = deflection_threshold*deflection_threshold;                  
        //int N1 = frame_window-1;
        for(int k=0; k<n; k++){
            List<Integer> list_nb = neigbors[k];            
            for(int i=0; i<list_nb.size(); i++){
                int j = list_nb.get(i);
                CU_Projection[] proj_f = new CU_Projection[frames];
                for(int f=0; f<frames; f++) proj_f[f] = proj_list[f][k].get(i);
                MovingCrossCorrelationNaN[] window_R = cross_correlation(proj_f, frame_window);
                MovingCrossCorrelationNaN window_Rp = window_R[0];
                MovingCrossCorrelationNaN window_Rs = window_R[1];
                
                MyPoint[] Rsp = new MyPoint[num_Rsp];
                for(int f=0; f<num_Rsp; f++){
                    //MyPoint R = new MyPoint(0,0);
                    boolean is_accept = false;
                    int N = window_Rp.countN[f];
                    
                    int N1 = N - 1;
                    double std_p1 = window_Rp.varA[f]/N1;
                    double std_p2 = window_Rp.varB[f]/N1;
                    double std_s1 = window_Rs.varA[f]/N1;
                    double std_s2 = window_Rs.varB[f]/N1;                    
                    if(Math.max(std_p1, std_p2)>Math.max(std_s1, std_s2))
                    {
                        if(check_contracting==CHECK_CONTRACTING || check_contracting==CHECK_STRECHING){
                            double avg_p1 = window_Rp.sumA[f]/N;
                            double avg_p2 = window_Rp.sumB[f]/N;
                            if(check_contracting==CHECK_CONTRACTING){
                                if(avg_p1>deflection_threshold && avg_p2<-deflection_threshold){
                                    is_accept = true;
                                }
                            }
                            else{
                                if(avg_p1<-deflection_threshold && avg_p2>deflection_threshold){
                                    is_accept = true;
                                }
                            }
                        }
                        else{
                            if(std_p1>threshold_d2 && std_p2>threshold_d2){
                                is_accept = true;
                            }
                        }
                    }
                    Rsp[f] = is_accept ? new MyPoint(window_Rs.R[f], window_Rp.R[f]) : new MyPoint(0,0); 
                }
                
                //search minimum rp
                List<Integer> minimas = new ArrayList();
                for(int f=0; f<num_Rsp; f++){
//                    double rs = window_Rs.R[f];
//                    double rp = window_Rp.R[f];
                    MyPoint R = Rsp[f];
                    double rs = R.x;
                    double rp = R.y;
                    boolean is_minimum = false;
                    if(Math.abs(rs)<threshold_rs && rp<threshold_rp){
                        is_minimum = true;
                        if(f-1>=0 && Rsp[f-1].y<=rp) is_minimum = false;
                        else if(f+1<num_Rsp && Rsp[f+1].y<=rp) is_minimum = false;
//                        if(f-1>=0 && window_Rp.R[f-1]<rp) is_minimum = false;
//                        else if(f+1<num_Rsp && window_Rp.R[f+1]<rp) is_minimum = false;
//                        if(is_minimum){
//                            is_minimum = false;
//                        }
                    }
                    
                    if(is_minimum) minimas.add(f);                    
                }
                
                int num_minimas = minimas.size();
                for(int m=0; m<num_minimas; m++)
                {                        
                    int f_start = minimas.get(m);          
                    int f_end = f_start+frame_window-1;                    
                    MyPoint R = Rsp[f_start];//new MyPoint();
                    
                    if(merging_overlaps){
                        while((m+1)<num_minimas && minimas.get(m+1)<f_end){                        
                            int f = minimas.get(m+1);
                            int ff_end = f+frame_window-1;

                            MyPoint RR = cross_correlation(proj_f, f_start, ff_end);
                            double rs = RR.x;
                            double rp = RR.y;
                            if(Math.abs(rs)<threshold_rs && rp<threshold_rp){
                                R = RR;                            
                                f_end = ff_end;
                                m++;                            
                            }
                            else break;
                        }
                    }
                    
                    double rs = R.x;
                    double rp = R.y;
                    ContractionUnit cu = new ContractionUnit(k,j,f_start,f_end,rs,rp);
                    cu_list.add(cu);
                }
            }
            IJ.showProgress(k, n);
        }
        IJ.showStatus("hunting done");
        return cu_list;
    }
    
    public List<ContractionUnit> hunt(double[]gx, double[]gy, double catch_radius, double[][] dx, double[][] dy, double pixel_size, 
            double threshold_rs, double threshold_rp, double threshold_deflection,int frame_window, boolean merging){        
        List[] neigbors = search_neighbors(gx, gy, catch_radius);
        int n = neigbors.length;
        MyPoint[] anchors = new MyPoint[n];
        for(int k=0; k<n; k++) anchors[k] = new MyPoint(gx[k], gy[k]);
        List[] distance = get_distance(neigbors, anchors);                        
        int frames = dx[0].length;
        IJ.showStatus("projecting");
        List<CU_Projection>[][] proj_list = new List[frames][];     
        for(int f=0; f<frames; f++){
            MyPoint[] centers = new MyPoint[n];
            for(int k=0; k<n; k++) centers[k] = new MyPoint(dx[k][f]/pixel_size, dy[k][f]/pixel_size);                
            proj_list[f] = project(neigbors, distance, anchors, centers);
            IJ.showProgress(f, frames);
        }
        
        IJ.showStatus("hunting contractile unit");
        List<ContractionUnit> cu_list = new ArrayList();  
        //frame_window = 500;
        if(frame_window<20 || frame_window>frames) frame_window = frames;  
        //int radius_f = frame_window/2;
        int num_Rsp = frames-frame_window+1;
        double deflection_threshold = threshold_deflection/pixel_size;//20/pixel_size;
        double threshold_d2 = deflection_threshold*deflection_threshold;                  
        //int N1 = frame_window-1;
        for(int k=0; k<n; k++){
            List<Integer> list_nb = neigbors[k];            
            for(int i=0; i<list_nb.size(); i++){
                int j = list_nb.get(i);
                CU_Projection[] proj_f = new CU_Projection[frames];
                for(int f=0; f<frames; f++) proj_f[f] = proj_list[f][k].get(i);
                MovingCrossCorrelationNaN[] window_R = cross_correlation(proj_f, frame_window);
                MovingCrossCorrelationNaN window_Rp = window_R[0];
                MovingCrossCorrelationNaN window_Rs = window_R[1];
                
                MyPoint[] Rsp = new MyPoint[num_Rsp];
                for(int f=0; f<num_Rsp; f++){
                    int N = window_Rp.countN[f];
                    int N1 = N - 1;
                    double std_p1 = window_Rp.varA[f]/N1;
                    double std_p2 = window_Rp.varB[f]/N1;
                    double std_s1 = window_Rs.varA[f]/N1;
                    double std_s2 = window_Rs.varB[f]/N1;
                    boolean is_accept = false;
                    if(Math.min(std_p1, std_p2)>Math.max(std_s1, std_s2)){
                        if(check_contracting==CHECK_CONTRACTING || check_contracting==CHECK_STRECHING){
                            double avg_p1 = window_Rp.sumA[f]/N;
                            double avg_p2 = window_Rp.sumB[f]/N;
                            if(check_contracting==CHECK_CONTRACTING){
                                if(avg_p1>deflection_threshold && avg_p2<-deflection_threshold){
                                    is_accept = true;
                                }
                            }
                            else{
                                if(avg_p1<-deflection_threshold && avg_p2>deflection_threshold){
                                    is_accept = true;
                                }
                            }
                        }
                        else{
                            if(std_p1>threshold_d2 && std_p2>threshold_d2){
                                is_accept = true;
                            }
                        }
                    }
                    Rsp[f] = is_accept ? new MyPoint(window_Rs.R[f], window_Rp.R[f]) : new MyPoint(0,0); 
                }
                
                //search minimum rp
                List<Integer> minimas = new ArrayList();
                for(int f=0; f<num_Rsp; f++){
//                    double rs = window_Rs.R[f];
//                    double rp = window_Rp.R[f];
                    MyPoint R = Rsp[f];
                    double rs = R.x;
                    double rp = R.y;
                    boolean is_minimum = false;
                    if(Math.abs(rs)<threshold_rs && rp<threshold_rp){
                        is_minimum = true;
                        if(f-1>=0 && Rsp[f-1].y<=rp) is_minimum = false;
                        else if(f+1<num_Rsp && Rsp[f+1].y<=rp) is_minimum = false;
//                        if(f-1>=0 && window_Rp.R[f-1]<rp) is_minimum = false;
//                        else if(f+1<num_Rsp && window_Rp.R[f+1]<rp) is_minimum = false;
//                        if(is_minimum){
//                            is_minimum = false;
//                        }
                    }
                    
                    if(is_minimum) minimas.add(f);                    
                }
                
                int num_minimas = minimas.size();
                for(int m=0; m<num_minimas; m++)
                {                        
                    int f_start = minimas.get(m);          
                    int f_end = f_start+frame_window-1;                    
                    MyPoint R = Rsp[f_start];//new MyPoint();
                    
                    if(merging){
                        while((m+1)<num_minimas && minimas.get(m+1)<f_end){                        
                            int f = minimas.get(m+1);
                            int ff_end = f+frame_window-1;

                            MyPoint RR = cross_correlation(proj_f, f_start, ff_end);
                            double rs = RR.x;
                            double rp = RR.y;
                            if(Math.abs(rs)<threshold_rs && rp<threshold_rp){
                                R = RR;                            
                                f_end = ff_end;
                                m++;                            
                            }
                            else break;
                        }
                    }
                    
                    double rs = R.x;
                    double rp = R.y;
                    ContractionUnit cu = new ContractionUnit(k,j,f_start,f_end,rs,rp);
                    cu_list.add(cu);
                }
            }
            IJ.showProgress(k, n);
        }
        IJ.showStatus("hunting done");
        return cu_list;
    }
    
    public static List<ContractionUnit> hunt_slow(double[]gx, double[]gy, double catch_radius, double[][] dx, double[][] dy, double pixel_size, 
            double threshold_rs, double threshold_rp, double threshold_deflection,int frame_window, boolean merging){        
        List[] neigbors = search_neighbors(gx, gy, catch_radius);
        int n = neigbors.length;
        MyPoint[] anchors = new MyPoint[n];
        for(int k=0; k<n; k++) anchors[k] = new MyPoint(gx[k], gy[k]);
        List[] distance = get_distance(neigbors, anchors);                        
        int frames = dx[0].length;
        IJ.showStatus("projecting");
        List<CU_Projection>[][] proj_list = new List[frames][];     
        for(int f=0; f<frames; f++){
            MyPoint[] centers = new MyPoint[n];
            for(int k=0; k<n; k++) centers[k] = new MyPoint(dx[k][f]/pixel_size, dy[k][f]/pixel_size);                
            proj_list[f] = project(neigbors, distance, anchors, centers);
            IJ.showProgress(f, frames);
        }
        
        IJ.showStatus("hunting contractile unit");
        List<ContractionUnit> cu_list = new ArrayList();  
        //frame_window = 500;
        if(frame_window<20 || frame_window>frames) frame_window = frames;  
        //int radius_f = frame_window/2;
        int num_Rsp = frames-frame_window+1;
        double deflection_threshold = threshold_deflection/pixel_size;//20/pixel_size;
                                
        for(int k=0; k<n; k++){
            List<Integer> list_nb = neigbors[k];            
            for(int i=0; i<list_nb.size(); i++){
                int j = list_nb.get(i);
                CU_Projection[] proj_f = new CU_Projection[frames];
                for(int f=0; f<frames; f++) proj_f[f] = proj_list[f][k].get(i);                
                MyPoint[] Rsp = new MyPoint[num_Rsp];
                for(int f=0; f<num_Rsp; f++){
                    int f_end = f+frame_window-1;
                    MyPoint length_perpd = get_length_perpendicular(proj_f, f, f_end);
                    MyPoint R = new MyPoint(0,0);
                    if(length_perpd.x>deflection_threshold && length_perpd.y>deflection_threshold){
                        MyPoint length_ortho = get_length_orthogonal(proj_f, f, f_end);
                        if(Math.max(length_perpd.x, length_perpd.y)>Math.min(length_ortho.x, length_ortho.y))
                            R = cross_correlation(proj_f, f, f_end);
                    }
                    Rsp[f] = R; 
                }
                //search minimum rp
                List<Integer> minimas = new ArrayList();
                for(int f=0; f<num_Rsp; f++){
                    int f_end = f+frame_window-1;
                    MyPoint R = Rsp[f];
                    double rs = R.x;
                    double rp = R.y;
                    boolean is_minimum = false;
                    if(Math.abs(rs)<threshold_rs && rp<threshold_rp){
                        is_minimum = true;
                        if(f-1>=0 && Rsp[f-1].y<=rp) is_minimum = false;
                        else if(f+1<num_Rsp && Rsp[f+1].y<=rp) is_minimum = false;
//                        if(is_minimum){
//                            //if(percentage_contracting(proj_f, f, f_end)<80) is_minimum = false;
//                            //else
////                            {
////                                MyPoint length = get_length(proj_f, f, f_end);
////                                if(length.x<deflection_threshold || length.y<deflection_threshold){
////                                    is_minimum = false;
////                                }                            
////                            }
//                        }
                    }
                    
                    if(is_minimum) minimas.add(f);                    
                }
                
                int num_minimas = minimas.size();
                for(int m=0; m<num_minimas; m++)
                {                        
                    int f_start = minimas.get(m);          
                    int f_end = f_start+frame_window-1;                    
                    MyPoint R = Rsp[f_start]; 
                    
                    if(merging){
                        while((m+1)<num_minimas && minimas.get(m+1)<f_end){                        
                            int f = minimas.get(m+1);
                            int ff_end = f+frame_window-1;

                            MyPoint RR = cross_correlation(proj_f, f_start, ff_end);
                            double rs = RR.x;
                            double rp = RR.y;
                            if(Math.abs(rs)<threshold_rs && rp<threshold_rp){
                                R = RR;                            
                                f_end = ff_end;
                                m++;                            
                            }
                            else break;
                        }
                    }
                    
                    double rs = R.x;
                    double rp = R.y;
                    ContractionUnit cu = new ContractionUnit(k,j,f_start,f_end,rs,rp);
                    cu_list.add(cu);
                }
            }
            IJ.showProgress(k, n);
        }
        IJ.showStatus("hunting done");
        return cu_list;
    }
    
//    public static List<ContractionUnit> hunt(double[]gx, double[]gy, double catch_radius, double[][] cx, double[][] cy, double pixel_size, 
//            double threshold_rs, double threshold_rp, double threshold_deflection,int frame_window, boolean merging){        
//        List[] neigbors = search_neighbors(gx, gy, catch_radius);
//        int n = neigbors.length;
//        MyPoint[] anchors = new MyPoint[n];
//        for(int k=0; k<n; k++) anchors[k] = new MyPoint(gx[k], gy[k]);
//        List[] distance = get_distance(neigbors, anchors);                        
//        int frames = cx[0].length;
//        IJ.showStatus("projecting");
//        List<CU_Projection>[][] proj_list = new List[frames][];     
//        for(int f=0; f<frames; f++){
//            MyPoint[] centers = new MyPoint[n];
//            for(int k=0; k<n; k++) centers[k] = new MyPoint(cx[k][f]/pixel_size, cy[k][f]/pixel_size);                
//            proj_list[f] = project(neigbors, distance, anchors, centers);
//            IJ.showProgress(f, frames);
//        }
//        
//        IJ.showStatus("hunting contractile unit");
//        List<ContractionUnit> cu_list = new ArrayList();  
//        //frame_window = 500;
//        if(frame_window<20 || frame_window>frames) frame_window = frames;  
//        //int radius_f = frame_window/2;
//        int num_Rsp = frames-frame_window+1;
//        for(int k=0; k<n; k++){
//            List<Integer> list_nb = neigbors[k];            
//            for(int i=0; i<list_nb.size(); i++){
//                int j = list_nb.get(i);
//                CU_Projection[] proj_f = new CU_Projection[frames];
//                for(int f=0; f<frames; f++) proj_f[f] = proj_list[f][k].get(i);
//                MyPoint[] Rsp = new MyPoint[num_Rsp];
//                for(int f=0; f<num_Rsp; f++){
//                    int f_end = f+frame_window-1;
//                    MyPoint R = cross_correlation(proj_f, f, f_end);
//                    Rsp[f] = R;
//                }
//                //search minimum rp
//                List<Integer> minimas = new ArrayList();
//                for(int f=0; f<num_Rsp; f++){
//                    int f_end = f+frame_window-1;
//                    MyPoint R = Rsp[f];
//                    double rs = R.x;
//                    double rp = R.y;
//                    boolean is_minimum = false;
//                    if(Math.abs(rs)<threshold_rs && rp<threshold_rp){
//                        is_minimum = true;
//                        if(f-1>=0 && Rsp[f-1].y<rp) is_minimum = false;
//                        else if(f+1<num_Rsp && Rsp[f+1].y<rp) is_minimum = false;
//                        
//                        if(is_minimum){
//                            //if(percentage_contracting(proj_f, f, f_end)<80) is_minimum = false;
//                            //else
//                            {
//                                double deflection_threshold = threshold_deflection/pixel_size;//20/pixel_size;
//                                MyPoint length = get_length(proj_f, f, f_end);
//                                if(length.x<deflection_threshold || length.y<deflection_threshold){
//                                    is_minimum = false;
//                                }                            
//                            }
//                        }
//                    }
//                    
//                    if(is_minimum) minimas.add(f);                    
//                }
//                
//                int num_minimas = minimas.size();
//                for(int m=0; m<num_minimas; m++)
//                {                        
//                    int f_start = minimas.get(m);          
//                    int f_end = f_start+frame_window-1;
//                    
//                    MyPoint R = Rsp[f_start];                                        
//                    while((m+1)<num_minimas && minimas.get(m+1)<f_end){                        
//                        int f = minimas.get(m+1);
//                        int ff_end = f+frame_window-1;
//                        
//                        MyPoint RR = cross_correlation(proj_f, f_start, ff_end);
//                        double rs = RR.x;
//                        double rp = RR.y;
//                        if(Math.abs(rs)<threshold_rs && rp<threshold_rp){
//                            R = RR;                            
//                            f_end = ff_end;
//                            m++;                            
//                        }
//                        else break;
//                    }
//                    
//                    double rs = R.x;
//                    double rp = R.y;
//                    ContractionUnit cu = new ContractionUnit(k,j,f_start,f_end,rs,rp);
//                    cu_list.add(cu);
//                }
//            }
//            IJ.showProgress(k, n);
//        }
//        IJ.showStatus("hunting done");
//        return cu_list;
//    }
    
//    public static List<ContractionUnit> hunt_frames(double[]gx, double[]gy, double catch_radius, double[][] cx, double[][] cy, double pixel_size, double threshold_rs, double threshold_rp){        
//        List[] neigbors = search_neighbors(gx, gy, catch_radius);
//        int n = neigbors.length;
//        MyPoint[] anchors = new MyPoint[n];
//        for(int k=0; k<n; k++) anchors[k] = new MyPoint(gx[k], gy[k]);
//        List[] distance = get_distance(neigbors, anchors);                        
//        int frames = cx[0].length;
//        List<ContractionUnit> cu_list = new ArrayList();  
//        for(int f=0; f<frames; f++){
//            List<ContractionUnit> cu_list_f = hunt_one_frame(neigbors, anchors, distance, cx, cy, f, pixel_size, threshold_rs, threshold_rp);
//            int num = cu_list_f.size();
//            for(int i=0; i<num; i++) cu_list.add(cu_list_f.get(i));
//        }
////        IJ.showStatus("projecting");
////        MyPoint[] centers = new MyPoint[n];
////        for(int k=0; k<n; k++) centers[k] = new MyPoint(cx[k][frame]/pixel_size, cy[k][frame]/pixel_size);                
////        List<CU_Projection>[] proj_list = project(neigbors, distance, anchors, centers);
////        
////        IJ.showStatus("hunting contractile unit");
////        List<ContractionUnit> cu_list = new ArrayList();  
////        
////        for(int k=0; k<n; k++){
////            List<Integer> list_nb = neigbors[k];            
////            for(int i=0; i<list_nb.size(); i++){
////                int j = list_nb.get(i);                
////                CU_Projection proj_f = proj_list[k].get(i);
////                
////                double rs = Math.abs(proj_f.proj_s1) + Math.abs(proj_f.proj_s2);
////                double rp = Math.abs(proj_f.proj_p1) + Math.abs(proj_f.proj_p2);
////                boolean is_minimum = false;
////                if(Math.abs(proj_f.cosa1)>threshold_rs && Math.abs(proj_f.cosa2)>threshold_rs && 
////                        proj_f.proj_p1*proj_f.proj_p2<0 && 
////                        Math.abs(proj_f.proj_p1)>threshold_rp && Math.abs(proj_f.proj_p2)>threshold_rp){
////                    is_minimum = true;   
////                }
////
////                if(is_minimum){
////                    ContractionUnit cu = new ContractionUnit(k,j,frame,frame,rs,rp);
////                    cu_list.add(cu);
////                }
////            }
////            IJ.showProgress(k, n);
////        }
////        IJ.showStatus("hunting done");
//        return cu_list;
//    }
//    
//    public static List<ContractionUnit> hunt_one_frame(List[] neigbors, MyPoint[] anchors, List[] distance, double[][] cx, double[][]cy, int frame, double pixel_size, double threshold_rs, double threshold_rp){        
////        List[] neigbors = search_neighbors(gx, gy, catch_radius);
//        int n = neigbors.length;
////        MyPoint[] anchors = new MyPoint[n];
////        for(int k=0; k<n; k++) anchors[k] = new MyPoint(gx[k], gy[k]);
////        List[] distance = get_distance(neigbors, anchors);                        
//        
//        IJ.showStatus("projecting");
//        MyPoint[] centers = new MyPoint[n];
//        for(int k=0; k<n; k++) centers[k] = new MyPoint(cx[k][frame]/pixel_size, cy[k][frame]/pixel_size);                
//        List<CU_Projection>[] proj_list = project(neigbors, distance, anchors, centers);
//        
//        IJ.showStatus("hunting contractile unit");
//        List<ContractionUnit> cu_list = new ArrayList();  
//        
//        for(int k=0; k<n; k++){
//            List<Integer> list_nb = neigbors[k];            
//            for(int i=0; i<list_nb.size(); i++){
//                int j = list_nb.get(i);                
//                CU_Projection proj_f = proj_list[k].get(i);
//                
//                double rs = Math.abs(proj_f.proj_s1) + Math.abs(proj_f.proj_s2);
//                double rp = Math.abs(proj_f.proj_p1) + Math.abs(proj_f.proj_p2);
//                boolean is_minimum = false;
//                if(Math.abs(proj_f.cosa1)>threshold_rs && Math.abs(proj_f.cosa2)>threshold_rs && 
//                        proj_f.proj_p1*proj_f.proj_p2<0 && 
//                        Math.abs(proj_f.proj_p1)>threshold_rp && Math.abs(proj_f.proj_p2)>threshold_rp){
//                    is_minimum = true;   
//                }
//
//                if(is_minimum){
//                    ContractionUnit cu = new ContractionUnit(k,j,frame,frame,rs,rp);
//                    cu_list.add(cu);
//                }
//            }
//            IJ.showProgress(k, n);
//        }
//        IJ.showStatus("hunting done");
//        return cu_list;
//    }
}
