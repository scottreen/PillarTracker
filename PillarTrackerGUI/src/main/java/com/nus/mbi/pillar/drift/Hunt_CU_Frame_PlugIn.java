package com.nus.mbi.pillar.drift;

import com.nus.mbi.pillar.detection.ContractionUnit;
import static com.nus.mbi.pillar.drift.Hunt_CU_PlugIn.get_distance;
import static com.nus.mbi.pillar.drift.Hunt_CU_PlugIn.project;
import ij.IJ;
import ij.gui.GenericDialog;
import java.util.ArrayList;
import java.util.List;
import com.nus.mbi.pillar.stat.MyPoint;
import static com.nus.mbi.pillar.drift.Hunt_CU_PlugIn.search_neighbors;

/**
 *
 * @author xiaochun
 */
public class Hunt_CU_Frame_PlugIn extends Hunt_CU_PlugIn{
    public double threshold_deflection = 50;
    public double threshold_angle_diff = 10;    
    public boolean process_all_frames = true;
    
    public double getCosThreshold(){
        return Math.cos(Math.PI*threshold_angle_diff/180);
    }

    public Hunt_CU_Frame_PlugIn(double catch_radius) {
        super(catch_radius);
    }
    
    public Hunt_CU_Frame_PlugIn(double catch_radius, double threshold_deflection, double threshold_angle_diff) {
        this.catch_radius = catch_radius;
        this.threshold_deflection = threshold_deflection;
        this.threshold_angle_diff = threshold_angle_diff;
    }
    
    @Override
    public boolean showDialog(int nframes){
        GenericDialog gd = new GenericDialog("Hunt Contractile Unit(CU)");
        gd.addNumericField("catch_radius:", catch_radius, 2, 10, "pixel");
        gd.addNumericField("Deflection_Threshold:", threshold_deflection, 2, 10, "nm");
        gd.addNumericField("Angle_Difference_Threshold:", threshold_angle_diff, 2, 10, "degree");
        gd.addCheckbox("show CU table?", show_table);
        gd.addCheckbox("process all frames? " + nframes, process_all_frames);
        
        gd.showDialog();
        if (gd.wasCanceled()) return false;
        
        catch_radius = (double)gd.getNextNumber();            
        threshold_deflection = (double)gd.getNextNumber();
        threshold_angle_diff = (double)gd.getNextNumber();
        show_table = gd.getNextBoolean();
        process_all_frames = gd.getNextBoolean();
        if (catch_radius<=0) {IJ.showMessage("the catch radius must be greater than zero");  return false;}
        if (threshold_deflection<0) {IJ.showMessage("deflection threshold must be greater than zero"); return false;}
        if (Math.abs(threshold_angle_diff)>90) {IJ.showMessage("angle difference threshold must be smaller than 90"); return false;}
        
        return true;
    }
    
    public static List<ContractionUnit> hunt(double[]gx, double[]gy, double catch_radius, double[][] cx, double[][] cy, double pixel_size, double threshold_deflections, double threshold_cosa){        
        List[] neigbors = search_neighbors(gx, gy, catch_radius);
        int n = neigbors.length;
        MyPoint[] anchors = new MyPoint[n];
        for(int k=0; k<n; k++) anchors[k] = new MyPoint(gx[k], gy[k]);
        List[] distance = get_distance(neigbors, anchors);                        
        int frames = cx[0].length;
        List<ContractionUnit> cu_list = new ArrayList();  
        for(int f=0; f<frames; f++){
            List<ContractionUnit> cu_list_f = hunt_one_frame(neigbors, anchors, distance, cx, cy, f, pixel_size, threshold_deflections, threshold_cosa);
            int num = cu_list_f.size();
            for(int i=0; i<num; i++) cu_list.add(cu_list_f.get(i));
            IJ.showProgress(f, frames);
        }
        IJ.showStatus("hunting done");
        return cu_list;
    }
    
    public static List<ContractionUnit> hunt(double[]gx, double[]gy, double catch_radius, int frame, double[][] cx, double[][] cy, double pixel_size, double threshold_deflections, double threshold_cosa){        
        List[] neigbors = search_neighbors(gx, gy, catch_radius);
        int n = neigbors.length;
        MyPoint[] anchors = new MyPoint[n];
        for(int k=0; k<n; k++) anchors[k] = new MyPoint(gx[k], gy[k]);
        List[] distance = get_distance(neigbors, anchors);                        
        int frames = cx[0].length;
        List<ContractionUnit> cu_list = new ArrayList();  
        if(frame<0 || frame>=frames) return cu_list;
        
        List<ContractionUnit> cu_list_f = hunt_one_frame(neigbors, anchors, distance, cx, cy, frame, pixel_size, threshold_deflections, threshold_cosa);
        int num = cu_list_f.size();
        for(int i=0; i<num; i++) cu_list.add(cu_list_f.get(i));
        IJ.showStatus("hunting done");
        return cu_list;
    }
   
    public static List<ContractionUnit> hunt_one_frame(List[] neigbors, MyPoint[] anchors, List[] distance, double[][] cx, double[][]cy, int frame, double pixel_size, double threshold_deflection, double threshold_cosa){        
        int n = neigbors.length;

        IJ.showStatus("projecting");
        MyPoint[] centers = new MyPoint[n];
        for(int k=0; k<n; k++) centers[k] = new MyPoint(cx[k][frame]/pixel_size, cy[k][frame]/pixel_size);                
        List<CU_Projection>[] proj_list = project(neigbors, distance, anchors, centers);
        
        IJ.showStatus("hunting contractile unit");
        List<ContractionUnit> cu_list = new ArrayList();  
        
        threshold_deflection = threshold_deflection/pixel_size;
        for(int k=0; k<n; k++){
            List<Integer> list_nb = neigbors[k];            
            for(int i=0; i<list_nb.size(); i++){
                int j = list_nb.get(i);                
                CU_Projection proj_f = proj_list[k].get(i);
                
                double rs = Math.abs(proj_f.proj_s1) + Math.abs(proj_f.proj_s2);
                double rp = Math.abs(proj_f.proj_p1) + Math.abs(proj_f.proj_p2);
                boolean is_minimum = false;
                if(Math.abs(proj_f.cosa1)>threshold_cosa && Math.abs(proj_f.cosa2)>threshold_cosa && 
                        proj_f.proj_p1*proj_f.proj_p2<0 && 
                        Math.abs(proj_f.proj_p1)>threshold_deflection && Math.abs(proj_f.proj_p2)>threshold_deflection){
                    is_minimum = true;   
                }

                if(is_minimum){
                    ContractionUnit cu = new ContractionUnit(k,j,frame,frame,rs,rp);
                    cu_list.add(cu);
                }
            }
            IJ.showProgress(k, n);
        }
        IJ.showStatus("hunting done");
        return cu_list;
    }
}
