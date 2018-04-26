package com.nus.mbi.pillar.grid;

import ij.IJ;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import static com.nus.mbi.pillar.stat.BasicStatisitic.*;
import com.nus.mbi.pillar.stat.PloyFitter;

/**
 *
 * @author xiaochun
 */
public class GridCreation {
    private double[] xpoints;
    private double[] ypoints;
    private double[] xround;
    private double[] yround;
    private int[] ix = null;
    private int[] iy = null;
    private boolean[] is_labeled = null;
    private int ix_off;
    private int iy_off;
    private double[][] map;
    private int map_width = 0;
    private int map_height = 0;
    private boolean is_grid_labelled = false;
    private boolean is_grid_estimated = false;
    private double estimated_seperation;
    private double estimated_grid_oblique;
    
    private double[] dx = null;
    private double[] dy = null;
    private boolean[] active_d = null;
    private boolean is_grid_created = false;
    
    public int getMapWidth(){
        return map_width;
    }
    
    public int getMapHeight(){
        return map_height;
    }
    
    public double[] getDX(){
        return dx;
    }
    
    public double[] getDY(){
        return dy;
    }
    
    public int[] getIX(){
        return ix;
    }
    
    public int[] getIY(){
        return iy;
    }
    
    public boolean[] IsLabeled(){
        return is_labeled;
    }
    
    public boolean IsGridCreated(){
        return is_grid_created;
    }
    
    public boolean IsGridLabelled(){
        return is_grid_labelled;
    }
    
    public boolean IsGridEstimated(){
        return is_grid_estimated;
    }
    
    public double GetEstimatedSeperation(){
        return estimated_seperation;
    }
    
    public double GetEstimatedGridOblique(){
        return estimated_grid_oblique;
    }
    
    public GridCreation(double[] x, double[] y){
        xpoints = x;
        ypoints = y;    
        is_grid_labelled = false;
        is_grid_estimated = false;
        is_grid_created = false;
    }
    
    public GridCreation(double[] x, double[] y, int[] ix, int[] iy, double[] dx, double[] dy, boolean[] dxy_active, boolean[] label_flag){
        this.xpoints = x;
        this.ypoints = y;            
        this.ix = ix;
        this.iy = iy;
        this.dx = dx;
        this.dy = dy;
        //this.active_d = dxy_active;
        int npoints = x.length;
        is_labeled = new boolean[npoints];
        active_d = new boolean[npoints];
        for(int i=0; i<npoints; i++){
            is_labeled[i] = label_flag[i];
            active_d[i] = dxy_active[i];
        }     
        create_map();
        is_grid_labelled = true;        
        is_grid_estimated = false;
        is_grid_created = false;
    }
    
    public boolean[] getActive(){
        return active_d;
    }
    
    public void create_map(){
        double max_x = max(xpoints);
        double max_y = max(ypoints);
        map_width = (int)Math.ceil(max_x);
        map_height = (int)Math.ceil(max_y);
        map = new double[map_height][map_width];
        for(int i=0; i<map_width; i++){
            for(int j=0; j<map_height; j++){
                map[j][i] = 0;
            }
        }
        int npillars = xpoints.length;
        xround = new double[npillars];
        yround = new double[npillars];
                
        for(int i=0; i<npillars; i++){
            int rxx = (int)Math.round(xpoints[i]);
            int ryy = (int)Math.round(ypoints[i]);
            xround[i] = rxx;
            yround[i] = ryy;
            if(rxx>=0 && ryy>=0 && rxx<map_width && ryy<map_height)
                map[ryy][rxx] = i+1;
        }
    }
    
    public double[][] search_best_pattern(double[][] map_image, double[][] search_area, double[] rx, double[] ry, double[] nx, double[] ny, double score_wall){
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
    
    public double[] search_neighbor(double x, double y, double[][] map_image, int map_w, int map_h, double[][] area, int area_size){               
        double min_dis = -1;
        double min_p = -1;
        int xc = (int)Math.round(x);
        int yc = (int)Math.round(y);
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
    
    public double[] search_neighbor(double x, double y, double[] cx, double[] cy){               
        double min_dis = -1;
        double min_p = -1;
        int n = cx.length;
        for(int k=0; k<n; k++){           
           double dx = cx[k] - x;
           double dy = cy[k] - y;
           double dis = dx*dx+dy*dy;
           if(min_dis<0 || min_dis>dis){
                min_dis = dis;
                min_p = k;
           }
        }   
        
        double[] ret = {min_p, min_dis};
        return ret;
    }
    
    public static double estimate_seperation(double[] cx, double[] cy){        
        int n = cx.length;        
        double avg_dis = 0;        
        int num = 0;   
        List<Double> list = new ArrayList<Double>();
        for(int i=0; i<n; i++){
            double min_dis = -1;
            for(int k=i+1; k<n; k++){           
               double dx = cx[k] - cx[i];
               double dy = cy[k] - cy[i];
               double dis = dx*dx+dy*dy;
               if(min_dis<0 || min_dis>dis) min_dis = dis;               
            }
            if(min_dis>0){
                list.add(min_dis);
                avg_dis+=min_dis;
                num++;
            }
        }
        
        if(num<=0) return Double.NaN;
        
        Collections.sort(list);
        double s = Math.sqrt(list.get(num/2)); 
        IJ.log("estimate grid: s="+s + " num="+ num);
        return s;
    }
    
    public double[] estimate_square_grid_properties(double seperation, double tol){        
        double[] sa = estimate_square_grid_properties(xpoints, ypoints, seperation,tol);
        if(sa==null){
            is_grid_estimated = false;
            return null;
        }
        
        estimated_seperation = sa[0];
        estimated_grid_oblique = sa[1];
        is_grid_estimated = true;
        return sa;
    }
    
    public double[] estimate_square_grid_properties(){        
        if(!is_grid_labelled) return null;
        return estimate_square_grid_properties(xpoints, ypoints, ix, iy, is_labeled);        
    }
    
    public int getFisrtlabeled(){
        if(!is_grid_labelled) return -1;
        int index = -1;
        double dis = Double.MAX_VALUE;
        for(int i=0; i<is_labeled.length; i++){
            if(is_labeled[i]){
                double dx = xpoints[i];
                double dy = ypoints[i];
                double d2 = dx*dx+dy*dy;
                if(d2<dis){
                    index = i;
                    dis = d2;
                }
            }
        }
        return index;
    }
    
    public double[] getFisrtlabeledXY(){
        int i = getFisrtlabeled();        
        if(i<0) return null;
        double[] xy = {xpoints[i], ypoints[i]};
        return xy;
    }
    
    // the grid model without pinchusion distortion:
    // x[i][j] = x0 + i*dx - j*dy;
    // y[i][j] = y0 + i*dy + j*dx;    
    public static double[] estimate_square_grid_properties(double[] cx, double[] cy, int[] ix, int[] iy, boolean[] is_labeled){
        int npillars = is_labeled.length;
        int ixmin = find_min(ix, is_labeled);        
        int iymin = find_min(iy, is_labeled);
        if(ixmin==Integer.MAX_VALUE || iymin==Integer.MAX_VALUE)  return null;
        
        double xoff = ixmin - 1;
        double sigmaI  = 0;
        double sigmaI2 = 0;
        double sigmaX  = 0;
        double sigmaXI = 0;
        
        double yoff = iymin - 1;
        double sigmaJ  = 0;
        double sigmaJ2 = 0;
        double sigmaY  = 0;
        double sigmaYJ = 0;
        
        double sigmaYI = 0;
        double sigmaXJ = 0;
        
        int num = 0;
        for(int k=0; k<npillars; k++){
            if(is_labeled[k]){
                double i=ix[k]-xoff;
                double x=cx[k];
                sigmaI  += i;
                sigmaI2 += i*i;
                sigmaX  += x;
                sigmaXI += x*i;
                
                double j=iy[k]-yoff;
                double y=cy[k];
                sigmaJ  += j;
                sigmaJ2 += j*j;
                sigmaY  += y;
                sigmaYJ += y*j;
                
                sigmaYI += y*i;
                sigmaXJ += x*j;
                
                num++;
            }
        }      
        double sigmaIJ2 = sigmaI2+sigmaJ2;        
        double det = num*sigmaIJ2 - (sigmaI*sigmaI + sigmaJ*sigmaJ);
        if(det == 0) return null;
        
        double sigmaXIYJ = sigmaXI+sigmaYJ;
        double sigmaYIXJ = sigmaYI-sigmaXJ;
        double x0 = sigmaX*sigmaIJ2 - (sigmaI*sigmaXIYJ - sigmaJ*sigmaYIXJ);
        double y0 = sigmaY*sigmaIJ2 - (sigmaI*sigmaYIXJ + sigmaJ*sigmaXIYJ);
        double dx = num*sigmaXIYJ - (sigmaI*sigmaX + sigmaJ*sigmaY);
        double dy = num*sigmaYIXJ - (sigmaI*sigmaY - sigmaJ*sigmaX);
        x0 = x0/det;
        y0 = y0/det;
        dx = dx/det;
        dy = dy/det;
        /*
        // the following code has the problem, only works on the points along the line.
        double detI = num*sigmaI2 - sigmaI*sigmaI;
        double detJ = num*sigmaJ2 - sigmaJ*sigmaJ;        
        if(detI==0 || detJ==0) return null;
        
        double x0   = (sigmaX*sigmaI2 - sigmaI*sigmaXI)/detI;
        double dx   = (num*sigmaXI - sigmaI*sigmaX)/detI;
        
        double y0   = (sigmaY*sigmaJ2 - sigmaJ*sigmaYJ)/detJ;
        double dy   = (num*sigmaYJ - sigmaJ*sigmaY)/detJ;
        */
        double s = Math.sqrt(dx*dx+dy*dy);
        double a = Math.atan2(-dy, dx)*180/Math.PI;
        IJ.log("estimate after grid creation: x0=" + x0 + " y0="+y0 + " s="+s + " a="+a);        
        double[] sa = {x0, y0, dx, dy};
        return sa;
    }
    
    public static double[] estimate_line_properties(double[] cx ,double[] cy, double ir2, double or2){                 
        int n = cx.length;        
        double avg_dis = 0;        
        int num = 0;
        double avg_dx = 0;
        double avg_dy = 0;
        int num_dxy = 0;
        for(int i=0; i<n; i++){
            for(int k=i+1; k<n; k++){           
               double dx = cx[k] - cx[i];
               double dy = cy[k] - cy[i];
               if(dx>=0){
                        avg_dx += dx;
                        avg_dy += dy;
                        num_dxy++;
               }
               double dis = dx*dx+dy*dy;               
               if(dis>=ir2 && dis<=or2){
                    avg_dis += dis;                    
                    num++;     
               }
            }
        }
        
        if(num<=0 || num_dxy<=0) return null;
        
        double s = Math.sqrt(avg_dis/num);
        double a = Math.atan2(-avg_dy, avg_dx)*180/Math.PI;         
        double[] ret = {s, a};        
        return ret;
    }
    
    
    public static double[] esimate_grid_properties(double[] cx, double[] cy, int[] ix, int[] iy, double seperation, double tol){        
        double r = Math.abs(seperation);
        double t = Math.abs(tol);
        double out_radius = r+t;        
        double or2 = out_radius*out_radius;
        double in_radius = r-t;
        double ir2 = in_radius*in_radius; 
        
        //int npillars = is_labeled.length;
        int ixmin = find_min(ix);        
        int iymin = find_min(iy);
        if(ixmin==Integer.MAX_VALUE || iymin==Integer.MAX_VALUE)  return null;
        
        int ixmax = find_max(ix);        
        int iymax = find_max(iy);
        if(ixmax==Integer.MIN_VALUE || iymax==Integer.MIN_VALUE)  return null;
        
        int xdim = ixmax - ixmin + 1;
        int ydim = iymax - iymin + 1;
        
        double avg_s = 0;
        double avg_a = 0;
        int num = 0;
        for(int ixx = ixmin; ixx<=ixmax; ixx++){
            List<Integer> list = getSelectedIndex(ixx, ix);
            if(list.size()>3){
                double[][] xy = getSelectedPoints(list, cx, cy);
                double[] sa = estimate_line_properties(xy[0], xy[1], ir2, or2);                
                if(sa!=null){                    
                    avg_s += sa[0];
                    avg_a += sa[1];
                    num++;
                }
            }            
        }
        
        if(num<=0) return null;
        double sx = avg_s/num;
        double ax = avg_a/num; 
        IJ.log("estimate grid in x axis: s="+sx + " num="+ num + " a="+ax);
        
        avg_s = 0;
        avg_a = 0;
        num = 0;
        for(int iyy = iymin; iyy<=iymax; iyy++){
            List<Integer> list = getSelectedIndex(iyy, iy);
            if(list.size()>3){
                double[][] xy = getSelectedPoints(list, cx, cy);                
                double[] sa = estimate_line_properties(xy[0], xy[1], ir2, or2);
                
                if(sa!=null){                    
                    avg_s += sa[0];
                    avg_a += sa[1];
                    num++;
                }
            }            
        }
        
        if(num<=0) return null;
        double sy = avg_s/num;
        double ay = avg_a/num; 
        IJ.log("estimate grid in y axis: s="+sy + " num="+ num + " a="+ay);
        double ss = (sx+sy)/2.0;
        double oo = ay; //if(a>45) a=a-90;
        double aa = Math.abs(ax-ay);
        double[] ret = {ss, oo, aa};    
        IJ.log("estimate grid: s="+ss + " oblique="+oo + " angle="+aa);
        return ret;
    }
    
    
    public static double[] estimate_square_grid_properties(double[] cx ,double[] cy, double seperation, double tol){        
        double r = Math.abs(seperation);
        double t = Math.abs(tol);
        double out_radius = r+t;        
        double r2 = out_radius*out_radius;
        double in_radius = r-t;
        double ir2 = in_radius*in_radius;
        
        int n = cx.length;        
        double avg_dis = 0;        
        int num = 0;
        double avg_dx = 0;
        double avg_dy = 0;
        int num_dxy = 0;
        for(int i=0; i<n; i++){
            for(int k=i+1; k<n; k++){           
               double dx = cx[k] - cx[i];
               double dy = cy[k] - cy[i];
               double dis = dx*dx+dy*dy;
               if(dis>=ir2 && dis<=r2){
                   avg_dis += dis;
                   num++;
                   if(dx>=0 && dy<=0){
                        avg_dx += dx;
                        avg_dy += dy;
                        num_dxy ++;
                   }
               }        
               
            }
        }
        
        if(num<=0 || num_dxy<=0) return null;
        
        double s = Math.sqrt(avg_dis/num);
        double a = Math.atan2(-avg_dy, avg_dx)*180/Math.PI; 
        if(a>45) a=a-90;
        IJ.log("estimate grid: s="+s + " num="+ num + " a="+a+" num_dxy="+num_dxy);
        double[] ret = {s, a};        
        return ret;
    }
        
    public void label_grid(double seperation, double grid_oblique, double grid_angle){
        SearchAreaCreation area_creator = new SearchAreaCreation();        
        double[][] area = area_creator.create_search_area(seperation, grid_oblique, grid_angle);
        double[] nx = area_creator.getNX();
        double[] ny = area_creator.getNY();
        double r = area_creator.getOuterRadius();
        double score_wall = r*r;
        
        label_grid(area, nx, ny, score_wall);        
    }
    
    public static void paint(int k, int[] ix, int[] iy, boolean[] is_labeled, double[][]np){
       //% search the neighbours along x forward
        int p = (int)(np[0][k]-1);
        if( p>=0 && !is_labeled[p]){
                ix[p] = ix[k] + 1;
                iy[p] = iy[k];
                is_labeled[p] = true;
                paint(p, ix, iy, is_labeled, np);
        }
        //% search the neighbours along y forward
        p = (int)(np[1][k]-1);
        if( p>=0 && !is_labeled[p]){   
                ix[p] = ix[k] ;
                iy[p] = iy[k] + 1;
                is_labeled[p] = true;
                paint(p, ix, iy, is_labeled, np);
        }
        //% search the neighbours along x backward
        p = (int)(np[2][k]-1);
        if( p>=0 && !is_labeled[p]){   
                ix[p] = ix[k]  - 1;
                iy[p] = iy[k];
                is_labeled[p] = true;
                paint(p, ix, iy, is_labeled, np);
        }
        //% search the neighbours along y backward
        p = (int)(np[3][k]-1);
        if( p>=0 && !is_labeled[p]){   
                ix[p] = ix[k] ;
                iy[p] = iy[k] - 1;
                is_labeled[p] = true;
                paint(p, ix, iy, is_labeled, np);
        }
    }
    
    public void label_grid(double[][] area, double[] nx, double[] ny, double score_wall){        
        int npillars = xpoints.length;
        create_map();
        double[][] np= search_best_pattern(map, area, xround, yround, nx,ny, score_wall);
        
        double[] match_error = np[4];
        double min_err = match_error[0];
        int min_err_index = 0;
        for(int i=1; i<match_error.length; i++){
            if(match_error[i]<min_err){
                min_err = match_error[i];
                min_err_index = i;
            }
        }
        
        int init_pillar = 0;   
        if(match_error[init_pillar] > min_err) init_pillar = min_err_index;
        ix = new int[npillars];
        iy = new int[npillars];
        is_labeled = new boolean[npillars];
        boolean[] is_painted = new boolean[npillars];
        for(int i=0; i<npillars; i++){
            ix[i] = Integer.MAX_VALUE;
            iy[i] = Integer.MAX_VALUE;
            is_labeled[i] = false;
            is_painted[i] = false;
        }
        //% the first pillar will be seeded as 2D index (0,0)
        //%ix(init_pillar) = 0;
       // %iy(init_pillar) = 0;
        is_labeled[init_pillar] = true;
        ix[init_pillar] = 0;
        iy[init_pillar] = 0;
        //%is_painted(init_pillar) = false;
        paint(init_pillar, ix, iy, is_labeled, np);
        //%% label the neighbours
//        int num_new_labels = 1;
//        while(num_new_labels>0){
//            num_new_labels = 0;
//            for(int k=0; k<npillars; k++){
//                if(is_labeled[k] && !is_painted[k]){
//                    is_painted[k] = true;
//                    //% search the neighbours along x forward
//                    int p = (int)(np[0][k]-1);
//                    if( p>=0 && !is_labeled[p]){
//                            ix[p] = ix[k] + 1;
//                            iy[p] = iy[k];
//                            is_labeled[p] = true;
//                            num_new_labels = num_new_labels+1;
//                    }
//                    //% search the neighbours along y forward
//                    p = (int)(np[1][k]-1);
//                    if( p>=0 && !is_labeled[p]){   
//                            ix[p] = ix[k] ;
//                            iy[p] = iy[k] + 1;
//                            is_labeled[p] = true;
//                            num_new_labels = num_new_labels+1;
//                    }
//                    //% search the neighbours along x backward
//                    p = (int)(np[2][k]-1);
//                    if( p>=0 && !is_labeled[p]){   
//                            ix[p] = ix[k]  - 1;
//                            iy[p] = iy[k];
//                            is_labeled[p] = true;
//                            num_new_labels = num_new_labels+1;
//                    }
//                    //% search the neighbours along y backward
//                    p = (int)(np[3][k]-1);
//                    if( p>=0 && !is_labeled[p]){   
//                            ix[p] = ix[k] ;
//                            iy[p] = iy[k] - 1;
//                            is_labeled[p] = true;
//                            num_new_labels = num_new_labels+1;
//                    }                
//                }
//            }
//        }
        
        is_grid_labelled = true;
        //estimate_square_grid_properties();
        //estimate_square_grid_properties(xpoints,ypoints,ix,iy,is_labeled);
    }
    
    public List<Integer> getXLines(int ixx){
//        if(!is_grid_labelled) return null;
//        int npillars = is_labeled.length;
//        List<Integer> list = new ArrayList<Integer>();
//        
//        for(int i=0; i<npillars; i++){
//            if(is_labeled[i] && ix[i] == ixx) list.add(i);            
//        }
//        
//        return list;
        return getXLines(ixx, null);
    }
    
    public List<Integer> getYLines(int iyy){
//        if(!is_grid_labelled) return null;
//        int npillars = is_labeled.length;
//        List<Integer> list = new ArrayList<Integer>();
//        
//        for(int i=0; i<npillars; i++){
//            if(is_labeled[i] && iy[i] == iyy) list.add(i);            
//        }
//        
//        return list;
        return getYLines(iyy, null);
    }
    
    
    public List<Integer> getXLines(int ixx, boolean[] active){
        if(!is_grid_labelled) return null;
        int npillars = is_labeled.length;
        List<Integer> list = new ArrayList<Integer>();
        
        if(active==null){
            for(int i=0; i<npillars; i++){
                if(is_labeled[i] && ix[i] == ixx) list.add(i);            
            }
            return list;
        }       
        
        for(int i=0; i<npillars; i++){
            if(active[i]){
                if(is_labeled[i] && ix[i] == ixx) list.add(i);            
            }
        }
        
        return list;
    }
    
    public List<Integer> getYLines(int iyy, boolean[] active){
        if(!is_grid_labelled) return null;
        int npillars = is_labeled.length;
        List<Integer> list = new ArrayList<Integer>();
        
        if(active==null){
            for(int i=0; i<npillars; i++){
                if(is_labeled[i] && iy[i] == iyy) list.add(i);            
            }
            return list;
        }       
        
        for(int i=0; i<npillars; i++){
            if(active[i]){
                if(is_labeled[i] && iy[i] == iyy) list.add(i);            
            }
        }
        
        return list;
    }
    
    public double[][] getXYPoints(List<Integer> list){
        if(list==null || list.isEmpty()) return null;
        int n = list.size();
        double[][] xy = new double[2][n];
        for(int i=0; i<n; i++){
            int ii = list.get(i);
            xy[0][i] = xpoints[ii];
            xy[1][i] = ypoints[ii];
        }
        return xy;
    }
        
    public static List<Integer> getSelectedIndex(int value, int[] array){
        int len = array.length;
        List<Integer> list = new ArrayList<Integer>();
        
        for(int i=0; i<len; i++) {
            if(array[i] == value) list.add(i);            
        }
        
        return list;
    }
        
    public static double[][] getSelectedPoints(List<Integer> list, double[] x, double[] y){
        if(list==null || list.isEmpty()) return null;
        int n = list.size();
        double[][] xy = new double[2][n];
        for(int i=0; i<n; i++){
            int ii = list.get(i);
            xy[0][i] = x[ii];
            xy[1][i] = y[ii];
        }
        return xy;
    }
    
    public void create_grid_points(double cacth_radius){
        create_grid_points(cacth_radius, null);
    }
    
    private int get_num_points_line(boolean[] active, int[] index, int target){
        if(!is_grid_labelled) return 0;
        int npillars = active.length;
        int num = 0;
        for(int i=0; i<npillars; i++){               
            if(active[i] && index[i]==target) num++;
        }
        return num;
    }
    
    private int get_num_points_zero(){
        if(!is_grid_labelled) return 0;
        int npillars = is_labeled.length;
        int num = 0;
        for(int i=0; i<npillars; i++){               
            if(is_labeled[i] && ix[i]==0 && iy[i]==0) num++;
        }
        return num;
    }
    
    public void create_grid_points(double cacth_radius, double threshold_dxy){
        if(!is_grid_labelled) return;
        //if(!is_grid_created) return;
        int npillars = is_labeled.length;
        boolean[] active = new boolean[npillars];
        for(int i=0;i<npillars; i++) active[i] = is_labeled[i];
        int num_zeros = get_num_points_zero();
        if(num_zeros != 1){
            for(int i=0; i<npillars; i++){
                if(is_labeled[i] && ix[i]==0 && iy[i]==0){
                    ix[i]=Integer.MAX_VALUE;
                    iy[i]=Integer.MAX_VALUE;
                }
            }
        }
        double threshold2 = threshold_dxy*threshold_dxy;
        for(int i=0;i<npillars; i++){            
            if(active_d[i]){
                double ddx = dx[i];
                double ddy = dy[i];
                double dis = ddx*ddx + ddy*ddy;
                if(dis>threshold2){ 
                    active[i] = false;
                    boolean removing = false;
                    int num_x=get_num_points_line(active, ix, ix[i]);
                    if(num_x>5){
                        int num_y=get_num_points_line(active, iy, iy[i]);
                        if(num_y>5) removing = true;
                    }
                    active[i] = (!removing);
                }
            }
        }
        
        create_grid_points(cacth_radius, active);
    }
    
    public void create_grid_points(double cacth_radius, boolean[] active){
        is_grid_created = false;
        if(!is_grid_labelled) return;
        
        //int npillars = is_labeled.length;
        int ixmin = find_min(ix, is_labeled);        
        int iymin = find_min(iy, is_labeled);
        if(ixmin==Integer.MAX_VALUE || iymin==Integer.MAX_VALUE)  return;
        ix_off = ixmin;
        iy_off = iymin;
        
        int ixmax = find_max(ix, is_labeled);        
        int iymax = find_max(iy, is_labeled);
        if(ixmax==Integer.MIN_VALUE || iymax==Integer.MIN_VALUE)  return;        
    
        int xdim = ixmax - ixmin + 1;
        int ydim = iymax - iymin + 1;
        double[] minx = new double[xdim];
        double[] maxx = new double[xdim];                
        for(int ixx=ixmin; ixx<=ixmax; ixx++){
            List<Integer> list = getXLines(ixx, null);
            int index = ixx - ixmin;            
            double[][] xy = getXYPoints(list);
            if(xy!=null){
                minx[index] = min(xy[0]);
                maxx[index] = max(xy[0]);            
            }
        }
        double[] miny = new double[ydim];
        double[] maxy = new double[ydim];        
        for(int iyy=iymin; iyy<=iymax; iyy++){
            List<Integer> list = getYLines(iyy, null);
            int index = iyy - iymin;           
            double[][] xy = getXYPoints(list);
            if(xy!=null){
                miny[index] = min(xy[0]);
                maxy[index] = max(xy[0]);          
            }
        }
        
        double[][] px = new double[xdim][];
        int[] suc_px = new int[xdim];
        for(int ixx=ixmin; ixx<=ixmax; ixx++){
            List<Integer> list = getXLines(ixx, active);
            int index = ixx - ixmin;
            if(list.size()>3){
                double[][] xy = getXYPoints(list);
                PloyFitter fitter = new PloyFitter(xy[1], xy[0]);
                fitter.doFit();                
                px[index] = fitter.getParameters();
                suc_px[index] = fitter.getStatus();
                //IJ.log(fitter.getResultString());
                //minx[index] = min(xy[0]);
                //maxx[index] = max(xy[0]);
            }
            else suc_px[index] = -1;
        }
        
        double[][] py = new double[ydim][];
        int[] suc_py = new int[ydim];
        for(int iyy = iymin; iyy<=iymax; iyy++){
            List<Integer> list = getYLines(iyy, active);
            int index = iyy - iymin;
            if(list.size()>3){
                double[][] xy = getXYPoints(list);
                PloyFitter fitter = new PloyFitter(xy[0], xy[1]);
                fitter.doFit();                
                py[index] = fitter.getParameters();
                suc_py[index] = fitter.getStatus();
                //IJ.log(fitter.getResultString());
                //miny[index] = min(xy[0]);
                //maxy[index] = max(xy[0]);
            }
            else suc_py[index] = -1;
        }
        
        computeDeflection(px, suc_px, minx, maxx, py, suc_py, miny, maxy, ixmin, iymin, cacth_radius);
        
        is_grid_created = true;
    }
    
    private static double[] zero_solver(double[] px, double[] py, double start, double end, int map_width, int map_height){
        double[] gxy = zero_solver(px, py, start, end);
        if(gxy==null) return null;
        double gx = gxy[0];
        double gy = gxy[1];
        if(gx>=0 && gy>=0 && gx<map_width && gy<map_height) return gxy;
        return null;
    }
    
    private static double[] zero_solver(double[] px, double[] py, double start, double end){
        if(end<start) return null;
        
        double fys = PloyFitter.f(py, start);
        double fye = PloyFitter.f(py, end);
        double fxs = PloyFitter.f(px, fys);
        double fxe = PloyFitter.f(px, fye);
        double eval_s = fxs - start;
        double eval_e = fxe - end;
        
        if(eval_s*eval_e>0) return null;
        
        double mid = (start + end)/2.0;
        double fym = PloyFitter.f(py,mid);
        double fxm = PloyFitter.f(px,fym);
        double eval_m = fxm-mid;   
        
        double step=(end-start)/2.0; 
        double min_step = 1.0e-10;
        int k = 0;
        int k_max = 1000;
        while(step>min_step && eval_m!=0 && k<k_max){
            if(eval_m*eval_s>0){
                step = mid-start;
                start = mid;
                eval_s = eval_m;
            }
            else{
                step = end-mid;
                end = mid;
            }
            
            mid = (start + end)/2.0;
            fym = PloyFitter.f(py,mid);
            fxm = PloyFitter.f(px,fym);
            eval_m = fxm-mid; 
            k++;
        } 
        if(Math.abs(eval_m)>1.0e-5) return null;
        double gx = mid;
        double gy = fym;  
        /*
        double step = end - start;
        double min_step = 1.0e-6;
        while(step>min_step){
            start = fxe;
            end = fxs;
            step = end - start;
            fys = PloyFitter.f(py, start);
            fye = PloyFitter.f(py, end);
            fxs = PloyFitter.f(px, fys);
            fxe = PloyFitter.f(px, fye);
        } 
        
        if(Math.abs(fys-fye)>1.0e-3) return null;
        
        double gx = (start+end)/2.0;
        double gy = (fys+fye)/2.0;
        */
        double[] gxy = {gx, gy};
        return gxy;
    }
    
    private static EvaluateRangeResult evaluate_range(double scan_start, double scan_end, double[] pxi, double[] pyj){
        double fys = PloyFitter.f(pyj,scan_start);                        
        double fye = PloyFitter.f(pyj,scan_end);  

        double fxs = PloyFitter.f(pxi,fys);
        double fxe = PloyFitter.f(pxi,fye);              

        double eval_s = fxs-scan_start;
        double eval_e = fxe-scan_end;
        boolean suc = (eval_s*eval_e<0);
        
        EvaluateRangeResult rst = new EvaluateRangeResult(fxs, fxe, suc);
        return rst;
    }    
    
    private static IntersectionRangeResult check_intersection_range(double xxmin,double xxmax,double yymin,double yymax,double[] pxi,double[] pyj) {           
        double scan_start = Math.max(xxmin, yymin);
        double scan_end = Math.min(xxmax, yymax);                  
        EvaluateRangeResult fevl = evaluate_range(scan_start, scan_end, pxi, pyj);
        boolean suc = fevl.Sucessful();
        if(!suc){
            double xxlen = xxmax-xxmin;
            double yylen = yymax-yymin;            
            if(xxlen<yylen){
                scan_start = xxmin;
                scan_end = xxmax;
                fevl = evaluate_range(scan_start, scan_end, pxi, pyj);
                suc = fevl.Sucessful();
                if(!suc){
                    scan_start = yymin;
                    scan_end = yymax;
                    fevl = evaluate_range(scan_start, scan_end, pxi, pyj);            
                    suc = fevl.Sucessful();
                }            
            }
            else{
                scan_start = yymin;
                scan_end = yymax;
                fevl = evaluate_range(scan_start, scan_end, pxi, pyj);         
                suc = fevl.Sucessful();
                if(!suc){
                    scan_start = xxmin;
                    scan_end = xxmax;
                    fevl = evaluate_range(scan_start, scan_end, pxi, pyj);      
                    suc = fevl.Sucessful();
                }
            }

            if(!suc){   
                scan_start = Math.min(xxmin, yymin);
                scan_end = Math.max(xxmax, yymax);                    
                fevl = evaluate_range(scan_start, scan_end, pxi, pyj);   
                suc = fevl.Sucessful();
            }
        }
        
        IntersectionRangeResult rst = new IntersectionRangeResult(scan_start, scan_end, fevl.get_fevl_start(), fevl.get_fevl_end(), suc);
        return rst;
    }
    
    public void computeDeflection(double[][] px, int[] suc_px, double[] minx, double[] maxx, double[][] py, int[] suc_py, double[] miny, double[] maxy, int ix_min, int iy_min, double radius){
        int xdim = px.length;
        int ydim = py.length;
        
        double[][] grid_x = new double[ydim][xdim];
        double[][] grid_y = new double[ydim][xdim];
        boolean[][] found = new boolean[ydim][xdim];
        
        for(int i=0; i<xdim; i++){ 
                for(int j=0; j<ydim; j++){
                    grid_x[j][i] = Double.NaN;
                    grid_y[j][i] = Double.NaN;
                    found[j][i] = false;
                }
        }
        
        int num_found = 0;
        double margin = radius;
        for(int i=0; i<xdim; i++){ // %for i=20:20 %
            if(suc_px[i]==0){
                double[] pxi=px[i];
                double xxmin = minx[i] - margin;
                double xxmax = maxx[i] + margin;
                        
                for(int j=0; j<ydim; j++){ //%for j=40:40 %
                    if(suc_py[j]==0){
                        double[] pyj = py[j];
                        double yymin = miny[j] - margin;
                        double yymax = maxy[j] + margin;
                        
//                        double scan_start = Math.max(xxmin, yymin);
//                        double scan_end = Math.min(xxmax, yymax);
                        IntersectionRangeResult range = check_intersection_range(xxmin, xxmax, yymin, yymax,pxi,pyj);
                        boolean suc = range.Sucessful();
                        if(suc){
                            double scan_start = range.get_start();
                            double scan_end = range.get_end();
                            double[] gxy = zero_solver(pxi, pyj, scan_start, scan_end,map_width,map_height);
                            if(gxy!=null){
                                grid_x[j][i] = gxy[0];
                                grid_y[j][i] = gxy[1];
                                found[j][i] = true;
                                num_found++;
                            }
                        }
                    }
                }
            }
        }
        double[] gx = new double[num_found];
        double[] gy = new double[num_found];
        int[][] grid_map = new int[map_height][map_width];
        for(int i=0; i<map_width; i++) for(int j=0; j<map_height; j++) grid_map[j][i] = -1;
        int[][] grid_ij = new int[ydim][xdim];
        int num = 0;
        for(int i=0; i<xdim; i++){ 
            for(int j=0; j<ydim; j++){
                grid_ij[j][i] = -1;
                if(found[j][i]){                    
                    grid_ij[j][i] = num;
                    gx[num] = grid_x[j][i];
                    gy[num] = grid_y[j][i];            
                    int igx = (int)Math.round(gx[num]);
                    int igy = (int)Math.round(gy[num]);                    
                    if(igx>=0 && igy>=0 && igx<map_width && igy<map_height)
                        grid_map[igy][igx] = num;
                    
                    num = num + 1;
                }
            }
        }

        int npillars = is_labeled.length;
        dx = new double[npillars];
        dy = new double[npillars];
        active_d = new boolean[npillars];        
        boolean[] valid_grid = new boolean[num];
        for(int k=0; k<num; k++) valid_grid[k] = false;
        
        int num_active = 0;
        int num_labeld = 0;
        for(int k=0; k<npillars; k++){
            active_d[k] = false;   
            if(is_labeled[k]){
                num_labeld++;
                if(ix[k]==0 && iy[k]==0){}
                else{
                    int i = ix[k] - ix_min;
                    int j = iy[k] - iy_min;                
                    if(i>=0 && j>=0 && found[j][i]){
                        dx[k] = xpoints[k] - grid_x[j][i];
                        dy[k] = ypoints[k] - grid_y[j][i];
                        active_d[k] = true;   
                        int p = grid_ij[j][i];
                        valid_grid[p] = true;
                        num_active++;
                    }
                }
            }
        }
        
        if(num_active<npillars){            
            int r = Math.max(1,(int)Math.ceil(radius));    
            double max_dis = radius*radius;
            for(int k=0; k<npillars; k++){
                if(!active_d[k]){
                    int x = (int)Math.round(xpoints[k]);
                    int y = (int)Math.round(ypoints[k]); 
                    //% search the neighbours
                    double min_dis = Double.MAX_VALUE;//radius*radius*2;
                    int min_p = -1;    
                    for(int i=-r; i<=r; i++){
                        for(int j=-r; j<=r; j++){
                            int xx = i + x;
                            int yy = j + y;
                            if(xx>=0 && yy>=0 && xx<map_width && yy<map_height){
                                int p = grid_map[yy][xx];
                                if(p>=0 && !valid_grid[p]) {
                                    double ddx = xpoints[k] - gx[p];
                                    double ddy = ypoints[k] - gy[p];
                                    double dis = ddx*ddx + ddy*ddy;
                                    //int dis = i*i + j*j;
                                    if(max_dis>dis){                                
                                        if(min_dis>dis){
                                            min_dis = dis;
                                            min_p = p;                        
                                        }
                                    }
                                }
                            }
                        }
                    }

                    if (min_p>=0){
                        dx[k] = xpoints[k] - gx[min_p];
                        dy[k] = ypoints[k] - gy[min_p];
                        active_d[k] = true;
                        valid_grid[min_p] = true;
                    }
                    else{
                        dx[k] = Double.NaN;
                        dy[k] = Double.NaN;
                        active_d[k] = false;
                    }
                }
            }
            
        }
    }
    
    public void computeDeflection(double[][] px, int[] suc_px, double[] minx, double[] maxx, double[][] py, int[] suc_py, double[] miny, double[] maxy, double radius){
        int xdim = px.length;
        int ydim = py.length;
        
        double[][] grid_x = new double[ydim][xdim];
        double[][]  grid_y = new double[ydim][xdim];
        boolean[][] found = new boolean[ydim][xdim];
        
        for(int i=0; i<xdim; i++){ 
                for(int j=0; j<ydim; j++){
                    grid_x[j][i] = Double.NaN;
                    grid_y[j][i] = Double.NaN;
                    found[j][i] = false;
                }
        }
        
        int num_found = 0;
        double margin = radius;
        for(int i=0; i<xdim; i++){ // %for i=20:20 %
            if(suc_px[i]==0){
                double[] pxi=px[i];
                double xxmin = minx[i] - margin;
                double xxmax = maxx[i] + margin;
                        
                for(int j=0; j<ydim; j++){ //%for j=40:40 %
                    if(suc_py[j]==0){
                        double[] pyj = py[j];
                        double yymin = miny[j] - margin;
                        double yymax = maxy[j] + margin;
                        
                        double scan_start = Math.max(xxmin, yymin);
                        double scan_end = Math.min(xxmax, yymax);
                        
                        double[] gxy = zero_solver(pxi, pyj, scan_start, scan_end,map_width,map_height);
                        if(gxy!=null){
                            grid_x[j][i] = gxy[0];
                            grid_y[j][i] = gxy[1];
                            found[j][i] = true;
                            num_found++;
                        }
                    }
                }
            }
        }
        
        double[] gx = new double[num_found];
        double[] gy = new double[num_found];
        int[][] grid_map = new int[map_height][map_width];
        for(int i=0; i<map_width; i++) for(int j=0; j<map_height; j++) grid_map[j][i] = -1;
        
        int num = 0;
        for(int i=0; i<xdim; i++){ 
            for(int j=0; j<ydim; j++){
                if(found[j][i]){                    
                    gx[num] = grid_x[j][i];
                    gy[num] = grid_y[j][i];            
                    int igx = (int)Math.round(gx[num]);
                    int igy = (int)Math.round(gy[num]);                    
                    if(igx>=0 && igy>=0 && igx<map_width && igy<map_height)
                        grid_map[igy][igx] = num;
                    
                    num = num + 1;
                }
            }
        }

        int npillars = is_labeled.length;
        dx = new double[npillars];
        dy = new double[npillars];
        active_d = new boolean[npillars];        
        boolean[] valid_grid = new boolean[num];
        for(int k=0; k<num; k++) valid_grid[k] = false;
        
        int r = Math.max(1,(int)Math.ceil(radius));    
        double max_dis = radius*radius;
        for(int k=0; k<npillars; k++){
            int x = (int)Math.round(xpoints[k]);
            int y = (int)Math.round(ypoints[k]); 
            //% search the neighbours
            double min_dis = Double.MAX_VALUE;//radius*radius*2;
            int min_p = -1;    
            for(int i=-r; i<=r; i++){
                for(int j=-r; j<=r; j++){
                    int xx = i + x;
                    int yy = j + y;
                    if(xx>=0 && yy>=0 && xx<map_width && yy<map_height){
                        int p = grid_map[yy][xx];
                        if(p>=0 && !valid_grid[p]) {
                            double ddx = xpoints[k] - gx[p];
                            double ddy = ypoints[k] - gy[p];
                            double dis = ddx*ddx + ddy*ddy;
                            //int dis = i*i + j*j;
                            if(max_dis>dis){                                
                                if(min_dis>dis){
                                    min_dis = dis;
                                    min_p = p;                        
                                }
                            }
                        }
                    }
                }
            }
            
            if (min_p>=0){
                dx[k] = xpoints[k] - gx[min_p];
                dy[k] = ypoints[k] - gy[min_p];
                active_d[k] = true;
                valid_grid[min_p] = true;
            }
            else{
                dx[k] = Double.NaN;
                dy[k] = Double.NaN;
                active_d[k] = false;
            }
        }        
    }
    
    public static double[][] computeDeflections(double[] xc, double[] yc, double[] gx, double[] gy, int map_width, int map_height, double radius){
        int[][] grid_map = new int[map_height][map_width];
        for(int i=0; i<map_width; i++) for(int j=0; j<map_height; j++) grid_map[j][i] = -1;            
        int npoints_grid = gx.length;            
        for(int i=0; i<npoints_grid; i++){ 
            int igx = (int)Math.round(gx[i]);
            int igy = (int)Math.round(gy[i]);                    
            if(igx>=0 && igy>=0 && igx<map_width && igy<map_height)
                grid_map[igy][igx] = i;
        }

        int npillars = xc.length;
        double[][] dxy = new double[2][npillars];
        boolean[] valid_grid = new boolean[npoints_grid];
        for(int k=0; k<npoints_grid; k++) valid_grid[k] = false;

        int r = Math.max(1,(int)Math.ceil(radius));    
        double max_dis = radius*radius;
        for(int k=0; k<npillars; k++){
            int x = (int)Math.round(xc[k]);
            int y = (int)Math.round(yc[k]); 
            //% search the neighbours
            double min_dis = Double.MAX_VALUE;//radius*radius*2;
            int min_p = -1;    
            for(int i=-r; i<=r; i++){
                for(int j=-r; j<=r; j++){
                    int xx = i + x;
                    int yy = j + y;
                    if(xx>=0 && yy>=0 && xx<map_width && yy<map_height){
                        int p = grid_map[yy][xx];
                        if(p>=0 && !valid_grid[p]) {
                            double ddx = xc[k] - gx[p];
                            double ddy = yc[k] - gy[p];
                            double dis = ddx*ddx + ddy*ddy;
                            //int dis = i*i + j*j;
                            if(max_dis>dis){                                
                                if(min_dis>dis){
                                    min_dis = dis;
                                    min_p = p;                        
                                }
                            }
                        }
                    }
                }
            }

            if (min_p>=0){
                dxy[0][k] = xc[k] - gx[min_p];
                dxy[1][k] = yc[k] - gy[min_p];
                //active_d[k] = true;
                valid_grid[min_p] = true;
            }
            else{
                dxy[0][k] = Double.NaN;
                dxy[1][k] = Double.NaN;
                //active_d[k] = false;
            }
        }
        return dxy;
    }
}
