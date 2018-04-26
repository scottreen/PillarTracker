package com.nus.mbi.pillar.drift;

import com.nus.mbi.pillar.stat.MyPoint;

/**
 * vector projection: https://en.wikipedia.org/wiki/Vector_projection
 * @author xiaochun
 */
public class CU_Projection {
    public double proj_p1;
    public double proj_s1;
    public double proj_p2;
    public double proj_s2;
    public double cosa1;
    public double cosa2;
    public boolean isNaN = true;

    public CU_Projection() {
        
    }
    
    public CU_Projection(double proj_p1, double proj_s1, double proj_p2, double proj_s2, double cosa1, double cosa2) {
        this.proj_p1 = proj_p1;
        this.proj_s1 = proj_s1;
        this.proj_p2 = proj_p2;
        this.proj_s2 = proj_s2;
        this.cosa1 = cosa1;
        this.cosa2 = cosa2;
        isNaN = Double.isNaN(proj_p1) || Double.isNaN(proj_p2);
    }
    
    public static CU_Projection create(MyPoint g1, MyPoint g2, double dis, MyPoint c1, MyPoint c2){        
        double gx1 = g1.x;
        double gy1 = g1.y;
        double gx2 = g2.x;
        double gy2 = g2.y;
        double xc1 = c1.x;
        double yc1 = c1.y;
        double xc2 = c2.x;
        double yc2 = c2.y;
        
        double bx = gx2-gx1;
        double by = gy2-gy1;
        double b2 = dis;
        double ax1 = xc1-gx1;
        double ay1 = yc1-gy1;        
        
        double a1 = (ax1*bx + ay1*by)/b2;          
        double ab1 = a1/b2;
        double cx1 = ax1 - ab1*bx; 
        double cy1 = ay1 - ab1*by;
        double a2 = Math.sqrt(cx1*cx1+cy1*cy1);     
        double sign = cy1*bx - cx1*by;
        if(sign<0) a2 = -a2;
        double dis1 = ax1*ax1+ay1*ay1;
        double cos1 = dis1>0 ? a1/Math.sqrt(dis1) : 0;

        double ax2 = xc2-gx2;
        double ay2 = yc2-gy2;
        double aa1 = (ax2*bx + ay2*by)/b2;                  
        double ab2 = aa1/b2;
        double cx2 = ax2 - ab2*bx; 
        double cy2 = ay2 - ab2*by;
        double aa2 = Math.sqrt(cx2*cx2+cy2*cy2);     
        sign = cy2*bx - cx2*by;
        if(sign<0) aa2 = -aa2;
        double dis2 = ax2*ax2+ay2*ay2;
        double cos2 = dis2>0 ? aa1/Math.sqrt(dis2) : 0; 
        return new CU_Projection(a1,a2,aa1,aa2,cos1,cos2); 
    }
    
    public static CU_Projection create(MyPoint g, double dis, MyPoint d1, MyPoint d2){
        double bx = g.x;
        double by = g.y;
        double b2 = dis;
        double ax1 = d1.x;
        double ay1 = d1.y;     
        
        double a1 = (ax1*bx + ay1*by)/b2;          
        double ab1 = a1/b2;
        double cx1 = ax1 - ab1*bx; 
        double cy1 = ay1 - ab1*by;
        double a2 = Math.sqrt(cx1*cx1+cy1*cy1);     
        double sign = cy1*bx - cx1*by;
        if(sign<0) a2 = -a2;
        double dis1 = ax1*ax1+ay1*ay1;
        double cos1 = dis1>0 ? a1/Math.sqrt(dis1) : 0;

        double ax2 = d2.x;
        double ay2 = d2.y;
        double aa1 = (ax2*bx + ay2*by)/b2;                  
        double ab2 = aa1/b2;
        double cx2 = ax2 - ab2*bx; 
        double cy2 = ay2 - ab2*by;
        double aa2 = Math.sqrt(cx2*cx2+cy2*cy2);     
        sign = cy2*bx - cx2*by;
        if(sign<0) aa2 = -aa2;
        double dis2 = ax2*ax2+ay2*ay2;
        double cos2 = dis2>0 ? aa1/Math.sqrt(dis2) : 0; 
        return new CU_Projection(a1,a2,aa1,aa2,cos1,cos2);
    }
}
