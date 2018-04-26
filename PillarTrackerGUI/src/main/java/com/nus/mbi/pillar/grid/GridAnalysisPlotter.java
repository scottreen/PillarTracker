package com.nus.mbi.pillar.grid;

import ij.ImagePlus;
import ij.gui.Arrow;
import ij.gui.Line;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.PointRoi;
import ij.process.FloatPolygon;
import java.awt.Color;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.List;
import static com.nus.mbi.pillar.stat.BasicStatisitic.*;

/**
 *
 * @author xiaochun
 */
public class GridAnalysisPlotter {
    
    public int arrow_head_size = 5;
    public Color arrow_color = Color.cyan;
    public int cross_size = 3;
    public Color cross_color = Color.white;
    
    public Overlay getPoints(double[] cx, double[] cy, double zoom){
            int nrois = cx.length;            
            
            FloatPolygon poly = new FloatPolygon();
            for(int i=0; i<nrois; i++){
                double rcx = (cx[i]+0.5)*zoom;
                double rcy = (cy[i]+0.5)*zoom;
                poly.addPoint(rcx, rcy);
            }
            
            PointRoi roi = new PointRoi(poly);
            Overlay overlay = new Overlay(roi);  
            return overlay;
    }
    
    public Overlay getPoints(double[] cx, double[] cy){
        return getPoints(cx, cy, 1);
    }
        
    public Overlay getArrows(double[] x0, double[] y0, double[] x1, double[] y1, double zoom_image, double zoom_arrow){
        int nrois = x0.length;        
        Overlay overlay = new Overlay();            
       
        for(int i=0; i<nrois; i++){
            double xx0 = (x0[i]+0.5)*zoom_image;
            double yy0 = (y0[i]+0.5)*zoom_image;
            double xx1 = (x1[i]+0.5)*zoom_image;
            double yy1 = (y1[i]+0.5)*zoom_image;
            double ax1 = (xx1-xx0)*zoom_arrow + xx0;
            double ay1 = (yy1-yy0)*zoom_arrow + yy0;            
            Arrow roi = new Arrow(xx0, yy0, ax1, ay1);      
            overlay.add(roi);              
        }
        //labels_overlay = overlay;
        return overlay;
    }
    
    public Overlay getArrows(double[] x0, double[] y0, double[] x1, double[] y1, double zoom_arrow){
        return getArrows(x0, y0, x1, y1, 1, zoom_arrow);
    }
    
    public Overlay getArrows(double[] x0, double[] y0, double[] x1, double[] y1){
        return getArrows(x0, y0, x1, y1, 1, 1);
    }
    
    public Overlay getPointsArrows(double[] x0, double[] y0, double[] x1, double[] y1, double zoom_image, double zoom_arrow){
        int nrois = x0.length;        
        Overlay overlay = new Overlay();                    
        double arrow_size = arrow_head_size;
        double t = arrow_size*arrow_size;    
        double size = cross_size;
        for(int i=0; i<nrois; i++){
            double xx0 = (x0[i]+0.5)*zoom_image;
            double yy0 = (y0[i]+0.5)*zoom_image;
            double xx1 = (x1[i]+0.5)*zoom_image;
            double yy1 = (y1[i]+0.5)*zoom_image;
            
            double ax1 = (xx1-xx0)*zoom_arrow;
            double ay1 = (yy1-yy0)*zoom_arrow;
            if(ax1*ax1 + ay1*ay1 > t){
                Arrow arrow = new Arrow(xx0, yy0, ax1+xx0, ay1+yy0);  
                arrow.setHeadSize(arrow_head_size); 
                arrow.setStrokeWidth(1);
                //arrow.setStyle(Arrow.OPEN);
                arrow.setFillColor(arrow_color);
                arrow.setStyle("OPEN");
                overlay.add(arrow);
            }
            
            //PointRoi point = new PointRoi(xx0, yy0);
            //point.setPointType(2);
            //point.setSize(1);
            //overlay.add(point);     
            if(size>0){
                Line line1 = new Line(xx0-size,yy0-size,xx0+size,yy0+size);
                Line line2 = new Line(xx0-size,yy0+size,xx0+size,yy0-size);
                line1.setStrokeColor(cross_color);//setFillColor(Color.white);
                line2.setStrokeColor(cross_color);//line2.setFillColor(Color.white);
                overlay.add(line1);
                overlay.add(line2);
            }
        }
        //labels_overlay = overlay;
        return overlay;
    }
    
    public Overlay getPointsArrows(double[] x0, double[] y0, double[] x1, double[] y1, boolean[]flags, double zoom_image, double zoom_arrow){
        int nrois = x0.length;        
        Overlay overlay = new Overlay();                    
        double arrow_size = arrow_head_size;
        double t = arrow_size*arrow_size;    
        double size = cross_size;
        for(int i=0; i<nrois; i++){
            double xx0 = (x0[i]+0.5)*zoom_image;
            double yy0 = (y0[i]+0.5)*zoom_image;
            double xx1 = (x1[i]+0.5)*zoom_image;
            double yy1 = (y1[i]+0.5)*zoom_image;
            
            double ax1 = (xx1-xx0)*zoom_arrow;
            double ay1 = (yy1-yy0)*zoom_arrow;
            Color arrow_clr = flags[i] ? Color.RED : arrow_color;
            if(ax1*ax1 + ay1*ay1 > t){
                Arrow arrow = new Arrow(xx0, yy0, ax1+xx0, ay1+yy0);  
                arrow.setHeadSize(arrow_head_size); 
                arrow.setStrokeWidth(1);
                //arrow.setStyle(Arrow.OPEN);
                arrow.setFillColor(arrow_clr);
                arrow.setStyle("OPEN");
                overlay.add(arrow);
            }
            
            //PointRoi point = new PointRoi(xx0, yy0);
            //point.setPointType(2);
            //point.setSize(1);
            //overlay.add(point);
            if(size>0){        
                    Line line1 = new Line(xx0-size,yy0-size,xx0+size,yy0+size);
                    Line line2 = new Line(xx0-size,yy0+size,xx0+size,yy0-size);
                    line1.setStrokeColor(cross_color);//setFillColor(Color.white);
                    line2.setStrokeColor(cross_color);//line2.setFillColor(Color.white);
                    overlay.add(line1);
                    overlay.add(line2);
            }
        }
        //labels_overlay = overlay;
        return overlay;
    }
    
    public Overlay getPointsArrows(double[] x0, double[] y0, double[] x1, double[] y1, double zoom_arrow){
        return getPointsArrows(x0, y0, x1, y1, 1, zoom_arrow);
    }
    
    public Overlay getPointsArrows(double[] x0, double[] y0, double[] x1, double[] y1){
        return getPointsArrows(x0, y0, x1, y1, 1, 1);
    }
    
    public void plot_arrow_overlay(Overlay arrow, ImagePlus ip, Overlay labels){
            if(arrow!=null){                                            
                if(labels!=null){
                    Overlay overlay = arrow.duplicate();
                    int n = labels.size();
                    for(int i=0; i<n; i++) overlay.add(labels.get(i));

                    ip.setOverlay(overlay);
                }                
                else{
                    ip.setOverlay(arrow);
                }
                
                if(ip.getWindow()==null) ip.show();
            }
    }
    /*
    public Overlay getLabels(double[] x0, double[] y0,double[] dx, double[] dy, double zoom_image, boolean showdxy){
        int nrois = x0.length;       
        int font_size = 8;//(int)Math.ceil(5);
        Font font = new Font("Arial",Font.PLAIN, font_size);
        Overlay overlay = new Overlay();  
        for(int i=0; i<nrois; i++){
            double x = x0[i]*zoom_image;
            double y = y0[i]*zoom_image;
            if(Double.isNaN(x) || Double.isNaN(y)) continue;
            
            double dxx = dx[i];
            double dyy = dy[i];
            String text = showdxy ? String.format("(%.0f,%.0f)", dxx, dyy) : String.format("%d", (i+1));
            TextRoi label = new TextRoi(x, y, text, font);
            label.setStrokeColor(Color.yellow);
            overlay.add(label);                            
        }
        //labels_overlay = overlay;
        return overlay;
    }
    */
    public Plot get_plot_driftXY(double[] xaxis, double[] disX, double[] disY, String title)
    {
        float plotlimitmin = Math.min((float)min(disX), (float)min(disY));
        float plotlimitmax = Math.max((float)max(disX), (float)max(disY));
        //Plot plot = new Plot("ROI-" + (i+1),"slice", "distance(nm)", xaxis, dis[i]);
        Plot plot = new Plot(title,"frame #", "deflection(nm)");
        plot.setLimits(0,disX.length,plotlimitmin,plotlimitmax);

        plot.setColor(Color.RED);	
        plot.addPoints(xaxis, disX,Plot.LINE);			
        plot.draw();	

        plot.setColor(Color.BLUE);	
        plot.addPoints(xaxis, disY,Plot.LINE);		
        plot.draw();

        return plot;
    }
    
        
    public static void plotGrid(double x0, double y0, double grid_separtion, double grid_oblique, double grid_angle, ImagePlus ip){
        if(ip==null || !ip.isVisible()) return;
        Rectangle rect = new Rectangle(0,0,ip.getWidth(),ip.getHeight());         
        Overlay overlay = plotLines(x0,y0,grid_separtion,grid_oblique,grid_angle,rect,1.0);
        ip.setOverlay(overlay);
    }    

final static class GridIJ {
    private final double I;
    private final double J;

    public GridIJ(double i, double j) {
        this.I = i;
        this.J = j;
    }

    public double getI() {
        return I;
    }

    public double getJ() {
        return J;
    }
}

    
    public static GridIJ solve_Grid_IJ(double x0, double y0, double dx1, double dy1, double dx2, double dy2, double x, double y){
        double xx = x - x0;
        double yy = y - y0;
        
        double t = dx2*dy1 - dx1*dy2;
        double i = -(xx*dy2 - yy*dx2)/t;
        double j = (xx*dy1 - yy*dx1)/t;
        GridIJ ij = new GridIJ(i,j);
        return ij;      
    }
    
    public static Overlay plotLines(double x0, double y0, double grid_separtion, double grid_oblique, double grid_angle, Rectangle rect, double zoom){
        Overlay overlay = new Overlay(); 
        
        double a = -grid_oblique*Math.PI/180;
        double b = -(grid_oblique+grid_angle)*Math.PI/180;   
        double s = grid_separtion;
        double cosa = Math.cos(a);
        double sina = Math.sin(a);
        double cosb = Math.cos(b);
        double sinb = Math.sin(b);        
        double dx1 = s*cosa;
        double dy1 = s*sina;        
        double dx2 = s*cosb;
        double dy2 = s*sinb;       
       
        GridIJ[] ij = new GridIJ[4];
        ij[0] = solve_Grid_IJ(x0, y0, dx1, dy1, dx2, dy2, rect.x,            rect.y);
        ij[1] = solve_Grid_IJ(x0, y0, dx1, dy1, dx2, dy2, rect.x,            rect.y+rect.height);
        ij[2] = solve_Grid_IJ(x0, y0, dx1, dy1, dx2, dy2, rect.x+rect.width, rect.y);
        ij[3] = solve_Grid_IJ(x0, y0, dx1, dy1, dx2, dy2, rect.x+rect.width, rect.y+rect.height);
        
        double imin = ij[0].I;
        double jmin = ij[0].J;
        double imax = ij[0].I;
        double jmax = ij[0].J;
        
        for(int k=1; k<4; k++){
            double ii = ij[k].I;
            if(imin>ii) imin = ii;
            if(imax<ii) imax = ii;
            
            double jj = ij[k].J;
            if(jmin>jj) jmin = jj;
            if(jmax<jj) jmax = jj;
        }
                            
        int xoff = (int) Math.floor(imin);
        int yoff = (int) Math.floor(jmin);
        int xend = (int) Math.ceil(imax);
        int yend = (int) Math.ceil(jmax);
//        double dxy = Math.abs(s*Math.sin(grid_angle*Math.PI/180));//Math.max(Math.abs(dx), Math.abs(dy));        
//        int xdim = (int) Math.ceil(Math.abs(rect.width/dxy));
//        int ydim = (int) Math.ceil(Math.abs(rect.height/dxy));                
//        int xoff = (int) Math.floor((rect.x-x0)/dxy);
//        int yoff = (int) Math.floor((rect.y-y0)/dxy);
//        int xend = xoff + xdim;
//        int yend = yoff + ydim;
        
        for(int j=yoff; j<=yend; j++){
            double x1=x0 + xoff*dx1 + j*dx2;
            double y1=y0 + xoff*dy1 + j*dy2;
            double x2=x0 + xend*dx1 + j*dx2;
            double y2=y0 + xend*dy1 + j*dy2;
            Line line = new Line(x1*zoom, y1*zoom, x2*zoom, y2*zoom);
            overlay.add(line);            
        }
        
        for(int i=xoff; i<=xend; i++){
            double x1=x0 + i*dx1 + yoff*dx2;
            double y1=y0 + i*dy1 + yoff*dy2;
            double x2=x0 + i*dx1 + yend*dx2;
            double y2=y0 + i*dy1 + yend*dy2;
            Line line = new Line(x1*zoom, y1*zoom, x2*zoom, y2*zoom);            
            overlay.add(line);            
        }
        
        return overlay;
    }
    
    
    public static void plotGrid(double x0, double y0, double dx, double dy, ImagePlus ip){
        if(ip==null || !ip.isVisible()) return;
        Rectangle rect = new Rectangle(0,0,ip.getWidth(),ip.getHeight());
        Overlay overlay = plotGrid(x0,y0,dx,dy,rect,1.0);
        ip.setOverlay(overlay);
    }
    
    public static Overlay plotGrid(double x0, double y0, double dx, double dy, Rectangle rect, double zoom){
        Overlay overlay = new Overlay(); 
        
        double dxy = Math.max(Math.abs(dx), Math.abs(dy));
        int xdim = (int) Math.ceil(Math.abs(rect.width/dxy));
        int ydim = (int) Math.ceil(Math.abs(rect.height/dxy));
        int xoff = (int) Math.floor((rect.x-x0)/dxy);
        int yoff = (int) Math.floor((rect.y-y0)/dxy);
        
        int xend = xoff + xdim;
        int yend = yoff + ydim;
        for(int i=xoff; i<=xend; i++){
            double x1=x0+i*dx-yoff*dy;
            double y1=y0+yoff*dx+i*dy;
            double x2=x0+i*dx-yend*dy;
            double y2=y0+yend*dx+i*dy;
            Line line = new Line(x1*zoom, y1*zoom, x2*zoom, y2*zoom);
            //line.setStrokeColor(cross_color);
            overlay.add(line);            
        }
        
        for(int j=yoff; j<=yend; j++){
            double x1=x0+xoff*dx-j*dy;
            double y1=y0+j*dx+xoff*dy;
            double x2=x0+xend*dx-j*dy;
            double y2=y0+j*dx+xend*dy;
            Line line = new Line(x1*zoom, y1*zoom, x2*zoom, y2*zoom);
            //line.setStrokeColor(cross_color);
            overlay.add(line);            
        }
        
        return overlay;
    }
    
    public void addLines(List<GridLinePoints> points, Overlay overlay){
        int num = points.size();
        if(num>3){                
            points.sort((a, b) -> b.compareTo(a));
            for(int i=0; i<num-1; i++){
                GridLinePoints p1 = points.get(i);                    
                GridLinePoints p2 = points.get(i+1);                    

                Line line = new Line(p1.x, p1.y, p2.x, p2.y);
                line.setStrokeColor(cross_color);
                overlay.add(line);
            }
        }
    }
    
    public Overlay plotGrid(double[] cx, double[] cy, int[] ix, int[] iy, boolean[] is_labeled, double zoom){
        Overlay overlay = new Overlay();   
        
        int ixmin = find_min(ix, is_labeled);
        int ixmax = find_max(ix, is_labeled);
        int iymin = find_min(iy, is_labeled);
        int iymax = find_max(iy, is_labeled);

//        int xdim = ixmax - ixmin + 1;
//        int ydim = iymax - iymin + 1;
        
        int npillars = is_labeled.length;
        //%% polyfit along x axis
        for(int ixx=ixmin; ixx<=ixmax; ixx++){            
            List<GridLinePoints> points = new ArrayList();
            for(int i=0; i<npillars; i++){               
                if(is_labeled[i] && ix[i]==ixx) {
                    double x = (cx[i]+0.5)*zoom;
                    double y = (cy[i]+0.5)*zoom;
                    points.add(new GridLinePoints(x, y, ix[i], iy[i]));
                }                                   
            }
            addLines(points, overlay);   
        }

        //%% polyfit along y axis
        for(int iyy=iymin; iyy<=iymax; iyy++){//for iyy=iymin:iymax
            List<GridLinePoints> points = new ArrayList();
            for(int i=0; i<npillars; i++){               
                if(is_labeled[i] && iy[i]==iyy){
                    double x = (cx[i]+0.5)*zoom;
                    double y = (cy[i]+0.5)*zoom;
                    points.add(new GridLinePoints(x, y, ix[i], iy[i]));
                }
            }
            addLines(points, overlay);   
//            int num = points.size();
//            if(num>3){               
//                points.sort((a, b) -> b.compareTo(a));
//                for(int i=0; i<num-1; i++){
//                    GridLinePoints p1 = points.get(i);                    
//                    GridLinePoints p2 = points.get(i+1);                    
//                    
//                    Line line = new Line(p1.x, p1.y, p2.x, p2.y);
//                    line.setStrokeColor(cross_color);
//                    overlay.add(line);
//                }
//            }
        }
        
        return overlay;
    }
}
