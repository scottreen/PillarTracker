package com.nus.mbi.pillar.grid;

/**
 *
 * @author xiaochun
 */
public class SearchAreaCreation {
    private int radius;
    private int outer_radius;
    private double[][] nxy = null;  
    private double[][] area = null;
    
    public int getRadius(){
        return radius;
    }    
    public int getOuterRadius(){
        return outer_radius;
    }
    
    public double[][] getArea(){
        return area;
    }
    
    public double[] getNX(){
        return nxy[0];
    }
    
    public double[] getNY(){
        return nxy[1];
    }
    
    public double[][] create_search_area(double seperation, double grid_oblique, double grid_angle){
        double s = seperation;
        double ax = -grid_oblique*Math.PI/180.0;
        double xdx = s*Math.cos(ax);
        double xdy = s*Math.sin(ax);

        double ay = -(grid_angle+grid_oblique)*Math.PI/180.0;
        double ydx = s*Math.cos(ay);
        double ydy = s*Math.sin(ay);    

        double[] nx = {xdx,ydx,-xdx,-ydx};
        double[] ny = {xdy,ydy,-xdy,-ydy};
        
        nxy = new double[2][4];
        for(int k=0; k<4; k++){
            nxy[0][k] = nx[k];
            nxy[1][k] = ny[k];
        }
        double tol = 10; // use the tolerance to search the neighbors;    
        double r = Math.abs(s*Math.sin(tol*Math.PI/180)); 
        radius = (int) Math.ceil(r);
        double r2 = radius*radius;
        outer_radius = (int)Math.ceil(s+r);
        
        int p = 0;
        for(int i=-radius; i<=radius; i++){   
            int i2 = i*i;
            for(int j=-radius; j<=radius; j++){                        
                int j2 = j*j; 
                int d2 = i2 + j2;
                if(d2<=r2) p++;
            }
        }

        int area_size = p;
        area = new double[area_size][3];

        p = 0;
        for(int i=-radius; i<=radius; i++){   
            int i2 = i*i;
            for(int j=-radius; j<=radius; j++){                        
                int j2 = j*j; 
                int d2 = i2 + j2;
                if(d2<=r2){                                                   
                    area[p][0] = d2;
                    area[p][1] = i;
                    area[p][2] = j;
                    p++;      
                }
            }
        }     
        return area;
    }
//    
//    public static double[][] create_search_area(double s, double tol){
//        //double radius = (int) Math.ceil(r);
//        double r = Math.abs(s);
//        double t = Math.abs(tol);
//        int radius = (int)Math.ceil(r+t);        
//        double r2 = radius*radius;
//        int in_radius = (int)Math.floor(r-t);
//        double ir2 = in_radius*in_radius;
//        
//        int p = 0;
//        for(int i=-radius; i<=radius; i++){   
//            int i2 = i*i;
//            for(int j=-radius; j<=radius; j++){                        
//                int j2 = j*j; 
//                int d2 = i2 + j2;
//                if(d2<=r2 && d2>=ir2) p++;
//            }
//        }
//
//        int area_size = p;
//        double[][] area = new double[area_size][3];
//
//        p = 0;
//        for(int i=-radius; i<=radius; i++){   
//            int i2 = i*i;
//            for(int j=-radius; j<=radius; j++){                        
//                int j2 = j*j; 
//                int d2 = i2 + j2;
//                if(d2<=r2 && d2>=ir2){                                                   
//                    area[p][0] = d2;
//                    area[p][1] = i;
//                    area[p][2] = j;
//                    p++;      
//                }
//            }
//        }     
//        return area;
//    }
}
