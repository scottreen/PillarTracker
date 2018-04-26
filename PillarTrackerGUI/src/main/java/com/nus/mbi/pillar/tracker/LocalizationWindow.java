package com.nus.mbi.pillar.tracker;

/**
 *
 * @author xiaochun
 */
public class LocalizationWindow {
    public int xc; 
    public int yc;
    
    public double[] win;
    public double[] winx;
    public double[] winy;
    public int win_size;
    public int window_radius;

    public LocalizationWindow(int xc, int yc, double[] win, double[] winx, double[] winy, int win_size, int window_radius) {
        this.xc = xc;
        this.yc = yc;
        this.win = win;
        this.winx = winx;
        this.winy = winy;
        this.win_size = win_size;
        this.window_radius = window_radius;
    }
    
    public static LocalizationWindow create(int xc, int yc, double[] pixels, int width, int height, int kernel_width){						
        int window_r = kernel_width/2;
        int window_w = window_r*2 + 1;
        int window_s = window_w*window_w;
        double[] win = new double[window_s];
        double[] winx= new double[window_s];
        double[] winy= new double[window_s];  
        
        int k=0;
        for(int x=-window_r; x<=window_r; x++){
            for(int y=-window_r; y<=window_r; y++){
                int i=y+yc;
                int j=x+xc;
                if(i>=0 && i<height && j>=0 && j<width){					
                        winx[k] = x;
                        winy[k] = y;
                        double t = pixels[i*width+j];					
                        win[k] = t;
                        k++;
                }
            }
        }
        
        return new LocalizationWindow(xc, yc, win, winx, winy, k, window_r);
    }    
}
