package com.nus.mbi.pillar.solver;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 *
 * @author xiaochun
 */
public class MinimumErrorSolver {
    private double[] psf;
    private int[] psf_dx;
    private int[] psf_dy;
    //private int psf_w;
    //private int psf_h;
    private int psf_size;
    private int psf_halfw;
    private int psf_halfh;
    private double sigmaY;
    private double sigmaY2;
    private double sigmaY2N_sigmaYY;
    
    public void setPSF(double[] psf, int[] psf_dx, int[] psf_dy, int psf_halfw, int psf_halfh, int psf_size){
        this.psf = psf;
        this.psf_dx = psf_dx;
        this.psf_dy = psf_dy;
        this.psf_size = psf_size;
        this.psf_halfw = psf_halfw;
        this.psf_halfh = psf_halfh;
        sigmaY = 0.0;
        sigmaY2 = 0.0;

        for(int i=0; i<psf_size; i++) 
        {
            sigmaY += psf[i];
            sigmaY2 += psf[i] * psf[i];
        }
        
        sigmaY2N_sigmaYY = sigmaY2*psf_size-sigmaY*sigmaY;
    }
    
    public void setPSF(ImageProcessor ip_psf){        
        int psf_w = ip_psf.getWidth();
        int psf_h = ip_psf.getHeight();   
        //if(psf_w!=psf_h) return;
        
        int k_size = psf_w*psf_h;
        double[] kernel = new double[k_size];
        int[] kernel_dx = new int[k_size];
        int[] kernel_dy = new int[k_size];        
        int k = 0;
        int half_w = psf_w/2;
        int half_h = psf_h/2;
        for(int y=-half_h; y<=half_h; y++){
            int j = y + half_h;
            for(int x=-half_w; x<=half_w; x++){
                int i = x + half_w;            
                double g = ip_psf.getf(i, j);  
                kernel[k] = g;
                kernel_dx[k] = x;
                kernel_dy[k] = y;                                
                k++;
            }
        } 
        setPSF(kernel, kernel_dx, kernel_dy, psf_w/2, psf_h/2,k_size);
    }
    
    public ImagePlus getPSF(){
        if(psf_halfw<1 || psf_halfh<1) return null;
        if(psf==null) return null;
        if(psf.length<1) return null;
        
        int k_w = psf_halfw*2 + 1;
        int k_h = psf_halfh*2 + 1;
        FloatProcessor ip_psf = new FloatProcessor(k_w, k_h);
        for(int k=0; k<psf_size; k++){
            int x = psf_dx[k] + psf_halfw;
            int y = psf_dy[k] + psf_halfh;
            float p = (float)psf[k];
            ip_psf.setf(x, y, (float)p);
        }            
        
        ImagePlus ip = new ImagePlus("PSF_AMP", ip_psf);
        return ip;
    }
    
    public ImagePlus setGaussianKernel(int radius, double sigma, boolean dark){        
        int k_w = 2*radius + 1;
        int k_size = k_w*k_w;
        FloatProcessor ip_psf = new FloatProcessor(k_w, k_w);
        double[] kernel = new double[k_size];
        int[] kernel_dx = new int[k_size];
        int[] kernel_dy = new int[k_size];        
        double sigma2 = sigma*sigma*2;
        int k = 0;
        for(int y=-radius; y<=radius; y++){
            int j=y+radius;
            for(int x=-radius; x<=radius; x++){
                int i=x+radius;            
                double e = Math.exp(-(x*x+y*y)/sigma2);                
                double g = dark ? -e : e;
                ip_psf.setf(i, j, (float)g);
                kernel[k] = g;  
                kernel_dx[k] = x;
                kernel_dy[k] = y;                                
                k++;
            }
        }  
        setPSF(kernel, kernel_dx, kernel_dy, radius, radius, k_size);
        ImagePlus ip = new ImagePlus("Gaussian_PSF", ip_psf);
        return ip;
    }
     
    public ImagePlus MinimumError_AMP(ImagePlus image) {
        return MinimumError_AMP(image, 1);
    }
    
    public ImagePlus MinimumError_AMP(ImagePlus image, int num_threads) {
        int n = image.getStackSize();
        int w = image.getWidth();
        int h = image.getHeight();
        ImageStack is = new ImageStack(w,h);       
        
        // slice numbers start with 1 for historical reasons
        for (int i = 1; i <= n; i++){
            ImageProcessor img = MinimumError_AMP(image.getStack().getProcessor(i), num_threads);
            is.addSlice(img);
            IJ.showProgress(i, n);
        }
        
        ImagePlus ip = new ImagePlus(image.getShortTitle()+"_AMP", is);
        return ip;                
    }
    
    public ImageProcessor MinimumError_AMP(ImageProcessor ip)
    {
        return MinimumError_AMP(ip, 1);
    }
    
    public float[] MinimumError_AMP(float[] pixels, int width, int height, int num_threads)
    {
        int size = width*height;                
        double[]image = new double[size];
        for(int k=0; k<size; k++) image[k] = pixels[k];
        
        double[] amp_img = num_threads<2 ? MinimumError_AMP(image, width, height) : MinimumError_AMP(image, width, height, num_threads);        
        
        float[] amp_pixels = new float[size];
        for(int k=0; k<size; k++) amp_pixels[k] = (float)amp_img[k];
        return amp_pixels;
    }
    
    public ImageProcessor MinimumError_AMP(ImageProcessor ip, int num_threads)
    {
        int width = ip.getWidth();
        int height = ip.getHeight();        
        float[][] pixels = ip.getFloatArray();        
        int size = width*height;                
        double[]image = new double[size];
        
        for(int y=0; y<height; y++){
            int off = y*width;
            for(int x=0; x<width; x++){
                int p = off + x;
                image[p] = pixels[x][y];
            }
        }
        //for(int k=0; k<size; k++) image[k] = ip.getf(k);
        
        double[] amp_img = num_threads<2 ? MinimumError_AMP(image, width, height) : MinimumError_AMP(image, width, height, num_threads);        
        FloatProcessor fp = new FloatProcessor(width, height);
        float[] amp_pixels = (float[])fp.getPixels();
        for(int k=0; k<size; k++){            
            amp_pixels[k] = (float)amp_img[k];
            //fp.setf(k, dark_object?-p:p);
        }
        
        return fp;
    }
    
    public double[] MinimumError_AMP(double[] imagef, int width, int height)
    {
        int NN = imagef.length;
        double[] amp_map = new double[NN];
        
        for(int y=psf_halfh; y<height-psf_halfh; y++){
            int off = y*width;        
            for(int x=psf_halfw; x<width-psf_halfw; x++){                    
                //int p= x + y*width;
                int p = x + off;
                double sigmaXY = 0.0;
                double sigmaX = 0.0;
                for(int i=0; i<psf_size; i++){
                        int x1 = x + psf_dx[i];
                        int y1 = y + psf_dy[i];

                        double mi = imagef[x1 + y1*width];
                        double dy = psf[i];

                        sigmaX  += mi;//sigmaX + mi;                                    
                        sigmaXY += dy*mi;//sigmaXY + //psf[i]*mi;
                }
                double lxy = psf_size*sigmaXY-sigmaX*sigmaY;		
                double amp = lxy/sigmaY2N_sigmaYY;
                amp_map[p] = amp;							
            }
        }

        // do the margins only
        int xmargin=psf_halfw;
        int ymargin=psf_halfh;
        int kernelw = psf_halfw*2 + 1;
        int kernelh = psf_halfh*2 + 1;
        int kernel_size=kernelw*kernelh;
        double[] mask=new double[kernel_size];
        for(int x=0;x<width;x++) {
            for(int y=0;y<height;y++) {
                if((x<xmargin) || (x>=width-xmargin) || (y<ymargin) || (y>=height-ymargin)) {
                    int i=0;
                    for(int y1=y-ymargin; y1<=y+ymargin; y1++) {
                        for(int x1=x-xmargin; x1<=x+xmargin; x1++) {
                            //if((x1>=0) && (x1<width) && (y1>=0) && (y1<height)) mask[i]=slice[x1+width*y1]; else mask[i]=0;
                            int x2=x1;
                            if(x1<0) x2=0; else if (x1>=width) x2=width-1;
                            int y2=y1;
                            if(y1<0) y2=0; else if (y1>=height) y2=height-1;
                            mask[i]=imagef[x2+width*y2];
                            i++;
                        }
                    }
                    double sigmaX=0.;
                    for(i=0; i<kernel_size; i++) sigmaX=sigmaX+mask[i];
//                                    double sigmaX2=0.;
//                                    for(i=0; i<kernelw2; i++) sigmaX2=sigmaX2+mask[i]*mask[i];
                    // calculate the variable solver parameters
                    double sigmaXY=0.;
                    for(i=0; i<kernel_size; i++) sigmaXY=sigmaXY+psf[i]*mask[i];
                    double amp=(kernel_size*sigmaXY-sigmaX*sigmaY)/sigmaY2N_sigmaYY;

                    int sindex=x+y*width;			  	
                    amp_map[sindex]=amp;		  			
                }	  		
            }		
        }

        return amp_map;
    }
    
    public double[] MinimumError_AMP(double[] imagef, int width, int height, int threadCounts)
    {
        int NN = imagef.length;
        double[] amp_map = new double[NN];
        
        ExecutorService exec=Executors.newFixedThreadPool(threadCounts);  
        List<Callable<Integer>> callList=new ArrayList<>();  
        List<Integer> list = new ArrayList<>();
        for(int y=psf_halfh; y<height-psf_halfh; y++) list.add(y);
        int list_size = list.size();
        int len = list_size/threadCounts;
        if(len==0){
            threadCounts = list_size;
            len = 1;
        }
        
        for(int i=0; i<threadCounts; i++){
            final List<Integer> sub_list;
            if(i==threadCounts-1) sub_list = list.subList(i*len, list_size);
            else sub_list=list.subList(i*len, len*(i+1)>list.size()?list.size():len*(i+1));
            
            //采用匿名内部类实现  
            callList.add((Callable<Integer>) () -> {
                for (Integer y : sub_list) {
                    int off = y*width;
                    for (int x = psf_halfw; x<width-psf_halfw; x++) {
                        //int p= x + y*width;
                        int p = x + off;
                        double sigmaXY = 0.0;
                        double sigmaX = 0.0;
                        for (int i1 = 0; i1 < psf_size; i1++) {
                            int x1 = x + psf_dx[i1];
                            int y1 = y + psf_dy[i1];
                            double mi = imagef[x1 + y1*width];
                            double dy = psf[i1];
                            sigmaX  += mi;//sigmaX + mi;                                    
                            sigmaXY += dy*mi;//sigmaXY + //psf[i]*mi;
                        }
                        double lxy = psf_size*sigmaXY-sigmaX*sigmaY;
                        double amp = lxy/sigmaY2N_sigmaYY;
                        amp_map[p] = amp;  
                    }
                }
                return 1;
            });  
        }
        
        int sum = 0;
        try{
            List<Future<Integer>> futureList=exec.invokeAll(callList);  
//            for(Future<Integer> future:futureList){  
//                sum+=future.get();  
//            }
        }
        catch(Exception ex){
            
        }
//        exec.shutdown();
        
        // do the margins only
        int xmargin=psf_halfw;
        int ymargin=psf_halfh;
        int kernelw = psf_halfw*2 + 1;
        int kernelh = psf_halfh*2 + 1;
        int kernel_size=kernelw*kernelh;
        double[] mask=new double[kernel_size];
        for(int x=0;x<width;x++) {
            for(int y=0;y<height;y++) {
                if((x<xmargin) || (x>=width-xmargin) || (y<ymargin) || (y>=height-ymargin)) {
                    int i=0;
                    for(int y1=y-ymargin; y1<=y+ymargin; y1++) {
                        for(int x1=x-xmargin; x1<=x+xmargin; x1++) {
                            //if((x1>=0) && (x1<width) && (y1>=0) && (y1<height)) mask[i]=slice[x1+width*y1]; else mask[i]=0;
                            int x2=x1;
                            if(x1<0) x2=0; else if (x1>=width) x2=width-1;
                            int y2=y1;
                            if(y1<0) y2=0; else if (y1>=height) y2=height-1;
                            mask[i]=imagef[x2+width*y2];
                            i++;
                        }
                    }
                    double sigmaX=0.;
                    for(i=0; i<kernel_size; i++) sigmaX=sigmaX+mask[i];
//                                    double sigmaX2=0.;
//                                    for(i=0; i<kernelw2; i++) sigmaX2=sigmaX2+mask[i]*mask[i];
                    // calculate the variable solver parameters
                    double sigmaXY=0.;
                    for(i=0; i<kernel_size; i++) sigmaXY=sigmaXY+psf[i]*mask[i];
                    double amp=(kernel_size*sigmaXY-sigmaX*sigmaY)/sigmaY2N_sigmaYY;

                    int sindex=x+y*width;			  	
                    amp_map[sindex]=amp;		  			
                }	  		
            }		
        }

        return amp_map;
    }
    
    private static double[] calc_sigma_psf(double[] psf, int size)
    {
            double[] array = new double[3];

            double sigmaY = 0.0;
            double sigmaY2 = 0.0;

            for(int i=0; i<size; i++) 
            {
                    sigmaY += psf[i];
                    sigmaY2 += psf[i] * psf[i];
            }

            double a_scale = size/(sigmaY2*size-sigmaY*sigmaY);

            array[0] = sigmaY;
            array[1] = sigmaY2;
            array[2] = a_scale;

            return array;
    }
    
    public static void MinimumError(double[] imagef, double[] amp_map, double[] offset_map, double[] err_map, int width, int height, double[] psf, int[] psf_dx, int[] psf_dy, int psf_size, int psf_halfw)
	{
		double[] array = calc_sigma_psf(psf, psf_size);
		double sigmaY = array[0];
		double sigmaY2 = array[1];
		double a_scale = array[2];
		
		//IJ.log("sigmaY=" + sigmaY);
		//IJ.log("sigmaY2=" + sigmaY2);
		//IJ.log("a_scale=" + a_scale);
		
		for(int x=psf_halfw; x<width-psf_halfw; x++)
		{
			for(int y=psf_halfw; y<height-psf_halfw; y++)
			{
				int p= x + y*width;
				
				double sigmaXY = 0.0;
				double sigmaX = 0.0;
				double sigmaX2 = 0.0;
				
				for(int i=0; i<psf_size; i++)
				{
					int x1 = x + psf_dx[i];
					int y1 = y + psf_dy[i];
		
					double mi = imagef[x1 + y1*width];
					double dy = psf[i];
		
					sigmaX  += mi;//sigmaX + mi;
					sigmaX2 += mi*mi;//sigmaX2 + mi*mi;
					sigmaXY += dy*mi;//sigmaXY + //psf[i]*mi;
				}
	
				double lxy = sigmaXY-sigmaX*sigmaY/psf_size;		
				double amp = lxy*a_scale;//sigmaXY/sigmaY2;//sigmaXY;//
				double o = (sigmaX - amp*sigmaY)/psf_size;
				double err = 0;
				
				double absdet = Math.abs(a_scale);
				if(absdet<1e-20){			
					amp = 0;
					o= 0;
					err = sigmaX2;			
				}				
				else{
					 err = amp*amp*sigmaY2 + psf_size*o*o + sigmaX2 + 2*amp*o*sigmaY - 2*amp*sigmaXY - 2*o*sigmaX;
				}
							
				amp_map[p] = amp;
				err_map[p] = err;
				offset_map[p] = o;			
			}
		}
	
		// do the margins only
		int margin=psf_halfw;
		int kernelw = psf_halfw*2 + 1;
		int kernelw2=kernelw*kernelw;
		double[] mask=new double[kernelw*kernelw];
		for(int x=0;x<width;x++) {
			for(int y=0;y<height;y++) {
				if((x<margin) || (x>=width-margin) || (y<margin) || (y>=height-margin)) {
					int i=0;
	  				for(int y1=y-margin; y1<=y+margin; y1++) {
			  			for(int x1=x-margin; x1<=x+margin; x1++) {
		  					//if((x1>=0) && (x1<width) && (y1>=0) && (y1<height)) mask[i]=slice[x1+width*y1]; else mask[i]=0;
		  					int x2=x1;
		  					if(x1<0) x2=0; else if (x1>=width) x2=width-1;
		  					int y2=y1;
		  					if(y1<0) y2=0; else if (y1>=height) y2=height-1;
		  					mask[i]=imagef[x2+width*y2];
		  					i++;
		  				}
		  			}
		  			double sigmaX=0.;
		  			for(i=0; i<kernelw2; i++) sigmaX=sigmaX+mask[i];
		  			double sigmaX2=0.;
		  			for(i=0; i<kernelw2; i++) sigmaX2=sigmaX2+mask[i]*mask[i];
		  			// calculate the variable solver parameters
				  	double sigmaXY=0.;
				  	for(i=0; i<kernelw2; i++) sigmaXY=sigmaXY+psf[i]*mask[i];
				  	double amp=(sigmaXY-(1./(double)kernelw2)*sigmaX*sigmaY)*a_scale;
				  	double offset  =(sigmaX-amp*sigmaY)/(double)kernelw2;
				  	double err=amp*amp*sigmaY2+(double)kernelw2*offset*offset+sigmaX2+2.*amp*offset*sigmaY-2.*amp*sigmaXY-2.*offset*sigmaX;
				  	double absdet = Math.abs(a_scale);
				  	if(absdet<1e-20){			
						amp = 0;
						offset= 0;
						err = sigmaX2;			
					}				
				  	int sindex=x+y*width;			  	
		  			amp_map[sindex]=amp;
		  			offset_map[sindex] = offset;			  	
				  	err_map[sindex]=err;
		  		}	  		
			}		
		}
	
		array = null;
	}
    
    public static void MinimumError_AMP(double[] imagef, double[] amp_map, int width, int height, double[] psf, int[] psf_dx, int[] psf_dy, int psf_size, int psf_halfw)
	{
		double[] array = calc_sigma_psf(psf, psf_size);
		double sigmaY = array[0];
		double sigmaY2 = array[1];
		double a_scale = array[2];

		for(int x=psf_halfw; x<width-psf_halfw; x++)
		{
			for(int y=psf_halfw; y<height-psf_halfw; y++)
			{
				int p= x + y*width;
				
				double sigmaXY = 0.0;
				double sigmaX = 0.0;
				double sigmaX2 = 0.0;
				
				for(int i=0; i<psf_size; i++)
				{
					int x1 = x + psf_dx[i];
					int y1 = y + psf_dy[i];
		
					double mi = imagef[x1 + y1*width];
					double dy = psf[i];
		
					sigmaX  += mi;//sigmaX + mi;
					sigmaX2 += mi*mi;//sigmaX2 + mi*mi;
					sigmaXY += dy*mi;//sigmaXY + //psf[i]*mi;
				}
	
				double lxy = sigmaXY-sigmaX*sigmaY/psf_size;		
				double amp = lxy*a_scale;//sigmaXY/sigmaY2;//sigmaXY;//
						
				amp_map[p] = amp;							
			}
		}
	
		// do the margins only
		int margin=psf_halfw;
		int kernelw = psf_halfw*2 + 1;
		int kernelw2=kernelw*kernelw;
		double[] mask=new double[kernelw*kernelw];
		for(int x=0;x<width;x++) {
			for(int y=0;y<height;y++) {
				if((x<margin) || (x>=width-margin) || (y<margin) || (y>=height-margin)) {
					int i=0;
	  				for(int y1=y-margin; y1<=y+margin; y1++) {
			  			for(int x1=x-margin; x1<=x+margin; x1++) {
		  					//if((x1>=0) && (x1<width) && (y1>=0) && (y1<height)) mask[i]=slice[x1+width*y1]; else mask[i]=0;
		  					int x2=x1;
		  					if(x1<0) x2=0; else if (x1>=width) x2=width-1;
		  					int y2=y1;
		  					if(y1<0) y2=0; else if (y1>=height) y2=height-1;
		  					mask[i]=imagef[x2+width*y2];
		  					i++;
		  				}
		  			}
		  			double sigmaX=0.;
		  			for(i=0; i<kernelw2; i++) sigmaX=sigmaX+mask[i];
		  			double sigmaX2=0.;
		  			for(i=0; i<kernelw2; i++) sigmaX2=sigmaX2+mask[i]*mask[i];
		  			// calculate the variable solver parameters
				  	double sigmaXY=0.;
				  	for(i=0; i<kernelw2; i++) sigmaXY=sigmaXY+psf[i]*mask[i];
				  	double amp=(sigmaXY-(1./(double)kernelw2)*sigmaX*sigmaY)*a_scale;
				  				
				  	int sindex=x+y*width;			  	
		  			amp_map[sindex]=amp;		  			
		  		}	  		
			}		
		}
	
		array = null;
	}
    
       /*direct solver for a and b*/
    public static double direct_solver(double[] win, double[] win_x, double[] win_y, int win_size, double sigma2, double[] v)
    {	
        double x0 = v[2];
        double y0 = v[3];
        //IJ.log("x0=" + x0 + "y0=" + y0 + "\n");

        double amp = 0; 
        double bkg = 0;
        double sigmaY = 0.0;
        double sigmaY2 = 0.0;		
        double sigmaXY = 0.0;		
        double sigmaX = 0;
        double sigmaX2 = 0;


        double min_bkg = win[0];
        for(int i=0; i<win_size; i++)
        {			
            if(min_bkg>win[i]) min_bkg = win[i];
            sigmaX += win[i];
            sigmaX2 += win[i]*win[i];			

            double dx = win_x[i]-x0;
            double dy = win_y[i]-y0;
            double psfi = Math.exp(-(dx*dx+dy*dy)/sigma2);			
            sigmaXY += psfi*win[i];
            sigmaY +=  psfi;
            sigmaY2 += psfi*psfi;		
        }	

        //return 0;		
        double err = sigmaX2;
        double det = sigmaY2*win_size-sigmaY*sigmaY;	
        if(det>1e-20 || det<-1e-20){
                amp = (sigmaXY*win_size-sigmaX*sigmaY)/det;
                bkg = (sigmaX - amp*sigmaY)/win_size;
                err = amp*amp*sigmaY2 + win_size*bkg*bkg + sigmaX2 + 2*amp*bkg*sigmaY - 2*amp*sigmaXY - 2*bkg*sigmaX;
        }			

//        double m10 = 0;
//        double m01 = 0;
//        double m00 = 0;
//        if(bkg<min_bkg && bkg>=0) min_bkg = bkg;
//
//        for(int i=0; i<win_size; i++)
//        {			
//                double w = win[i] - min_bkg;
//                m10 += win_x[i]*w;
//                m01 += win_y[i]*w;
//                m00 += w;
//        }
//
//        v[0] = amp;
//        v[1] = bkg;
//        v[2] = m10/m00;
//        v[3] = m01/m00;
        v[0] = amp;
        v[1] = bkg;
        v[2] = 0;
        v[3] = 0;
        //err = err/2;
        return err;
    }
    
    //load window, dx0, dy0
    public static int load_window(double[] pixels, int width, int height, int xc, int yc, int window_r, double[] win, double[] winx, double[] winy)
    {
        //IJ.log("x0=" + x0 + "y0=" + y0 + "\n");
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
        
        return k;
    }
    
    public static double[] estimate_amp_bkg(int xc, int yc, double[] pixels, int width, int height, int kernel_width, double gauss_sigma)
    {						
        int window_r = kernel_width/2;
        int window_w = window_r*2 + 1;
        int window_s = window_w*window_w;
        double[] win = new double[window_s];
        double[] winx= new double[window_s];
        double[] winy= new double[window_s];            
        int win_size = load_window(pixels,width,height,xc,yc,window_r, win, winx, winy);            
        double sigma2 = gauss_sigma*gauss_sigma*2;

        double[] par_seed = {0, 0, 0, 0};
        direct_solver(win, winx, winy, win_size, sigma2, par_seed);
                   
        return par_seed;		
    }
}
