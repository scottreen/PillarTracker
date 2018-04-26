package com.nus.mbi.pillar.detection;

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
public class CrossCorrelation {
    private double[] kernel;    
    private int[] kernel_dx;
    private int[] kernel_dy;
    
    private double[][] image;
    private int img_width;
    private int img_height;
    
    private int kernel_halfw;
    private int kernel_halfh; 
    private double avg_kernel;
    private double std_kernel;
            
    public void setImage(ImageProcessor ip) {        
        this.img_width = ip.getWidth();
        this.img_height = ip.getHeight();
        image = new double[img_width][img_height];
        float[][] pixels = ip.getFloatArray();  
        for(int x=0; x<img_width; x++){
            for(int y=0; y<img_height; y++){
                //image[x][y] = ip.getf(x, y);
                image[x][y] = pixels[x][y];
            }
        }
    }
    
    public void setImage(float[] pixels, int width, int height) {        
        this.img_width = width;
        this.img_height = height;
        image = new double[img_width][img_height];        
        for(int y=0; y<img_height; y++){
            int off = y*width;
            for(int x=0; x<img_width; x++){
                image[x][y] = pixels[off+x];
            }
        }
    }
       
    public void setPSFKernel(ImageProcessor ip_psf){        
        int psf_w = ip_psf.getWidth();
        int psf_h = ip_psf.getHeight();        
        int k_size = psf_w*psf_h;
        kernel = new double[k_size];
        kernel_dx = new int[k_size];
        kernel_dy = new int[k_size];
        avg_kernel = 0;
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
                avg_kernel += g;
                k++;
            }
        } 
        
        avg_kernel /= k_size;
        std_kernel = 0;
        for(k=0; k<k_size; k++){            
            double d = kernel[k] - avg_kernel;
            std_kernel += (d*d);
        }
        std_kernel = Math.sqrt(std_kernel/(k_size-1));
        
        kernel_halfw = half_w;
        kernel_halfh = half_h;
    }
    
    public ImagePlus getPSF(){
        if(kernel_halfw<1 || kernel_halfh<1) return null;
        if(kernel==null) return null;
        if(kernel.length<1) return null;
        
        int k_w = kernel_halfw*2 + 1;
        int k_h = kernel_halfh*2 + 1;
        FloatProcessor ip_psf = new FloatProcessor(k_w, k_h);
        int kernel_size = kernel.length;
        for(int k=0; k<kernel_size; k++){
            int x = kernel_dx[k] + kernel_halfw;
            int y = kernel_dy[k] + kernel_halfh;
            float p = (float)kernel[k];
            ip_psf.setf(x, y, (float)p);
        }            
        ImagePlus ip = new ImagePlus("PSF_CC", ip_psf);
        return ip;
    }
    
    public FloatProcessor getKernal(){
        if(kernel_halfw<1 || kernel_halfh<1) return null;
        if(kernel==null) return null;
        if(kernel.length<1) return null;
        
        int k_w = kernel_halfw*2 + 1;
        int k_h = kernel_halfh*2 + 1;
        FloatProcessor ip_psf = new FloatProcessor(k_w, k_h);
        int kernel_size = kernel.length;
        for(int k=0; k<kernel_size; k++){
            int x = kernel_dx[k] + kernel_halfw;
            int y = kernel_dy[k] + kernel_halfh;
            float p = (float)kernel[k];
            ip_psf.setf(x, y, (float)p);
        }            
        
        return ip_psf;
    }
    
    public ImagePlus setGaussianKernel(int radius, double sigma, boolean dark){
        int k_w = 2*radius + 1;
        int k_size = k_w*k_w;
        FloatProcessor ip_psf = new FloatProcessor(k_w, k_w);
        kernel = new double[k_size];
        kernel_dx = new int[k_size];
        kernel_dy = new int[k_size];
        avg_kernel = 0;
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
                avg_kernel += g;
                k++;
            }
        } 
        
        avg_kernel /= k_size;
        std_kernel = 0;
        for(k=0; k<k_size; k++){            
            double d = kernel[k] - avg_kernel;
            std_kernel += (d*d);
        }
        std_kernel = Math.sqrt(std_kernel/(k_size-1));
        
        kernel_halfw = radius;
        kernel_halfh = radius;
        
        ImagePlus ip = new ImagePlus("Gaussian_PSF", ip_psf);
        return ip;
    }
    
    public double getNormalizedCrossCorrelation(int x, int y){
        int k_size = kernel.length;
        
        double[] img = new double[k_size];
        double avg_img = 0;
        for(int k=0; k<k_size; k++){
            int dx = kernel_dx[k];
            int dy = kernel_dy[k];
            int xx = x + dx;
            int yy = y + dy;
            double I = image[xx][yy];
            img[k] = I;
            avg_img += I;
        }
        
        avg_img /= k_size;
        double std_img = 0;
        double sum_ft = 0;
        for(int k=0; k<k_size; k++){            
            double df = img[k] - avg_img;
            std_img += (df*df);
            
            double dt = kernel[k]-avg_kernel;
            sum_ft += df*dt;            
        }
        std_img = Math.sqrt(std_img/(k_size-1));
        
        double cross_cor = sum_ft/(k_size*std_img*std_kernel);
        return cross_cor;
    }
    
    public double getEdgeNormalizedCrossCorrelation(int x, int y){
        int k_size = kernel.length;
        
        double[] img = new double[k_size];
        int[] index = new int[k_size];
        int counter = 0;
        double avg_img = 0;
        double avg_psf = 0;
        for(int k=0; k<k_size; k++){
            int dx = kernel_dx[k];
            int dy = kernel_dy[k];
            int xx = x + dx;
            int yy = y + dy;            
            if(xx>=0 && xx<img_width && yy>=0 && yy<img_height){
                double I = image[xx][yy];
                img[counter] = I;
                index[counter] = k;
                counter++;
                avg_img += I;
                avg_psf += kernel[k];
            }
        }
        
        int num = counter;
        if(num<1) return 0;
        
        avg_img /= num;
        avg_psf /= num;
        double std_img = 0;
        double std_psf = 0;
        double sum_ft = 0;
        for(int k=0; k<num; k++){                        
            double df = img[k] - avg_img;
            std_img += (df*df);
            int p = index[k];            
            double dt = kernel[p]-avg_psf;
            std_psf += (dt*dt);
            sum_ft += df*dt;            
        }
        std_img = Math.sqrt(std_img/(num-1));
        std_psf = Math.sqrt(std_psf/(num-1));
        
        double cross_cor = sum_ft/(num*std_img*std_psf);
        return cross_cor;
    }
    
//    public double[] process(double[] img, int width, int height){
//        int size = width*height;        
//        this.img_width = width;
//        this.img_height = height;
//        image = new double[img_width][img_height];
//        for(int y=0; y<img_height; y++){
//                int off = y*width;
//                for(int x=0; x<img_width; x++) image[x][y] = img[off + x];
//        }
//        
//        double[][] img_cc = process();
//        double[] processed_img = new double[size];
//        for(int y=0; y<img_height; y++){
//                int off = y*width;
//                for(int x=0; x<img_width; x++) processed_img[off+x] = img_cc[x][y];
//        }
//        
//        return processed_img;
//    }    
    
    public double[][] process(){
        int xmin = kernel_halfw;
        int ymin = kernel_halfh;
        int xmax = img_width - 1 - kernel_halfw;
        int ymax = img_height - 1 - kernel_halfh;
        double[][] processed_img = new double[img_width][img_height];
        boolean[][] mask = new boolean[img_width][img_height];
        for(int x=0; x<img_width; x++){
            for(int y=0; y<img_height; y++){
                processed_img[x][y] = 0;
                mask[x][y] = true;
            }
        }
        
        for(int y=ymin; y<=ymax; y++){
            for(int x=xmin; x<=xmax; x++){            
                processed_img[x][y] = getNormalizedCrossCorrelation(x,y);
                mask[x][y] = false;
            }
        }
        
        for(int y=0; y<img_height; y++){
            for(int x=0; x<img_width; x++) {            
                if(mask[x][y]) processed_img[x][y] = getEdgeNormalizedCrossCorrelation(x,y);
            }
        }
        
        return processed_img;
    }    
    
    public double[][] process(int threadCounts){
        int xmin = kernel_halfw;
        int ymin = kernel_halfh;
        int xmax = img_width - 1 - kernel_halfw;
        int ymax = img_height - 1 - kernel_halfh;
        double[][] processed_img = new double[img_width][img_height];
        boolean[][] mask = new boolean[img_width][img_height];
        for(int x=0; x<img_width; x++){
            for(int y=0; y<img_height; y++){
                processed_img[x][y] = 0;
                mask[x][y] = true;
            }
        }
        
        ExecutorService exec=Executors.newFixedThreadPool(threadCounts);  
        List<Callable<Integer>> callList=new ArrayList<>();  
        List<Integer> list = new ArrayList<>();
        for(int y=ymin; y<=ymax; y++) list.add(y);
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
                    for(int x=xmin; x<=xmax; x++){            
                        processed_img[x][y] = getNormalizedCrossCorrelation(x,y);
                        mask[x][y] = false;
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
        
//        for(int y=ymin; y<=ymax; y++){
//            for(int x=xmin; x<=xmax; x++){            
//                processed_img[x][y] = getNormalizedCrossCorrelation(x,y);
//                mask[x][y] = false;
//            }
//        }
        
        for(int y=0; y<img_height; y++){
            for(int x=0; x<img_width; x++) {            
                if(mask[x][y]) processed_img[x][y] = getEdgeNormalizedCrossCorrelation(x,y);
            }
        }
        
        return processed_img;
    }    
    
    public ImageProcessor process(ImageProcessor ip){
        return process(ip, 1);
//        setImage(ip);
//        double[][] processed_img = process();
//        FloatProcessor fp = new FloatProcessor(img_width, img_height);
//        float[] cc_pixels = (float[])fp.getPixels();
//        for(int y=0; y<img_height; y++){
//            int off = y*img_width;
//            for(int x=0; x<img_width; x++){            
//                //fp.setf(x, y, (float)processed_img[x][y]);
//                cc_pixels[off + x] = (float)processed_img[x][y];
//            }
//        }
//        
//        return fp;
    }
    
    public ImageProcessor process(ImageProcessor ip, int num_threads){
        setImage(ip);
        double[][] processed_img = num_threads>1 ? process(num_threads) : process();
        FloatProcessor fp = new FloatProcessor(img_width, img_height);
        float[] cc_pixels = (float[])fp.getPixels();
        for(int y=0; y<img_height; y++){
            int off = y*img_width;
            for(int x=0; x<img_width; x++){            
                //fp.setf(x, y, (float)processed_img[x][y]);
                cc_pixels[off + x] = (float)processed_img[x][y];
            }
        }
        
        return fp;
    }
    
    
    public float[] process(float[] image, int width, int height, int num_threads){
        setImage(image,width, height);
        double[][] processed_img = num_threads>1 ? process(num_threads) : process();
        float[] cc_pixels = new float[width*height];
        for(int y=0; y<height; y++){
            int off = y*width;
            for(int x=0; x<width; x++){            
                cc_pixels[off + x] = (float)processed_img[x][y];
            }
        }
        
        return cc_pixels;
    }
    
    public ImagePlus process(ImagePlus image) {
        return process(image, 1);
//        int n = image.getStackSize();
//        int w = image.getWidth();
//        int h = image.getHeight();
//        ImageStack is = new ImageStack(w,h);       
//        
//        // slice numbers start with 1 for historical reasons
//        for (int i = 1; i <= n; i++){
//            ImageProcessor img = process(image.getStack().getProcessor(i));
//            is.addSlice(img);
//            IJ.showProgress(i, n);
//        }
//        
//        ImagePlus ip = new ImagePlus(image.getShortTitle()+"_CC", is);
//        return ip;                
    }
    
    public ImagePlus process(ImagePlus image, int num_threads) {
        int n = image.getStackSize();
        int w = image.getWidth();
        int h = image.getHeight();
        ImageStack is = new ImageStack(w,h);       
        
        // slice numbers start with 1 for historical reasons
        for (int i = 1; i <= n; i++){
            ImageProcessor img = process(image.getStack().getProcessor(i), num_threads);
            is.addSlice(img);
            IJ.showProgress(i, n);
        }
        
        ImagePlus ip = new ImagePlus(image.getShortTitle()+"_CC", is);
        return ip;                
    }
}
