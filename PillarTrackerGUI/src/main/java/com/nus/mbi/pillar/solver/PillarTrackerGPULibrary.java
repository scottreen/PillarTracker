package com.nus.mbi.pillar.solver;

import com.sun.jna.Native;
import com.sun.jna.NativeLibrary;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.io.File;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.OutputStream;

/**
 *
 * @author xiaochun
 */
public class PillarTrackerGPULibrary {    
    protected static final String resoucefolder = "/lib/";
    protected static final String libray_name = "MinimumErrorGPU";
    //protected static final String librayfnameMacOSX = "libMinimumErrorGPUGSD.dylib";
    //protected static final String cudafname64MacOSX = "libcudart.dylib";

    protected static final String librayfname = "MinimumErrorGPU.dll";
    protected static final String cudafname64 = "cudart64_80.dll";

    protected static String libraypath = librayfname;
    private static int maximum_kernel_size = 1000;
    public GPUCrossCorrelationLibrary gpu_instance = null;
    private float[] psf_kernel = null;
    private int[] psf_kernel_dx = null;
    private int[] psf_kernel_dy = null;
    int psf_kernel_size = 0;
    
    public static boolean check_psf_size(int psf_w, int psf_h){
        int size = psf_w*psf_h+3;
        return (size<maximum_kernel_size);
    }
    
    public void copyFile(String source, String dest) throws Exception
    {
        InputStream stream = getClass().getResourceAsStream(source);
        //in the "jar tree" been the jar the root of the tree"
        if (stream == null) {
            Exception ex = new Exception("the inputstream is null.");
            throw ex;
        }	   
        OutputStream resStreamOut = null;	    
        byte[] buffer = new byte[4096];
        try {
            resStreamOut = new FileOutputStream(new File(dest));
            if(resStreamOut == null) IJ.log("resStreamOut is null");
            int readBytes;
            long totalBytes=0; 
            do{
                readBytes = stream.read(buffer);	            
                if(readBytes>0){
                    resStreamOut.write(buffer, 0, readBytes);	            
                    totalBytes += readBytes;
                }
            }while (readBytes>0);
            IJ.log("copy " + totalBytes + " Bytes");	
            resStreamOut.close();
            stream.close();

        } catch (Exception e1) {	    	
            // TODO Auto-generated catch block
            //e1.printStackTrace();
            IJ.log(e1.getMessage());
            if(resStreamOut!=null) resStreamOut.close();
            if(stream!=null) stream.close();
            throw e1;
        }
    }

    private boolean checkLibrays(String newfname, String fname)
    {
            boolean suc = true;
            try{
                    File newDir = new File(newfname);
                    if(!newDir.exists()){			
                            IJ.log("libray is missing");
                            IJ.log("Extracting \"" + fname + " from the local");
                            IJ.log("->" + newfname);
                            copyFile(fname,newfname);
                            IJ.log("Extracting is successful");
                    }
            }
            catch (Exception e) {	            
                suc = false;
                String msg = e.getMessage();
                if (msg==null || msg.equals("")) msg = "" + e;	
                IJ.log(msg + "\n Extracting the dll files is not successful");	            
            }
            return suc;
    }
    
    public boolean overwriteLibrays(String newfname, String fname)
    {
            boolean suc = true;
            try{
                    File newDir = new File(newfname);
                    //if(!newDir.exists())
                    {			
                            //IJ.log("libray is missing");
                            IJ.log("Extracting \"" + fname + " from the local");
                            IJ.log("->" + newfname);
                            copyFile(fname,newfname);
                            IJ.log("Extracting is successful");
                    }
            }
            catch (Exception e) {	            
                suc = false;
                String msg = e.getMessage();
                if (msg==null || msg.equals("")) msg = "" + e;	
                IJ.log(msg + "\n Extracting the dll files is not successful");	            
            }
            return suc;
    }

    public void updateLibrary(){
        if(!IJ.is64Bit() || IJ.isLinux() || IJ.isMacOSX()) {
                IJ.showMessage("GPU accelaration Only support on 64-bit Windows OS");
                return;				
        }

        String lib_fname = librayfname;
        String cu_fname = cudafname64;                        

        File classpathRoot =new File(IJ.getDir("imagej"));

        String rootfolder = classpathRoot.getAbsolutePath();
        NativeLibrary.addSearchPath(libray_name,rootfolder);                        
        System.setProperty("jna.library.path",  rootfolder);
        System.setProperty("java.library.path", rootfolder);
        //IJ.log("the library folder:" + rootfolder);

        String newlibrayfname = rootfolder + "/" + lib_fname;
        boolean suc_load_dll = overwriteLibrays(newlibrayfname,resoucefolder+lib_fname);        
        if(suc_load_dll){
                String newcudafname64 = rootfolder + "/" + cu_fname;
                suc_load_dll = overwriteLibrays(newcudafname64,resoucefolder+cu_fname);               
        }
    }
    
    public boolean prepare(){        
        if(!IJ.is64Bit() || IJ.isLinux() || IJ.isMacOSX()) {
                IJ.showMessage("GPU accelaration Only support on 64-bit Windows OS");
                return false;				
        }

        String lib_fname = librayfname;
        String cu_fname = cudafname64;                        
//			if(IJ.isMacOSX()){
//				lib_fname = librayfnameMacOSX;
//				cu_fname = cudafname64MacOSX;
//			}
//			
        File classpathRoot = new File(this.getClass().getProtectionDomain().getCodeSource().getLocation().getPath());
        if(classpathRoot.isDirectory()==false) classpathRoot=new File(IJ.getDir("imagej"));
        
        String rootfolder = classpathRoot.getAbsolutePath();
        NativeLibrary.addSearchPath(libray_name,rootfolder);                        
        System.setProperty("jna.library.path",  rootfolder);
        System.setProperty("java.library.path", rootfolder);
        //IJ.log("the library folder:" + rootfolder);

        String newlibrayfname = rootfolder + "/" + lib_fname;
        boolean suc_load_dll = checkLibrays(newlibrayfname,resoucefolder+lib_fname);
        //boolean suc_load_dll = checkLibrays(newlibrayfname,lib_fname);
        if(suc_load_dll){
                String newcudafname64 = rootfolder + "/" + cu_fname;
                suc_load_dll = checkLibrays(newcudafname64,resoucefolder+cu_fname);
                //suc_load_dll = checkLibrays(newcudafname64,cu_fname);
        }

        if(suc_load_dll){
//                            libraypath = newlibrayfname;
            IJ.log("the library folder:" + rootfolder);
            //IJ.log("JNA version:" + Native.VERSION + "  Native Verison:" + Native.VERSION_NATIVE);
            if(gpu_instance==null)
                gpu_instance = (GPUCrossCorrelationLibrary)Native.loadLibrary(libray_name, GPUCrossCorrelationLibrary.class);
            float cuda_sm = gpu_instance.checkGPU();       
            if(cuda_sm>0) IJ.log("GPU-cuda compute capablity is " + cuda_sm);
            else{
                    IJ.showMessage("PillarTracker","GPU is not avaliable, error code=" + cuda_sm);
                    return false;
            }
        }
        else{
            IJ.log("The GPU library can not be loaded!");
            IJ.showMessage("PillarTracker","The GPU library can not be loaded");
            return false;
        }
        
        return true;
    }
    
    public void setKernel(ImageProcessor psf){
        int psf_w = psf.getWidth();
        int psf_h = psf.getHeight();        
        int k_size = psf_w*psf_h;
        psf_kernel_size = k_size;
        psf_kernel = new float[k_size];
        psf_kernel_dx = new int[k_size];
        psf_kernel_dy = new int[k_size];        
        int k = 0;
        int half_w = psf_w/2;
        int half_h = psf_h/2;
        for(int y=-half_h; y<=half_h; y++){
            int j = y + half_h;
            for(int x=-half_w; x<=half_w; x++){
                int i = x + half_w;            
                float g = psf.getf(i, j);  
                psf_kernel[k] = g;
                psf_kernel_dx[k] = x;
                psf_kernel_dy[k] = y;                                
                k++;
            }
        } 
        
    }
    
    public FloatProcessor process(ImageProcessor ip, ImageProcessor psf, boolean use_normalized){
        float[][] slice=ip.getFloatArray();
        
        int width = ip.getWidth();
        int height = ip.getHeight();
        int size = width*height;
        float[]image = new float[size];        
        for(int y=0; y<height; y++){
            int off = y*width;
            for(int x=0; x<width; x++){
                int p = off + x;
                image[p] = slice[x][y];
            }
        }
        int psf_w = psf.getWidth();
        int psf_h = psf.getHeight();        
        int k_size = psf_w*psf_h;
        float[] kernel = new float[k_size];
        int[] kernel_dx = new int[k_size];
        int[] kernel_dy = new int[k_size];        
        int k = 0;
        int half_w = psf_w/2;
        int half_h = psf_h/2;
        for(int y=-half_h; y<=half_h; y++){
            int j = y + half_h;
            for(int x=-half_w; x<=half_w; x++){
                int i = x + half_w;            
                float g = psf.getf(i, j);  
                kernel[k] = g;
                kernel_dx[k] = x;
                kernel_dy[k] = y;                                
                k++;
            }
        } 
        
        FloatProcessor fp = new FloatProcessor(width, height);
        float[] amp_pixels = (float[])fp.getPixels();        
        if(use_normalized)
            gpu_instance.CrossCorrelationGPU(image,amp_pixels, width, height, kernel, kernel_dx, kernel_dy,k_size);
        else
            gpu_instance.MinimumErrorAMPGPU(image,amp_pixels,width, height, kernel, kernel_dx, kernel_dy,k_size);
        
        return fp;
    }
    
    public ImagePlus process(ImagePlus image, ImageProcessor psf, boolean use_normalized) {
        int n = image.getStackSize();
        int w = image.getWidth();
        int h = image.getHeight();
        ImageStack is = new ImageStack(w,h);       
        
        // slice numbers start with 1 for historical reasons
        for (int i = 1; i <= n; i++){
            ImageProcessor img = process(image.getStack().getProcessor(i), psf, use_normalized);
            is.addSlice(img);
            IJ.showProgress(i, n);
        }
        String title = use_normalized ? "_CC": "_AMP";
        ImagePlus ip = new ImagePlus(image.getShortTitle()+title, is);
        return ip;                
    }
    
    public FloatProcessor process(ImageProcessor ip, boolean use_normalized){
        float[][] slice=ip.getFloatArray();        
        int width = ip.getWidth();
        int height = ip.getHeight();
        int size = width*height;
        float[]image = new float[size];        
        for(int y=0; y<height; y++){
            int off = y*width;
            for(int x=0; x<width; x++){
                int p = off + x;
                image[p] = slice[x][y];
            }
        }
        
        FloatProcessor fp = new FloatProcessor(width, height);
        float[] amp_pixels = (float[])fp.getPixels();        
        if(use_normalized)
            gpu_instance.CrossCorrelationGPU(image,amp_pixels, width, height, psf_kernel, psf_kernel_dx, psf_kernel_dy,psf_kernel_size);
        else
            gpu_instance.MinimumErrorAMPGPU(image,amp_pixels,width, height, psf_kernel, psf_kernel_dx, psf_kernel_dy,psf_kernel_size);
        
        return fp;
    }
    
    public float[] process(float[] image, int width, int height, boolean use_normalized){
        int size = width*height;
        float[] amp_pixels = new float[size];
        if(use_normalized)
            gpu_instance.CrossCorrelationGPU(image,amp_pixels, width, height, psf_kernel, psf_kernel_dx, psf_kernel_dy,psf_kernel_size);
        else
            gpu_instance.MinimumErrorAMPGPU(image,amp_pixels,width, height, psf_kernel, psf_kernel_dx, psf_kernel_dy,psf_kernel_size);
        
        return amp_pixels;
    }
    
    public ImagePlus process(ImagePlus image, boolean use_normalized){
        int n = image.getStackSize();
        int w = image.getWidth();
        int h = image.getHeight();
        ImageStack is = new ImageStack(w,h);       
        
        // slice numbers start with 1 for historical reasons
        for (int i = 1; i <= n; i++){
            ImageProcessor img = process(image.getStack().getProcessor(i), use_normalized);
            is.addSlice(img);
            IJ.showProgress(i, n);
        }
        
        String title = use_normalized ? "_CC": "_AMP";
        ImagePlus ip = new ImagePlus(image.getShortTitle()+title, is);
        return ip;                        
    }
}
