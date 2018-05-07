package com.nus.mbi.pillar.detection;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import java.util.ArrayList;
import com.nus.mbi.pillar.solver.MinimumErrorSolver;
import com.nus.mbi.pillar.solver.PillarTrackerGPULibrary;

/**
 *
 * @author xiaochun
 */
public class CrossCorreclation_Plugin implements PlugIn{
    private ImagePlus experiment;
    private ImagePlus psf_ip;
    private boolean use_gaussian_psf = false;
    private int kernel_radius = 4;
    private double sigma = 1.74;
    private boolean dark = true;
    public  boolean keep_current_kernel = true;
    private boolean mode_cc = true;
    private boolean gpu_on = IJ.isWindows();
    private boolean is_process_current = true;
    private boolean check_once = false;
    MinimumErrorSolver amp_solver = new MinimumErrorSolver();
    CrossCorrelation cc = new CrossCorrelation();
    public PillarTrackerGPULibrary gpu_lib = new PillarTrackerGPULibrary();
    //private boolean isGPULibLoaded = false;
    private int num_threads = 1;
    public void setup(ImagePlus ip, boolean is_cc_mode, boolean is_single){
        experiment = ip;
        mode_cc = is_cc_mode;
        is_process_current = is_single;
    }
    
    public void setMode(boolean is_cc_mode){
        mode_cc = is_cc_mode;
    }
    
    public void setGPUOn(boolean use_gpu){
        gpu_on = use_gpu;
    }
    
    public boolean isGPUOn(){
        return gpu_on;
    }
    
    public void setToProcessSingleFrame(boolean is_single){
        is_process_current = is_single;
    }
    
    public void setCustomPSF(ImagePlus psf){
        use_gaussian_psf = false;
        psf_ip = psf;
        
        if(mode_cc) cc.setPSFKernel(psf_ip.getProcessor());  
        else amp_solver.setPSF(psf_ip.getProcessor());  
        
        if(gpu_on && psf_ip!=null){
            gpu_on=PillarTrackerGPULibrary.check_psf_size(psf_ip.getWidth(), psf_ip.getHeight());
            if(gpu_on) gpu_lib.setKernel(psf_ip.getProcessor());
            else IJ.log("the kernel size is not suitable for the GPU computation,use cpu instead");                            
        }
    }
    
    public void setGaussianPSF(double sigma, int radius, boolean dark){
        use_gaussian_psf = true;
        if(radius>0) kernel_radius = radius;
        if(sigma>0) this.sigma = sigma;
        this.dark = dark;
        
        if(mode_cc) psf_ip = cc.setGaussianKernel(kernel_radius, sigma, dark);
        else psf_ip = amp_solver.setGaussianKernel(kernel_radius, sigma, dark);    
        
        if(gpu_on){
            gpu_on=PillarTrackerGPULibrary.check_psf_size(psf_ip.getWidth(), psf_ip.getHeight());
            if(gpu_on) gpu_lib.setKernel(psf_ip.getProcessor());
            else IJ.log("the kernel size is not suitable for the GPU computation,use cpu instead");                            
        }
    }
    
    public void setGaussianParameters(double sigma, int radius, boolean dark){
        if(radius>0) kernel_radius = radius;
        if(sigma>0) this.sigma = sigma;
        this.dark = dark;
    }
    
    public void setNumThreads(int num_threads){
        if(num_threads>0) this.num_threads = num_threads;
    }
    
    public boolean getMode(){
        return mode_cc;
    }
    
    public ImagePlus getPSF(){
        psf_ip = mode_cc? cc.getPSF():amp_solver.getPSF();
        return psf_ip;
    }
    
    public boolean isUsingGaussianPSF(){
        return use_gaussian_psf;
    }
    
    public boolean isDarkObject(){
        return dark;
    }
    
    public double getGaussianSigma(){
        return sigma;
    }
    
    public int getGaussianRaidus(){
        return kernel_radius;
    }
    
    /**
    * This method gets called by ImageJ / Fiji.
    *
    * @param arg can be specified in plugins.config
    * @see ij.plugin.PlugIn#run(java.lang.String)
    */
    @Override
    public void run(String arg) {        
        check_once = arg.contains("once");
        boolean consistent;	
        do {
                consistent=true;
                if (!showDialog()) return;					
                if (kernel_radius<0) {IJ.showMessage("the radius of kernel must greater than zero"); consistent=false;}
                else if (sigma<0) {IJ.showMessage("the gaussian sigma must greater than zero"); consistent=false;}
        } while (!consistent && !check_once);		
        if(!consistent) return;
        
        //isGPULibLoaded = false;
        if(gpu_on){            
            boolean suc = gpu_lib.prepare();
            //isGPULibLoaded = suc;
            if(!suc){
                gpu_on = false;
                return;
            }
        }

        setKernel();
                
        
        ImagePlus ip;
        if(gpu_on){
            //gpu_lib.setKernel(psf_ip.getProcessor());
            if(is_process_current){ 
                ImageProcessor img = gpu_lib.process(experiment.getProcessor(), mode_cc);
                String title = mode_cc ? "_CC": "_AMP";
                ip = new ImagePlus(experiment.getShortTitle()+title, img);
               
            }
            else ip = gpu_lib.process(experiment, mode_cc);            
            ip.show();
            return; 
        }
        
        if(!mode_cc){
            if(is_process_current){                
                ImageProcessor img = amp_solver.MinimumError_AMP(experiment.getProcessor(), num_threads);
                ip = new ImagePlus(experiment.getShortTitle()+"_AMP", img);
            }
            else ip = amp_solver.MinimumError_AMP(experiment, num_threads);               
        }
        else{            
            if(is_process_current){
               ImageProcessor img = cc.process(experiment.getProcessor(), num_threads);
               ip = new ImagePlus(experiment.getShortTitle()+"_CC", img);
            }
            else ip = cc.process(experiment, num_threads);            
        }
        
        ip.show();
    }
    
    boolean showDialog() {
	
	int[] wList = WindowManager.getIDList();	
	if (wList==null || wList.length<1) {
		IJ.showMessage("Cross Correlation", "one image is required");	
		return false;
	}

	int wlistlen = wList.length;
	
	ArrayList<String> titles_list = new ArrayList<String>();
        ArrayList<String> titles_psf_list = new ArrayList<String>();
	for (int i=0; i<wlistlen; i++) {
		ImagePlus imp = WindowManager.getImage(wList[i]);
		if(imp!=null){
                    String title = imp.getTitle();
                    titles_list.add(title);
                    if(imp.getStackSize()==1) titles_psf_list.add(title);                    
                }
	}

	if(psf_ip==null && titles_psf_list.isEmpty()){
		IJ.showMessage("Cross Correlation", "one kernel (pillar profile) image is required");	
		return false;
	}
	
	String titles[] = titles_list.toArray(new String[0]);
        String gd_title = "Image Enhancement";        
	GenericDialog gd = new GenericDialog(gd_title);
	gd.addChoice("ONLY one channel time series:", titles, titles[0]);
        
        String use_gaussian_psf_str = "Use Gaussian Kernel";
        if(!check_once) titles_psf_list.add(use_gaussian_psf_str); 
        String use_current_kernel_str = "Use Current Kernel";        
        if(psf_ip!=null)titles_psf_list.add(use_current_kernel_str);  
        String titles_psf[] = titles_psf_list.toArray(new String[0]);        
        int titles_len = titles_psf.length;
        gd.addChoice("kernel (pillar profile):", titles_psf, titles_psf[titles_len-1]);
        if(psf_ip!=null) gd.addImage(psf_ip);
        
        gd.addCheckbox("processing_one frame?", is_process_current);
        gd.addCheckbox("normalized_cross_correlation? ", mode_cc);
        gd.addCheckbox("use GPU acceleration? ", gpu_on);
        if(!check_once){
            gd.addMessage("------------settings if using Gaussian Kernel-------------");            
            gd.addNumericField("   kernel radius:", kernel_radius, 0);
            gd.addNumericField("   Gaussian sigma:", sigma, 2);
            gd.addCheckbox("dark object?", dark);
        }
        gd.showDialog();
	if (gd.wasCanceled()) return false;
	String choice_name1 = gd.getNextChoice();
	experiment = WindowManager.getImage(choice_name1);
        
        String choice_name2 = gd.getNextChoice();
        keep_current_kernel = choice_name2.equalsIgnoreCase(use_current_kernel_str);
        if(!keep_current_kernel){
            use_gaussian_psf = choice_name2.equalsIgnoreCase(use_gaussian_psf_str);
            if(!use_gaussian_psf) psf_ip = WindowManager.getImage(choice_name2);
        }
        is_process_current = gd.getNextBoolean();
        mode_cc = gd.getNextBoolean();  
        gpu_on = gd.getNextBoolean();  
	if(!check_once){ 
            int radius = (int) gd.getNextNumber();            
            double s = (double)gd.getNextNumber();
            if(use_gaussian_psf){
                kernel_radius = radius;
                sigma = s;
                dark = gd.getNextBoolean();        
            }
        }
	return true;
    }
    
    public ImagePlus process(){
        if(experiment==null) return null;
        return process(experiment);
    }
    
    public ImagePlus process(ImagePlus experiment){        
        ImagePlus ip;
        if(gpu_on){
            gpu_lib.setKernel(psf_ip.getProcessor());
            if(is_process_current){ 
                ImageProcessor img = gpu_lib.process(experiment.getProcessor(), mode_cc);
                String title = mode_cc ? "_CC": "_AMP";
                ip = new ImagePlus(experiment.getShortTitle()+title, img);
               
            }
            else ip = gpu_lib.process(experiment, mode_cc);            
            //ip.show();
            return ip; 
        }
        
        if(!mode_cc){
            if(use_gaussian_psf) psf_ip = amp_solver.setGaussianKernel(kernel_radius, sigma, dark);
            else amp_solver.setPSF(psf_ip.getProcessor());
            if(is_process_current){
                ImageProcessor img = amp_solver.MinimumError_AMP(experiment.getProcessor(), num_threads);
                ip = new ImagePlus(experiment.getShortTitle()+"_AMP", img);
            }
            else ip = amp_solver.MinimumError_AMP(experiment, num_threads);                
            
            //ip.show();
        }
        else{
            //CrossCorrelation cc = new CrossCorrelation();
            if(use_gaussian_psf) psf_ip = cc.setGaussianKernel(kernel_radius, sigma, dark);
            else cc.setPSFKernel(psf_ip.getProcessor());
            
            if(is_process_current){
               ImageProcessor img = cc.process(experiment.getProcessor(), num_threads);
               ip = new ImagePlus(experiment.getShortTitle()+"_CC", img);
            }
            else ip = cc.process(experiment, num_threads);
            //ip.show();
        }
        return ip;
    }
    
    public void setKernel(){        
        if(use_gaussian_psf){
            if(mode_cc) psf_ip = cc.setGaussianKernel(kernel_radius, sigma, dark);
            else psf_ip = amp_solver.setGaussianKernel(kernel_radius, sigma, dark);
        }
        else if(psf_ip!=null){
            if(mode_cc) cc.setPSFKernel(psf_ip.getProcessor());  
            else amp_solver.setPSF(psf_ip.getProcessor());  
        }

        if(gpu_on && psf_ip!=null){
            gpu_on=PillarTrackerGPULibrary.check_psf_size(psf_ip.getWidth(), psf_ip.getHeight());
            if(gpu_on) gpu_lib.setKernel(psf_ip.getProcessor());
            else IJ.log("the kernel size is not suitable for the GPU computation,use cpu instead");                            
        }
    }
    
    public ImageProcessor process(ImageProcessor ip){
        if(gpu_on){
            ImageProcessor img = gpu_lib.process(ip, mode_cc);
            return img; 
        }
        
        ImageProcessor img = mode_cc ? cc.process(ip, num_threads) : amp_solver.MinimumError_AMP(ip, num_threads);        
        return img;
    }
    
    public float[] process(float[] pixels, int width, int height){
        if(gpu_on){
            float[] img = gpu_lib.process(pixels, width, height, mode_cc);
            return img; 
        }
        
        float[] img = mode_cc ? cc.process(pixels,width, height,num_threads) : amp_solver.MinimumError_AMP(pixels, width, height, num_threads);        
        return img;
    }
    
    public boolean showSetKernelDialog() {	
	int[] wList = WindowManager.getIDList();
        int wlistlen = 0;
        ArrayList<String> titles_list = new ArrayList<String>();
        if(wList!=null){
            wlistlen = wList.length;	            
            for (int i=0; i<wlistlen; i++) {
                    ImagePlus imp = WindowManager.getImage(wList[i]);
                    if(imp!=null && imp.getStackSize()<2) titles_list.add(imp.getTitle());
            }
        }
	
        String use_gaussian_psf_str = "Use Gaussian Kernel";
        String use_current_kernel_str = "Use Current Kernel";
        titles_list.add(use_gaussian_psf_str);
        if(psf_ip!=null)titles_list.add(use_current_kernel_str);        
        String titles_psf[] = titles_list.toArray(new String[0]);        
                
	GenericDialog gd = new GenericDialog("Image Enhancement");       
        int titles_len = titles_psf.length;
        gd.addChoice("kernel (pillar profile):", titles_psf, titles_psf[titles_len-1]);  
        if(psf_ip!=null) gd.addImage(psf_ip);
        //gd.addCheckbox("use current kernel?", keep_current_kernel);  
        gd.addCheckbox("normalized_cross_correlation? ", mode_cc);
        gd.addCheckbox("use GPU acceleration? ", gpu_on);
        gd.addMessage("------------settings if using Gaussian Kernel-------------");
	//gd.addCheckbox("using Gaussian?", use_gaussian_psf);
	gd.addNumericField("   kernel radius:", kernel_radius, 0, 10, "pixel");
        gd.addNumericField("   Gaussian sigma:", sigma, 2, 10, "pixel");
        gd.addCheckbox("dark object?", dark);
	gd.showDialog();
	if (gd.wasCanceled()) return false;
	
        String choice_name2 = gd.getNextChoice();
        keep_current_kernel = choice_name2.equalsIgnoreCase(use_current_kernel_str);
        if(!keep_current_kernel){
            use_gaussian_psf = choice_name2.equalsIgnoreCase(use_gaussian_psf_str);
            if(!use_gaussian_psf)psf_ip = WindowManager.getImage(choice_name2);
        }
        else use_gaussian_psf = false;
        //keep_current_kernel = gd.getNextBoolean(); 
        mode_cc = gd.getNextBoolean();        
        boolean use_gpu = gd.getNextBoolean();        
        
        if(use_gaussian_psf){ 
            int radius = (int) gd.getNextNumber();        
            double s = (double)gd.getNextNumber();  
            if(radius>0) kernel_radius = radius;
            else{
                IJ.showMessage("the radius of kernel must greater than zero");                
                return false;
            }
            if(s>0) sigma = s;
            else{
                IJ.showMessage("the gaussian sigma must greater than zero");
                return false;
            }
            dark = gd.getNextBoolean();                        
        }     
        
        if(use_gpu){                        
            boolean suc = gpu_lib.prepare();            
            if(!suc){
                use_gpu = false;            
                IJ.log("GPU computation is not avaliable");
                //return false;
            }
        }       
        
        gpu_on = use_gpu;                    
        setKernel();        
	return true;
    }

}
