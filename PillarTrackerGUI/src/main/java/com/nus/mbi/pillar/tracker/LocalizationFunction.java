package com.nus.mbi.pillar.tracker;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import com.nus.mbi.pillar.solver.LevmarSolver;
/**
 *
 * @author xiaochun
 */
class LocalizationFunction{
        
        public static double[] metric_centroid_safe(double[] pixels, int width, int height, int xc, int yc, int radius)
        {
            int k = radius;
            if(xc<k || xc>=width-k) return null;
            if(yc<k || yc>=height-k) return null;
            
            double M00, M10, M01;
            M00 = M10 = M01 = 0;
                   
            int r2 = radius*radius;
            int num = 0;
            for (int i=-k; i<=k; i++)
            {
                int y = i + yc;
                int off = y*width;
                int i2 = i*i;
                for (int j=-k; j<=k; j++)
                {
                    if(j*j+i2<=r2){
                        int x = j + xc;
                        double t = pixels[off + x];
                        if(t>0){
                            M00 += t;
                            M10 += x*t;
                            M01 += y*t;		            	   
                            num++;
                        }		   
                    }
                }
            }

            if(num<1) return null;

            double xm = M10 / M00;
            double ym = M01 / M00;

            return new double[] {xm, ym};		    
        }
        
        private static boolean check_boundray( double[] metric_mlxy, int xc, int yc, double max_drift, double[] pixels, int width, int height){
            boolean accept = false;
            if(metric_mlxy!=null){
                double mcx = metric_mlxy[0];
                double mcy = metric_mlxy[1];
                double dx = Math.abs(mcx-xc);
                double dy = Math.abs(mcy-yc); 
                if(dx<max_drift && dy<max_drift){
                    xc = (int)Math.round(mcx);
                    yc = (int)Math.round(mcy);
                    if(xc>=0 && xc<width && yc>=0 && yc<height){
                        if(pixels[yc*width+xc]>0){
                            accept=true;
                        }
                    }
                }
            }
            return accept;
        }
        
        public static double[] metric_centroid_safe(double[] pixels, int width, int height, int xc, int yc, int radius, int times, double max_drift){            
            int xx = xc;
            int yy = yc;
            double[] metric_mlxy = metric_centroid_safe(pixels, width, height, xc,yc, radius);            
            boolean suc = check_boundray(metric_mlxy, xc, yc, max_drift, pixels, width, height);
            if(!suc) return null;
            
            double[] mlxy = metric_mlxy;
            
            int k = 1;
            while(k<=times && metric_mlxy!=null){                
                k++;                                
                double mcx = metric_mlxy[0];            
                double mcy = metric_mlxy[1];  
                double dx = Math.abs(mcx-xc);
                double dy = Math.abs(mcy-yc);   
                if(dx>1 || dy>1){
                    xc = (int)Math.round(mcx);
                    yc = (int)Math.round(mcy);
                    metric_mlxy = metric_centroid_safe(pixels, width, height, xc,yc, radius);  
                    suc = check_boundray(metric_mlxy, xx, yy, max_drift, pixels, width, height);
                    if(suc) mlxy = metric_mlxy;
                    else metric_mlxy = null;                    
                }                    
                else break;           
            }          
            
            return mlxy;
        }
        
        public static double[][] metric_centroid_safe(double[] pixels, int width, int height, double[] xc, double[] yc, int radius)
        {
            int num = xc.length;
            double[][] metric_mlxy = new double[num][];
            for(int i=0; i<num; i++){
                if((Double.isNaN(xc[i])) || Double.isNaN(yc[i])){
                    metric_mlxy[i] = null;
                }
                else{
                    int xx = (int)Math.round(xc[i]);
                    int yy = (int)Math.round(yc[i]);
                    metric_mlxy[i] = metric_centroid_safe(pixels, width, height, xx, yy, radius);                
                }
            }
            return metric_mlxy;
        }
       
        public static double[][] metric_centroid_safe(int num_threads, double[] pixels, int width, int height, double[] xc, double[] yc, int radius, int times, double max_drift)
        {
            if(xc==null || xc.length<1) return null; 
            int num = xc.length;
            double[][] metric_mlxy = new double[num][];            
            if(num_threads<2){
                for(int k=0; k<num; k++){
                    if((Double.isNaN(xc[k])) || Double.isNaN(yc[k])){
                        metric_mlxy[k] = null;
                    }
                    else{
                        int xx = (int)Math.round(xc[k]);
                        int yy = (int)Math.round(yc[k]);
                        metric_mlxy[k] = metric_centroid_safe(pixels, width, height, xx, yy, radius, times, max_drift);
                    }
                }
                return metric_mlxy;
            }                        
            
            ExecutorService exec=Executors.newFixedThreadPool(num_threads);  
            List<Callable<Integer>> callList=new ArrayList<>();  
            List<Integer> list = new ArrayList<>();
            for(int i=0; i<num; i++) list.add(i);
            int list_size = list.size();
            int len = list_size/num_threads;
            if(len==0){
                num_threads = list_size;
                len = 1;
            }

            for(int i=0; i<num_threads; i++){
                final List<Integer> sub_list;
                if(i==num_threads-1) sub_list = list.subList(i*len, list_size);
                else sub_list=list.subList(i*len, len*(i+1)>list.size()?list.size():len*(i+1));

                //  
                callList.add((Callable<Integer>) () -> {
                    for(Integer k : sub_list) {
                        if((Double.isNaN(xc[k])) || Double.isNaN(yc[k])){
                            metric_mlxy[k] = null;
                        }
                        else{
                            int xx = (int)Math.round(xc[k]);
                            int yy = (int)Math.round(yc[k]);
                            metric_mlxy[k] = metric_centroid_safe(pixels, width, height, xx, yy, radius, times, max_drift);
                        }
                    }
                    return 1;
                });
            }

            int sum = 0;
            try{
                List<Future<Integer>> futureList=exec.invokeAll(callList);  
//                for(Future<Integer> future:futureList){  
//                    sum+=future.get();  
//                }
            }
            catch(Exception ex){ }            
            //exec.shutdown();
            
            
//            for(int i=0; i<num; i++){
//                int xx = (int)Math.round(xc[i]);
//                int yy = (int)Math.round(yc[i]);
//                metric_mlxy[i] = metric_centroid_safe(pixels, width, height, xx, yy, radius, times, max_drift);
//            }                
            return metric_mlxy;
        }
        
        public static double[][] metric_centroid_safe(double[] pixels, int width, int height, int[] xc, int[] yc, int radius)
        {
            int num = xc.length;
            double[][] metric_mlxy = new double[num][];
            for(int i=0; i<num; i++) metric_mlxy[i] = metric_centroid_safe(pixels, width, height, xc[i], yc[i], radius);                
            return metric_mlxy;
        }
        
        public static double[][] metric_centroid_safe(double[] pixels, int width, int height, int[] xc, int[] yc, int radius, int times, double max_drift)
        {
            int num = xc.length;
            double[][] metric_mlxy = new double[num][];
            for(int i=0; i<num; i++) metric_mlxy[i] = metric_centroid_safe(pixels, width, height, xc[i], yc[i], radius, times, max_drift);            
            return metric_mlxy;
        }
        
                public static double[] localization(int x, int y, double[] pixels, int width, int height, int kernel_w, double g_sigmax, double g_sigmay, int box_constrian_R, boolean dark_object){             
			double[] lm_para = {0,0,0,0};
                        double[] mlxy = null;
                        
                        LevmarSolver lms = new LevmarSolver();
                        boolean suc = lms.estimate_gaussian_fit(x,y,pixels,width,height,kernel_w,g_sigmax,box_constrian_R,dark_object,lm_para);
                        //IJ.log("optimizing lm: suc=" + suc + "	a=" + lm_para[0] + "	b=" + lm_para[1] + "	x0=" + lm_para[2] + "	y0=" + lm_para[3]);			
                        if(suc){  
                                mlxy = new double[2];
				mlxy[0] = x + lm_para[2];
				mlxy[1] = y + lm_para[3];
                        }
			return mlxy;	
		}
        
                public static double[] localization(LocalizationWindow window, double g_sigmax, double g_sigmay, int box_constrian_R, boolean dark_object){	
			double[] lm_para = {0,0,0,0};
                        double[] mlxy = null;
                        
                        LevmarSolver lms = new LevmarSolver();
                        boolean suc = lms.estimate_gaussian_fit(window,g_sigmax,box_constrian_R,dark_object,lm_para);
                        //IJ.log("optimizing lm: suc=" + suc + "	a=" + lm_para[0] + "	b=" + lm_para[1] + "	x0=" + lm_para[2] + "	y0=" + lm_para[3]);			
                        if(suc){                                        
                                mlxy = new double[2];
				mlxy[0] = window.xc + lm_para[2];
				mlxy[1] = window.yc + lm_para[3];
                        }
                        
			return mlxy;	
		}
                
                public static double[][] ML_localization(double[] pixels, int width, int height, double[] xc, double[] yc, int kernel_w, double g_sigmax, double g_sigmay, int box_constrian_R, boolean dark_object)
                {
                    double[][] mlxy = ML_localization(1, pixels, width, height, xc, yc, kernel_w, g_sigmax, g_sigmay, box_constrian_R, dark_object);                    
                    return mlxy;
                }
                
                public static double[][] ML_localization(int num_threads, double[] pixels, int width, int height, double[] xc, double[] yc, int kernel_w, double g_sigmax, double g_sigmay, int box_constrian_R, boolean dark_object)
		{			
                    if(xc==null || xc.length<1) return null; 
                    int num = xc.length;
                    
                    if(num_threads<2){                        
                        double[][] mlxy = new double[num][];
                        for(int i=0; i<num; i++){                        
                            if(Double.isNaN(xc[i]) || Double.isNaN(yc[i])){
                                mlxy[i] = null;
                            }
                            else{
                                int xx = (int)Math.round(xc[i]);
                                int yy = (int)Math.round(yc[i]);   
                                mlxy[i] = localization(xx, yy, pixels, width, height, kernel_w, g_sigmax, g_sigmay, box_constrian_R, dark_object);
                            }
                        }                
                        return mlxy;
                    }                       
                                       
                    double[][] mlxy = new double[num][];
                    
                    ExecutorService exec=Executors.newFixedThreadPool(num_threads);  
                    List<Callable<Integer>> callList=new ArrayList<>();  
                    List<Integer> list = new ArrayList<>();
                    for(int i=0; i<num; i++) list.add(i);
                    int list_size = list.size();
                    int len = list_size/num_threads;
                    if(len==0){
                        num_threads = list_size;
                        len = 1;
                    }

                    for(int i=0; i<num_threads; i++){
                        final List<Integer> sub_list;
                        if(i==num_threads-1) sub_list = list.subList(i*len, list_size);
                        else sub_list=list.subList(i*len, len*(i+1)>list.size()?list.size():len*(i+1));

                        //  
                        callList.add((Callable<Integer>) () -> {
                            for(Integer k : sub_list) {
                                if(Double.isNaN(xc[k]) || Double.isNaN(yc[k])){
                                    mlxy[k] = null;
                                }
                                else{
                                    int xx = (int)Math.round(xc[k]);
                                    int yy = (int)Math.round(yc[k]);  
                                    LocalizationWindow window = LocalizationWindow.create(xx, yy, pixels, width, height, kernel_w);
                                    mlxy[k] = localization(window, g_sigmax, g_sigmay, box_constrian_R, dark_object);
//                                    mlxy[k] = localization(xx, yy, pixels, width, height, kernel_w, g_sigmax, g_sigmay, box_constrian_R, dark_object);
                                }                                
                                
                            }
                            return 1;
                        });
                    }
                    
                    int sum = 0;
                    try{
                        List<Future<Integer>> futureList=exec.invokeAll(callList);  
//                        for(Future<Integer> future:futureList){  
//                            sum+=future.get();  
//                        }
                    }
                    catch(Exception ex){ }
//                    exec.shutdown();
                    
//                    for(int i=0; i<num; i++){                        
//                        mlxy[i] = localization(x[i], y[i], pixels, width, height, kernel_w, g_sigmax, g_sigmay, box_constrian_R, dark_object);
//                    }                
                    return mlxy;
                } 
                
                public static double[][] ML_localization(int num_threads, int[] x, int[] y, double[] pixels, int width, int height, int kernel_w, double g_sigmax, double g_sigmay, int box_constrian_R, boolean dark_object)
		{			
                    if(x==null || x.length<1) return null; 
                    
                    if(num_threads<2){
                        int num = x.length;
                        double[][] mlxy = new double[num][];
                        for(int i=0; i<num; i++){                        
                            mlxy[i] = localization(x[i], y[i], pixels, width, height, kernel_w, g_sigmax, g_sigmay, box_constrian_R, dark_object);
                        }                
                        return mlxy;
                    }                        
                                       
                    int num = x.length;                    
                    double[][] mlxy = new double[num][];
//                    LocalizationWindow[] windows = new LocalizationWindow[num];
//                    for(int i=0; i<num; i++) windows[i] = LocalizationWindow.create(x[i], y[i], pixels, width, height, kernel_w);
                                    
                    ExecutorService exec=Executors.newFixedThreadPool(num_threads);  
                    List<Callable<Integer>> callList=new ArrayList<>();  
                    List<Integer> list = new ArrayList<>();
                    for(int i=0; i<num; i++) list.add(i);
                    int list_size = list.size();
                    int len = list_size/num_threads;
                    if(len==0){
                        num_threads = list_size;
                        len = 1;
                    }

                    for(int i=0; i<num_threads; i++){
                        final List<Integer> sub_list;
                        if(i==num_threads-1) sub_list = list.subList(i*len, list_size);
                        else sub_list=list.subList(i*len, len*(i+1)>list.size()?list.size():len*(i+1));

                        //  
                        callList.add((Callable<Integer>) () -> {
                            for(Integer k : sub_list) {
                                LocalizationWindow window = LocalizationWindow.create(x[k], y[k], pixels, width, height, kernel_w);                                    
                                mlxy[k] = localization(window, g_sigmax, g_sigmay, box_constrian_R, dark_object);
                                //mlxy[k] = localization(x[k], y[k], pixels, width, height, kernel_w, g_sigmax, g_sigmay, box_constrian_R, dark_object);
                            }
                            return 1;
                        });
                    }
                    
                    int sum = 0;
                    try{
                        List<Future<Integer>> futureList=exec.invokeAll(callList);  
//                        for(Future<Integer> future:futureList){  
//                            sum+=future.get();  
//                        }
                    }
                    catch(Exception ex){ }
//                    exec.shutdown();

//                    for(int i=0; i<num; i++){                        
//                        mlxy[i] = localization(x[i], y[i], pixels, width, height, kernel_w, g_sigmax, g_sigmay, box_constrian_R, dark_object);
//                    }                
                    return mlxy;
                } 
                
                public static double[][] ML_localization(int num_threads, LocalizationWindowQueue queue, double g_sigmax, double g_sigmay, int box_constrian_R, boolean dark_object)
		{			
                    if(queue==null || queue.getCount()<1) return null; 
                    
                    int num = queue.getCount();
                    LocalizationWindow[] windows = queue.getWindows();
                    if(num_threads<2){                        
                        double[][] mlxy = new double[num][];
                        for(int i=0; i<num; i++){                        
                            mlxy[i] = localization(windows[i], g_sigmax, g_sigmay, box_constrian_R, dark_object);
                        }                
                        return mlxy;
                    }                        
                                
                    double[][] mlxy = new double[num][];
                                    
                    ExecutorService exec=Executors.newFixedThreadPool(num_threads);  
                    List<Callable<Integer>> callList=new ArrayList<>();  
                    List<Integer> list = new ArrayList<>();
                    for(int i=0; i<num; i++) list.add(i);
                    int list_size = list.size();
                    int len = list_size/num_threads;
                    if(len==0){
                        num_threads = list_size;
                        len = 1;
                    }

                    for(int i=0; i<num_threads; i++){
                        final List<Integer> sub_list;
                        if(i==num_threads-1) sub_list = list.subList(i*len, list_size);
                        else sub_list=list.subList(i*len, len*(i+1)>list.size()?list.size():len*(i+1));

                        //  
                        callList.add((Callable<Integer>) () -> {
                            for(Integer k : sub_list) {
                                if(windows[k]!=null) mlxy[k] = localization(windows[k], g_sigmax, g_sigmay, box_constrian_R, dark_object);
                                else mlxy[k] = null;
                            }
                            return 1;
                        });
                    }
                    
                    try{
                        exec.invokeAll(callList);  
                    }
                    catch(Exception ex){ }
               
                    return mlxy;
                } 

}

