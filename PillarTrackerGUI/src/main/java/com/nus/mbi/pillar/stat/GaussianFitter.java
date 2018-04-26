package com.nus.mbi.pillar.stat;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.ProfilePlot;
import ij.measure.CurveFitter;
import ij.plugin.frame.Fitter;

/**
 *
 * @author xiaochun
 */
public class GaussianFitter {
    public static double[] do_fit(){               
        int n = WindowManager.getImageCount();
        if(n<1){
            IJ.showMessage("there is no image opened!");
            return null;
        }
        ImagePlus ip = IJ.getImage();
        ProfilePlot profile = new ProfilePlot(ip);
        double[] dataY = profile.getProfile();
        if(dataY==null) return null;
        int len = dataY.length;
        if(len<=0) return null;
        
        double[] dataX = new double[len];
        double miny, maxy;
        miny = maxy = dataY[0];            
        for(int i=1; i<len; i++){                
            double y = dataY[i];                                
            if(miny>y) miny = y;
            if(maxy<y) maxy = y;

        }

        if(miny<0){
            for(int i=0; i<len; i++) dataY[i] -= miny;                 
            maxy -= miny;
            miny = 0;                
        }

        double m0=0;
        double m1=0; 
        for(int i=0; i<len; i++){
            dataX[i] = i;
            double y = dataY[i];
            m0 += y;
            m1 += i*y;                       
        }          

        double halfy = (maxy-miny)/2 + miny;
        double halfy_index1 = 0;
        double halfy_index2 = len-1;
        boolean forward_pos = dataY[0]<halfy;
        boolean backward_pos = dataY[len-1]<halfy;
        for(int i=1; i<len; i++){                
            double y = dataY[i];                
            if(y>halfy && forward_pos){
                halfy_index1 = i;
                break;
            }
            else if(y<halfy && !forward_pos){
                halfy_index1 = i;
                break;
            }
        }

        for(int i=len-2; i>=0; i--){                
            double y = dataY[i];                
            if(y>halfy && backward_pos){
                halfy_index2 = i;
                break;
            }
            else if(y<halfy && !backward_pos){
                halfy_index2 = i;
                break;
            }
        }

        //compute Gaussian sigma from FWHM: https://en.wikipedia.org/wiki/Full_width_at_half_maximum
        double fwhm = halfy_index2 - halfy_index1;
        boolean suc_sigma = fwhm>0;
        double sigma = suc_sigma? fwhm/2.355 : len/2.0;

        double centroid = m1/m0;
        int center = (int)Math.round(centroid);
        boolean dark = dataY[center]*2 < (miny + maxy);

        double a = dark ? maxy : miny;
        double b = dark ? miny : maxy;
        double c = centroid;
        double d = sigma;
        double[] init = {a, b, c, d};
        IJ.log("Initial guess:\n dark?=" + dark + "\n a=" + a + "\n b=" + b + "\n c=" + c + "\n d=" + d);
        CurveFitter cf = new CurveFitter(dataX, dataY);
        cf.setInitialParameters(init);            
        cf.doFit(CurveFitter.GAUSSIAN, false);   
        
        //double[] rst = .scaled() ? null : cf.getParams();
        double[] rst = cf.getParams();
        Fitter.plot(cf);
        IJ.log(cf.getResultString());        
        //Plot plot = profile.getPlot();
        
        return rst;
    }
    
    public static double[] do_fit(double[] dataY){
        int len = dataY.length;
        double[] dataX = new double[len];
        double miny, maxy;
        miny = maxy = dataY[0];            
        for(int i=1; i<len; i++){                
            double y = dataY[i];                                
            if(miny>y) miny = y;
            if(maxy<y) maxy = y;
        }

        if(miny<0){
            for(int i=0; i<len; i++) dataY[i] -= miny;                 
            maxy -= miny;
            miny = 0;                
        }

        double m0=0;
        double m1=0; 
        for(int i=0; i<len; i++){
            dataX[i] = i;
            double y = dataY[i];
            m0 += y;
            m1 += i*y;                       
        }          

        double halfy = (maxy-miny)/2 + miny;
        double halfy_index1 = 0;
        double halfy_index2 = len-1;
        boolean forward_pos = dataY[0]<halfy;
        boolean backward_pos = dataY[len-1]<halfy;
        for(int i=1; i<len; i++){                
            double y = dataY[i];                
            if(y>halfy && forward_pos){
                halfy_index1 = i;
                break;
            }
            else if(y<halfy && !forward_pos){
                halfy_index1 = i;
                break;
            }
        }

        for(int i=len-2; i>=0; i--){                
            double y = dataY[i];                
            if(y>halfy && backward_pos){
                halfy_index2 = i;
                break;
            }
            else if(y<halfy && !backward_pos){
                halfy_index2 = i;
                break;
            }
        }

        //compute Gaussian sigma from FWHM: https://en.wikipedia.org/wiki/Full_width_at_half_maximum
        double fwhm = halfy_index2 - halfy_index1;
        boolean suc_sigma = fwhm>0;
        double sigma = suc_sigma? fwhm/2.355 : len/2.0;

        double centroid = m1/m0;
        int center = (int)Math.round(centroid);
        boolean dark = dataY[center]*2 < (miny + maxy);

        double a = dark ? maxy : miny;
        double b = dark ? miny : maxy;
        double c = centroid;
        double d = sigma;
        double[] init = {a, b, c, d};
        IJ.log("Initial guess:\n dark?=" + dark + "\n a=" + a + "\n b=" + b + "\n c=" + c + "\n d=" + d);
        CurveFitter cf = new CurveFitter(dataX, dataY);
        cf.setInitialParameters(init);            
        cf.doFit(CurveFitter.GAUSSIAN, false);   
        
        //double[] rst = .scaled() ? null : cf.getParams();
        double[] rst = cf.getParams();
        
        return rst;
    }
}
