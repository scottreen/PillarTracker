package com.nus.mbi.pillar.stat;

/**
 *
 * @author xiaochun
 */
public class BasicStatisitic {
    	
    public static double max(double[] values) {
            int len = values.length;
            double max = Double.NaN;//values[0];
            for(int k=0; k<len; k++) {                	
                double v = values[k];
                    if(!Double.isNaN(v)){
                    if(Double.isNaN(max))  max =v;
                    else if(v > max) max =v;
                    }
            }
            return max;
    }

    public static double min(double[] values) {
                    int len = values.length;
            double min = Double.NaN;//values[0];
            for(int k=0; k<len; k++){		                	

                    double v = values[k];
                    if(!Double.isNaN(v)){                	
                    if(Double.isNaN(min)) min = v;
                    else if(v < min) min = v;
                    }
            }
            return min;
    } 

    //support nan double
    public static double avg(double[] values) {
                    int len = values.length;
                    int num = 0;
                double avg = 0;//values[0];
                for(int k=0; k<len; k++) {                	
                    double v = values[k];
                    if(!Double.isNaN(v)){
                    avg += v;
                    num++;
                    }
                }

                avg = num>0? avg/num : Double.NaN;
                return avg;
    } 

    public static double avg(double[] values, int start_index, int end_index) {
                    int len = values.length;
                    if(start_index<0) start_index = 0;
                    if(end_index>len) end_index = len;
                    int num = 0;
                double avg = 0;//values[0];
                for(int k=start_index; k<end_index; k++) {                	
                    double v = values[k];
                    if(!Double.isNaN(v)){
                    avg += v;
                    num++;
                    }
                }

                avg = num>0? avg/num : Double.NaN;
                return avg;
    } 
    
    public static double mean_sqr(double[] values)
    {
                    int len = values.length;                    
                    int num = 0;
                    double var = 0;
                    for(int k=0; k<len; k++){
                            double v = values[k];
                            if(!Double.isNaN(v)){				                                    
                                    var += v*v;
                                    num++;
                            }
                    }

                    var = num>0? var/num : Double.NaN;
                    return var;	
    }

    public static double var(double[] values)
    {
                    int len = values.length;
                    double a = avg(values);
                    int num = 0;
                    double var = 0;
                    for(int k=0; k<len; k++){
                            double v = values[k];
                            if(!Double.isNaN(v)){				
                                    double da = v-a;
                                    var += da*da;
                                    num++;
                            }
                    }

                    var = num>0? var/num : Double.NaN;
                    return var;	
    }
    
    public static double var(double[] values, boolean[] flags)
    {
                    int len = values.length;
                    double a = avg(values);
                    int num = 0;
                    double var = 0;
                    for(int k=0; k<len; k++){
                            double v = values[k];
                            if(flags[k] && !Double.isNaN(v)){				
                                    double da = v-a;
                                    var += da*da;
                                    num++;
                            }
                    }

                    var = num>0? var/num : Double.NaN;
                    return var;	
    }

    public static double avg(double[] values, boolean[] flags) {
            int len = values.length;	
            int size = 0;
            double avg = 0;	
            for(int i=0; i<len; i++){
                    double v = values[i];
                    if(flags[i] && !Double.isNaN(v)){
                            avg += v;
                            size++;
                    }
            }
            avg = size>0 ? avg/size : Double.NaN;	
            return avg;
    }     
    
    public static boolean IsNotNaN(double a, double b, double c, double d){
            if(Double.isNaN(a)) return false;
            if(Double.isNaN(b)) return false;
            if(Double.isNaN(c)) return false;
            if(Double.isNaN(d)) return false;
            return true;
    }
    
    
    public static int find_min(int[] array, boolean[] active){
        int min=Integer.MAX_VALUE;
        int n = array.length;
        for(int i=0; i<n; i++){
            int a = array[i];
            if(a!=Integer.MAX_VALUE && a!=Integer.MIN_VALUE){
                if(active[i]){
                    if(min>a) min = a;                    
                }
            }
        }        
        return min;        
    }
    
    public static int find_max(int[] array, boolean[] active){
        int max=Integer.MIN_VALUE;
        int n = array.length;
        for(int i=0; i<n; i++){
            int a = array[i];
            if(a!=Integer.MAX_VALUE && a!=Integer.MIN_VALUE){
                if(active[i]){
                    if(max<a) max = a;                    
                }
            }
        }        
        return max;        
    }
    
    public static int find_min(int[] array){
        int min=Integer.MAX_VALUE;
        int n = array.length;
        for(int i=0; i<n; i++){            
            int a = array[i];
            if(a!=Integer.MAX_VALUE && a!=Integer.MIN_VALUE){
                if(min>a) min = a;
            }            
        }        
        return min;        
    }
    
    public static int find_max(int[] array){
        int max=Integer.MIN_VALUE;
        int n = array.length;
        for(int i=0; i<n; i++){            
            int a = array[i];
            if(a!=Integer.MAX_VALUE && a!=Integer.MIN_VALUE){
                if(max<a) max = a;
            }            
        }        
        return max;        
    }
}
