package com.nus.mbi.pillar.drift;

/**
 *
 * @author xiaochun
 */
public class VarianceReduction {
    double[][] data;
    boolean[] flags;
    
    double[] sum_n;
    double[] sum_m;
    double[] sum_diff2;
    double sum_nm;
    int count_m; 
    
    double[] var;
    
    private final int nrois;
    private final int sliceN;
    
    public double[] getVAR(){
        return var;
    }
    
    public double[] getMeanM(){
        //int nframes = sum_m.length;
        double[] mean = new double[sliceN];
        for(int i=0; i<sliceN; i++) mean[i] = sum_m[i]/count_m;
        return mean;
    }    
    
    public int getCountM(){
        return count_m;
    }
    
    public VarianceReduction(double[][] data, boolean[] flags) {
        this.data = data;
        this.flags = flags;
        
        nrois = data.length;
        sliceN = data[0].length;
    }
    
    public void start_var(){
        sum_n = sumN();        
        sum_m = sumM();
        sum_diff2 = sumDiff2();
        sum_nm = sumNM();
        count_m= countM();
        
        compute_var();
    }
    
    public void delete(int k){
        //int nrois = flags.length;
        if(k<0 || k>=nrois) return;
        if(false==flags[k]) return;
        
        flags[k] = false;
        count_m--;
        sum_nm -= sum_n[k];
        
        //int sliceN = data[0].length;
        for(int i=0; i<sliceN; i++) sum_m[i] -= data[k][i];                
        sum_diff2 = sumDiff2();
        
        compute_var();
    }
    
    private void compute_var(){
        //int nrois = data.length;
        //int sliceN = data[0].length;
        if(var==null) var = new double[nrois];
        int M = count_m;        
        int N = sliceN;
        int N2 = N*N;
        double mean_nm = sum_nm/M;        
        for(int j=0; j<nrois; j++) if(flags[j]){
            double sum_diff = sum_n[j]-mean_nm;
            double E1 = sum_diff2[j]/N;
            double E2 = sum_diff*sum_diff/N2;
            var[j] = E1-E2;
        }
        //return var;
    }
    
    private int countM(){
        //int nrois = flags.length;			
        int count = 0;
        for(int j=0; j<nrois; j++){            
            if(flags[j]) count++;
        }
        return count;
    }
    
    private double[] sumN(){
        //int nrois = data.length;			
        //int sliceN = data[0].length;		
        double[] sum = new double[nrois];
        for(int j=0; j<nrois; j++){
            sum[j] = 0;
            if(flags[j]) for(int i=0; i<sliceN; i++) sum[j] += data[j][i];            
        }
        return sum;
    }
    
    private double[] sumN2(){
        //int nrois = data.length;			
        //int sliceN = data[0].length;		
        double[] sum = new double[nrois];
        for(int j=0; j<nrois; j++){
            sum[j] = 0;
            if(flags[j]) for(int i=0; i<sliceN; i++){
                double t = data[j][i];
                sum[j] += (t*t);
            }            
        }
        return sum;
    }
    
    private double[] sumM(){
        //int nrois = data.length;			
        //int sliceN = data[0].length;		
        double[] sum = new double[sliceN];
        for(int i=0; i<sliceN; i++){
            sum[i] = 0;
            for(int j=0; j<nrois; j++) if(flags[j]) sum[i] += data[j][i];            
        }
        return sum;
    }
    
    private double sumNM(){
        //int nrois = data.length;			
        //int sliceN = data[0].length;		
        double sum = 0;
        for(int j=0; j<nrois; j++){
            //if(flags[j]) for(int i=0; i<sliceN; i++) sum += data[j][i];            
            if(flags[j]) sum += sum_n[j];            
        }
        return sum;
    }
    
    private double[] sumDiff2(){
        //int nrois = data.length;			
        //int sliceN = data[0].length;		
        double[] sum = new double[nrois];
        for(int j=0; j<nrois; j++){
            sum[j]=0;
            if(flags[j]){
                for(int i=0; i<sliceN; i++){
                    double diff = data[j][i]-sum_m[i]/count_m;
                    sum[j] += (diff*diff);
                }
            }            
        }
        return sum;
    }  
}
