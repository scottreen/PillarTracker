package com.nus.mbi.pillar.drift;

/**
 *
 * @author xiaochun
 */
public class VarianceMatrix {
    double[][] data;
    boolean[] flags;
    
    private final int nrois;
    private final int sliceN;
    
    double[] sum_n;
    double[] sum_n2;
    double[] var;    
    public double[] getVAR(){
        return var;
    }
    
    public VarianceMatrix(double[][] data, boolean[] flags) {
        this.data = data;
        this.flags = flags;
        nrois = data.length;
        sliceN = data[0].length;
    }
    
    public void compute_var(){           
        if(var==null) var = new double[nrois];
        if(sum_n==null) sum_n = new double[nrois];
        if(sum_n2==null) sum_n2 = new double[nrois];
        
        int N = sliceN;
        int N2 = N*N;
              
        for(int j=0; j<nrois; j++) if(flags[j]){
            double sumN = 0;            
            double sumN2= 0;            
            for(int i=0; i<sliceN; i++){
                double t = data[j][i];
                sumN += t;
                sumN2 += (t*t);
            }     

            sum_n[j] = sumN;
            sum_n2[j]= sumN2;
            var[j] = sumN2/N-sumN*sumN/N2;
        }       
    }
}
