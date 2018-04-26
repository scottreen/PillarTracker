package com.nus.mbi.pillar.solver;

/**
 *
 * @author xiaochun
 */
public class MovingCrossCorrelation {
    public double[] dataA;
    public double[] dataB;
    public int k;
    
    public double[] sumA;
    public double[] sumB;
    public double[] sumA2;
    public double[] sumB2;
    public double[] sumAB;
    public double[] varA;
    public double[] varB;
    public double[] varAB;
    
    public double[] R;

    public MovingCrossCorrelation(double[] data1, double[] data2, int k) {
        this.dataA = data1;
        this.dataB = data2;
        this.k = k;
    }
    
    public double[] moving_cross_correlation(){
        int n = dataA.length;
        if(k>n) return null;
        
        sumA  = moving_sum(dataA,k);
        sumB  = moving_sum(dataB,k);
        sumA2 = moving_sum2(dataA,k);
        sumB2 = moving_sum2(dataB,k);
        sumAB = moving_dot(dataA,dataB,k);
        
        int counter = k;  
        int m = sumA.length;
        varA   = new double[m];
        varB   = new double[m];
        varAB  = new double[m];
        R      = new double[m];
        
        for(int i=0; i<m; i++){
            varA[i]  = sumA2[i] - sumA[i]*sumA[i]/counter;
            varB[i]  = sumB2[i] - sumB[i]*sumB[i]/counter;
            varAB[i] = sumAB[i] - sumA[i]*sumB[i]/counter;						
            R[i] = varAB[i]*(counter-1)/(counter*Math.sqrt(varA[i]*varB[i]));        
        }

        return R;
    }
    
    public static double cross_correlation(double[] a, double[]b){
        int n = a.length;
        double sumI = 0;
        double sumK = 0;
        double sumI2 = 0;
        double sumK2 = 0;
        double sumIK = 0;
        int counter = n;

        for(int i=0; i<n; i++){
            double mi = a[i];
            double ki = b[i];
            sumI += mi;
            sumK += ki;
            sumIK += mi*ki;
            sumI2 += mi*mi;
            sumK2 += ki*ki;
        }

        double std_img = sumI2 - sumI*sumI/counter;
        double std_psf = sumK2 - sumK*sumK/counter;
        double sum_img_psf = sumIK-sumI*sumK/counter;						
        double cc = sum_img_psf*(counter-1)/(counter*Math.sqrt(std_img*std_psf));
        return cc;    
    }
    
    public static double[] moving_cross_correlation(double[] a, double[]b, int k){
        int n = a.length;
        if(k>n) return null;
        
        double[] sumI  = moving_sum(a,k);
        double[] sumK  = moving_sum(b,k);
        double[] sumI2 = moving_sum2(a,k);
        double[] sumK2 = moving_sum2(b,k);
        double[] sumIK = moving_dot(a,b,k);
        
        int counter = k;  
        int m = sumI.length;
        double[] cc = new double[m];
        for(int i=0; i<m; i++){
            double std_img = sumI2[i] - sumI[i]*sumI[i]/counter;
            double std_psf = sumK2[i] - sumK[i]*sumK[i]/counter;
            double sum_img_psf = sumIK[i]-sumI[i]*sumK[i]/counter;						
            cc[i] = sum_img_psf*(counter-1)/(counter*Math.sqrt(std_img*std_psf));        
        }

        return cc;
    }
    
    public static double[] moving_sum(double[] data, int k){
        int n = data.length;
        if(k>n) return null;
        int m = n-k+1;
        double[] rst = new double[m];
        double sum = 0;
        for(int i=0; i<k; i++) sum += data[i];
        rst[0] = sum;
        for(int i=1; i<m; i++){
            sum = sum-data[i-1]+data[k+i-1];
            rst[i] = sum;
        }
        return rst;
    }
    
    public static double[] moving_avg(double[] data, int k){
        double[] sum = moving_sum(data,k);
        if(sum==null) return null;        
        for(int i=0; i<sum.length; i++) sum[i] /= k;        
        return sum;
    }
    
    public static double[] moving_sum2(double[] data, int k){
        int n = data.length;
        double[] data2 = new double[n];
        for(int i=0; i<n; i++) data2[i] = data[i]*data[i];
        return moving_sum(data2,k);
    }
    
    public static double[] moving_dot(double[] data1, double[] data2, int k){
        int n = data1.length;
        double[] data = new double[n];
        for(int i=0; i<n; i++) data[i] = data1[i]*data2[i];
        return moving_sum(data,k);
    }
    
    
    
}
