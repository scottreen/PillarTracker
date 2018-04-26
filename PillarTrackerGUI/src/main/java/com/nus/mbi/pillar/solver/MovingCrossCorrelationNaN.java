package com.nus.mbi.pillar.solver;

/**
 *
 * @author xiaochun
 */
public class MovingCrossCorrelationNaN extends MovingCrossCorrelation{
    public boolean[] nan_flags;
    public int[] countN;

    public MovingCrossCorrelationNaN(boolean[] nan_flags, double[] data1, double[] data2, int k) {
        super(data1, data2, k);
        this.nan_flags = nan_flags;
    }
    
    @Override
    public double[] moving_cross_correlation(){
        int n = dataA.length;
        if(k>n) return null;
        
        sumA  = moving_sum(dataA,nan_flags,k);
        sumB  = moving_sum(dataB,nan_flags,k);
        sumA2 = moving_sum2(dataA,nan_flags,k);
        sumB2 = moving_sum2(dataB,nan_flags,k);
        sumAB = moving_dot(dataA,dataB,nan_flags,k);        
        countN= moving_count(nan_flags,k);  
        
        int m = sumA.length;
        varA   = new double[m];
        varB   = new double[m];
        varAB  = new double[m];
        R      = new double[m];
        
        for(int i=0; i<m; i++){
            int counter = countN[i];
            if(counter>0){
                varA[i]  = sumA2[i] - sumA[i]*sumA[i]/counter;
                varB[i]  = sumB2[i] - sumB[i]*sumB[i]/counter;
                varAB[i] = sumAB[i] - sumA[i]*sumB[i]/counter;						
                R[i] = varAB[i]*(counter-1)/(counter*Math.sqrt(varA[i]*varB[i]));        
                //R[i] = varAB[i]/(Math.sqrt(varA[i]*varB[i]));        
            }
            else{
                varA[i]  = 0;
                varB[i]  = 0;
                varAB[i] = 0;						
                R[i] = 0;  
            }
        }

        return R;
    }
    
    public CrossCorrelationNaN get(int i){        
        double[] win_data1 = new double[k];
        double[] win_data2 = new double[k];
        boolean[] win_nan_flags = new boolean[k];
        for(int c=0; c<k; c++){
            int shift = i+c;
            win_data1[c] = dataA[shift];
            win_data2[c] = dataB[shift];
            win_nan_flags[c] = nan_flags[shift];
        }
        
        CrossCorrelationNaN cc = new CrossCorrelationNaN(win_data1, win_data2, win_nan_flags);
        cc.sumA  = sumA[i];
        cc.sumB  = sumB[i];
        cc.sumA2 = sumA2[i];
        cc.sumB2 = sumA2[i];
        cc.sumAB = sumAB[i];       
        cc.countN= countN[i];
        cc.varA = varA[i];
        cc.varB = varB[i];
        cc.varAB = varAB[i];
        cc.R = R[i];
        return cc;
    }
    
    public static double cross_correlation(double[] a, double[]b, boolean[] nan_flag){
        int n = a.length;
        double sumI = 0;
        double sumK = 0;
        double sumI2 = 0;
        double sumK2 = 0;
        double sumIK = 0;
        int counter = 0;

        for(int i=0; i<n; i++){
            if(!nan_flag[i]){
                double mi = a[i];
                double ki = b[i];
                sumI += mi;
                sumK += ki;
                sumIK += mi*ki;
                sumI2 += mi*mi;
                sumK2 += ki*ki;
                counter++;
            }
        }
        
        double cc = 0;
        if(counter>0){
            double std_img = sumI2 - sumI*sumI/counter;
            double std_psf = sumK2 - sumK*sumK/counter;
            double sum_img_psf = sumIK-sumI*sumK/counter;						
            cc = sum_img_psf*(counter-1)/(counter*Math.sqrt(std_img*std_psf));
            //cc = sum_img_psf/Math.sqrt(std_img*std_psf);
        }
        
        return cc;    
    }
    
    public static int[] moving_count(boolean[] nan_flags, int k){
        int n = nan_flags.length;
        if(k>n) return null;
        int m = n-k+1;
        int[] rst = new int[m];
        int sum = 0;
        for(int i=0; i<k; i++){
            if(!nan_flags[i]) sum++;
        }
        rst[0] = sum;
        for(int i=1; i<m; i++){
            int first = i-1;
            int last  = first + k;
            if(!nan_flags[first]) sum--;
            if(!nan_flags[last])  sum++;
            rst[i] = sum;
        }
        return rst;
    }
    
    public static double[] moving_sum(double[] data, boolean[] nan_flags, int k){
        int n = data.length;
        if(k>n) return null;
        int m = n-k+1;
        double[] rst = new double[m];
        double sum = 0;
        for(int i=0; i<k; i++){
            if(!nan_flags[i]) sum += data[i];
        }
        rst[0] = sum;
        for(int i=1; i<m; i++){
            int first = i-1;
            int last  = first + k;
            if(!nan_flags[first]) sum -= data[first];
            if(!nan_flags[last])  sum += data[last];
            rst[i] = sum;
        }
        return rst;
    }
    
    public static double[] moving_sum2(double[] data, boolean[] nan_flags, int k){
        int n = data.length;
        double[] data2 = new double[n];
        for(int i=0; i<n; i++) data2[i] = data[i]*data[i];
        return moving_sum(data2,nan_flags,k);
    }
    
    public static double[] moving_dot(double[] data1, double[] data2, boolean[] nan_flags, int k){
        int n = data1.length;
        double[] data = new double[n];
        for(int i=0; i<n; i++) data[i] = data1[i]*data2[i];
        return moving_sum(data,nan_flags,k);
    }
}
