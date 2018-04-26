package com.nus.mbi.pillar.solver;

/**
 *
 * @author xiaochun
 */
public class CrossCorrelationNaN {
    public double[] dataA;
    public double[] dataB;
    public boolean[] nan_flags;
    
    public double sumA;
    public double sumB;
    public double sumA2;
    public double sumB2;
    public double sumAB;
    public double varA;
    public double varB;
    public double varAB;
    public int countN;
    
    public double R;

    public CrossCorrelationNaN(double[] dataA, double[] dataB, boolean[] nan_flags) {
        this.dataA = dataA;
        this.dataB = dataB;
        this.nan_flags = nan_flags;
    }
    
    public double compute_R_from_sum(){
        int counter = countN;
        if(counter>0){
            varA  = sumA2 - sumA*sumA/counter;
            varB  = sumB2 - sumB*sumB/counter;
            varAB = sumAB - sumA*sumB/counter;						                
            R = varAB*(counter-1)/(counter*Math.sqrt(varA*varB));        
        }
        else{
            varA  = 0;
            varB  = 0;
            varAB = 0;						
            R = 0;  
        }
        return R;
    }
    
    public CrossCorrelationNaN merge_overlaps(CrossCorrelationNaN cc, int shift){
        return merge_overlaps(this, cc, shift);
    }
    
    //some bugs not fixed here.. 
    public static CrossCorrelationNaN merge_overlaps(CrossCorrelationNaN cc1, CrossCorrelationNaN cc2, int shift){
        int n1 = cc1.dataA.length;
        int n2 = cc2.dataA.length;
        
        if(shift<1) return cc1;
        if(shift>n1) shift = n1;
        
        int n = n2+shift;
        double[] win_data1 = new double[n];
        double[] win_data2 = new double[n];    
        boolean[] nan_flags = new boolean[n];        
        for(int i=0; i<shift; i++){
            win_data1[i] = cc1.dataA[i];
            win_data2[i] = cc1.dataB[i];
            nan_flags[i] = cc1.nan_flags[i];
        }            
        for(int i=0; i<n2; i++){
            int off = shift+i;
            win_data1[off] = cc2.dataA[i];
            win_data2[off] = cc2.dataB[i];
            nan_flags[off] = cc2.nan_flags[i];
        }
        
        CrossCorrelationNaN cc = new CrossCorrelationNaN(win_data1, win_data2, nan_flags);
        cc.sumA  = cc1.sumA + cc2.sumA;
        cc.sumB  = cc1.sumB + cc2.sumB;
        cc.sumA2 = cc1.sumA2 + cc2.sumA2;
        cc.sumB2 = cc1.sumB2 + cc2.sumB2;
        cc.sumAB = cc1.sumAB + cc2.sumAB;
        cc.countN= cc1.countN + cc2.countN; 
        
        for(int i=shift; i<n1; i++){
            if(!cc1.nan_flags[i]){
                double a = cc1.dataA[i];
                double b = cc1.dataB[i];
                cc.sumA  -= a;
                cc.sumB  -= b;
                cc.sumA2 -= (a*a);
                cc.sumB2 -= (b*b);
                cc.sumAB -= (a*b);
                cc.countN--;
            }
        }  
        cc.compute_R_from_sum();                    
        return cc;
    }
}
