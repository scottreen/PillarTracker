package com.nus.mbi.pillar.stat;

import ij.measure.CurveFitter;

/**
 *
 * @author xiaochun
 */
public class PloyFitter {
    private CurveFitter cf = null;
    public static int type = CurveFitter.POLY2;
    
    public PloyFitter(double[] x, double[] y){
        cf = new CurveFitter(x,y);
    }
    
    public void doFit(){
        cf.doFit(type);
    }
    
    public double[] getParameters(){
        return cf.getParams();
    }
    
    public int getStatus(){
        return cf.getStatus();
    }
    
    public String getResultString(){
        return cf.getResultString();
    }
    
    public static double f(double[] p, double x){
        return CurveFitter.f(type, p, x);
    }
}
