package com.nus.mbi.pillar.grid;

/**
 *
 * @author xiaochun
 */
public class EvaluateRangeResult {
    private double fs;
    private double fe;
    private boolean suc;
    
    public EvaluateRangeResult(double fevl_start, double fevl_end, boolean suc){
        this.fs = fevl_start;
        this.fe = fevl_end;
        this.suc = suc;
    }
    
    public double get_fevl_start(){
        return fs;
    }
    
    public double get_fevl_end(){
        return fe;
    }
       
    public boolean Sucessful(){
        return suc;
    }
}
