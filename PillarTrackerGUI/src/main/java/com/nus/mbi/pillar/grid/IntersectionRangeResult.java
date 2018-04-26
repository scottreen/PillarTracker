package com.nus.mbi.pillar.grid;

/**
 *
 * @author xiaochun
 */
public class IntersectionRangeResult {
    private double start;
    private double end;
    private double fevl_start;
    private double fevl_end;
    private boolean suc;
    
    public IntersectionRangeResult(double start, double end, double fevl_start, double fevl_end, boolean suc){
        this.start = start;
        this.end = end;        
        this.fevl_start = fevl_start;
        this.fevl_end = fevl_end;
        this.suc = suc;
    }
        
    public double get_start(){
        return start;
    }
    
    public double get_end(){
        return end;
    }
    
    public double get_fevl_start(){
        return fevl_start;
    }
    
    public double get_fevl_end(){
        return fevl_end;
    }
       
    public boolean Sucessful(){
        return suc;
    }
}
