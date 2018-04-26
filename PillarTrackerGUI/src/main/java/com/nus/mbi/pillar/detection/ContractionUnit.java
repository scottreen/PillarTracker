package com.nus.mbi.pillar.detection;

import ij.IJ;
import ij.measure.ResultsTable;
import java.util.List;

/**
 *
 * @author xiaochun
 */
public class ContractionUnit {
    public int i1;
    public int i2;
    public int f1;
    public int f2;
    public double rs;
    public double rp;    
    
    public static int SIZE = Integer.SIZE*4 + Double.SIZE*2;
    public static int BYTES = Integer.BYTES*4 + Double.BYTES*2;
    
    public ContractionUnit(int i1, int i2, int f1, int f2, double rs, double rp) {
        this.i1 = i1;
        this.i2 = i2;
        this.f1 = f1;
        this.f2 = f2;
        this.rs = rs;
        this.rp = rp;
    }
    
    public static void print_list(List<ContractionUnit> cu_list){
        int n = cu_list.size();
        if(n>0){
            IJ.log("index1\t index2\t start_frame\t end_frame\t R_perpendicular\t R_orthogonal");
            for(int i=0; i<n; i++){
                ContractionUnit cu = cu_list.get(i);
                IJ.log((cu.i1+1)+"\t "+(cu.i2+1)+"\t "+(cu.f1+1)+"\t "+(cu.f2+1)+"\t "+cu.rp+"\t "+cu.rs);
            }
            IJ.log("List of Pillar Index");
            IJ.log("index1 index2");
            for(int i=0; i<n; i++){
                ContractionUnit cu = cu_list.get(i);
                IJ.log((cu.i1+1)+" "+(cu.i2+1));
            }
        }
        else IJ.log("the contractile unit(CU) list is empty");
    }
    
    public static void show_table_list(List<ContractionUnit> cu_list){
        int n = cu_list.size();
        if(n>0){
            ResultsTable rt = new ResultsTable();        
            rt.setNaNEmptyCells(true);
            rt.showRowNumbers(false);
            //rt.incrementCounter();
            for(int i=0; i<n; i++){
                ContractionUnit cu = cu_list.get(i);
                rt.incrementCounter();                         
                rt.addValue("Index 1",          cu.i1+1);
                rt.addValue("Index 2",          cu.i2+1);
                rt.addValue("Start Frame",      cu.f1+1);
                rt.addValue("End Frame",        cu.f2+1);
                rt.addValue("R Perpendicular",  cu.rp);  
                rt.addValue("R Orthogonal",     cu.rs);                               
            }
            rt.show("Contractile Units");            
        }
    }
    
    public static int[] get_pillar_index_list(List<ContractionUnit> cu_list, int length){
        int n = cu_list.size();
        if(n<1) return null;
        boolean[] flag = new boolean[length];
        for(int i=0; i<length; i++) flag[i] = false;
        
        for(int i=0; i<n; i++){
            ContractionUnit cu = cu_list.get(i);
            int i1=cu.i1;
            int i2=cu.i2;
            flag[i1] = true;
            flag[i2] = true;
        }      
        
        int num = 0;
        for(int i=0; i<length; i++){
            if(flag[i]) num++;
        }
        if(num<1) return null;
        
        int[] list = new int[num];
        num = 0;
        //IJ.log("List of Pillar Index");
        for(int i=0; i<length; i++){
            if(flag[i] && i>=0){                
                //IJ.log(""+(i+1));
                list[num] = i;                
                num++;
            }
        }
        return list;
    }
    
    public static int[] get_pillar_index_list(List<ContractionUnit> cu_list, int frame, int length){
        int n = cu_list.size();
        if(n<1) return null;
        boolean[] flag = new boolean[length];
        for(int i=0; i<length; i++) flag[i] = false;
        
        for(int i=0; i<n; i++){
            ContractionUnit cu = cu_list.get(i);
            if(frame>=cu.f1 && frame<=cu.f2){
                int i1=cu.i1;
                int i2=cu.i2;
                flag[i1] = true;
                flag[i2] = true;
            }
        }      
        
        int num = 0;
        for(int i=0; i<length; i++){
            if(flag[i]) num++;
        }
        if(num<1) return null;
        
        int[] list = new int[num];
        num = 0;
        //IJ.log("List of Pillar Index");
        for(int i=0; i<length; i++){
            if(flag[i] && i>=0){                
                //IJ.log(""+(i+1));
                list[num] = i;                
                num++;
            }
        }
        return list;
    }
}
