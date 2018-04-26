package com.nus.mbi.pillar.drift;

import ij.IJ;
import com.nus.mbi.pillar.solver.Matrix_Util;
import com.nus.mbi.pillar.stat.BasicStatisitic;

/**
 *
 * @author xiaochun
 */
public class DriftCorrection {
    
    public static void search_slient_pillars(double[][]cx, double[][]cy, boolean[] best_flags){
                    int nrois = cx.length;				
                    int sliceN = cx[0].length;			
                    boolean[] flags = new boolean[nrois];
                    //best_flags = new boolean[nrois];
                    int nrois_checked = 0; 		
                    for(int k=0; k<nrois; k++){			
                            flags[k] = true;			
                            for(int i=0; i<sliceN; i++){	
                                    if(flags[k]){
                                            if(Double.isNaN(cx[k][i]) || Double.isNaN(cy[k][i])) flags[k] = false;
                                    }
                                    else break;
                            }			
                            best_flags[k] = flags[k];
                            if(flags[k]) nrois_checked++;
                    }
                    IJ.log("total pillars:" + nrois + "\r\nremaining pillars:" + nrois_checked);
                    
                    boolean[] flags_x = new boolean[nrois];
                    boolean[] flags_y = new boolean[nrois];
                    for(int k=0; k<nrois; k++){
                        flags_x[k]=flags[k];
                        flags_y[k]=flags[k];
                    }
                    
                    VarianceMatrix varX = new VarianceMatrix(cx, flags_x);
                    VarianceMatrix varY = new VarianceMatrix(cy, flags_y);
                    varX.compute_var();
                    varY.compute_var();
                    
                    double[] varsum_array = new double[nrois];
                    double min_varsum = -1;		
                    int min_varsum_k = -1;		
                    int nleft = nrois_checked;
                    boolean find_turning_point = false;
                    int checktimes = 3;
                    for(int deletetimes=0; deletetimes<nrois_checked-1; deletetimes++){
                                    double varsum = 0;					
                                    double maxvarxy = -1;
                                    int maxvarxy_k = -1;
                                    for(int k=0; k<nrois; k++) if(flags[k]){
                                            double varxy = varX.getVAR()[k] + varY.getVAR()[k];
                                            varsum += varxy;
                                            if(maxvarxy<0 || maxvarxy<varxy){
                                                    maxvarxy = varxy;
                                                    maxvarxy_k = k;
                                            }
                                    }

                                    varsum = varsum/nleft/nleft;				
                                    varsum_array[nleft-1] = varsum;				

                                    if(!find_turning_point){
                                            int nondecreasing_times = 0;				
                                            for(int j=0; j<checktimes; j++){
                                                    int current_k = nleft-1+j;
                                                    int last_k = current_k+1;
                                                    if(last_k<nrois_checked){						
                                                            if(varsum_array[current_k]>=varsum_array[last_k]) nondecreasing_times++;													
                                                            else break;
                                                    }
                                            }			

                                            if(nondecreasing_times == checktimes){ //found the first turning point				
                                                    min_varsum = varsum;
                                                    min_varsum_k = nleft;                                                    	
                                                    find_turning_point = true;															
                                            }					
                                    }

                                    if(find_turning_point){
                                            IJ.showStatus("found the first turning point at" + nleft);
                                            IJ.showProgress(1,1);
                                            break;		
                                    }
                                    
                                    flags[maxvarxy_k] = false; //deleting
                                    nleft--;	
                                    IJ.showStatus("searching the reference pillars " + deletetimes + "/" + (nrois_checked-1));		
                                    IJ.showProgress(deletetimes+1,nrois_checked-1);
                    }		
                    
                    if(find_turning_point && nleft>0){
                        for(int k=0; k<nrois; k++) best_flags[k] = flags[k];
                        IJ.log("minimum std is " + Math.sqrt(min_varsum) + " at " + min_varsum_k);			
                    }
                    else{
                        IJ.log("can't found the reference pillars" + " at " + min_varsum_k);	
                    }
                    
    }
    
    public static double[] get_driftXY_reduced_std(double[][]cx, double[][]cy, boolean[] best_flags, double[]varsum_array){
                    int nrois = cx.length;				
                    int sliceN = cx[0].length;			
                    boolean[] flags = new boolean[nrois];
                    //best_flags = new boolean[nrois];
                    int nrois_checked = 0; 		
                    for(int k=0; k<nrois; k++){			
                            flags[k] = true;			
                            for(int i=0; i<sliceN; i++){	
                                    if(flags[k]){
                                            if(Double.isNaN(cx[k][i]) || Double.isNaN(cy[k][i])) flags[k] = false;
                                    }
                                    else break;
                            }			
                            best_flags[k] = flags[k];
                            if(flags[k]) nrois_checked++;
                    }
                    IJ.log("total pillars:" + nrois + "\r\nremaining pillars:" + nrois_checked);
                    
                    boolean[] flags_x = new boolean[nrois];
                    boolean[] flags_y = new boolean[nrois];
                    for(int k=0; k<nrois; k++){
                        flags_x[k]=flags[k];
                        flags_y[k]=flags[k];
                    }
                    
                    VarianceReduction varX = new VarianceReduction(cx, flags_x);
                    VarianceReduction varY = new VarianceReduction(cy, flags_y);
                    varX.start_var();
                    varY.start_var();
                    
                    double min_varsum = -1;		
                    int min_varsum_k = -1;		
                    int nleft = nrois_checked;
                    boolean find_turning_point = false;
                    int checktimes = 3;
                    for(int deletetimes=0; deletetimes<nrois_checked-1; deletetimes++){
                                    double varsum = 0;					
                                    double maxvarxy = -1;
                                    int maxvarxy_k = -1;
                                    for(int k=0; k<nrois; k++) if(flags[k]){
                                            double varxy = varX.getVAR()[k] + varY.getVAR()[k];
                                            varsum += varxy;
                                            if(maxvarxy<0 || maxvarxy<varxy){
                                                    maxvarxy = varxy;
                                                    maxvarxy_k = k;
                                            }
                                    }

                                    varsum = varsum/nleft/nleft;				
                                    varsum_array[nleft-1] = varsum;				

                                    if(!find_turning_point){
                                            int nondecreasing_times = 0;				
                                            for(int j=0; j<checktimes; j++){
                                                    int current_k = nleft-1+j;
                                                    int last_k = current_k+1;
                                                    if(last_k<nrois_checked){						
                                                            if(varsum_array[current_k]>=varsum_array[last_k]) nondecreasing_times++;													
                                                            else break;
                                                    }
                                            }			

                                            if(nondecreasing_times == checktimes){ //found the first turning point				
                                                    min_varsum = varsum;
                                                    min_varsum_k = nleft;
                                                    for(int k=0; k<nrois; k++) best_flags[k] = flags[k];	
                                                    find_turning_point = true;															
                                            }					
                                    }

                                    if(find_turning_point){
                                            IJ.showStatus("found the first turning point at" + nleft);
                                            IJ.showProgress(1,1);
                                            break;		
                                    }
                                    
                                    varX.delete(maxvarxy_k);
                                    varY.delete(maxvarxy_k);
                                    flags[maxvarxy_k] = false; //deleting
                                    nleft--;	
                                    IJ.showStatus("searching the reference pillars " + deletetimes + "/" + (nrois_checked-1));		
                                    IJ.showProgress(deletetimes+1,nrois_checked-1);
                    }		
                    IJ.log("minimum std is " + Math.sqrt(min_varsum) + " at " + min_varsum_k);			
                    
                    double[] driftXY = null;
                    if(find_turning_point && nleft>0){
                        driftXY = new double[2*sliceN];
                        double[] mx = varX.getMeanM();
                        double[] my = varY.getMeanM();
                        for(int i=0; i<sliceN; i++){	
                                double dx = mx[0] - mx[i];
                                double dy = my[0] - my[i];				
                                driftXY[2*i]   = dx;
                                driftXY[2*i+1] = dy;
                        }
                    }                    
                    
                    return driftXY;
    }
    
    public static double[] get_driftXY_minimum_std(double[][]cx, double[][]cy, boolean[] best_flags, double[]varsum_array){
                    int nrois = cx.length;				
                    int sliceN = cx[0].length;			
                    boolean[] flags = new boolean[nrois];
                    //best_flags = new boolean[nrois];
                    int nrois_checked = 0; 		
                    for(int k=0; k<nrois; k++){			
                            flags[k] = true;			
                            for(int i=0; i<sliceN; i++){	
                                    if(flags[k]){
                                            if(Double.isNaN(cx[k][i]) || Double.isNaN(cy[k][i])) flags[k] = false;
                                    }
                                    else break;
                            }			
                            best_flags[k] = flags[k];
                            if(flags[k]) nrois_checked++;
                    }
                    IJ.log("total pillars:" + nrois + "\r\nremaining pillars:" + nrois_checked);
                    
                    boolean[] flags_x = new boolean[nrois];
                    boolean[] flags_y = new boolean[nrois];
                    for(int k=0; k<nrois; k++){
                        flags_x[k]=flags[k];
                        flags_y[k]=flags[k];
                    }
                    VarianceReduction varX = new VarianceReduction(cx, flags_x);
                    VarianceReduction varY = new VarianceReduction(cy, flags_y);
                    varX.start_var();
                    varY.start_var();
                    
                    double min_varsum = -1;		
                    int min_varsum_k = -1;		
                    int nleft = nrois_checked;
                    for(int deletetimes=0; deletetimes<nrois_checked-1; deletetimes++){
                                    
                                    double varsum = 0;					
                                    double maxvarxy = -1;
                                    int maxvarxy_k = -1;
                                    for(int k=0; k<nrois; k++) if(flags[k]){
                                            double varxy = varX.getVAR()[k] + varY.getVAR()[k];
                                            varsum += varxy;
                                            if(maxvarxy<0 || maxvarxy<varxy){
                                                    maxvarxy = varxy;
                                                    maxvarxy_k = k;
                                            }
                                    }

                                    varsum = varsum/nleft/nleft;
                                    if(min_varsum<0 || min_varsum>varsum){
                                            min_varsum = varsum;
                                            min_varsum_k = nleft;
                                            for(int k=0; k<nrois; k++) best_flags[k] = flags[k];			
                                    }				
                                    
                                    varsum_array[nleft-1] = varsum;
                                    
                                    varX.delete(maxvarxy_k);
                                    varY.delete(maxvarxy_k);
                                    flags[maxvarxy_k] = false;
                                    nleft--;	
                                    IJ.showStatus("searching the reference pillars " + deletetimes + "/" + (nrois_checked-1));		
                                    IJ.showProgress(deletetimes+1,nrois_checked-1);
                    }		
                    IJ.log("minimum std is " + Math.sqrt(min_varsum) + " at " + min_varsum_k);			

                    double[] driftXY = null;
                    if(nleft>0){
                        driftXY = new double[2*sliceN];
                        double[] mx = varX.getMeanM();
                        double[] my = varY.getMeanM();
                        for(int i=0; i<sliceN; i++){	
                                double dx = mx[0] - mx[i];
                                double dy = my[0] - my[i];				
                                driftXY[2*i]   = dx;
                                driftXY[2*i+1] = dy;
                        }
                    }                    
                    
                    return driftXY;
    }
    
    public static double[] get_driftXY_reduced_std_slow(double[][]cx, double[][]cy, boolean[] best_flags, double[]varsum_array){
                    int nrois = cx.length;				
                    int sliceN = cx[0].length;			
                    boolean[] flags = new boolean[nrois];
                    //best_flags = new boolean[nrois];
                    int nrois_checked = 0; 		
                    for(int k=0; k<nrois; k++){			
                            flags[k] = true;			
                            for(int i=0; i<sliceN; i++){	
                                    if(flags[k]){
                                            if(Double.isNaN(cx[k][i]) || Double.isNaN(cy[k][i])) flags[k] = false;
                                    }
                                    else break;
                            }			
                            best_flags[k] = flags[k];
                            if(flags[k]) nrois_checked++;
                    }
                    IJ.log("total pillars:" + nrois + "\r\nremaining pillars:" + nrois_checked);

                    //varsum_array = new double[nrois_checked];
                    double[] mx = new double[sliceN];			
                    double[] my = new double[sliceN];
                    //double[] var = new double[nrois];	
                    double[][] cxt = Matrix_Util.transpose(cx);		
                    double[][] cyt = Matrix_Util.transpose(cy);
                    double[][] diffx = new double[nrois][sliceN];
                    double[][] diffy = new double[nrois][sliceN];

                    double min_varsum = -1;		
                    int min_varsum_k = -1;		
                    int nleft = nrois_checked;
                    boolean find_turning_point = false;
                    int checktimes = 3;
                    for(int deletetimes=0; deletetimes<nrois_checked-1; deletetimes++){
                                    for(int i=0; i<sliceN; i++){	
                                            mx[i] = BasicStatisitic.avg(cxt[i],flags);
                                            my[i] = BasicStatisitic.avg(cyt[i],flags);

                                            for(int k=0; k<nrois; k++) if(flags[k]){
                                                    diffx[k][i] = cxt[i][k] - mx[i];
                                                    diffy[k][i] = cyt[i][k] - my[i];
                                            }					
                                    }

                                    double varsum = 0;					
                                    double maxvarxy = -1;
                                    int maxvarxy_k = -1;
                                    for(int k=0; k<nrois; k++) if(flags[k]){
                                            double varxy = BasicStatisitic.var(diffx[k]) + BasicStatisitic.var(diffy[k]);
                                            varsum += varxy;
                                            if(maxvarxy<0 || maxvarxy<varxy){
                                                    maxvarxy = varxy;
                                                    maxvarxy_k = k;
                                            }
                                    }

                                    varsum = varsum/nleft/nleft;				
                                    varsum_array[nleft-1] = varsum;				

                                    if(!find_turning_point){
                                            int nondecreasing_times = 0;				
                                            for(int j=0; j<checktimes; j++){
                                                    int current_k = nleft-1+j;
                                                    int last_k = current_k+1;
                                                    if(last_k<nrois_checked){						
                                                            if(varsum_array[current_k]>=varsum_array[last_k]) nondecreasing_times++;													
                                                            else break;
                                                    }
                                            }			

                                            if(nondecreasing_times == checktimes){ //found the first turning point				
                                                    min_varsum = varsum;
                                                    min_varsum_k = nleft;
                                                    for(int k=0; k<nrois; k++) best_flags[k] = flags[k];	
                                                    find_turning_point = true;															
                                            }					
                                    }

                                    if(find_turning_point){
                                            IJ.showStatus("found the first turning point at" + nleft);
                                            IJ.showProgress(1,1);
                                            break;		
                                    }

                                    flags[maxvarxy_k] = false; //deleting
                                    nleft--;	
                                    IJ.showStatus("searching the reference pillars " + deletetimes + "/" + (nrois_checked-1));		
                                    IJ.showProgress(deletetimes+1,nrois_checked-1);
                    }		
                    IJ.log("minimum std is " + Math.sqrt(min_varsum) + " at " + min_varsum_k);			

                    double[] driftXY = new double[2*sliceN];
                    for(int i=0; i<sliceN; i++){	
                            mx[i] = BasicStatisitic.avg(cxt[i],best_flags);
                            my[i] = BasicStatisitic.avg(cyt[i],best_flags);				
                            double dx = mx[0] - mx[i];
                            double dy = my[0] - my[i];				
                            driftXY[2*i]   = dx;
                            driftXY[2*i+1] = dy;
                    }

                    return driftXY;
    }

    public static double[] get_driftXY_minimum_std_slow(double[][]cx, double[][]cy, boolean[] best_flags, double[]varsum_array){
                    int nrois = cx.length;				
                    int sliceN = cx[0].length;			
                    boolean[] flags = new boolean[nrois];
                    //best_flags = new boolean[nrois];
                    int nrois_checked = 0; 		
                    for(int k=0; k<nrois; k++){			
                            flags[k] = true;			
                            for(int i=0; i<sliceN; i++){	
                                    if(flags[k]){
                                            if(Double.isNaN(cx[k][i]) || Double.isNaN(cy[k][i])) flags[k] = false;
                                    }
                                    else break;
                            }			
                            best_flags[k] = flags[k];
                            if(flags[k]) nrois_checked++;
                    }
                    IJ.log("total pillars:" + nrois + "\r\nremaining pillars:" + nrois_checked);

                    //varsum_array = new double[nrois_checked];
                    double[] mx = new double[sliceN];			
                    double[] my = new double[sliceN];
                    //double[] var = new double[nrois];	
                    double[][] cxt = Matrix_Util.transpose(cx);		
                    double[][] cyt = Matrix_Util.transpose(cy);
                    double[][] diffx = new double[nrois][sliceN];
                    double[][] diffy = new double[nrois][sliceN];

                    double min_varsum = -1;		
                    int min_varsum_k = -1;		
                    int nleft = nrois_checked;
                    for(int deletetimes=0; deletetimes<nrois_checked-1; deletetimes++){
                                    for(int i=0; i<sliceN; i++){	
                                            mx[i] = BasicStatisitic.avg(cxt[i],flags);
                                            my[i] = BasicStatisitic.avg(cyt[i],flags);

                                            for(int k=0; k<nrois; k++) if(flags[k]){
                                                    diffx[k][i] = cxt[i][k] - mx[i];
                                                    diffy[k][i] = cyt[i][k] - my[i];
                                            }					
                                    }

                                    double varsum = 0;					
                                    double maxvarxy = -1;
                                    int maxvarxy_k = -1;
                                    for(int k=0; k<nrois; k++) if(flags[k]){
                                            double varxy = BasicStatisitic.var(diffx[k]) + BasicStatisitic.var(diffy[k]);
                                            varsum += varxy;
                                            if(maxvarxy<0 || maxvarxy<varxy){
                                                    maxvarxy = varxy;
                                                    maxvarxy_k = k;
                                            }
                                    }

                                    varsum = varsum/nleft/nleft;
                                    if(min_varsum<0 || min_varsum>varsum){
                                            min_varsum = varsum;
                                            min_varsum_k = nleft;
                                            for(int k=0; k<nrois; k++) best_flags[k] = flags[k];			
                                    }				

                                    varsum_array[nleft-1] = varsum;
                                    flags[maxvarxy_k] = false;
                                    nleft--;	
                                    IJ.showStatus("searching the reference pillars " + deletetimes + "/" + (nrois_checked-1));		
                                    IJ.showProgress(deletetimes+1,nrois_checked-1);
                    }		
                    IJ.log("minimum std is " + Math.sqrt(min_varsum) + " at " + min_varsum_k);			

                    double[] driftXY = new double[2*sliceN];
                    for(int i=0; i<sliceN; i++){	
                            mx[i] = BasicStatisitic.avg(cxt[i],best_flags);
                            my[i] = BasicStatisitic.avg(cyt[i],best_flags);				
                            double dx = mx[0] - mx[i];
                            double dy = my[0] - my[i];				
                            driftXY[2*i]   = dx;
                            driftXY[2*i+1] = dy;
                    }

                    return driftXY;
    }
 
}
