package com.nus.mbi.pillar.detection;

/**
 *
 * @author xiaochun
 */
public class MaskFilterDialogResults {
    public boolean dialog_return = false;
    public int radius_start=1;
    public int radius_end=1;
    public int radius_step=1;

    public MaskFilterDialogResults() {
    }

    
    
    public void setRaiuds(int start, int end, int step) {
        this.radius_start = start;
        this.radius_end = end;
        this.radius_step = step;
    }
    
}
