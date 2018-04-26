package com.nus.mbi.pillar.tracker;

/**
 *
 * @author xiaochun
 */
public class LocalizationWindowQueue {
    public static final int MAX_SIZE = 102400;
    
    LocalizationWindow[] windows;
    int[] pillar_id;
    int[] frame_id;
    
    int count = 0;
    int size = 0;

    public LocalizationWindow[] getWindows() {
        return windows;
    }

    public int[] getPillarID() {
        return pillar_id;
    }

    public int[] getFrameID() {
        return frame_id;
    }
    
    public int getCount() {
        return count;
    }

    public int getSize() {
        return size;
    }

    public LocalizationWindowQueue() {
        windows = new LocalizationWindow[MAX_SIZE];
        pillar_id = new int[MAX_SIZE];
        frame_id = new int[MAX_SIZE];
        this.size = MAX_SIZE;
        this.count = 0;
    }
    
    public void LocalizationQueue(int size){
        if(size<=0 || size>MAX_SIZE) size = MAX_SIZE;
        windows = new LocalizationWindow[size];
        pillar_id = new int[size];
        frame_id = new int[size];
        this.size = size;
        this.count = 0;
    }
    
    public void add(LocalizationWindow win, int ipillar, int iframe){
        windows[count] = win;
        pillar_id[count] = ipillar;
        frame_id[count] = iframe;
        count++;
    }
    
    public boolean isEmpty(){
        return count==0;
    }
    
    public boolean isFull(){
        return count==size;
    }
    
    public void reset(){
        count = 0;
    }
}
