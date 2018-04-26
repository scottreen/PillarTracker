package com.nus.mbi.pillar.detection;

import ij.process.FHT;
import ij.process.ImageProcessor;

/**
 *
 * @author xiaochun
 */
public class FDResult {
    public FHT fft_img;
    public FHT fft_mask;                
    public ImageProcessor ifft_img;

    public FDResult(FHT fft_img, FHT fft_mask, ImageProcessor ifft_img) {
        this.fft_img = fft_img;
        this.fft_mask = fft_mask;
        this.ifft_img = ifft_img;
    }
    
    
}
