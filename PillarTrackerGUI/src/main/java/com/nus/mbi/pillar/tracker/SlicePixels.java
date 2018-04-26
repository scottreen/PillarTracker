package com.nus.mbi.pillar.tracker;

import com.nus.mbi.pillar.detection.FDResult;

/**
 *
 * @author xiaochun
 */
public class SlicePixels {
    double[] pixels;
    double[] ref_pixels;
    FDResult fft_data;
    FDResult ref_fft_data;
    

    public SlicePixels(double[] pixels, double[] ref_pixels, FDResult fft_data, FDResult ref_fft_data) {
        this.pixels = pixels;
        this.ref_pixels = ref_pixels;
        this.fft_data = fft_data;
        this.ref_fft_data = ref_fft_data;
    }
}
