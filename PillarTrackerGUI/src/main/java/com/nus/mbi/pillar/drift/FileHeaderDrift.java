package com.nus.mbi.pillar.drift;

/**
 *
 * @author xiaochun
 */
public class FileHeaderDrift {
    //public int file_version;
    public int npillars;
    public int nframes;
    public double lattice;
    public double diameter;
    public double spacing;
    public double oblique;
    public double grid_angle;
    public double sigma_PSF;
    public double catch_radius;
    public int kernel_w;
    public int box_constrian_R;
    //public int num_start_points;
    public boolean dark_pillars;
    public boolean apply_mean_rank_filter;    
    //public boolean use_minimum_std;
    public boolean use_metric_CG;
    public boolean use_enhancer;
    public boolean use_fft;    
}
