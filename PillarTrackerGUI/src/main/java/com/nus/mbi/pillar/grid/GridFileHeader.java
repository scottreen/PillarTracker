package com.nus.mbi.pillar.grid;

/**
 *
 * @author xiaochun
 */
public class GridFileHeader {
    public double lattice;
    public double diameter;
    public double spacing;
    public double grid_oblique;
    public double grid_angle;
    public double deflection_threshold;

    public GridFileHeader(double lattice, double diameter, double spacing, double grid_oblique, double grid_angle, double threshold) {
        this.lattice = lattice;
        this.diameter = diameter;
        this.spacing = spacing;
        this.grid_oblique = grid_oblique;
        this.grid_angle = grid_angle;
        this.deflection_threshold = threshold;
    }
    
    public int ByteSize(){
        return Double.SIZE*6;
    }
}
