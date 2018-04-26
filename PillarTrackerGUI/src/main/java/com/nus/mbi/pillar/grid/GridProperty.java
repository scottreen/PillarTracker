package com.nus.mbi.pillar.grid;

/**
 *
 * @author xiaochun
 */
public class GridProperty {
    public double spacing;
    public double oblique;
    public double grid_angle = 90;

    public GridProperty() {
    }

    public GridProperty(double spacing, double oblique) {
        this.spacing = spacing;
        this.oblique = oblique;
        this.grid_angle = 90;
    }

    public GridProperty(double spacing, double oblique, double grid_angle) {
        this.spacing = spacing;
        this.oblique = oblique;
        this.grid_angle = grid_angle;
    }
    
    
}
