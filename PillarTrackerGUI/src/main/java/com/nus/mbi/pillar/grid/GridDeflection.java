package com.nus.mbi.pillar.grid;

/**
 *
 * @author xiaochun
 */
public class GridDeflection {
    public double x;
    public double y;
    public double dx;
    public double dy;
    public boolean dxy_active;

    public GridDeflection(double x, double y, double dx, double dy, boolean dxy_active) {
        this.x = x;
        this.y = y;
        this.dx = dx;
        this.dy = dy;
        this.dxy_active = dxy_active;
    }
}
