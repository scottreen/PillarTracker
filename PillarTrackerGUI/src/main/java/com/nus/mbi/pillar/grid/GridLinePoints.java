package com.nus.mbi.pillar.grid;

/**
 *
 * @author xiaochun
 */
public class GridLinePoints implements Comparable<GridLinePoints> {
    public double x;
    public double y;
    public int ix;
    public int iy;

    public GridLinePoints(double x, double y, int ix, int iy) {
        this.x = x;
        this.y = y;
        this.ix = ix;
        this.iy = iy;
    }

    @Override
    public int compareTo(GridLinePoints t) {
        if(t.ix==ix)
            return t.iy>iy ? 1 : -1;
        else if(t.iy==iy)
            return t.ix>ix ? 1 : -1;
        else return -1;
    }
}
