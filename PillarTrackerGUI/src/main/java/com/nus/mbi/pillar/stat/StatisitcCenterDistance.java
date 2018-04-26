package com.nus.mbi.pillar.stat;

/**
 *
 * @author xiaochun
 */
public class StatisitcCenterDistance {
    public double avg_x, avg_y, avg_dis, max_dis, min_dis, std_dis;

    public StatisitcCenterDistance(double avg_raw_x, double avg_raw_y, double avg_dis, double max_dis, double min_dis, double std_dis) {
        this.avg_x = avg_raw_x;
        this.avg_y = avg_raw_y;
        this.avg_dis = avg_dis;
        this.max_dis = max_dis;
        this.min_dis = min_dis;
        this.std_dis = std_dis;
    }
    
}
