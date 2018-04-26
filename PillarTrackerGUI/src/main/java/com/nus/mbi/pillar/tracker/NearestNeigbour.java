package com.nus.mbi.pillar.tracker;

import java.util.List;

/**
 *
 * @author nus
 */
public class NearestNeigbour {
    public List<Integer> neigbours;
    public int distance;

    public NearestNeigbour(List<Integer> neigbours, int min_distance) {
        this.neigbours = neigbours;
        this.distance = min_distance;
    }
}
