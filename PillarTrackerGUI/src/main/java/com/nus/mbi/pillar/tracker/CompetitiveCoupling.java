package com.nus.mbi.pillar.tracker;


import java.awt.Point;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author xiaochun
 */
public class CompetitiveCoupling {
    private final int width;
    private final int height;
    private final double catch_radius;
    private final Point[] points1;
    private final Point[] points2;
    
    int[] map1;
    int[] map2;
    
    int[] forward;
    int[] backward;
    
    private int num_found = 0;
    int num_total_found = 0;

    public int[] getMap1() {
        return map1;
    }

    public int[] getMap2() {
        return map2;
    }

    public int[] getForward() {
        return forward;
    }

    public int[] getBackward() {
        return backward;
    }

    public int getNumTotalFound() {
        return num_total_found;
    }

    public CompetitiveCoupling(Point[] points1, Point[] points2, int width, int height, double catch_radius) {
        this.width = width;
        this.height = height;
        this.catch_radius = catch_radius;
        this.points1 = points1;
        this.points2 = points2;
        
        map1 = create_map(points1, width, height);
        map2 = create_map(points2, width, height);
        
        int n1=points1.length;
        forward = new int[n1];
        for(int i=0; i<n1; i++) forward[i] = -1;
        int n2=points2.length;
        backward = new int[n2];
        for(int i=0; i<n2; i++) backward[i] = -1;
    }
    
    public static int[] create_map(Point[] p, int width, int height){
        int n = p.length;
        int size = width*height;
        int[] map = new int[size];
        for(int i=0;i<size; i++) map[i]=-1;
        for(int i=0;i<n;i++){
            int x = p[i].x;
            int y = p[i].y;
            map[y*width+x] = i;
        }
        return map;
    }
    
    public void coupling(){
        int num_iter_found = 0;
        do{
            num_iter_found = 0;
            for(int i=0; i<points1.length; i++){
                if(forward[i]<0){
                    num_found = 0;
                    coupling(i,points1,points2,map1,map2,forward,backward,catch_radius);
                    num_iter_found += num_found;
                }
            }
            num_total_found += num_iter_found;
        }
        while(num_iter_found>0);
    }
    
    public NearestNeigbour coupling(int i, Point[] p1, Point[] p2, int[] m1, int[] m2, int[] forward, int[] backward, double radius){
        int x = p1[i].x;
        int y = p1[i].y;
        NearestNeigbour nb1 = search_nearest_neigbors(x,y,m2,width,height,radius);
        int num = nb1.neigbours.size();
        if(num>0){
            double shrink_radius = Math.sqrt(nb1.distance);
            for(int k=0; k<num; k++){
                int j=nb1.neigbours.get(k);
                NearestNeigbour nb2 = coupling(j,p2,p1,m2,m1,backward,forward,shrink_radius);
                if(nb2.distance>=nb1.distance){
                    forward[i] = j;
                    backward[j] = i;
                    m1[i] = -1;
                    m2[j] = -1;
                    num_found++;
                    break;
                }
            }
        }
        
        return nb1;
    }
    
    public static NearestNeigbour search_nearest_neigbors(int x, int y, int[] map, int width, int height, double catch_radius){
        int r = (int)Math.ceil(catch_radius);
        List<Integer> neigbors = new ArrayList();
        int min_dis = Integer.MAX_VALUE;
        for(int i=-r; i<=r; i++){
            int yy = y+i;
            int off = yy*width;
            for(int j=-r; j<=r; j++){
                int xx = x+j;
                
                int p = map[off+xx];
                if(p>=0){
                    int dis = i*i+j*j;
                    if(dis<min_dis){
                        neigbors.clear();
                        neigbors.add(p);
                        min_dis = dis;
                    }
                    else if(dis==min_dis){
                        neigbors.add(p);
                    }
                }
            }
        }
        
        return new NearestNeigbour(neigbors, min_dis);
    }
}
