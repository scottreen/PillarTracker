package com.nus.mbi.pillar.drift;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.HistogramWindow;
import ij.gui.Line;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.PointRoi;
import ij.gui.TextRoi;
import ij.process.ColorProcessor;
import ij.process.FloatPolygon;
import ij.process.FloatProcessor;
import ij.process.LUT;
import java.awt.Color;
import java.awt.Font;
import static com.nus.mbi.pillar.stat.BasicStatisitic.*;

/**
 *
 * @author xiaochun
 */
public class DriftAnalysisPlotter {
    int shown_pillar_index = 1;    
    boolean show_points_in_energy_map = false;
    double zoom_in_tracks = 15;
    double zoom_in_image = 6;
    boolean energy_map_raw = false;

    boolean plotting_energy = false;
    boolean plotting_histogram = false;
    boolean plotting_stdxy = false;
    boolean plotting_driftxy = false;
    boolean plotting_track = false;
    boolean plotting_super_track = false;

    boolean plotting_total_track = true;
    boolean plotting_label_track = false;

    boolean plotting_flattenRGB = false;
    
    public static boolean SHOW_TRACK_HEADER = true;

    static public int[] color_lut = {
                                     -16776964,-16776964,-16776964,-16513796,-16513796,-16513796,-16250628,-16250628,-15987460,-15987460,-15987464,-15724296,-15724296,-15461128,
                                     -15461128,-15461128,-15197960,-15197960,-14934792,-14934792,-14934796,-14671628,-14671628,-14408460,-14408460,-14408460,-14145292,-14145292,
                                     -13882124,-13882124,-13882128,-13618960,-13618960,-13355792,-13355792,-13355792,-13092624,-13092624,-12829456,-12829460,-12829460,-12566292,
                                     -12566292,-12303124,-12303124,-12303124,-12039956,-12039956,-11776788,-11776792,-11776792,-11513624,-11513624,-11250456,-11250456,-11250456,
                                     -10987288,-10987288,-10724120,-10724124,-10724124,-10460956,-10460956,-10460956,-10197788,-10197788,-9934620,-9934620,-9934624,-9671456,
                                     -9671456,-9408288,-9408288,-9408288,-9145120,-9145120,-8881952,-8881952,-8881956,-8618788,-8618788,-8355620,-8355620,-8355620,-8092452,
                                     -8092452,-7829284,-7829284,-7829288,-7566120,-7566120,-7302952,-7302952,-7302952,-7039784,-7039784,-6776616,-6776620,-6776620,-6513452,
                                     -6513452,-6250284,-6250284,-6250284,-5987116,-5987116,-5723948,-5723952,-5723952,-5460784,-5460784,-5197616,-5197616,-5197616,-4934448,
                                     -4934448,-4671280,-4671284,-4671284,-4408116,-4408116,-4144948,-4144948,-4144948,-3881780,-3881780,-3618616,-3618616,-3618616,-3619644,
                                     -3619644,-3619644,-3620672,-3620672,-3621700,-3621700,-3359556,-3360584,-3360584,-3361612,-3361612,-3361612,-3362640,-3362640,-3363668,
                                     -3363668,-3101524,-3102552,-3102552,-3102552,-3103580,-3103580,-3104608,-3104608,-3104608,-3105636,-2843492,-2844520,-2844520,-2844520,
                                     -2845548,-2845548,-2846576,-2846576,-2846576,-2847604,-2585460,-2586488,-2586488,-2586488,-2587516,-2587516,-2587516,-2588544,-2588544,
                                     -2589572,-2327428,-2327428,-2328456,-2328456,-2329484,-2329484,-2329484,-2330512,-2330512,-2331540,-2069396,-2069396,-2070424,-2070424,
                                     -2071452,-2071452,-2071452,-2072480,-2072480,-1810336,-1811364,-1811364,-1812392,-1812392,-1812392,-1813420,-1813420,-1814448,-1814448,
                                     -1552304,-1553332,-1553332,-1554360,-1554360,-1554360,-1555388,-1555388,-1555388,-1556416,-1294272,-1295300,-1295300,-1295300,-1296328,
                                     -1296328,-1297356,-1297356,-1297356,-1298384,-1036240,-1037268,-1037268,-1037268,-1038296,-1038296,-1039324,-1039324,-1039324,-1040352,
                                     -778208,-778208,-779236,-779236,-780264,-780264,-780264,-781292,-781292,-782320,-520176,-520176,-521204,-521204,-522232,-522232,-522232,-523260,-523260,-262144,-262144
                                     }; //phase lut
                                     
    private Overlay labels_overlay = null;
    private Overlay tracks_overlay = null;
    
	public static Plot get_plot_driftXY(double[] xaxis, double[] disX, double[] disY, String title)
	{
                return plot3(xaxis, disX, disY, title, "frame #", "displacement(nm) refer to 1st frame");
	}
        
        public static Plot plot3(double[] xaxis, double[] disX, double[] disY, String title, String xlabel, String ylabel)
	{
                int n = disX.length;
		double[] dis = new double[n];
                for(int i=0; i<n; i++) dis[i] = Math.sqrt(disX[i]*disX[i]+disY[i]*disY[i]);                            
                float plotlimitmin = Math.min((float)min(disX), (float)min(disY));
		float plotlimitmax = Math.max((float)max(disX), (float)max(disY));
                plotlimitmax = Math.max((float)max(dis), plotlimitmax);
		//Plot plot = new Plot("ROI-" + (i+1),"slice", "distance(nm)", xaxis, dis[i]);
		Plot plot = new Plot(title, xlabel, ylabel);
		plot.setLimits(0,disX.length,plotlimitmin,plotlimitmax);
		
		plot.setColor(Color.RED);	
		plot.addPoints(xaxis, disX,Plot.LINE);			
		plot.draw();	
		
		plot.setColor(Color.BLUE);	
		plot.addPoints(xaxis, disY,Plot.LINE);		
		plot.draw();
                
                plot.setColor(Color.CYAN);	
		plot.addPoints(xaxis, dis,Plot.LINE);		
		plot.draw();
	
		return plot;
	}
        
        public static Plot plot2(double[] xaxis, double[] disX, double[] disY, String title, String xlabel, String ylabel)
	{
                                          
                float plotlimitmin = Math.min((float)min(disX), (float)min(disY));
		float plotlimitmax = Math.max((float)max(disX), (float)max(disY));
                
		Plot plot = new Plot(title, xlabel, ylabel);
		plot.setLimits(0,disX.length,plotlimitmin,plotlimitmax);
		
		plot.setColor(Color.RED);	
		plot.addPoints(xaxis, disX,Plot.LINE);			
		plot.draw();	
		
		plot.setColor(Color.BLUE);	
		plot.addPoints(xaxis, disY,Plot.LINE);		
		plot.draw();
                	
		return plot;
	}
        
        public void show_histogram(double[] data, double bin_width, String title){
		int len = data.length;
		double hist_min = min(data);
		double hist_max = max(data);
		int bins = (int)Math.ceil((hist_max-hist_min)/bin_width);
		if(bins>256) bins = 256;
                else if(bins<32) bins = 32;
                FloatProcessor ip = new FloatProcessor(1,len);
		for(int i=0; i<len; i++) ip.setf(i,(float)data[i]);
		ImagePlus img = new ImagePlus(title, ip);
		
                HistogramWindow hist= new HistogramWindow(title,img,bins);                
		//hist.show();
	}       
        
        public void show_jump_length_histogram(double[][] cx, double[][] cy, double bin_width, String title){
		int npillars = cx.length;                
                double[] data = getJumpLength(cx, cy, 0);
                int nframes = data.length;
                
		double hist_min = min(data);
		double hist_max = max(data);                
                FloatProcessor ip = new FloatProcessor(nframes,npillars);		
                for(int f=0; f<nframes; f++) ip.setf(f,0,(float)data[f]);
                
                for(int i=1; i<npillars; i++){
                    data = getJumpLength(cx, cy, i);
                    hist_min = Math.min(hist_min, min(data));
                    hist_max = Math.max(hist_max, max(data));
                    for(int f=0; f<nframes; f++) ip.setf(f,i,(float)data[f]);
                }
                
		int bins = (int)Math.ceil((hist_max-hist_min)/bin_width);
		if(bins>256) bins = 256;
                else if(bins<32) bins = 32;
		ImagePlus img = new ImagePlus(title, ip);
		
                HistogramWindow hist= new HistogramWindow(title,img,bins);                
		//hist.show();
	} 
        
        public double[] getJumpLength(double[][] cx, double[][] cy, int ipillar){             
            int nframes = cx[0].length;
            double [] jump = new double[nframes-1];
            for (int f = 0; f < nframes-1; f++){                 			
                double dx = cx[ipillar][f+1] - cx[ipillar][f];
                double dy = cy[ipillar][f+1] - cy[ipillar][f];
                jump[f] = Math.sqrt(dx*dx+dy*dy);
            }
            return jump;
        }
        
        public Overlay getLabels(double[][] cx, double[][] cy, double pixel_size){
            Overlay overlay = new Overlay();  
            if(cx==null) return overlay;            
            int nrois = cx.length;
            if(nrois<1) return overlay;
            int sliceN = cx[0].length;
            if(sliceN<1) return overlay;
            //int nrois = cx.length;
            Font font = new Font("Arial",Font.PLAIN, 3);
            //Overlay overlay = new Overlay();            
            for(int i=0; i<nrois; i++){
                double rcx = cx[i][0]/pixel_size;//avg(cx[i])/pixel_size;			
                double rcy = cy[i][0]/pixel_size;//avg(cy[i])/pixel_size;						

                TextRoi label = new TextRoi(rcx, rcy, ""+(i+1), font);
                label.setStrokeColor(Color.yellow);
                overlay.add(label);                
            }
            labels_overlay = overlay;
            return overlay;        
        }
        
        public Overlay getLabels(double[][] cx, double[][] cy, int[] index, double pixel_size){            
            //if(index==null) return getLabels(cx, cy, pixel_size);   
            Overlay overlay = new Overlay();  
            if(index==null) return overlay;
            if(cx==null) return overlay;            
            int nrois = cx.length;
            if(nrois<1) return overlay;
            int sliceN = cx[0].length;
            if(sliceN<1) return overlay;
            int n = index.length;
            //int nrois = cx.length;
            Font font = new Font("Arial",Font.PLAIN, 3);            
            for(int j=0; j<n; j++){
                int i = index[j];
                    if(i>=0 && i<nrois){
                    double rcx = cx[i][0]/pixel_size;//avg(cx[i])/pixel_size;			
                    double rcy = cy[i][0]/pixel_size;//avg(cy[i])/pixel_size;						

                    TextRoi label = new TextRoi(rcx, rcy, ""+(i+1), font);
                    label.setStrokeColor(Color.yellow);
                    overlay.add(label); 
                }
            }
            labels_overlay = overlay;
            
            return overlay;
        }
        
        public Overlay getLabels(double[][] cx, double[][] cy, double pixel_size, double image_zoom){
            Overlay overlay = new Overlay();            
            if(cx==null) return overlay;            
            int nrois = cx.length;
            if(nrois<1) return overlay;
            int sliceN = cx[0].length;
            if(sliceN<1) return overlay;
            //int nrois = cx.length;
            
            int font_size = (int)Math.ceil(3*image_zoom);
            Font font = new Font("Arial",Font.PLAIN, font_size);            
            for(int i=0; i<nrois; i++){
                double rcx = cx[i][0]/pixel_size;//avg(cx[i])/pixel_size;			
                double rcy = cy[i][0]/pixel_size;//avg(cy[i])/pixel_size;						
                double zoom_rcx = image_zoom*rcx;
                double zoom_rcy = image_zoom*rcy;
                TextRoi label = new TextRoi(zoom_rcx, zoom_rcy, ""+(i+1), font);
                label.setStrokeColor(Color.yellow);
                overlay.add(label);                
            }
            return overlay;
        }       
        
        public Overlay getTracks(double[][] cx, double[][] cy, double pixel_size, double vector_zoom){
            return getTracks(cx,cy,pixel_size,vector_zoom, color_lut);
        }
        
        public Overlay getTracks(double[][] cx, double[][] cy, double pixel_size, double vector_zoom, int[] color_table){
            Overlay overlay = new Overlay();
            if(cx==null) return overlay;            
            int nrois = cx.length;
            if(nrois<1) return overlay;
            int sliceN = cx[0].length;
            if(sliceN<1) return overlay;
            return getTracks(cx,cy,pixel_size,0, sliceN-1, vector_zoom, color_lut);
	}
        
        public Overlay getTracks(double[][] cx, double[][] cy, double pixel_size,int start_frame, int end_frame, double vector_zoom){
            return getTracks(cx,cy,pixel_size,start_frame, end_frame, vector_zoom, color_lut);
        }
        
        public Overlay getTracks(double[][] cx, double[][] cy, double pixel_size, int start_frame, int end_frame, double vector_zoom, int[] color_table){                                   
            Overlay overlay = new Overlay();
            if(cx==null) return overlay;            
            int nrois = cx.length;
            if(nrois<1) return overlay;
            int sliceN = cx[0].length;
            if(sliceN<1) return overlay;
            
            if(start_frame<0) start_frame = 0;
            if(end_frame>sliceN-1) end_frame = sliceN-1;
            if(start_frame>end_frame){
                tracks_overlay = overlay;
                return overlay;
            }
            
            FloatPolygon poly = new FloatPolygon();
            
            int nframes = end_frame-start_frame+1;
            int color_length = color_table.length;
            for(int i=0; i<nrois; i++){
                    double rcx = cx[i][0]/pixel_size + 0.5;//avg(cx[i])/pixel_size;			
                    double rcy = cy[i][0]/pixel_size + 0.5;//avg(cy[i])/pixel_size;
                    
                    for(int p=start_frame; p<=end_frame; p++){							 											 	
                            double mlx = cx[i][p]/pixel_size + 0.5;
                            double mly = cy[i][p]/pixel_size + 0.5;
                            double xx = (mlx-rcx)*vector_zoom + rcx;
                            double yy = (mly-rcy)*vector_zoom + rcy;
                            //if(p<sliceN-1)
                            if(p<end_frame)
                            {
                                    double mlx1 = cx[i][p+1]/pixel_size + 0.5;
                                    double mly1 = cy[i][p+1]/pixel_size + 0.5;				
                                    double xx1 = (mlx1-rcx)*vector_zoom + rcx;
                                    double yy1 = (mly1-rcy)*vector_zoom + rcy;
                                    if(IsNotNaN(xx,yy,xx1,yy1))
                                    {				
                                            Line line = new Line(xx,yy,xx1,yy1);					
                                            //int indexcolor = p*color_length/sliceN; //normalized within(0~256);
                                            int indexcolor = (p-start_frame)*color_length/nframes; //normalized within(0~256);
                                            line.setStrokeColor(new Color(color_table[indexcolor]));						
                                            overlay.add(line);
                                    }
                            }
                            else if(!Double.isNaN(xx)){
                                poly.addPoint(xx, yy);
                            }
                    }
            }
            
            if(SHOW_TRACK_HEADER){
                //poly.addPoint(1, 1);
                PointRoi points = new PointRoi(poly);
                points.setPointType(2);
                points.setStrokeColor(Color.RED);
                points.setSize(2);
                overlay.add(points);
            }
            
            tracks_overlay = overlay;
            return overlay;
	}
        
        public Overlay getTracks(double[][] cx, double[][] cy, double pixel_size, double vector_zoom, LUT lut){
		int nrois = cx.length;
		int sliceN = cx[0].length;
		Overlay overlay = new Overlay();		
		for(int i=0; i<nrois; i++){
			double rcx = cx[i][0]/pixel_size + 0.5;//avg(cx[i])/pixel_size;			
			double rcy = cy[i][0]/pixel_size + 0.5;//avg(cy[i])/pixel_size;
			
			for(int p=0; p<sliceN; p++){							 											 	
				double mlx = cx[i][p]/pixel_size + 0.5;
				double mly = cy[i][p]/pixel_size + 0.5;
				double xx = (mlx-rcx)*vector_zoom + rcx;
				double yy = (mly-rcy)*vector_zoom + rcy;
				if(p<sliceN-1)
				{
					double mlx1 = cx[i][p+1]/pixel_size + 0.5;
					double mly1 = cy[i][p+1]/pixel_size + 0.5;				
					double xx1 = (mlx1-rcx)*vector_zoom + rcx;
					double yy1 = (mly1-rcy)*vector_zoom + rcy;
					if(IsNotNaN(xx,yy,xx1,yy1))
					{				
						Line line = new Line(xx,yy,xx1,yy1);					
						int indexcolor = p*256/sliceN; //normalized within(0~256);
						line.setStrokeColor(new Color(lut.getRGB(indexcolor)));						
						overlay.add(line);
					}
				}			
			}
		}
                tracks_overlay = overlay;
		return overlay;
	}
                
        public ImagePlus plot_flatten_tracks(double[][] cx, double[][] cy, double pixel_size, int width, int height, double vector_zoom, double image_zoom, String title){
            return plot_flatten_tracks( cx,  cy,  pixel_size,  width,  height,  vector_zoom,  image_zoom, color_lut, title);
        }
        
        public void plot_flatten_tracks(double[][] cx, double[][] cy, double pixel_size, double vector_zoom, double image_zoom, ImagePlus tracks_ip){
                int nrois = cx.length;
		int sliceN = cx[0].length;
		int zoom_width = tracks_ip.getWidth();
		int zoom_height = tracks_ip.getHeight();
                int zoom_size = zoom_width*zoom_height;
                
		ColorProcessor tracks_overlays_image = (ColorProcessor)tracks_ip.getProcessor(); 
                for(int i=0; i<zoom_width; i++)
                    for(int j=0; j<zoom_height; j++){
                        tracks_overlays_image.set(i, j, 0);
                    }
                
		int color_length = color_lut.length;
		for(int i=0; i<nrois; i++){
			double rcx = cx[i][0]/pixel_size + 0.5;//avg(cx[i])/pixel_size;			
			double rcy = cy[i][0]/pixel_size + 0.5;//avg(cy[i])/pixel_size;
			
			double zoom_rcx = image_zoom*rcx;
			double zoom_rcy = image_zoom*rcy;			
			
			for(int p=0; p<sliceN; p++){							 											 	
				double mlx = cx[i][p]/pixel_size + 0.5;
				double mly = cy[i][p]/pixel_size + 0.5;
				double xx = (mlx-rcx)*vector_zoom + zoom_rcx;
				double yy = (mly-rcy)*vector_zoom + zoom_rcy;
				//int ixx = (int)Math.round(xx);
				//int iyy = (int)Math.round(yy);
				//int indexcolor = p*color_length/sliceN; //normalized within(0~256);
				//int color = color_table[indexcolor];
				//tracks_overlays_image.putPixel(ixx, iyy, color);		
				if(p<sliceN-1)
				{
					double mlx1 = cx[i][p+1]/pixel_size + 0.5;
					double mly1 = cy[i][p+1]/pixel_size + 0.5;				
					double xx1 = (mlx1-rcx)*vector_zoom + zoom_rcx;
					double yy1 = (mly1-rcy)*vector_zoom + zoom_rcy;
					if(IsNotNaN(xx,yy,xx1,yy1))
					{				
						int indexcolor = p*color_length/sliceN; //normalized within(0~256);
						int color = color_lut[indexcolor];
						tracks_overlays_image.setColor(color);
						
						int ixx = (int)Math.round(xx);
						int iyy = (int)Math.round(yy);
						int ixx1 = (int)Math.round(xx1);
						int iyy1 = (int)Math.round(yy1);	
						
						tracks_overlays_image.drawLine(ixx, iyy, ixx1, iyy1);		
						//line.setStrokeColor(new Color(color_table[indexcolor]));						
						//overlay.add(line);
					}
				}						
			}		
		}
                
                tracks_ip.repaintWindow();
        }
        
	public ImagePlus plot_flatten_tracks(double[][] cx, double[][] cy, double pixel_size, int width, int height, double vector_zoom, double image_zoom, int[] color_table, String title){
		int nrois = cx.length;
		int sliceN = cx[0].length;
		
		int zoom_width = (int)Math.ceil(image_zoom*width);
		int zoom_height = (int)Math.ceil(image_zoom*height);
		
		ColorProcessor tracks_overlays_image = new ColorProcessor(zoom_width, zoom_height); 
		ImagePlus tracks_overlays_plus = new ImagePlus(title + vector_zoom, tracks_overlays_image);
		
		int color_length = color_table.length;
		for(int i=0; i<nrois; i++){
			double rcx = cx[i][0]/pixel_size + 0.5;//avg(cx[i])/pixel_size;			
			double rcy = cy[i][0]/pixel_size + 0.5;//avg(cy[i])/pixel_size;
			
			double zoom_rcx = image_zoom*rcx;
			double zoom_rcy = image_zoom*rcy;			
			
			for(int p=0; p<sliceN; p++){							 											 	
				double mlx = cx[i][p]/pixel_size + 0.5;
				double mly = cy[i][p]/pixel_size + 0.5;
				double xx = (mlx-rcx)*vector_zoom + zoom_rcx;
				double yy = (mly-rcy)*vector_zoom + zoom_rcy;
					
				if(p<sliceN-1)
				{
					double mlx1 = cx[i][p+1]/pixel_size + 0.5;
					double mly1 = cy[i][p+1]/pixel_size + 0.5;				
					double xx1 = (mlx1-rcx)*vector_zoom + zoom_rcx;
					double yy1 = (mly1-rcy)*vector_zoom + zoom_rcy;
					if(IsNotNaN(xx,yy,xx1,yy1))
					{				
						int indexcolor = p*color_length/sliceN; //normalized within(0~256);
						int color = color_table[indexcolor];
						tracks_overlays_image.setColor(color);
						
						int ixx = (int)Math.round(xx);
						int iyy = (int)Math.round(yy);
						int ixx1 = (int)Math.round(xx1);
						int iyy1 = (int)Math.round(yy1);	
						
						tracks_overlays_image.drawLine(ixx, iyy, ixx1, iyy1);		
						//line.setStrokeColor(new Color(color_table[indexcolor]));						
						//overlay.add(line);
					}
				}						
			}
		}

		//tracks_overlays_plus.setOverlay(overlay);
		tracks_overlays_plus.show();
		return tracks_overlays_plus;
	}
	
	public ImagePlus create_empty_image(int width, int height, String title){
		ColorProcessor tracks_overlays_image = new ColorProcessor(width, height); 
		ImagePlus tracks_overlays_plus = new ImagePlus(title, tracks_overlays_image);				
		return tracks_overlays_plus;
	}
        
        public void plot_tracks_overlay(int width, int height, String title, boolean show_label){
            if(tracks_overlay!=null){
                ImagePlus tracks_overlays_plus = create_empty_image(width, height, title);
                plot_tracks_overlay(tracks_overlays_plus, show_label);
            }
        }
        
//        public void empty_tracks_overply(ImagePlus ip){
//            ip.setOverlay(new Overlay());
//        }
        
        public void plot_tracks_overlay(ImagePlus ip, boolean show_label){
            //if(tracks_overlay!=null)
            {                            
                if(show_label){
                    if(labels_overlay!=null){
                        Overlay overlay = tracks_overlay.duplicate();
                        int n = labels_overlay.size();
                        for(int i=0; i<n; i++) overlay.add(labels_overlay.get(i));
                        
                        ip.setOverlay(overlay);
                        
                    }
                }
                else{
                    ip.setOverlay(tracks_overlay);
                }
                
                if(ip.getWindow()==null) ip.show();
            }
        }
        
        public void plot_flatten_tracks(ImagePlus ip, Overlay labels){            
            ip.setOverlay(labels);                
        }
        
        public ImagePlus create_empty_map(double[][] cx, double[][] cy, double pixel_size, int margin, String title){
            int nrois = cx.length;	
            double xmax = 0;
            double ymax = 0;
            for(int i=0; i<nrois; i++){
                    double x = max(cx[i]);
                    double y = max(cy[i]);
                    if(!Double.isNaN(x) && !Double.isNaN(y)){
                            if(x>xmax) xmax = x;
                            if(y>ymax) ymax = y;
                    }
            }

            int map_img_width = (int)Math.ceil(xmax/pixel_size + margin);
            int map_img_height = (int)Math.ceil(ymax/pixel_size + margin);
            IJ.log("xmax="+ xmax + " ymax=" + ymax + " map_img_width="+ map_img_width + " map_img_height=" + map_img_height);
            return create_empty_image(map_img_width, map_img_height, title);
        }
              
}
