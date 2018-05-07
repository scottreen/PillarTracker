package com.nus.mbi.pillar.tracker.gui;
/*
update history:
2017-09-14 (V1.1.4):
    1).hunt CU using sliding window to calculate cross correlation
    2).faster create movie. 
    3).draw dragon tail for tracks.
    4).support the batch processing for pillar tracking FD
    5).save and load the cu file
    6).faster drift correction for searching the silent pillars
    7).invert the minima in direct solver for the lev-mar algorithm
    8).bugs fixed
    9).relaxiation degree.
    10).muliti-thread acceleration x1.5~2
    11).add the rank mean filter to find maximas.
*/

import ij.IJ;
import ij.ImageJ;
import ij.plugin.PlugIn;

public class Pillar_Tracker_Plugin implements PlugIn {
    
    public static final String PLUGIN_NAME = "PillarTracker";
    public static final String VERSION = version();
    
    private static String version() {
        String version = null;
        final Package pack = Pillar_Tracker_Plugin.class.getPackage();        
        if (pack != null) version = pack.getImplementationVersion();        
        return version == null ? "DEVELOPMENT" : version;        
    }

    @Override
    public void run(String arg) {
        //IJ.log("Running " + PLUGIN_NAME + " version:" + VERSION);
        MainAppFrame app = new MainAppFrame();
        if(arg.contains("debug")) app.setPlatformImplicitExit(true);
        
        //app.dispatchEvent(new WindowEvent(app, WindowEvent.WINDOW_CLOSING));
        app.setTitle(PLUGIN_NAME + " Version:" + VERSION);               
        //app.setTitle(PLUGIN_NAME);               
        // Launch JavaFX interface                    
        app.init();         
    }
    
    public static void main(String[] args) throws Exception {
        Class<?> clazz = Pillar_Tracker_Plugin.class;
        //debug = true;
        //lunch ImageJ
        new ImageJ();        
        //IJ.log(clazz.getName()); 
        IJ.runPlugIn(clazz.getName(), "debug");    
    }
    
}
