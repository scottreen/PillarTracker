package com.nus.mbi.pillar.tracker.gui;

import ij.plugin.PlugIn;
import java.awt.Desktop;
import java.net.URI;

/**
 *
 * @author xiaochun
 */
public class Pillar_Tracker_Documentation implements PlugIn {
    String url_path = "https://drive.google.com/file/d/0B3hxvkn3VvhCVWJsdUN3eDUyNkk";
    //"https://drive.google.com/open?id=0B3hxvkn3VvhCVWJsdUN3eDUyNkk";
    
    public void run(String arg) {
        if (Desktop.isDesktopSupported()) {
                    try {
                    URI uri = new URI(url_path);
                    Desktop.getDesktop().browse(uri);
                } catch (Exception ex) {	        
            }
	}     
    }    
}
