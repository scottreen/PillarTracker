package com.nus.mbi.pillar.tracker.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;
import java.awt.Image;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.Robot;

/**
 *
 * @author xiaochun
 */
public class MyCaptureImage {
     public static ImagePlus captureImage(ImagePlus imp, int wait_ms) {
        //ImagePlus imp = IJ.getImage();
        if (imp==null) {
            //IJ.noImage();
            return null;
        }
        ImageWindow win = imp.getWindow();
        if (win==null) return null;
        win.toFront();
        //IJ.wait(500);
        Point loc = win.getLocation();
        ImageCanvas ic = win.getCanvas();
        Rectangle bounds = ic.getBounds();
        loc.x += bounds.x;
        loc.y += bounds.y;
        Rectangle r = new Rectangle(loc.x, loc.y, bounds.width, bounds.height);
        ImagePlus imp2 = null;
        Image img = null;
        boolean wasHidden = ic.hideZoomIndicator(true);
        if(wait_ms>0) IJ.wait(wait_ms);
        try {
            Robot robot = new Robot();
            img = robot.createScreenCapture(r);
        } catch(Exception e) { }
        ic.hideZoomIndicator(wasHidden);
        if (img!=null) {
            String title = WindowManager.getUniqueName(imp.getTitle());
            imp2 = new ImagePlus(title, img);
        }
        return imp2;
    }
}
