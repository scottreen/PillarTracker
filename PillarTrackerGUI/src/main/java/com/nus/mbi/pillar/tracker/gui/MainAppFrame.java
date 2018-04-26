package com.nus.mbi.pillar.tracker.gui;

import java.awt.Dimension;
import java.awt.Toolkit;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.IOException;
import javafx.application.Platform;
import javafx.embed.swing.JFXPanel;
import javafx.fxml.FXMLLoader;
import javafx.scene.Scene;
import javafx.scene.layout.TilePane;
import javax.swing.JFrame;

/**
 * This class is called from the ImageJ plugin.
 *
 * @author Xiaochun Xu
 */
public class MainAppFrame extends JFrame {

    private JFXPanel fxPanel;
    private boolean implicitExit = false;
    //private Stage primaryStage;
    /*
    public MainAppFrame(ImageJ ij) {
        ij.context().inject(this);
        this.ij = ij;
    }
    */
    
    public void setPlatformImplicitExit(boolean implicit_exit){
        implicitExit = implicit_exit;
    }
    
    public MainAppFrame() {
        
    }
    
    /**
     * Create the JFXPanel that make the link between Swing (IJ) and JavaFX plugin.
     */
    public void init() {
        this.fxPanel = new JFXPanel();
        this.add(this.fxPanel);
        
        Platform.setImplicitExit(implicitExit);
        this.addWindowListener(new WindowAdapter(){
           public void windowClosing(WindowEvent e)
           {               
             remove(fxPanel);      
             fxPanel = null;
             //fxPanel.setVisible(false);   
             //dispose();  
           }
        });
        
        this.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);//DISPOSE_ON_CLOSE);//EXIT_ON_CLOSE);
        this.setVisible(true);
        
        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
	final int WIDTH = screenSize.width;
	final int HEIGHT = screenSize.height;
        this.setLocation(WIDTH/2, 0);
        
        //fxPanel.disdispatchEvent(new WindowEvent(this, WindowEvent.WINDOW_CLOSING));
        //
        //initFX(fxPanel);
        
        // The call to runLater() avoid a mix between JavaFX thread and Swing thread.
        Platform.runLater(new Runnable() {
            @Override
            public void run() {
                initFX(fxPanel);
            }
        });
                
    }

    public void initFX(JFXPanel fxPanel) {
        // Init the root layout
        try {
            FXMLLoader loader = new FXMLLoader();
            loader.setLocation(MainAppFrame.class.getResource("RootLayout.fxml"));
            TilePane rootLayout = (TilePane) loader.load();
            
            // Show the scene containing the root layout.
            Scene scene = new Scene(rootLayout);
            this.fxPanel.setScene(scene); 
            this.fxPanel.setVisible(true); 
                                  
            // Resize the JFrame to the JavaFX scene
            this.setSize((int) scene.getWidth(), (int) scene.getHeight());                    
            
            // Get the controller and add an ImageJ context to it.
            RootLayoutController controller = loader.getController();
            controller.setOwnerWindow(this);
            
        } catch (IOException e) {
            //e.printStackTrace();
            //IJ.log("unable to load the fxml file");
        }
    }
    
    
    public static void main(String[] args) {
        /*
        SwingUtilities.invokeLater(new Runnable() {
            @Override public void run() {
                abstractOptionsFx.dispose();
            }
        });
        */
        MainAppFrame app = new MainAppFrame();        
        app.init();
    }
    
}
