package com.nus.mbi.pillar.tracker.gui;

import com.nus.mbi.pillar.tracker.Compute_Drift_Phase_Correlation_Plugin;
import com.nus.mbi.pillar.tracker.pillar_detector_match_filter;
import com.nus.mbi.pillar.tracker.pillar_tracking;
import com.nus.mbi.pillar.detection.ContractionUnit;
import com.nus.mbi.pillar.detection.CrossCorreclation_Plugin;
import com.nus.mbi.pillar.detection.HighPassFFT;
import com.nus.mbi.pillar.tracker.PSF_Extraction_Plugin;
import com.nus.mbi.pillar.tracker.SetMaskFilterRadius_Plugin;
import com.nus.mbi.pillar.tracker.pillar_tracking_FD;
import com.nus.mbi.pillar.detection.Threshold_Multi_Points;
import com.nus.mbi.pillar.drift.ComputeDrift_Phase_Plugin;
import com.nus.mbi.pillar.drift.DriftAnalysisPlotter;
import com.nus.mbi.pillar.drift.DriftDataLoader;
import com.nus.mbi.pillar.drift.FileHeaderDrift;
import com.nus.mbi.pillar.drift.Hunt_CU_Frame_PlugIn;
import com.nus.mbi.pillar.drift.Hunt_CU_PlugIn;
import com.nus.mbi.pillar.grid.GridAnalysisPlotter;
import com.nus.mbi.pillar.grid.GridCreation;
import com.nus.mbi.pillar.grid.GridDataLoader;
import com.nus.mbi.pillar.grid.GridDeflection;
import com.nus.mbi.pillar.grid.GridFileHeader;
import com.nus.mbi.pillar.grid.GridProperty;
import com.nus.mbi.pillar.grid.Grid_Creation_Plugin;
import com.nus.mbi.pillar.grid.ReCreate_Grid_Plugin;
import com.nus.mbi.pillar.grid.SetGridDialogResult;

import fiji.util.gui.GenericDialogPlus;
import ij.IJ;

import java.net.URL;
import java.util.ResourceBundle;
import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.fxml.Initializable;
import javafx.scene.control.Button;

import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.io.FileInfo;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.filter.RankFilters;
import ij.process.FloatPolygon;
import ij.process.ImageProcessor;
import ij.util.ThreadUtil;
import java.awt.Color;
import java.awt.Desktop;
import java.awt.Window;

import java.io.File;
import java.net.URI;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.Optional;
import javafx.application.Platform;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.event.EventHandler;
import javafx.scene.control.Alert;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.control.ButtonType;
import javafx.scene.control.CheckBox;
import javafx.scene.control.CheckMenuItem;
import javafx.scene.control.ChoiceBox;
import javafx.scene.control.ColorPicker;
import javafx.scene.control.ComboBox;
import javafx.scene.control.ContentDisplay;
import javafx.scene.control.Hyperlink;
import javafx.scene.control.Label;
import javafx.scene.control.MenuItem;
import javafx.scene.control.SingleSelectionModel;
import javafx.scene.control.Slider;
import javafx.scene.control.Tab;
import javafx.scene.control.TabPane;
import javafx.scene.control.TextArea;
import javafx.scene.control.TextField;
import javafx.scene.control.TextFormatter;
import javafx.scene.control.TitledPane;
import javafx.scene.control.ToggleButton;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.scene.input.DragEvent;
import javafx.scene.input.Dragboard;
import javafx.scene.input.KeyCode;
import javafx.scene.input.KeyCodeCombination;
import javafx.scene.input.KeyCombination;
import javafx.scene.input.KeyEvent;
import javafx.scene.input.MouseButton;
import javafx.scene.input.MouseEvent;
import javafx.scene.input.TransferMode;
import javafx.scene.layout.AnchorPane;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.Priority;
import javafx.util.converter.DoubleStringConverter;
import javafx.util.converter.IntegerStringConverter;
import javax.swing.JFileChooser;
import javax.swing.SwingUtilities;
import javax.swing.filechooser.FileNameExtensionFilter;
import com.nus.mbi.pillar.solver.PillarTrackerGPULibrary;
import com.nus.mbi.pillar.stat.BasicStatisitic;
import com.nus.mbi.pillar.stat.GaussianFitter;
import com.nus.mbi.pillar.stat.StatisitcCenterDistance;
import ij.plugin.filter.MaximumFinder;
import ij.plugin.filter.PlugInFilterRunner;

/**
 * FXML Controller class
 *
 * @author Xiaochun Xu
 */
public class RootLayoutController implements Initializable {    

    @FXML private TabPane tabs_pane;
    @FXML private AnchorPane tab_drift;
    @FXML private AnchorPane tab_grid;
    KeyCombination ctrlP = KeyCodeCombination.keyCombination("Ctrl+P");
    KeyCombination ctrlT = KeyCodeCombination.keyCombination("Ctrl+T");
    KeyCombination ctrlD = KeyCodeCombination.keyCombination("Ctrl+D");
    KeyCombination ctrlG = KeyCodeCombination.keyCombination("Ctrl+G");       
    KeyCombination ctrlR = KeyCodeCombination.keyCombination("Ctrl+R");
    
    // properties
    private DoubleStringConverter double_converter = new DoubleStringConverter();
    private Properties settings = new Properties();
    @FXML private AnchorPane pane_properties;
    @FXML private Button button_preview_grid;
    @FXML private TextField textfiled_spacing;
    @FXML private TextField textfiled_grid_oblique;
    @FXML private TextField textfiled_grid_angle;
    @FXML private TextField textfiled_diameter;
    @FXML private TextField textfiled_gauss_sigma;
    @FXML private TextField textfiled_pixel_size;
    @FXML private CheckBox checkbox_dark_pillar;
    @FXML private CheckBox checkbox_pixel_unit;
    @FXML private CheckBox checkbox_square_grid;
    @FXML private CheckBox checkbox_apply_enhancement;
    @FXML private CheckBox checkbox_show_current_kernel;
    @FXML private CheckMenuItem meun_apply_rank_filter;
    @FXML private TextField textfiled_unit_nm;
    @FXML private TextField textfiled_unit_pixel;
    
    @FXML private ToggleButton toggle_converter_back;
    @FXML private ToggleButton toggle_converter_forward;
    
    @FXML private Button button_FFT;
    @FXML private Button button_create_movie_drift;
    @FXML private MenuItem meun_export_grid_format;
    @FXML private Button button_pillar_detection;    
    @FXML private Button button_remove_overlays;
    @FXML private Button button_cross_correlation;
    @FXML private Button button_extract_PSF;
    @FXML private Button button_set_PSF;
    @FXML private Button button_show_current_kernel;
    //ImagePlus current_psf_ip = null;
    CrossCorreclation_Plugin image_enhancer_plugin = null;  
    HighPassFFT fft = null;
    //private Overlay FFT_Overlay = null;
    
    // tracking
    @FXML private ComboBox combo_images;
    @FXML private ChoiceBox choicebox_localization_algorithm;
    @FXML private Button button_do_tracking;
    @FXML private Button button_browse_save_track_file;
    @FXML private TextField textfiled_searchwindow;
    @FXML private TextField textfiled_constrain_radius;
    @FXML private TextField textfiled_max_drift;
    @FXML private TextField textfiled_num_threads;
    @FXML private TextField textfiled_output_fname;   
    @FXML private CheckBox checkbox_apply_pillar_recontruction;
    @FXML private CheckMenuItem checkmenu_show_fft;
    @FXML private CheckMenuItem checkmenu_show_recontruction;
    @FXML private CheckMenuItem checkmenu_drift_correction;
    @FXML private CheckMenuItem checkmenu_create_grid_poly_fit;
    
    private List<Integer> list_image_ids = null;
    private JFileChooser fileChooser_track_fname = new JFileChooser();
    private String tracker_save_fname;
    
    private int localization_algorithm = pillar_tracking.localization_algorithm_Levmar;
    
    // drift analysis     
    @FXML private Button button_drift;
    @FXML private Button button_std;
    @FXML private Button button_hist;
    @FXML private Button button_hunt_cu;
    @FXML private Button button_tracks_map;
    @FXML private Button button_load_file_drift;
    @FXML private TitledPane accordion_drift;    
    @FXML private Label label_drift_load_info; 
    //@FXML private Button button_get_pillars_list;
    @FXML private CheckBox checkbox_raw;
    @FXML private CheckBox checkbox_flatten;
    @FXML private CheckBox checkbox_show_labels;
    @FXML private CheckBox checkbox_still_pillars;
    //@FXML private CheckBox checkbox_show_CU;
    
    @FXML private Slider slider_zoomin_tracks;
    @FXML private Slider slider_zoomin_images;
    @FXML private TextField textfield_zoomin_tracks;
    @FXML private TextField textfield_zoomin_images;
    
    @FXML private TextField textfield_pillar_no_drift; 
    @FXML private Slider slider_pillar_no_drift;
    @FXML private TextField textfield_frame_window;
    @FXML private TextField textfield_frame_no_drift;
    @FXML private Slider slider_frame_no_drift;
    @FXML private Button button_defelction_drift;
    @FXML private CheckMenuItem checkbox_shownin;
    @FXML private CheckMenuItem checkbox_show_track_head; 
    @FXML private TextField textfile_drift_fname;
    @FXML private Button button_drift_file_browse;
    @FXML private Button button_plot_pair_seperation;
    @FXML private Button button_draw_list_pillars;
    private int[] selected_pillars_indexs = null;
    PairPillars pair_pillars = new PairPillars();
    private String list_pillar_string = "";
    private IntegerStringConverter int_converter = new IntegerStringConverter();
    //private TextFormatter int_formater = new TextFormatter(int_converter);
    
    private Hunt_CU_PlugIn cu_hunter = new Hunt_CU_PlugIn();
    private Hunt_CU_PlugIn cu_hunter_grid = new Hunt_CU_PlugIn();
    private ComputeDrift_Phase_Plugin cdpp = new ComputeDrift_Phase_Plugin();
    
    private double pixel_size = 1.0;
    private DriftDataLoader reader = new DriftDataLoader();
    private boolean load_suc_drift = false;
    private DriftAnalysisPlotter plotter = new DriftAnalysisPlotter();
    
    private String drift_fname = "";//"H:\\pillar\\2015-06-29-test\\75fps_cell2_dish2_large-output\\75fps_cell2_dish2_large_X11-pillar-tracks-big.bin";//"E:/pillar-tracks-big.bin"; ;
    private ImagePlus flatten_tracks_ip = null;
    private ImagePlus overlay_tracks_ip = null;
    private double zoom_in_tracks = 15;
    private double zoom_in_image = 6;
    
    private PlotWindow plot_window_dxy = null;
    private PlotWindow plot_window_dis = null;
    
    private boolean need_recreate_map_drift = true;
    private ImagePlus created_map_ip = null;
    private boolean show_drift_dxy = false;
    //private FileChooser drift_fileChooser = new FileChooser();
    private JFileChooser drift_fileChooser = new JFileChooser();
    
    // grid analysis
    @FXML private TitledPane accordion_grid;    
    @FXML private Label label_grid_load_info;     
    @FXML private Button button_load_file_grid;
    @FXML private Button button_deflection_map;
    @FXML private Slider slider_zoomin_deflection;
    @FXML private Slider slider_zoomin_image_grid;
    @FXML private TextField textfield_zoomin_deflection;
    @FXML private TextField textfield_zoomin_image_grid;
    @FXML private Slider slider_frame_no_grid;
    @FXML private TextField textfield_frame_no_grid;    
    @FXML private CheckMenuItem checkbox_shownin_grid; 
    @FXML private CheckBox checkbox_show_labels_grid;
    
    @FXML private ChoiceBox choice_crosshair;
    @FXML private ChoiceBox choice_arrows;
    @FXML private ChoiceBox choice_labels;    
    @FXML private ColorPicker color_crosshair;
    @FXML private ColorPicker color_arrows;
    @FXML private ColorPicker color_labels;
    
    @FXML private CheckBox checkbox_remove_large;
    @FXML private TextField textfield_large_threshold;
    
    @FXML private Button button_get_data_frame;
    @FXML private Button button_get_data_pillar;
    @FXML private TextField textfield_pillar_no_grid; 
    @FXML private Slider slider_pillar_no_grid;
    @FXML private CheckBox checkbox_show_linegrid;
    //@FXML private Button button_set_grid_properites;
    @FXML private TextField textfile_grid_fname;
    @FXML private Button button_grid_file_browse;
    @FXML private Button button_create_tiff_movie;
    @FXML private Button button_recreate_grid;
    @FXML private Button button_hunt_cu_grid;
    private double threshold_deflection_dxy = Double.NaN;
    private String grid_fname = "";//"H:\\pillar\\2015-06-29-test\\75fps_cell2_dish2_large-output\\75fps_cell2_dish2_large_X11-pillar-tracks-big.bin.dxy";
    private double zoom_in_deflection = 15;
    private double zoom_in_image_grid = 1;    
    private GridDataLoader grid_reader = new GridDataLoader();
    private GridAnalysisPlotter grid_plotter = new GridAnalysisPlotter();
    private GridCreation grid_recreator = null;
    private boolean load_suc_grid = false;
    private boolean need_recreate_map_grid = true;
    private int current_frame_index = -1;
    private ImagePlus overlay_deflections_ip = null;
    private ImagePlus created_deflection_map_ip = null;
    private boolean show_grid_dxy = false;
    private PlotWindow plot_window_dxy_grid = null;
    private PlotWindow plot_window_dis_grid = null;
    
    private int[] candidate_grid = null;
    private int[] contraction_unit_list = null;
    List<Integer> cross_size = Arrays.asList(0,1,2,3,4,5,6,7,8);
    List<Integer> arrow_size = Arrays.asList(1,2,3,4,5,6,7,8,9,10);
    List<Integer> label_size = Arrays.asList(1,2,3,4,5,6,7,8,9,10,11,12);
    
    //about        
    @FXML private Hyperlink link_mbi;
    @FXML private Hyperlink link_update_gpu_lib;
    @FXML private Hyperlink link_copyright; 
    @FXML private Hyperlink link_version;
    
    //private FileChooser grid_fileChooser = new FileChooser();
    private JFileChooser grid_fileChooser = new JFileChooser();
    //private Window owner_stage = null;
    private Window owner_window = null;
    private String grid_list_pillar_string;
    private List<ContractionUnit> cu_list_grid;
    public void setOwnerWindow(Window s){
        owner_window = s;
    }
    
    @FXML private void shortcutKeyTabs(KeyEvent event){               
       //KeyCode key = event.getCode(); 
        SingleSelectionModel<Tab> tab = tabs_pane.getSelectionModel();
        int current_tab = tab.getSelectedIndex();

        if (ctrlP.match(event)) {
             tab.select(0);
         }
         else if(ctrlT.match(event)){
             tab.select(1);
         }
         else if(ctrlD.match(event)){
             tab.select(2);
         }
         else if(ctrlG.match(event)){
             tab.select(3);
         }
         else if (ctrlR.match(event)) {
             if(current_tab==3)
                 handleButtonLoadFileGrid(null);
             else if(current_tab==2)
                 handleButtonLoadFileDrift(null);
             else if(current_tab==1)
                 handleButtonDoTracking(null);
         }
    }
    
    @FXML private void shortcutKeyPanes(KeyEvent event){               
        KeyCode key = event.getCode(); 
        SingleSelectionModel<Tab> tab = tabs_pane.getSelectionModel();
        int current_tab = tab.getSelectedIndex();
        if (ctrlP.match(event)) {
             tab.select(0);
         }
         else if(ctrlT.match(event)){
             tab.select(1);
         }
         else if(ctrlD.match(event)){
             tab.select(2);
         }
         else if(ctrlG.match(event)){
             tab.select(3);
         }
         else if(key==KeyCode.RIGHT || key==KeyCode.DOWN){
             int num_tab = tabs_pane.getTabs().size();        
             tab.select((current_tab+1)%num_tab); 
         }else if(key==KeyCode.LEFT || key==KeyCode.UP){
             int num_tab = tabs_pane.getTabs().size();        
             tab.select((current_tab-1)%num_tab);            
         }else if (ctrlR.match(event)) {
             if(current_tab==3)
                 handleButtonLoadFileGrid(null);
             else if(current_tab==2)
                 handleButtonLoadFileDrift(null);
             else if(current_tab==1)
                 handleButtonDoTracking(null);
         }
    }
    
    // <editor-fold defaultstate="collapsed" desc="properties events handlding">
    @FXML private void handleToggleForward(ActionEvent event){
        String text = textfiled_unit_nm.getText();
        if(text.isEmpty()) return;
        double nm = double_converter.fromString(text);
        textfiled_unit_pixel.setText(String.format("%.3f", settings.convertNanoToPixel(nm)));
    }
    
    @FXML private void handleToggleBackward(ActionEvent event){
        String text = textfiled_unit_pixel.getText();
        if(text.isEmpty()) return;
        double pixel = double_converter.fromString(text);
        textfiled_unit_nm.setText(String.format("%.3f", settings.convertPixelToNano(pixel)));
    }
    
    @FXML private void handleButtonPreviewGrid(ActionEvent event){
        
        SwingUtilities.invokeLater(new Runnable() {
            @Override public void run() {
                showPrevewGridDialog();
            }
        });
        
    }
    
    @FXML private void handleButtonRemoveOverlays(ActionEvent event){
        int n = WindowManager.getImageCount();
        if(n<1){
            //IJ.showMessage("there is no image opened!");
            return;
        }
        
        IJ.getImage().setOverlay(null);
    }
    
    @FXML private void handleButtonFFT(ActionEvent event){
        int n = WindowManager.getImageCount();
        if(n<1){
            IJ_showMessage("there is no image opened!");                    
            return;
        }
        
        ImagePlus ip = IJ.getImage();
        Calibration cal = ip.getLocalCalibration();
        if(cal!=null && cal.scaled()){
            IJ.log("Selected image for FFT has calibration information");
            IJ.log("    pixel width  = " + cal.pixelWidth + " " +cal.getUnits());
            IJ.log("    pixel height = " + cal.pixelHeight + " " +cal.getUnits());
        }
        if(fft==null) fft = new HighPassFFT();
        fft.run("");
//        
//        FFT fft = new FFT();
//        fft.run("");         
    }
    
    @FXML private void handleButtonHighPassFFT(ActionEvent event){
        int n = WindowManager.getImageCount();
        if(n<1){
            IJ_showMessage("there is no image opened!");                    
            return;
        }
        
        if(fft==null || fft.isEmptyMask()){
            IJ_showMessage("The mask has not been set yet!");                    
            return;
        }
                
        Thread worker = new Thread() {
            public void run() {		            	            
                 try {
                   // Something that takes a long time . . . in real life,
                   fft.process();  
                 } catch (Exception ex) {
                 }
               }
        };
        worker.start();
//        ImagePlus ip = IJ.getImage();        
//        
//        if(FFT_Overlay==null) return;
//        
////        FHT fht = new FHT(ip.getProcessor());
////        fht.setShowProgress(false);
////        fht.transform();
////        ImageProcessor fft_img = fht.getPowerSpectrum();
////        fft_img.setColor(255);
////        for(int i=0; i<FFT_Overlay.size(); i++) fft_img.fill(FFT_Overlay.get(i));            
////        FHT ifft = new FHT(fft_img, true);
////        ifft.inverseTransform();          
////        ip.setProcessor(ifft);
////        ip.updateAndDraw();
//        ImageStack stack = ip.getImageStack();
//        int nslice = stack.getSize();
//        ImageStack is = new ImageStack(ip.getWidth(), ip.getHeight());
//        int img_id = ip.getID();
//        FFT fft = new FFT();
//        for(int s=0; s<nslice; s++){            
//            IJ.selectWindow(img_id);
//            ip.setSlice(s+1);            
//            fft.run(""); 
//            ImagePlus fft_ip = IJ.getImage();
//            ImageProcessor img = fft_ip.getProcessor();
//            for(int i=0; i<FFT_Overlay.size(); i++){
//                img.setColor(255);
//                img.fill(FFT_Overlay.get(i));
//            }
//            fft.run("Inverse FFT"); 
//            ImagePlus ifft_ip = IJ.getImage();
//            is.addSlice(ifft_ip.getProcessor());
//            fft_ip.close();
//            ifft_ip.close();
//        }        
//        ImagePlus new_ip = new ImagePlus("High_Pass-" + ip.getShortTitle(), is);
//        new_ip.show();
    }
    
    @FXML private void handleButtonSetMaskHighPass(ActionEvent event){
        int n = WindowManager.getImageCount();
        if(n<1){
            IJ_showMessage("there is no image opened!");                    
            return;
        }
        
        ImagePlus imp = IJ.getImage();                
        if (!HighPassFFT.isInFrequencyDomain(imp)) {
            IJ_showMessage("Frequency domain image required");
            return;
        }
        
        if(fft==null) fft = new HighPassFFT();
//        if(fft.showHighPassDialog(imp)){    
//            fft.run("");
//            //FFT fft = new FFT();
//            //fft.run("Inverse FFT");  
//        }
        Thread worker = new Thread() {
            public void run() {		            	            
                 try {
                   // Something that takes a long time . . . in real life,
                    //if(fft==null) fft = new HighPassFFT();
                    if(fft.showHighPassDialog(imp)) fft.run(""); 
                 } catch (Exception ex) {
                 }
               }
        };

        worker.start(); // So we don't hold up the dispatch thread. 
        
    }

    
    @FXML private void handleButtonComputeDriftPhaseDiff(ActionEvent event){
        int n = WindowManager.getImageCount();
        if(n<1){
            IJ_showMessage("there is no image opened!");                    
            return;
        }
        
        //ImagePlus imp = IJ.getImage();                
        
        Thread worker = new Thread() {
            public void run() {		            	            
                 try {                    
                    cdpp.run(""); 
                 } catch (Exception ex) {
                 }
               }
        };

        worker.start(); // So we don't hold up the dispatch thread. 
        
    }

    
    private void pixel_unit_change(boolean is_pixel_unit){
        int n = WindowManager.getImageCount();
        if(n<1){
            //if(!is_pixel_unit) IJ.showMessage("there is no image opened!");
            return;
        }
        
        ImagePlus ip = IJ.getImage();
        Calibration cal = ip.getLocalCalibration();
        if(cal.scaled()){
            String unit = cal.getUnit();
            double v = cal.pixelWidth;
            if("um".equalsIgnoreCase(unit)) v = v * 1000;
            else if("micron".equalsIgnoreCase(unit)) v = v * 1000;
            else if("mm".equalsIgnoreCase(unit)) v = v*1e6;
            else if(!("nm".equalsIgnoreCase(unit))) v = Double.NaN;
            if(!Double.isNaN(v)) this.textfiled_pixel_size.setText(String.format(Locale.ENGLISH, "%.2f", v));
        }  
        
        Calibration global_cal = ImagePlus.getStaticGlobalCalibration();
        if(global_cal==null || global_cal.scaled()){
            Calibration cal_null = new Calibration();
            //cal_null.
            ip.setGlobalCalibration(cal_null);
        }    
        
        if(is_pixel_unit){                    
            //ip.setCalibration(new Calibration());                    
            ip.setIgnoreGlobalCalibration(false);
        }
        else{
            //ip.setCalibration(cal);
            ip.setIgnoreGlobalCalibration(true);
        }
        
        ip.updateAndRepaintWindow();
    }
    
    private double round90(double degree){
        int t = (int)Math.floor(degree/90);
        double r = degree - t*90;
        if(r>45) r = r-90;        
        return r;
    }
    
    private double round180(double degree){
        int t = (int)Math.floor(degree/180);
        double r = degree - t*180;
        if(r>90) r = 180 - r;        
        return r;
    }
    
    private boolean checkSettings(boolean show_message){
        String text = textfiled_spacing.getText();
        if(text.isEmpty()) {
            if(show_message) IJ_showMessage("the spacing is not valid!");
            return false;
        }else{
            double s = double_converter.fromString(text);
            if(s<0){
                if(show_message) IJ_showMessage("The spacing must be larger than 0!");
                return false;
            }
            settings.spacing = s;
        }
        
        text = textfiled_diameter.getText();
        if(text.isEmpty()) {
            if(show_message) IJ_showMessage("the gaussian simga is not valid!");
            return false;
        }else{
            double s = double_converter.fromString(text);
            if(s<0){
                if(show_message) IJ_showMessage("The gaussian simga must be larger than 0!");
                return false;
            }
            settings.diameter = s;
        }          
        
        text = textfiled_gauss_sigma.getText();
        if(text.isEmpty()) {
            if(show_message) IJ_showMessage("the gaussian simga is not valid!");
            return false;
        }else{
            double s = double_converter.fromString(text);
            if(s<0){
                if(show_message) IJ_showMessage("The gaussian simga must be larger than 0!");
                return false;
            }
            settings.gauss_sigma = s;
        }             
        
        text = textfiled_grid_oblique.getText();
        if(text.isEmpty()) {
            if(show_message) IJ_showMessage("the grid oblique is not valid!");
            return false;
        }else{
            double s = double_converter.fromString(text);
            s = round90(s);
            if(s>45 || s<-45){
                if(show_message) IJ_showMessage("The grid_oblique must be between (-45, 45)!");
                return false;
            }
            settings.grid_oblique = s;
        }    
        
        text = textfiled_grid_angle.getText();
        if(text.isEmpty()) {
            if(show_message) IJ_showMessage("the grid angle is not valid!");
            return false;
        }else{
            double s = double_converter.fromString(text);
            s = round180(s);
            if(s<=0 || s>90){
                if(show_message) IJ_showMessage("The angle must be between (0, 90]!");
                return false;
            }
            settings.grid_angle = s;
        }    
        
        settings.dark_pillar = this.checkbox_dark_pillar.isSelected();
        
        return true;
    }
    
    private void IJ_showMessage(String msg){
        SwingUtilities.invokeLater(new Runnable() {
            @Override public void run() {                    
                IJ.showMessage(msg);
            }
        });   
    }
    
    private boolean checkSettings(){ 
        return checkSettings(true); 
    }
    
    private void setDefaultParametersForTracking(boolean check_old_value){
       
        int max_drift = (int)Math.floor(settings.spacing*0.5);
        if(max_drift>0){            
            boolean need = true;
            if(check_old_value){
                String text = textfiled_max_drift.getText();
                if(!text.isEmpty()){
                    double s = int_converter.fromString(text);
                    if(s>0 && s<settings.spacing*0.5) need = false;                
                }
            }
            if(need) textfiled_max_drift.setText(int_converter.toString(max_drift));
        }
        
        int constrain_radius = (int)Math.round(settings.gauss_sigma*1.5);
        if(constrain_radius>settings.spacing*0.5) constrain_radius = (int)Math.floor(settings.spacing*0.5);
        if(constrain_radius>0){
            boolean need = true;
            if(check_old_value){
                String text = textfiled_constrain_radius.getText();

                if(!text.isEmpty()){
                    int s = int_converter.fromString(text);
                    if(s>0 && s<settings.spacing*0.5) need = false;                
                }
            }
            if(need) textfiled_constrain_radius.setText(int_converter.toString(constrain_radius));
        }
        
//        int search_window = (int)Math.round(settings.gauss_sigma*1.5)*2+1;
//        if(localization_algorithm != localization_algorithm_CG) search_window = (int)Math.floor(settings.spacing*0.5-0.5)*2+1;
        int search_window = (int)Math.floor(settings.spacing*0.5-0.5)*2+1;
        if(search_window>0){
            boolean need = true;
            if(check_old_value){
                String text = textfiled_searchwindow.getText();                
                if(!text.isEmpty()){
                    int s = int_converter.fromString(text);
                    if(s>0 && s<settings.spacing) need = false;                
                }
            }
            
            if(need) textfiled_searchwindow.setText(int_converter.toString(search_window));            
        }
        
    }
    
    @FXML private void handleMenuResetParameters(){
        //if(!checkSettings(false)) return;
        setDefaultParametersForTracking(false);
    }
    
    @FXML private void handleButtonDetection(ActionEvent event){        
        int n = WindowManager.getImageCount();
        if(n<1){
            IJ_showMessage("there is no image opened!");
            return;
        }
        
        ImagePlus ip = IJ.getImage();
        if(ip==null || !ip.isVisible()){
            IJ_showMessage("there is no image opened to draw on!");
            return;
        }
        
        if(!checkSettings()) return;
        setDefaultParametersForTracking(false);
        
        boolean apply_enhancer = this.checkbox_apply_enhancement.selectedProperty().get();             
        
        pillar_detector_match_filter detector = new pillar_detector_match_filter();
        detector.setShowDialog(false);
        detector.setup(ip,settings.spacing, settings.gauss_sigma, settings.grid_oblique, settings.grid_angle, settings.dark_pillar);
        //detector.run("");
        
        ImageProcessor img = ip.getProcessor();
        
        boolean apply_rank_filter = this.meun_apply_rank_filter.isSelected();   
        if(apply_rank_filter){
            RankFilters mean_filter = new RankFilters();
            ImageProcessor slice_filter = img.duplicate();
            mean_filter.rank(slice_filter, 1, RankFilters.MEAN);
            img = slice_filter;
        }
        
        boolean use_minimum_error = true;
        if(apply_enhancer && this.image_enhancer_plugin!=null){
            //ImageProcessor img_enhanced = image_enhancer_plugin.process(img);
            img = image_enhancer_plugin.process(img);
            detector.setDarkPillar(false);
            use_minimum_error = false;
            //detector.process(img_enhanced, false);            
        }//else detector.process(img, true);
        detector.process(img, use_minimum_error);
        
//        int num_detect_max = detector.getNumLocalMaximas();        
//        int num_local_max = detector.getNumLocalMaximas(img, 0);//detector.getNumDetectedMaximas();
//        IJ.log("local maximas:"+num_local_max+"-"+num_detect_max);
//        if(num_local_max>num_detect_max*2){
//            //IJ.log("local maximas are overwhelming:"+num_local_max+"-"+num_detect_max);
//            LocalMaximaThreashold lmt = new LocalMaximaThreashold(img, 0, detector);
//            Thread worker = new Thread() {
//                public void run() {		            	            
//                     try {
//                       // Something that takes a long time . . . in real life,
//                       boolean suc = lmt.showDialog();
//                       if(suc){
//                           IJ.log("set the threshold="+lmt.threshold);
//                       }
//                     } catch (Exception ex) {
//                     }
//                   }
//                 };
//
//            worker.start(); // So we don't hold up the dispatch thread.
//        }
    }
    
    @FXML private void handleButtonShowCurrentKernel(ActionEvent event){        
        if(this.image_enhancer_plugin!=null){
            ImagePlus current_psf_ip = image_enhancer_plugin.getPSF();             
            if(current_psf_ip!=null) current_psf_ip.show();
        }
    }
    
    @FXML private void handleButtonExtactPSF(ActionEvent event){        
        int n = WindowManager.getImageCount();
        if(n<1){
            IJ_showMessage("there is no image opened!");
            return;
        }
        
        ImagePlus ip = IJ.getImage();
        if(ip==null || !ip.isVisible()){
            IJ_showMessage("there is no image opened to draw on!");
            return;
        }
        
//        Roi point_roi= ip.getRoi();			
//        if(point_roi==null || point_roi.isArea() || point_roi.isLine()){
//            IJ_showMessage("there is no point selections");
//            return;
//        }
        
//        String text = textfiled_spacing.getText();
//        if(!text.isEmpty()){
//            double s = double_converter.fromString(text);
//            if(s>0) settings.spacing = s;
//        }
        checkSettings(false);
               
        PSF_Extraction_Plugin psf_extractor = new PSF_Extraction_Plugin();
        int kernel_radius = (int)Math.floor(settings.spacing/2-0.5);
        psf_extractor.setup(ip, kernel_radius);
        psf_extractor.dark_pillar = settings.dark_pillar;
        psf_extractor.sigma = settings.gauss_sigma;
        Thread worker = new Thread() {
            public void run() {		            	            
                 try {
                   // Something that takes a long time . . . in real life,
                   psf_extractor.run("once");
                 } catch (Exception ex) {
                 }
               }
             };

        worker.start(); // So we don't hold up the dispatch thread.         
    }
    
    @FXML private void handleMenuComputeDrift(ActionEvent event){        
        int n = WindowManager.getImageCount();
        if(n<1){
            IJ_showMessage("there is no image opened!");
            return;
        }
        
        ImagePlus ip = IJ.getImage();
        if(ip==null || !ip.isVisible()){
            IJ_showMessage("there is no image opened to draw on!");
            return;
        }
 
        checkSettings(false);
               
        Compute_Drift_Phase_Correlation_Plugin plugin = new Compute_Drift_Phase_Correlation_Plugin();
        int kernel_radius = (int)Math.floor(settings.spacing/2-0.5);
        plugin.setup(ip, kernel_radius);
        //plugin.dark_pillar = settings.dark_pillar;
        plugin.sigma = settings.gauss_sigma;
        Thread worker = new Thread() {
            public void run() {		            	            
                 try {
                   // Something that takes a long time . . . in real life,
                   plugin.run("once");
                 } catch (Exception ex) {
                 }
               }
             };

        worker.start(); // So we don't hold up the dispatch thread.         
    }
    
    @FXML private void handleMenuFindMaxima(ActionEvent event){
        int n = WindowManager.getImageCount();
        if(n<1) return;
        
        ImagePlus ip = IJ.getImage();
        if(ip==null || !ip.isVisible()) return;
        
        Thread worker = new Thread() {
            public void run() {		            	            
                 try {                   
                   IJ.run("Find Maxima...");                   
                 } catch (Exception ex) {
                 }
               }
             };

        worker.start(); // So we don't hold up the dispatch thread.         
    }
    
    @FXML private void handleMenuThresholdPoints(ActionEvent event){
        int n = WindowManager.getImageCount();
        if(n<1) return;
        
        ImagePlus ip = IJ.getImage();
        if(ip==null || !ip.isVisible()) return;
        
        Threshold_Multi_Points thresholder = new Threshold_Multi_Points();        
        Thread worker = new Thread() {
            public void run() {		            	            
                 try {
                   // Something that takes a long time . . . in real life,
                   thresholder.showDialog(ip);
                 } catch (Exception ex) {
                 }
               }
             };

        worker.start(); // So we don't hold up the dispatch thread.         
    }
            
    @FXML private void handleButtonCrossCorrelation(ActionEvent event){        
        int n = WindowManager.getImageCount();
        if(n<1){
            IJ_showMessage("there is no image opened!");
            return;
        }
        
        ImagePlus ip = IJ.getImage();
        if(ip==null || !ip.isVisible()){
            IJ_showMessage("there is no image opened!");
            return;
        }
        
        String text = textfiled_spacing.getText();
        if(!text.isEmpty()){
            double s = double_converter.fromString(text);
            if(s>0) settings.spacing = s;
        }
        
        text = textfiled_gauss_sigma.getText();
        if(!text.isEmpty()){
            double s = double_converter.fromString(text);
            if(s>0) settings.gauss_sigma = s;
        }
        
        settings.dark_pillar = this.checkbox_dark_pillar.isSelected();
        
        int num_threads = 1;
        text = textfiled_num_threads.getText();
        if(!text.isEmpty()) {
            int s = int_converter.fromString(text);
            if(s>0) num_threads = s;
        }             
        
        CrossCorreclation_Plugin image_enhancer = new CrossCorreclation_Plugin();
        if(image_enhancer_plugin==null){
            image_enhancer.setup(ip, true, true);
            image_enhancer.setNumThreads(num_threads);
            int kernel_radius = (int)Math.floor(settings.spacing/2-0.5);
            image_enhancer.setGaussianParameters(settings.gauss_sigma, kernel_radius, settings.dark_pillar);            
        }
        else{
            image_enhancer.setup(ip, image_enhancer_plugin.getMode() ,true);
            image_enhancer.gpu_lib = image_enhancer_plugin.gpu_lib;
            image_enhancer.setGPUOn(image_enhancer_plugin.isGPUOn());
            image_enhancer.setNumThreads(num_threads);
            image_enhancer.setCustomPSF(image_enhancer_plugin.getPSF());
            //boolean use_gaussian = image_enhancer_plugin.isUsingGaussianPSF();
            //if(use_gaussian) image_enhancer.setGaussianPSF(image_enhancer_plugin.getGaussianSigma(), image_enhancer_plugin.getGaussianRaidus(), image_enhancer_plugin.isDarkObject());
            //else image_enhancer.setCustomPSF(image_enhancer_plugin.getPSF());
        }
                
        Thread worker = new Thread() {
            public void run() {		            	            
                 try {
                   // Something that takes a long time . . . in real life,
                   image_enhancer.run("once");                   
                 } catch (Exception ex) {
                 }
               }
             };

        worker.start(); // So we don't hold up the dispatch thread. 
    }
    
    private void setKernelforEnhancement(boolean update_checkbox_apply){        
        String text = textfiled_spacing.getText();
        if(!text.isEmpty()){
            double s = double_converter.fromString(text);
            if(s>0) settings.spacing = s;
        }
        
        text = textfiled_gauss_sigma.getText();
        if(!text.isEmpty()){
            double s = double_converter.fromString(text);
            if(s>0) settings.gauss_sigma = s;
        }
        
        settings.dark_pillar = this.checkbox_dark_pillar.isSelected();
        int n = WindowManager.getImageCount();
        ImagePlus ip = null;
        if(n>0) ip = IJ.getImage();        
        
        int num_threads = 1;
        text = textfiled_num_threads.getText();
        if(!text.isEmpty()) {
            int s = int_converter.fromString(text);
            if(s>0) num_threads = s;
        }   
        
        CrossCorreclation_Plugin image_enhancer = new CrossCorreclation_Plugin();
        if(image_enhancer_plugin==null){
            image_enhancer.setup(ip, true, true);
            image_enhancer.setNumThreads(num_threads);
            int kernel_radius = (int)Math.floor(settings.spacing/2-0.5);
            image_enhancer.setGaussianParameters(settings.gauss_sigma, kernel_radius, settings.dark_pillar);
        }
        else{
            image_enhancer.setup(ip, image_enhancer_plugin.getMode() ,true);
            image_enhancer.gpu_lib = image_enhancer_plugin.gpu_lib;
            image_enhancer.setGPUOn(image_enhancer_plugin.isGPUOn());
            image_enhancer.setNumThreads(num_threads);
            //image_enhancer.keep_current_kernel = image_enhancer_plugin.keep_current_kernel;
            boolean use_gaussian = image_enhancer_plugin.isUsingGaussianPSF();
            if(use_gaussian) image_enhancer.setGaussianPSF(image_enhancer_plugin.getGaussianSigma(), image_enhancer_plugin.getGaussianRaidus(), image_enhancer_plugin.isDarkObject());
            else{
                int kernel_radius = (int)Math.floor(settings.spacing/2-0.5);
                image_enhancer.setGaussianParameters(settings.gauss_sigma, kernel_radius, settings.dark_pillar);
            }
            image_enhancer.setCustomPSF(image_enhancer_plugin.getPSF());
        }
        
        
        SwingUtilities.invokeLater(new Runnable() {
            @Override public void run() {                 
                boolean suc = image_enhancer.showSetKernelDialog();
                if(suc){                       
//                    if(image_enhancer.keep_current_kernel) {
//                        ImagePlus processed_ip = image_enhancer.process();
//                        if(processed_ip!=null) processed_ip.show();                        
//                    }
                    image_enhancer_plugin=image_enhancer;   
//                    boolean show = checkbox_show_current_kernel.selectedProperty().get();
//                    if(show && image_enhancer_plugin.isUsingGaussianPSF()){         
//                        ImagePlus current_psf_ip = image_enhancer_plugin.getPSF();                    
//                        current_psf_ip.show();                                 
//                    }
                }
                else if(update_checkbox_apply){
                    checkbox_apply_enhancement.selectedProperty().set(false);
                    //checkbox_show_current_kernel.selectedProperty().set(false);
                }
            }
        }); 
    }
    
    @FXML private void handleButtonSetKenel(ActionEvent event){        
        setKernelforEnhancement(false);
    }
    
    private void setMaskFilterRadius(){
        if(fft!=null && fft.hasOffCenterPoints()){
            SetMaskFilterRadius_Plugin smfr = new SetMaskFilterRadius_Plugin();
            smfr.imp = null;
            int n = WindowManager.getImageCount();
            if(n>0){
                ImagePlus ip = IJ.getImage();
                if(HighPassFFT.isInFrequencyDomain(ip)) smfr.imp = ip;
            }
            smfr.low_radius = fft.off_center_radius;
            smfr.larger_mask_radius = fft.end_radius_off_center;
            smfr.step_radius = fft.step_radius_off_center;
            smfr.points = fft.points;                                      
            SwingUtilities.invokeLater(new Runnable() {
                @Override public void run() {                 
                    boolean suc = smfr.showFFTDialog();
                    if(suc){                       
                       //if(smfr.larger_mask_radius>smfr.low_radius) 
                       fft.end_radius_off_center = smfr.larger_mask_radius;
                       fft.step_radius_off_center = smfr.step_radius;
                    }
                    else{
                        checkbox_apply_pillar_recontruction.selectedProperty().set(false);
                        
                    }
                }
            }); 
        }
    }
    
    // </editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="tracking">
    @FXML private void handleButtonBrowseTrackSaveFile(ActionEvent event){
        
        SwingUtilities.invokeLater(new Runnable() {
            @Override public void run() {
                File file = new File(textfiled_output_fname.getText());
                //fileChooser_track_fname.setCurrentDirectory(file);
                fileChooser_track_fname.setSelectedFile(file);
                int retval= fileChooser_track_fname.showSaveDialog(owner_window);        
                if (retval == JFileChooser.APPROVE_OPTION){
                    textfiled_output_fname.setText(fileChooser_track_fname.getSelectedFile().getAbsolutePath());            
                    //handleButtonLoadFileDrift(event);
                }      
            }
        });
              
    }
    
    private List<String> getImageList(){
        int[] wList = WindowManager.getIDList();	
	if (wList==null || wList.length<1) {
		//IJ.showMessage("Pillar Localization", "1 channel time series are required");
		return null;
	}

	int wlistlen = wList.length;
	
        list_image_ids = new ArrayList<Integer>();
	ArrayList<String> titles_list = new ArrayList<String>();
	for (int i=0; i<wlistlen; i++) {
                int id = wList[i];
		ImagePlus imp = WindowManager.getImage(id);
		if(imp!=null && imp.getNDimensions()<=3){
                //if(imp!=null){
                    titles_list.add("["+id+"]:"+imp.getTitle());
                    list_image_ids.add(wList[i]);
                }
		//if(imp!=null && imp.getNDimensions()<=3) titles_list.add(imp.getTitle());
	}

        return titles_list;
    }
    
    private void updateComboImages(){        
        this.combo_images.getItems().clear();
        
        List<String> title_list = getImageList();
        if(title_list==null || title_list.isEmpty()) return;
        
        ObservableList<String> options = FXCollections.observableArrayList(title_list);                
        this.combo_images.getItems().addAll(options);
    }
    
    private ImagePlus getSelectedImage(boolean show_msg){
        int n = WindowManager.getImageCount();
        if(n<1){
            if(show_msg)IJ_showMessage("there is no image opened!");
            return null;
        }
        int index = this.combo_images.getSelectionModel().getSelectedIndex();
        if(index<0 || index>=list_image_ids.size()){
            if(show_msg) IJ_showMessage("there is no image selected!");
            return null;
        }
        
        int id = list_image_ids.get(index);
        ImagePlus ip = WindowManager.getImage(id);
        if(ip==null || !ip.isVisible()){
            if(show_msg) IJ_showMessage("the selected image may have been closed!");
            return null;
        }
        return ip;
    }
    
    
    @FXML private void handleMenuLoadContractileUnits(ActionEvent event){  
        
    }
    
    @FXML private void handleMenuSaveCurrentSettings(ActionEvent event){                 
        
        if(!checkSettings()) return;
        
        int search_window = 0;
        String text = textfiled_searchwindow.getText();
        if(text.isEmpty()) {
            IJ_showMessage("the search_window is not valid!");
            return;
        }else{
            int s = int_converter.fromString(text);
            if(s<0){
                IJ_showMessage("The search_window must be larger than 0!");
                return;
            }
            search_window = s;
        }             
                
        double max_drift = 0;
        text = textfiled_max_drift.getText();
        if(text.isEmpty()) {
            IJ_showMessage("the max_drift is not valid!");
            return;
        }else{
            double s = double_converter.fromString(text);
            if(s<=0 || s>settings.spacing){
                IJ_showMessage("The max_drift must be larger than 0 and small than spacing!");
                return;
            }
            max_drift = s;
        }             
        
        int constrain_radius = 0;
        text = textfiled_constrain_radius.getText();
        if(text.isEmpty()) {
            IJ_showMessage("the constrain_radius is not valid!");
            return;
        }else{
            int s = int_converter.fromString(text);
            if(s<=0 || s>settings.spacing){
                IJ_showMessage("The constrain_radius must be larger than 0 and small than spacing!");
                return;
            }
            constrain_radius = s;
        }             
        
        boolean apply_reconstruction = this.checkbox_apply_pillar_recontruction.selectedProperty().get();        
        if(apply_reconstruction){
            if(fft==null || !fft.hasOffCenterPoints()){
                IJ_showMessage("the mask filter is not set yet!");
                return;
            }
        }
        
        int num_threads = 1;
        text = textfiled_num_threads.getText();
        if(!text.isEmpty()) {
            int s = int_converter.fromString(text);
            if(s>0) num_threads = s;
        }             
        
        ImagePlus ip = getSelectedImage(false);  
        String defaultPath = "";
        if(ip!=null){
            FileInfo fi = ip.getOriginalFileInfo();
            if(fi != null && fi.directory != null && fi.fileName!=null){
                String fname = ip.getShortTitle();
                if(fname==null) fname = fi.fileName;
                String info = fi.directory + fname;
                defaultPath = info + ".profile";
            }
        }        
        
        boolean apply_enhancer = checkbox_apply_enhancement.selectedProperty().get();
        pillar_tracking_FD tracker_fd = new pillar_tracking_FD();
        pillar_tracking tracker = new pillar_tracking();  
        if(apply_reconstruction){
            tracker_fd.fft = fft;             
            set_parameters_tracker(tracker_fd, ip, search_window, constrain_radius, max_drift, num_threads, apply_enhancer);
        }
        else{
            set_parameters_tracker(tracker, ip, search_window, constrain_radius, max_drift, num_threads, apply_enhancer);
        }                    
        
        SaveFileChooserDialog sfcd = new SaveFileChooserDialog(defaultPath);        
        Thread worker = new Thread() {
            
            public void run(){	
                try {  
                    if(sfcd.showDialog()){
                        String output_filename = sfcd.save_fname;
                        if(!check_writable_file(output_filename)){            
                            IJ.showMessage("can not write the data to output file!");
                            return;
                        } 
                        if(apply_reconstruction) tracker_fd.save_current_settings(output_filename);
                        else tracker.save_current_settings(output_filename);
                    }
                }
                catch(Exception ex){}
            }
        };                
        worker.start();                  
    }
    
    @FXML private void handleButtonDoTracking(ActionEvent event){        
        ImagePlus ip = getSelectedImage(true);
        if(ip==null) return;       
                
        String output_filename = "";
        String text = textfiled_output_fname.getText();
        if(text.isEmpty()) {
            IJ_showMessage("the output filename is void!");
            return;
        }
        else{            
            if(check_writable_file(text)){
                output_filename = text;
            }
            else{
                IJ_showMessage("can not write the data to output file!");
                return;
            }
        }        
        
        if(!checkSettings()) return;
        
        int search_window = 0;
        text = textfiled_searchwindow.getText();
        if(text.isEmpty()) {
            IJ_showMessage("the search_window is not valid!");
            return;
        }else{
            int s = int_converter.fromString(text);
            if(s<0){
                IJ_showMessage("The search_window must be larger than 0!");
                return;
            }
            search_window = s;
        }             
                
        double max_drift = 0;
        text = textfiled_max_drift.getText();
        if(text.isEmpty()) {
            IJ_showMessage("the max_drift is not valid!");
            return;
        }else{
            double s = double_converter.fromString(text);
            if(s<=0 || s>settings.spacing){
                IJ_showMessage("The max_drift must be larger than 0 and small than spacing!");
                return;
            }
            max_drift = s;
        }             
        
        int constrain_radius = 0;
        text = textfiled_constrain_radius.getText();
        if(text.isEmpty()) {
            IJ_showMessage("the constrain_radius is not valid!");
            return;
        }else{
            int s = int_converter.fromString(text);
            if(s<=0 || s>settings.spacing){
                IJ_showMessage("The constrain_radius must be larger than 0 and small than spacing!");
                return;
            }
            constrain_radius = s;
        }             
        
        boolean apply_reconstruction = this.checkbox_apply_pillar_recontruction.selectedProperty().get();        
        if(apply_reconstruction){
            if(fft==null || !fft.hasOffCenterPoints()){
                IJ_showMessage("the mask filter is not set yet!");
                return;
            }
        }
        
        int num_threads = 1;
        text = textfiled_num_threads.getText();
        if(!text.isEmpty()) {
            int s = int_converter.fromString(text);
            if(s>0) num_threads = s;
        }             
                
        boolean apply_enhancer = this.checkbox_apply_enhancement.selectedProperty().get();
        boolean apply_rank_filter = this.meun_apply_rank_filter.isSelected();                
        
        tracker_save_fname = output_filename;
        pillar_tracking tracker = new pillar_tracking();    
        pillar_tracking_FD tracker_fd = new pillar_tracking_FD();
        if(apply_reconstruction){
            tracker_fd.apply_mean_rank_filter = apply_rank_filter;
            tracker_fd.fft = fft; 
            tracker_fd.show_FFT = checkmenu_show_fft.isSelected();
            tracker_fd.show_reconstruction = checkmenu_show_recontruction.isSelected();            
            set_parameters_tracker(tracker_fd, ip, search_window, constrain_radius, max_drift, num_threads, apply_enhancer);
        }
        else{
            tracker.apply_mean_rank_filter = apply_rank_filter;
            tracker.apply_drift_creation = this.checkmenu_drift_correction.isSelected();
            set_parameters_tracker(tracker, ip, search_window, constrain_radius, max_drift, num_threads, apply_enhancer);            
        }        
        
        boolean creating_grid = checkmenu_create_grid_poly_fit.isSelected();
        Grid_Creation_Plugin creator = new Grid_Creation_Plugin();       
        creator.setup(tracker_save_fname, settings.pixel_size, settings.spacing, settings.grid_oblique, settings.grid_angle, settings.diameter);  
        creator.setShowDialog(false);
        
        Thread worker = new Thread() {
        public void run() {		            	            
             try {               
               if(apply_reconstruction) tracker_fd.run("checkonce");               
               else{
                   tracker.run("checkonce");
                   if(creating_grid) creator.run("");
               }
             } catch (Exception ex) {
             }

             // Report the result using invokeLater().
             SwingUtilities.invokeLater(new Runnable() {
               public void run() { 
                   String text = creator.getInputFilename();
                   if(text!=null) textfile_drift_fname.setText(text); 
                   if(!apply_reconstruction) text = creator.getOutputFilename();
                   if(text!=null) textfile_grid_fname.setText(text);
               }
             });
           }
         };

        worker.start(); // So we don't hold up the dispatch thread.
        //text = textfile_drift_fname.getText();
        //if(text==null)
        //textfile_drift_fname.setText(tracker_save_fname); 
    }
    
    private void set_parameters_tracker(pillar_tracking tracker, ImagePlus ip, int search_window, int constrain_radius, double max_drift, int num_threads, boolean apply_enhancer){
        tracker.setShowDialog(false);
        tracker.lattice = settings.pixel_size;
        tracker.diameter = settings.diameter;
        tracker.setup(ip, settings.gauss_sigma, search_window, constrain_radius, settings.dark_pillar, settings.spacing,  settings.grid_oblique, max_drift, num_threads, tracker_save_fname);  
        tracker.setGridAngle(settings.grid_angle);
        tracker.setAlgorithm(localization_algorithm);
        if(apply_enhancer) tracker.setEnhancementSolver(image_enhancer_plugin);        
        else tracker.setEnhancementSolver(null);
    }
    
    private boolean check_writable_file(String fname){
        boolean can_write = false;
        try{
            Path file = new File(fname).toPath();
            if(Files.exists(file)){
                boolean iw = Files.isWritable(file);
                if(iw) can_write = true;
            }
            else if(Files.notExists(file)){
                Files.createFile(file);
                can_write = true;
            }            
        }
        catch(Exception e){
            IJ.log("file is not writable->" + fname);
        }
        
        return can_write;
    }
    
    // </editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="drift analysis events handlding">
    @FXML private void handleButtonDeflectionDrift(ActionEvent event){
        updateDriftPlotWindow();
    }
    
    @FXML private void handleButtonBrowseDriftFile(ActionEvent event){
        
        SwingUtilities.invokeLater(new Runnable() {
            @Override public void run() {
                File file = new File(textfile_drift_fname.getText());        
                //drift_fileChooser.setCurrentDirectory(file);
                drift_fileChooser.setSelectedFile(file);
                int retval= drift_fileChooser.showOpenDialog(owner_window);        
                if (retval == JFileChooser.APPROVE_OPTION){
                    
                    Platform.runLater(new Runnable() {
                        @Override
                        public void run() {
                            textfile_drift_fname.setText(drift_fileChooser.getSelectedFile().getAbsolutePath());            
                            handleButtonLoadFileDrift(event);
                        }
                    });
                    
                }       
            }
        });
             
    }
    
    @FXML private void handleButtonBrowseGridFile(ActionEvent event){
        
        SwingUtilities.invokeLater(new Runnable() {
            @Override public void run() {
                File file = new File(textfile_grid_fname.getText());        
                //grid_fileChooser.setCurrentDirectory(file);
                grid_fileChooser.setSelectedFile(file);
                int retval= grid_fileChooser.showOpenDialog(owner_window);
                if (retval == JFileChooser.APPROVE_OPTION){
                   
                    Platform.runLater(new Runnable() {
                        @Override
                        public void run() {
                            textfile_grid_fname.setText(grid_fileChooser.getSelectedFile().getAbsolutePath());
                            handleButtonLoadFileGrid(event);
                        }
                    });
                }      
            }
        });
        
    }
    
    private boolean showPairSettings(){
        GenericDialog gd = new GenericDialog("Settings for Drawing a Pair of Pillars");
        int n = reader.npillars;    
        int frames = reader.nframes;
        
        gd.addNumericField("        pillar index 1:",  pair_pillars.p1, 0, 10, "(1-" + n + ")");
        gd.addNumericField("        pillar index 2:",  pair_pillars.p2, 0, 10, "(1-" + n + ")");		
        gd.addNumericField("        start frame:",  pair_pillars.start_frame, 0, 10, "(1-" + frames + ")");
        gd.addNumericField("        end frame:",  pair_pillars.end_frame, 0, 10, "<"+(frames+1));
        
        gd.showDialog();
        if (gd.wasCanceled()) return false;
        
        int p1 = (int)gd.getNextNumber();        
        int p2 = (int)gd.getNextNumber();
        int fs = (int)gd.getNextNumber();        
        int fe = (int)gd.getNextNumber();        
          
        if(p1<=0 || p1>n){
            this.IJ_showMessage("pillar index 1 is not in the valid range (1-" + n + ")" );
            return false;
        }
        
        if(p2<=0 || p2>n){
            this.IJ_showMessage("pillar index 2 is not in the valid range (1-" + n + ")" );
            return false;
        }
        
        if(p1==p2){
            this.IJ_showMessage("pillar index 1 and 2 is identical");
            return false;
        }
        
        
        if(fs<=0 || fs>frames){
            this.IJ_showMessage("start frame is not in the valid range (1-" + frames + ")" );
            return false;
        }
        
        if(fe<=0) fe = frames;
        
        if(fs>=fe){
            this.IJ_showMessage("start frame should be small than end frame");
            return false;
        }
        
        double[] dis = reader.getPairDistance(p1-1, p2-1);
        if(dis==null) return false;
        
        double[] plot_x = reader.getFrameSeries();                        
        Plot ploter_dis = new Plot("distance between pillar "+ p1 + " and " + p2, "frame #", "distance(nm)",plot_x,dis);
        ploter_dis.draw();
        ploter_dis.setColor(Color.CYAN);
        ploter_dis.setLineWidth(2);
        double[] xylimits = ploter_dis.getLimits();        
        ploter_dis.drawLine(fs, xylimits[2], fs, xylimits[3]);
        ploter_dis.drawLine(fe, xylimits[2], fe, xylimits[3]);
        ploter_dis.show(); 
        
        //plot the deflection information if existed                    
        double[] dis1 = new double[reader.nframes];
        double[] dis2 = new double[reader.nframes]; 
        if(reader.DX!=null){
            for(int f=0; f<reader.nframes; f++){                        
                double dx = reader.DX[p1-1][f];
                double dy = reader.DY[p1-1][f];
                dis1[f] = Math.sqrt(dx*dx+dy*dy);             

                dx = reader.DX[p2-1][f];
                dy = reader.DY[p2-1][f];
                dis2[f] = Math.sqrt(dx*dx+dy*dy);                                 
            }
        }
        else{
            if(reader.correted_data){
                dis1 = reader.getDeflections(p1-1)[2];
                dis2 = reader.getDeflections(p2-1)[2];            
            }
            else{
                dis1 = reader.getDeflectionsRAW(p1-1)[2];
                dis2 = reader.getDeflectionsRAW(p2-1)[2];            
            }
        }
        String title = "Deflections of a Pair of Pillars";
        Plot ploter_def = DriftAnalysisPlotter.plot2(plot_x, dis1, dis2, title ,"frame #", "displacement(nm)");
        ploter_def.setLegend("pillar ID:"+ p1 + "\npillar ID:" + p2,Plot.AUTO_POSITION|Plot.LEGEND_TRANSPARENT);                    
        xylimits = ploter_def.getLimits();
        ploter_def.setColor(Color.CYAN);
        ploter_def.setLineWidth(2);
        ploter_def.drawLine(fs, xylimits[2], fs, xylimits[3]);
        ploter_def.drawLine(fe, xylimits[2], fe, xylimits[3]);        
        ploter_def.show();	
        
        int[] indexs = new int[]{p1-1,p2-1};
        selected_pillars_indexs = indexs;
        draw_list_pillars(indexs, fs-1, fe-1);
        
        pair_pillars.p1 = p1;
        pair_pillars.p2 = p2;
        pair_pillars.start_frame = fs;
        pair_pillars.end_frame = fe;
        return true;
    }
    
    @FXML private void handleButtonPlotPairSepartion(ActionEvent event){
        if(!load_suc_drift) return;
        SwingUtilities.invokeLater(new Runnable() {
            @Override public void run() {
                showPairSettings();
            }
        });
    }
    
    @FXML private void handleButtonLoadFileDrift(ActionEvent event){        
        //accordion_drift.setDisable(true);
        //accordion_drift.setExpanded(false);
        //accordion_drift.setExpandedPane(null);
        String fpath = textfile_drift_fname.getText();    
        List<ContractionUnit> cu = Hunt_CU_PlugIn.load(fpath);
        if(cu!=null) {
            if(load_suc_drift){
                ContractionUnit.show_table_list(cu);
                selected_pillars_indexs = ContractionUnit.get_pillar_index_list(cu, reader.npillars);
                draw_list_pillars(selected_pillars_indexs);
            }
            return;
        }
        
        double current_pixel_size = settings.pixel_size;
        boolean unknown_pixle_size = Double.isNaN(current_pixel_size);
        pixel_size = unknown_pixle_size ? 1.0 : current_pixel_size;
        load_suc_drift = reader.load_tracks(fpath, pixel_size);                
        accordion_drift.setDisable(!load_suc_drift);  
        accordion_drift.setExpanded(load_suc_drift);
        if(load_suc_drift){
            //ObservableList<TitledPane> panes = accordion_drift.getPanes();
            //accordion_drift.setExpandedPane(panes.get(0));
            //accordion_drift.setExpanded(true);
            checkbox_raw.setDisable(!reader.correted_data);
            button_drift.setDisable(!reader.correted_data);
            meun_export_grid_format.setDisable(reader.getFileVersion()==pillar_tracking_FD.fileversion2);
            if(!reader.correted_data) checkbox_raw.setSelected(true);
            checkbox_still_pillars.setDisable(!reader.correted_data);
            if(!reader.correted_data) checkbox_still_pillars.setSelected(false);
            
            pixel_size = reader.pixel_size;
            String label = "Frames:" + reader.nframes + "  Pillars:" + reader.npillars + "  Pixel Size:" + reader.pixel_size + "nm";            
            if(reader.file_header==null || reader.file_header.lattice!=pixel_size) label +=  "(Unknown)";           
            
            label_drift_load_info.setText(label);
            
            selected_pillars_indexs = null;
            cu_hunter = new Hunt_CU_PlugIn();
            pair_pillars = new PairPillars();   
            
            overlay_tracks_ip = null;            
            slider_pillar_no_drift.maxProperty().set(reader.npillars);
            slider_frame_no_drift.maxProperty().set(reader.nframes);
            slider_frame_no_drift.setValue(reader.nframes);
            //int frame_window = int_converter.fromString(textfield_frame_window.getText());
            //if(frame_window>0 && frame_window<=reader.nframes){
            //    slider_frame_no_drift.setValue(frame_window);
            //}
            //File drift_file = new File(fpath);		            
            drift_fname = fpath; //drift_file.getName();
            need_recreate_map_drift = true;
        }
    }
    
    @FXML
    private void handleButtonDriftAction(ActionEvent event) {
        if(load_suc_drift && reader.correted_data){  
            //DriftAnalysisPlotter plotter = new DriftAnalysisPlotter();
            double[] xaxis = reader.getFrameSeries();            
            Plot ploter_driftXY = DriftAnalysisPlotter.get_plot_driftXY(xaxis, reader.driftX, reader.driftY, "drift x&y");
            ploter_driftXY.setLegend("x-drift\ny-drift\ndrift",Plot.AUTO_POSITION|Plot.LEGEND_TRANSPARENT);
            ploter_driftXY.show();	
        }
    }
    
    @FXML
    private void handleButtonStd(ActionEvent event){
        if(!load_suc_drift) return;
        
        int npillars = this.reader.npillars;
        int nframes = this.reader.nframes;
        if(npillars<1 || nframes<2) return;
        
        int[] indexs = this.selected_pillars_indexs;
        
        if(indexs==null || indexs.length<=0) {
            indexs = new int[npillars];
            for(int i=0; i<npillars; i++) indexs[i] = i;
        }
        this.showStatisticTable(indexs, "Results");
    }
    
    @FXML private void handleButtonHist(ActionEvent event){
        if(load_suc_drift && reader.correted_data){
            double[] plot_std = reader.getStdXY();            
            plotter.show_histogram(plot_std, 0.5, "Histogram of Std X&Y");	
            //plotter.show_jump_length_histogram(reader.correct_tracksX, reader.correct_tracksY, 0.5, "Histogram of Jump Length");            
        }
    }
    
    private void hunt_CU(){
        if(!load_suc_drift) return;
        if(!reader.correted_data) return;
        
        double[][] cx = reader.correct_tracksX;
        double[][] cy = reader.correct_tracksY;
        int npillars = reader.npillars;
        int frames = reader.nframes;
        double[] gx = new double[npillars];
        double[] gy = new double[npillars];
        for(int i=0; i<npillars; i++){
            //gx[i] = cx[i][0];
            //gy[i] = cx[i][0];
            int num = 0;
            gx[i] = 0;
            gy[i] = 0;
            for(int f=0; f<frames; f++){
                double x = cx[i][f]/reader.pixel_size;
                double y = cy[i][f]/reader.pixel_size;
                if(!Double.isNaN(x)){
                    gx[i] += x;
                    gy[i] += y;
                    num++;
                }
            }
            if(num>0){
                gx[i] /= num;
                gy[i] /= num;                    
            }
            else{
                gx[i] = Double.NaN;
                gy[i] = Double.NaN;                    
            }
        }
        double catch_radius = cu_hunter.catch_radius;
        if(Double.isNaN(catch_radius) || catch_radius<=0){
            double s = GridCreation.estimate_seperation(gx, gy); 
            if(Double.isNaN(s)) s = settings.spacing;
            catch_radius = s*1.25;
        }
        Hunt_CU_PlugIn hunter = new Hunt_CU_PlugIn(catch_radius, cu_hunter.threshold_rs, cu_hunter.threshold_rp, cu_hunter.frame_window);
        hunter.threshold_deflection = cu_hunter.threshold_deflection;
        hunter.merging_overlaps = cu_hunter.merging_overlaps;
        hunter.output_fname = reader.getFileName() + ".hcu";
        if(hunter.showDialog(frames)){
            if(cu_hunter!=null) cu_hunter.copy_parameters(hunter);   
            else cu_hunter = hunter;
            
            Hunt_CU_PlugIn.using_deflections = (reader.DX!=null && reader.getFileVersion()>=pillar_tracking_FD.fileversion2);
            if(Hunt_CU_PlugIn.using_deflections){
                cx = reader.DX;
                cy = reader.DY;
            }
            //List<ContractionUnit> cu_list = Hunt_CU_PlugIn.hunt(gx, gy, hunter.catch_radius, cx, cy, reader.pixel_size, hunter.threshold_rs, hunter.threshold_rp,hunter.threshold_deflection,hunter.frame_window, hunter.merging_overlaps);
            if(catch_radius!= cu_hunter.catch_radius) cu_hunter.setToReProjecting();            
            List<ContractionUnit> cu_list = cu_hunter.hunt(gx, gy, cx, cy, reader.pixel_size);
            
            if(check_writable_file(hunter.output_fname)) hunter.save(hunter.output_fname, cu_list);
            
            ContractionUnit.print_list(cu_list);
            if(hunter.show_table) ContractionUnit.show_table_list(cu_list);
            
            int[] list = ContractionUnit.get_pillar_index_list(cu_list, npillars);
            selected_pillars_indexs = list;
            
//            int frame_window = int_converter.fromString(textfield_frame_window.getText());
//            if(frame_window>0){
//                int start_f = (int)this.slider_frame_no_drift.getValue()-1;
//                int end_f = start_f + frame_window;
//                draw_list_pillars(selected_pillars_indexs,start_f,end_f);
//                return;
//            }
            
            draw_list_pillars(list);
        }
    }
    
    @FXML private void handleButtonHuntCU(ActionEvent event){        
        Thread worker = new Thread() {
                public void run() {		            	            
                         try {
                             hunt_CU();
                         } catch (Exception ex) { }
                }
        };

        worker.start(); // So we don't hold up the dispatch thread.
    }
    
    private void draw_highlight_cu(int frame){
        if(cu_list_grid!=null)            
            contraction_unit_list = ContractionUnit.get_pillar_index_list(cu_list_grid,frame,grid_reader.cens);
        if(contraction_unit_list==null){
            showLabelsArrows(candidate_grid, frame);
            return;
        }        
        if(candidate_grid!=null && candidate_grid.length>0){
            int num_candidate = candidate_grid.length;
            boolean[] contraction_unit_flag = new boolean[num_candidate];
            for(int i=0; i<num_candidate; i++){
                contraction_unit_flag[i] = false;
                for(int j=0; j<contraction_unit_list.length; j++){
                    if(contraction_unit_list[j]==candidate_grid[i]){
                        contraction_unit_flag[i] = true;
                        break;
                    }
                }
            }
            showLabelsArrows(candidate_grid, contraction_unit_flag, frame);                
        }
        else showLabelsArrows(contraction_unit_list, frame);    
    }
    
    private void hunt_CU_grid(boolean individual_frame){
        if(!load_suc_grid) return;
        if(grid_reader.CX==null) return;
        
        double[][] cx = grid_reader.CX;
        double[][] cy = grid_reader.CY;
        int npillars = grid_reader.cens;
        int frames = grid_reader.frames;
        double[][] gxy = grid_reader.getGridXY(0);
        if(gxy==null) return;
        
        double[] gx = gxy[0];
        double[] gy = gxy[1];
        double catch_radius = cu_hunter_grid.catch_radius;
        if(Double.isNaN(catch_radius) || catch_radius<=0){
            double s = grid_reader.spacing; 
            if(Double.isNaN(s)) s = GridCreation.estimate_seperation(gx, gy);         
            catch_radius = s*Math.sqrt(2);
        }
        
        if(!individual_frame){
            Hunt_CU_PlugIn hunter = new Hunt_CU_PlugIn(catch_radius, cu_hunter_grid.threshold_rs, cu_hunter_grid.threshold_rp, cu_hunter_grid.frame_window);
            hunter.threshold_deflection = cu_hunter_grid.threshold_deflection;
            hunter.merging_overlaps = cu_hunter_grid.merging_overlaps;
            hunter.output_fname = grid_reader.getFileName() + ".hcu";
            if(hunter.showDialog(frames)){
                if(cu_hunter_grid!=null) cu_hunter_grid.copy_parameters(hunter);                 
                else cu_hunter_grid = hunter;
                
                Hunt_CU_PlugIn.using_deflections = (grid_reader.DX!=null && grid_reader.getFileVersion()>=pillar_tracking_FD.fileversion2);
                if(Hunt_CU_PlugIn.using_deflections){
                    cx = grid_reader.DX;
                    cy = grid_reader.DY;
                }
                //cu_list_grid = Hunt_CU_PlugIn.hunt(gx, gy, hunter.catch_radius, cx, cy, grid_reader.lattice, hunter.threshold_rs, hunter.threshold_rp, hunter.threshold_deflection,hunter.frame_window, hunter.merging_overlaps);                                
                if(catch_radius != cu_hunter_grid.catch_radius) cu_hunter_grid.setToReProjecting();
                cu_list_grid = cu_hunter_grid.hunt(gx, gy, cx, cy, grid_reader.lattice);
                
                if(check_writable_file(hunter.output_fname)) hunter.save(hunter.output_fname, cu_list_grid);
                
                ContractionUnit.print_list(cu_list_grid);                
                if(hunter.show_table) ContractionUnit.show_table_list(cu_list_grid);
                
                //contraction_unit_list = ContractionUnit.get_pillar_index_list(cu_list_grid,f,npillars);                
                int f = (int)this.slider_frame_no_grid.getValue()-1;
                draw_highlight_cu(f);
            }
        }
        else{
            Hunt_CU_Frame_PlugIn hunter = new Hunt_CU_Frame_PlugIn(catch_radius);
            if(cu_hunter_grid instanceof  Hunt_CU_Frame_PlugIn){
                Hunt_CU_Frame_PlugIn cu_hunter = (Hunt_CU_Frame_PlugIn)cu_hunter_grid;
                hunter = new Hunt_CU_Frame_PlugIn(catch_radius, cu_hunter.threshold_deflection, cu_hunter.threshold_angle_diff);
            }            
            if(hunter.showDialog(frames)){
                cu_hunter_grid = hunter;      
                
                Hunt_CU_PlugIn.using_deflections = (grid_reader.DX!=null && grid_reader.getFileVersion()>=pillar_tracking_FD.fileversion2);
                if(Hunt_CU_PlugIn.using_deflections){
                    cx = grid_reader.DX;
                    cy = grid_reader.DY;
                }
                
                int f = (int)this.slider_frame_no_grid.getValue()-1;
                
                if(hunter.process_all_frames)
                    cu_list_grid = Hunt_CU_Frame_PlugIn.hunt(gx, gy, hunter.catch_radius, cx, cy, grid_reader.lattice, hunter.threshold_deflection, hunter.getCosThreshold());                
                else
                    cu_list_grid = Hunt_CU_Frame_PlugIn.hunt(gx, gy, hunter.catch_radius, f, cx, cy, grid_reader.lattice, hunter.threshold_deflection, hunter.getCosThreshold());
                
                ContractionUnit.print_list(cu_list_grid);
                if(hunter.show_table) ContractionUnit.show_table_list(cu_list_grid);
                draw_highlight_cu(f);
            }
        }
    }
    
    @FXML private void handleButtonHuntCUGrid(ActionEvent event){        
        Thread worker = new Thread() {
                public void run() {		            	            
                         try {
                             hunt_CU_grid(false);
                         } catch (Exception ex) { }
                }
        };

        worker.start(); // So we don't hold up the dispatch thread.
    }
    
    @FXML private void handleMenuHuntCUFrameGrid(ActionEvent event){        
        Thread worker = new Thread() {
                public void run() {		            	            
                         try {
                             hunt_CU_grid(true);
                         } catch (Exception ex) { }
                }
        };

        worker.start(); // So we don't hold up the dispatch thread.
    }
    
    
    @FXML private void handleButtonDrawListPillars(ActionEvent event){
        if(!load_suc_drift) return;
        Alert alert = new Alert(AlertType.CONFIRMATION);
        alert.setTitle("Pillar Tacker Version:" + Pillar_Tracker_Plugin.VERSION);
        alert.setHeaderText("Draw a List of Pillars");
        alert.setContentText("Input A List of Pillar Index,support separators(comma,colon,dash,tab,space and line)");

        TextArea textArea = new TextArea();
        textArea.setEditable(true);
        textArea.setWrapText(true);
        if(list_pillar_string!=null && !list_pillar_string.isEmpty())
            textArea.setText(list_pillar_string);

        textArea.setMaxWidth(Double.MAX_VALUE);
        textArea.setMaxHeight(Double.MAX_VALUE);
        GridPane.setVgrow(textArea, Priority.ALWAYS);
        GridPane.setHgrow(textArea, Priority.ALWAYS);

        GridPane expContent = new GridPane();
        expContent.setMaxWidth(Double.MAX_VALUE);        
        expContent.add(textArea, 0, 1);        

        // Set expandable Exception into the dialog pane.
        alert.getDialogPane().setExpandableContent(expContent);
        alert.getDialogPane().setExpanded(true);
        //alert.showAndWait();
        Optional<ButtonType> result = alert.showAndWait();
        if (result.get() == ButtonType.OK){
            // ... user chose OK
            list_pillar_string = textArea.getText();
            int[] list = get_int_list(list_pillar_string, reader.npillars+1);
            selected_pillars_indexs = list;
            if(overlay_tracks_ip==null) draw_list_pillars(list);
            else{
                Roi roi = overlay_tracks_ip.getRoi();            
                if(roi==null || !roi.isArea()) draw_list_pillars(list);
                else{
                    int[] new_list = reader.getSelectedPillars(roi, list);                    
                    draw_list_pillars(new_list);
                    selected_pillars_indexs = new_list;
                } 
            }            
        } else {
            // ... user chose CANCEL or closed the dialog
        }
    }
    
    @FXML private void handleButtonHighlightListPillars(ActionEvent event){
        if(!load_suc_grid) return;
        Alert alert = new Alert(AlertType.CONFIRMATION);
        alert.setTitle("Pillar Tacker Version:" + Pillar_Tracker_Plugin.VERSION);
        alert.setHeaderText("Draw a List of Pillars");
        alert.setContentText("Input A List of Pillar Index,support separators(comma,colon,dash,tab,space and line)");

        TextArea textArea = new TextArea();
        textArea.setEditable(true);
        textArea.setWrapText(true);
        if(grid_list_pillar_string!=null && !grid_list_pillar_string.isEmpty())
            textArea.setText(grid_list_pillar_string);

        textArea.setMaxWidth(Double.MAX_VALUE);
        textArea.setMaxHeight(Double.MAX_VALUE);
        GridPane.setVgrow(textArea, Priority.ALWAYS);
        GridPane.setHgrow(textArea, Priority.ALWAYS);

        GridPane expContent = new GridPane();
        expContent.setMaxWidth(Double.MAX_VALUE);        
        expContent.add(textArea, 0, 1);        

        // Set expandable Exception into the dialog pane.
        alert.getDialogPane().setExpandableContent(expContent);
        alert.getDialogPane().setExpanded(true);
        //alert.showAndWait();
        Optional<ButtonType> result = alert.showAndWait();
        if (result.get() == ButtonType.OK){
            // ... user chose OK
            int f = (int)this.slider_frame_no_grid.getValue()-1;            
            grid_list_pillar_string = textArea.getText();
            int[] list = get_int_list(grid_list_pillar_string, grid_reader.cens+1);            
            contraction_unit_list = list;
            this.cu_list_grid = null;
            draw_highlight_cu(f);
            
        } else {
            // ... user chose CANCEL or closed the dialog
        }
    }
    
    private boolean parse_string(String line, String sep, boolean[] flag){
        int length = flag.length;
        boolean conained = line.contains(sep);
        if(conained){
            try{
                for (String s : line.split(sep)){
                    int i = Integer.parseInt(s);
                    if(i>=0 && i<length) flag[i] = true;                       
                }
            }
            catch(Exception ex){
                
            }
        }
        return conained;
    }
    
    private void parse_string(String line, boolean[] flag){
        int length = flag.length;
        if(!parse_string(line, ",", flag)){
            if(!parse_string(line, " ", flag)){
                if(!parse_string(line, "-", flag)){
                    if(!parse_string(line, "\t", flag)){
                        if(line!=null && !line.isEmpty()) {
                            try{
                                int k = Integer.parseInt(line);
                                if(k>=0 && k<length) flag[k] = true;
                            }
                            catch(Exception ex){
                                
                            }
                        }
                    }
                }
            }
        }               
    }
    
    private int[] get_int_list(String text, int length){
        String line_sep = System.getProperty("line.separator");
        boolean[] flag = new boolean[length];
        for(int i=0; i<length; i++) flag[i] = false;
        
        if(text.contains(line_sep)){            
            for (String line : text.split(line_sep)){
                if(line==null || line.isEmpty()) continue;
                parse_string(line, flag);      
            }
        }
        else{
            if(text.contains("\n")){ 
                for (String line : text.split("\n")){
                    if(line==null || line.isEmpty()) continue;
                    parse_string(line, flag);        
                 }
            }
            else{
                parse_string(text, flag);         
            }
        }
        int num = 0;
        for(int i=0; i<length; i++){
            if(flag[i]) num++;
        }
        if(num<1) return null;
        
        int[] list = new int[num];
        num = 0;
        for(int i=0; i<length; i++){
            if(flag[i] && i>0){
                list[num] = i-1;
                num++;
            }
        }      
        
        return list;
    }
    
    private void draw_list_pillars(int[] list)
    {
        draw_list_pillars(list,0,reader.nframes-1);
    }
    
    private void draw_list_pillars(int[] list, int start_frame, int end_frame)
    {
        if(!load_suc_drift) return;
        if(list==null || list.length<1) return;
        
            //boolean flatten = checkbox_flatten.isSelected();
            boolean show=checkbox_show_labels.isSelected();
            boolean raw = checkbox_raw.isSelected();   
            
            File input_file = new File(drift_fname);            
            String title = input_file.getName();
            double[][] cx = raw ? reader.raw_tracksX : reader.correct_tracksX;
            double[][] cy = raw ? reader.raw_tracksY : reader.correct_tracksY;            
            
//            if(flatten){                
//                ImagePlus ip = plotter.create_empty_map(cx, cy, pixel_size, 5, title);
//                
//                cx = raw ? reader.getSubset(list,reader.raw_tracksX) : reader.getSubset(list,reader.correct_tracksX);
//                cy = raw ? reader.getSubset(list,reader.raw_tracksY) : reader.getSubset(list,reader.correct_tracksY);
//                Overlay labels = plotter.getLabels(cx, cy, pixel_size, zoom_in_image);
//                flatten_tracks_ip = plotter.plot_flatten_tracks(cx, cy, pixel_size, ip.getWidth(), ip.getHeight(),zoom_in_tracks*zoom_in_image, zoom_in_image, title);
//                if(show) flatten_tracks_ip.setOverlay(labels);
//            }
//            else
            { 
                boolean shownin = checkbox_shownin.isSelected();
                if(shownin && WindowManager.getImageCount()>0){
                    ImagePlus ip = WindowManager.getCurrentImage();//IJ.getImage();
                    if(ip!=null){
                            overlay_tracks_ip = ip;
                            need_recreate_map_drift = false;
                    }
                }
                else{
                    overlay_tracks_ip = created_map_ip;
                }
                
                if(overlay_tracks_ip==null || !overlay_tracks_ip.isVisible() || need_recreate_map_drift){
                    overlay_tracks_ip = created_map_ip = plotter.create_empty_map(cx, cy, pixel_size, 5, title + " x" + zoom_in_tracks);
                    need_recreate_map_drift = false;
                }
                
                if(show) plotter.getLabels(cx, cy, list, pixel_size);      
                
                cx = raw ? reader.getSubset(list,reader.raw_tracksX) : reader.getSubset(list,reader.correct_tracksX);
                cy = raw ? reader.getSubset(list,reader.raw_tracksY) : reader.getSubset(list,reader.correct_tracksY);                
                plotter.getTracks(cx, cy, pixel_size, start_frame, end_frame, zoom_in_tracks);
                plotter.plot_tracks_overlay(overlay_tracks_ip, show);    
            }
    }
    
    @FXML private void handleButtonTracksMap(ActionEvent event){
        if(load_suc_drift){
            boolean flatten = checkbox_flatten.isSelected();
            boolean show = checkbox_show_labels.isSelected();
            boolean raw = checkbox_raw.isSelected();   
            boolean still_pillars = checkbox_still_pillars.isSelected();   
            
            File input_file = new File(drift_fname);            
            String title = input_file.getName();
            double[][] cx = raw ? reader.raw_tracksX : reader.correct_tracksX;
            double[][] cy = raw ? reader.raw_tracksY : reader.correct_tracksY;
                
            if(flatten){
                ImagePlus ip = plotter.create_empty_map(cx, cy, pixel_size, 5, title);
                Overlay labels = plotter.getLabels(cx, cy, pixel_size, zoom_in_image);
                flatten_tracks_ip = plotter.plot_flatten_tracks(cx, cy, pixel_size, ip.getWidth(), ip.getHeight(),zoom_in_tracks*zoom_in_image, zoom_in_image, title);
                if(show) flatten_tracks_ip.setOverlay(labels);
            }
            else{ 
                boolean shownin = checkbox_shownin.isSelected();
                if(shownin && WindowManager.getImageCount()>0){
                    ImagePlus ip = WindowManager.getCurrentImage();//IJ.getImage();
                    if(ip!=null){
                            overlay_tracks_ip = ip;
                            need_recreate_map_drift = false;
                    }
                }
                else{
                    overlay_tracks_ip = created_map_ip;
                }
                
                if(overlay_tracks_ip==null || !overlay_tracks_ip.isVisible() || need_recreate_map_drift){
                    overlay_tracks_ip = created_map_ip = plotter.create_empty_map(cx, cy, pixel_size, 5, title + " x" + zoom_in_tracks);
                    need_recreate_map_drift = false;
                }
                //updateSliderZoominTracks();
                Roi roi = overlay_tracks_ip.getRoi();     
                int[] indexs = reader.getSelectedPillars(roi);                
                if(roi==null || !roi.isArea()){
                    int num_pillars = reader.npillars;
                    indexs = new int[num_pillars];
                    for(int i=0; i<num_pillars; i++) indexs[i] = i;                     
                } 
                if(!still_pillars) indexs = exculde_still_pillars(indexs);  
                selected_pillars_indexs = indexs;                
                
                int frame_window = int_converter.fromString(textfield_frame_window.getText());
                if(frame_window>0 && frame_window<reader.nframes){
                    int end_f = (int)this.slider_frame_no_drift.getValue()-1;
                    int start_f = end_f - frame_window;
                    draw_list_pillars(selected_pillars_indexs,start_f,end_f);
                    return;
                }
                
                if(indexs!=null && indexs.length>=0) { 
                    if(show) plotter.getLabels(cx, cy, indexs, pixel_size);  
                    cx = raw ? reader.getSubset(indexs,reader.raw_tracksX) : reader.getSubset(indexs,reader.correct_tracksX);
                    cy = raw ? reader.getSubset(indexs,reader.raw_tracksY) : reader.getSubset(indexs,reader.correct_tracksY);
                    
                    plotter.getTracks(cx, cy, pixel_size, zoom_in_tracks);
                    plotter.plot_tracks_overlay(overlay_tracks_ip, show);   
                }
                else{
                    overlay_tracks_ip.setOverlay(null);
                }
//                if(indexs==null || indexs.length<=0) {                    
//                    plotter.getTracks(cx, cy, pixel_size, zoom_in_tracks);
//                    plotter.plot_tracks_overlay(overlay_tracks_ip, show);
//                }
//                else{                    
//                    cx = raw ? reader.getSubset(indexs,reader.raw_tracksX) : reader.getSubset(indexs,reader.correct_tracksX);
//                    cy = raw ? reader.getSubset(indexs,reader.raw_tracksY) : reader.getSubset(indexs,reader.correct_tracksY);
//                    
//                    plotter.getTracks(cx, cy, pixel_size, zoom_in_tracks);
//                    plotter.plot_tracks_overlay(overlay_tracks_ip, show);   
//                }
                
            }
        }
    }
    
    private int[] exculde_still_pillars(int[] indexs){
        if(indexs==null || indexs.length<=0) return null;
        if(!load_suc_drift) return null;
        if(!reader.correted_data) return indexs;
        boolean[] flags = reader.flags_still_pillars;
        int numflags = flags.length;
        int num = 0;
        for(int i=0; i<indexs.length; i++){
            int k = indexs[i];
            if(k<numflags){
                if(!flags[k]) num++;
            }
        }
        if(num==0) return null;        
        int[] new_indexs = new int[num];
        num = 0;
        for(int i=0; i<indexs.length; i++){
            int k = indexs[i];
            if(k<numflags){
                if(!flags[k]){
                    new_indexs[num] = k;
                    num++;
                }
            }
        }
        return new_indexs;
    }
    
    @FXML private void handleExportButtonAction(ActionEvent event) {        
        //if(!checkSettings()) return;  
        checkSettings(false);        
        SwingUtilities.invokeLater(new Runnable() {
            @Override public void run() {
                Grid_Creation_Plugin creator = new Grid_Creation_Plugin();    
                if(reader.file_header!=null){
                    FileHeaderDrift header = reader.file_header;
                    creator.setup(drift_fname, header.lattice, header.spacing, header.oblique, header.grid_angle, header.diameter);   
                }
                else{ 
                    creator.setup(drift_fname, settings.pixel_size, settings.spacing, settings.grid_oblique, settings.grid_angle, settings.diameter);   
                }
                Thread worker = new Thread() {
                    public void run() {		            	            
                         try {
                           // Something that takes a long time . . . in real life,
                           creator.run("");
                           if(!creator.IsDialogCanceled()){
                                String text = creator.getOutputFilename();
                                if(text!=null) textfile_grid_fname.setText(text);
                            }
                         } catch (Exception ex) {
                         }
                    }
                };

                worker.start();
            }
        });        
    }
    
    @FXML private void handleButtonAction(ActionEvent event) {
        SwingUtilities.invokeLater(new Runnable() {
            @Override public void run() {
                double[] rst = GaussianFitter.do_fit();
                if(rst!=null && rst.length>3){
                    double a = rst[0];
                    double b = rst[1];
                    checkbox_dark_pillar.setSelected(a>b);            
                    double sigma = rst[3];
                    textfiled_gauss_sigma.textProperty().setValue(String.format(Locale.ENGLISH, "%.3f", sigma));
                    //textfiled_gauss_sigma.setText(String.format("%.3f", sigma));            
                }
            }
        });        
    }    
    
    private void showStatisticTable(int[] indexs, String table_name){        
        if(!load_suc_drift) return;        
        
        ResultsTable rt = new ResultsTable();        
        rt.setNaNEmptyCells(true);
        rt.showRowNumbers(true);       

        int n = indexs.length;
        int num = this.reader.npillars;
        boolean having_dxy = this.reader.DX!=null;
        for(int i=0; i<n; i++){
            IJ.showProgress(i, n);
            rt.incrementCounter();
            int c = indexs[i];
            if(c>=0 && c<num){
                rt.addValue("Pillar Index", c+1);
                StatisitcCenterDistance rst_raw = reader.getDistanceCenterRAW(c);
                double var_raw_x = BasicStatisitic.var(reader.raw_tracksX[c]);
                double var_raw_y = BasicStatisitic.var(reader.raw_tracksY[c]);
                double std_raw = Math.sqrt(var_raw_x + var_raw_y);
               
                rt.addValue("AVG_X_RAW", rst_raw.avg_x);
                rt.addValue("AVG_Y_RAW", rst_raw.avg_y);
                rt.addValue("STD_RAW", std_raw);
                rt.addValue("AVG_DIS_RAW", rst_raw.avg_dis);
                rt.addValue("MAX_DIS_RAW", rst_raw.max_dis);
                rt.addValue("MIN_DIS_RAW", rst_raw.min_dis);
                rt.addValue("STD_DIS_RAW", rst_raw.std_dis);
                
                if(reader.correted_data){
                    double var_x = BasicStatisitic.var(reader.correct_tracksX[c]);
                    double var_y = BasicStatisitic.var(reader.correct_tracksY[c]);
                    double std = Math.sqrt(var_x + var_y);
                    
                    StatisitcCenterDistance rst = reader.getDistanceCenter(c);
                    rt.addValue("AVG_X", rst.avg_x);
                    rt.addValue("AVG_Y", rst.avg_y);
                    rt.addValue("STD", std);
                    rt.addValue("AVG_DIS", rst.avg_dis);
                    rt.addValue("MAX_DIS", rst.max_dis);
                    rt.addValue("MIN_DIS", rst.min_dis);
                    rt.addValue("STD_DIS", rst.std_dis);
                }
                
                if(having_dxy){
                    double[] def = reader.getAbsDeflections(c)[2];                    
                    double avg_def = BasicStatisitic.avg(def);
                    double max_def = BasicStatisitic.max(def);
                    double min_def = BasicStatisitic.min(def);
                    double var_def = BasicStatisitic.var(def);
                    double std_def = Math.sqrt(var_def);
                    rt.addValue("AVG_Deflection", avg_def);
                    rt.addValue("MAX_Deflection", max_def);
                    rt.addValue("MIN_Deflection", min_def);
                    rt.addValue("STD_Deflection", std_def);
                }
            }                            
        }
        IJ.showProgress(1);
        rt.show(table_name);        
    }
    
    
    @FXML private void handleButtonGetPillarsList(ActionEvent event){
        if(!load_suc_drift) return;        
        int[] indexs = this.selected_pillars_indexs;
        if(indexs==null || indexs.length<=0) return;
        
        ResultsTable rt = new ResultsTable();        
        rt.setNaNEmptyCells(true);
        rt.showRowNumbers(true);       

        int n = indexs.length;
        int num = this.reader.npillars;
        boolean having_dxy = this.reader.DX!=null;
        int mid_index = this.reader.nframes/2;
        for(int i=0; i<n; i++){
            IJ.showProgress(i, n);
            rt.incrementCounter();
            int c = indexs[i];
            if(c>=0 && c<num){
                rt.addValue("Pillar Index", c+1);
                StatisitcCenterDistance rst_raw = reader.getDistanceCenterRAW(c);
                double var_raw_x = BasicStatisitic.var(reader.raw_tracksX[c]);
                double var_raw_y = BasicStatisitic.var(reader.raw_tracksY[c]);
                double std_raw = Math.sqrt(var_raw_x + var_raw_y);
                //double avg_raw_x = rst_raw.avg_raw_x;//BasicStatisitic.avg(reader.raw_tracksX[c]);
                //double avg_raw_y = rst_raw.avg_raw_y;//BasicStatisitic.avg(reader.raw_tracksY[c]);
                
                rt.addValue("AVG_X_RAW", rst_raw.avg_x);
                rt.addValue("AVG_Y_RAW", rst_raw.avg_y);
                //rt.addValue("VAR_RAW_X", var_raw_x);
                //rt.addValue("VAR_RAW_Y", var_raw_y);
                rt.addValue("STD_RAW", std_raw);
                rt.addValue("AVG_DIS_RAW", rst_raw.avg_dis);
                rt.addValue("MAX_DIS_RAW", rst_raw.max_dis);
                rt.addValue("MIN_DIS_RAW", rst_raw.min_dis);
                rt.addValue("STD_DIS_RAW", rst_raw.std_dis);
                
                if(reader.correted_data){
                    double var_x = BasicStatisitic.var(reader.correct_tracksX[c]);
                    double var_y = BasicStatisitic.var(reader.correct_tracksY[c]);
                    //double avg_x = BasicStatisitic.avg(reader.correct_tracksX[c]);
                    //double avg_y = BasicStatisitic.avg(reader.correct_tracksY[c]);
                    double std = Math.sqrt(var_x + var_y);
                    
                    StatisitcCenterDistance rst = reader.getDistanceCenter(c);
                    rt.addValue("AVG_X", rst.avg_x);
                    rt.addValue("AVG_Y", rst.avg_y);
                    //rt.addValue("VAR_X", var_x);
                    //rt.addValue("VAR_Y", var_y);
                    rt.addValue("STD", std);
                    rt.addValue("AVG_DIS", rst.avg_dis);
                    rt.addValue("MAX_DIS", rst.max_dis);
                    rt.addValue("MIN_DIS", rst.min_dis);
                    rt.addValue("STD_DIS", rst.std_dis);
                }
                
                if(having_dxy){
                    double[] def = reader.getAbsDeflections(c)[2];                    
                    double avg_def = BasicStatisitic.avg(def);
                    double max_def = BasicStatisitic.max(def);
                    double min_def = BasicStatisitic.min(def);
                    double var_def = BasicStatisitic.var(def);
                    double std_def = Math.sqrt(var_def);
                    //Arrays.sort(def);
                    //double med_def = def[mid_index];
                    rt.addValue("AVG_Deflection", avg_def);
                    //rt.addValue("MEDIAN_Deflection", med_def);
                    rt.addValue("MAX_Deflection", max_def);
                    rt.addValue("MIN_Deflection", min_def);
                    //rt.addValue("VAR_Deflection", var_def);
                    rt.addValue("STD_Deflection", std_def);
                }
            }                            
        }
        IJ.showProgress(1);
        rt.show("List of Pillars with Basic Statistic");        
    }
        
    private void plot_pillar_deflection(boolean type){
        boolean raw = checkbox_raw.isSelected();   
        int shown_pillar_index = int_converter.fromString(textfield_pillar_no_drift.getText()) - 1;
        if(shown_pillar_index>=reader.npillars) return;
        double[] xaxis = reader.getFrameSeries();
        double[][] dxy = raw ? reader.getDeflectionsRAW(shown_pillar_index) : reader.getDeflections(shown_pillar_index);
        String name = raw ? "pillar track#" + (shown_pillar_index+1) + "(raw)" : "pillar track#" + (shown_pillar_index+1) + "(drift corrected)";
        if(type){
            //Plot ploter_ipillar = DriftAnalysisPlotter.get_plot_driftXY(xaxis,dxy[0],dxy[1],name);
            Plot ploter_ipillar = DriftAnalysisPlotter.plot2(xaxis, dxy[0],dxy[1],name, "frame #", "displacement(nm) refer to 1st frame");
            ploter_ipillar.setLegend("deflection X Axis\ndeflection Y Axis",Plot.AUTO_POSITION|Plot.LEGEND_TRANSPARENT);
            ploter_ipillar.show();	
        }
        else{
            Plot plot_time = new Plot(name,"frame #", "displacement(nm)");				
            plot_time.addPoints(xaxis, dxy[2], Plot.LINE);			
            plot_time.draw();	
            plot_time.show();        
        }
    }
    
    private void updateDriftPlotWindow(){
        if(!load_suc_drift) return;
//        boolean raw = checkbox_raw.isSelected();   
        int shown_pillar_index = int_converter.fromString(textfield_pillar_no_drift.getText()) - 1;
        if(shown_pillar_index>=reader.npillars) return;
        
        if(overlay_tracks_ip!=null && overlay_tracks_ip.isVisible()){
            double[] xy = reader.getXY(shown_pillar_index);
            if(xy!=null){
                PointRoi roi = new PointRoi(xy[0]+0.5, xy[1]+0.5);
                overlay_tracks_ip.setRoi(roi);
                overlay_tracks_ip.updateAndDraw();
            }
        }
    }
    
    private void updateCheckboxShowLabels(){
        if(!load_suc_drift) return;
                
        //boolean flatten = checkbox_flatten.isSelected();
        boolean show = checkbox_show_labels.isSelected();
        boolean raw = checkbox_raw.isSelected();        
        
        if(flatten_tracks_ip!=null && flatten_tracks_ip.isVisible() ){            
            if(show){
                double[][] cx = raw ? reader.raw_tracksX : reader.correct_tracksX;
                double[][] cy = raw ? reader.raw_tracksY : reader.correct_tracksY;

                Overlay labels = plotter.getLabels(cx, cy, pixel_size, zoom_in_image);
                flatten_tracks_ip.setOverlay(labels);
            }
            else flatten_tracks_ip.setOverlay(null);
        }
        
        if(overlay_tracks_ip!=null && overlay_tracks_ip.isVisible()){
            if(show){
                double[][] cx = raw ? reader.raw_tracksX : reader.correct_tracksX;
                double[][] cy = raw ? reader.raw_tracksY : reader.correct_tracksY;                
                //Roi roi = overlay_tracks_ip.getRoi();            
                //int[] indexs = reader.getSelectedPillars(roi);
                int[] indexs = this.selected_pillars_indexs;
                plotter.getLabels(cx, cy, indexs, pixel_size);         
                //plotter.getLabels(cx, cy, pixel_size);                
            }
            
            plotter.plot_tracks_overlay(overlay_tracks_ip, show);
        }
        
    }
    
    private void updateSliderZoominTracks(){
       if(!load_suc_drift) return;        
        
       if(overlay_tracks_ip!=null && overlay_tracks_ip.isVisible()){
            //boolean flatten = checkbox_flatten.isSelected();
            boolean show = checkbox_show_labels.isSelected();
            boolean raw = checkbox_raw.isSelected(); 
            boolean still_pillars = checkbox_still_pillars.isSelected();  
                        
            boolean shownin = checkbox_shownin.isSelected();
            if(!shownin){
                File input_file = new File(drift_fname); 
                String title = input_file.getName() + " x" + zoom_in_tracks;
                overlay_tracks_ip.setTitle(title);
            }            
            
            int[] indexs = selected_pillars_indexs;

            int frame_window = int_converter.fromString(textfield_frame_window.getText());
            if(frame_window>0){
                int end_f = (int)this.slider_frame_no_drift.getValue()-1;
                int start_f = end_f - frame_window;
                draw_list_pillars(selected_pillars_indexs,start_f,end_f);
                return;
            }

            if(indexs!=null && indexs.length>0){
                if(show){
                    double[][] cx = raw ? reader.raw_tracksX : reader.correct_tracksX;
                    double[][] cy = raw ? reader.raw_tracksY : reader.correct_tracksY;
                    plotter.getLabels(cx, cy, indexs, pixel_size);
                } 
                
                double[][] cx = raw ? reader.getSubset(indexs,reader.raw_tracksX) : reader.getSubset(indexs,reader.correct_tracksX);
                double[][] cy = raw ? reader.getSubset(indexs,reader.raw_tracksY) : reader.getSubset(indexs,reader.correct_tracksY);
                plotter.getTracks(cx, cy, pixel_size, zoom_in_tracks); 
            }
                        
            plotter.plot_tracks_overlay(overlay_tracks_ip, show); 
        }
    }
    // </editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="grid analysis event handing">    
    private void updateDeflectionMap(){
        if(!load_suc_grid) return;
        int shown_frame_index = int_converter.fromString(textfield_frame_no_grid.getText()) - 1;
        if(shown_frame_index>=grid_reader.frames) return;
        if(overlay_deflections_ip!=null && overlay_deflections_ip.isVisible()){
            Roi roi = overlay_deflections_ip.getRoi();
            if(roi==null){
                roi = new Roi(0,0,overlay_deflections_ip.getWidth(), overlay_deflections_ip.getHeight());
            }

            //current_frame_index = shown_frame_index;
            candidate_grid = grid_reader.getSelections(roi, zoom_in_image_grid);   
            if(this.checkbox_remove_large.isSelected()){
                String text = this.textfield_large_threshold.getText();
                if(text!=null && !text.isEmpty()){
                    double threshold = double_converter.fromString(text);
                    if(threshold>0) candidate_grid = grid_reader.remove_large_deflection(candidate_grid, shown_frame_index, threshold);
                }
            }
            //showLabelsArrows(candidate_grid, shown_frame_index);
            draw_highlight_cu(shown_frame_index);
           
        }
    }
    
    private void showLabelsArrows(int[] candidate, boolean[] cu_flag, int shown_frame_index){        
        if(shown_frame_index<0 || shown_frame_index>=grid_reader.frames) return;
        if(overlay_deflections_ip==null || !overlay_deflections_ip.isVisible()) return;
        double[][] xy = grid_reader.getXY(candidate, shown_frame_index);
        if(xy!=null){
            double[][] gxy = grid_reader.getGrids(candidate, shown_frame_index);
            if(gxy!=null){
                if(shown_frame_index<overlay_deflections_ip.getNSlices()){
                    overlay_deflections_ip.setSlice(shown_frame_index+1);
                    overlay_deflections_ip.updateAndRepaintWindow();
                }
                else if(shown_frame_index<overlay_deflections_ip.getNFrames()){
                    overlay_deflections_ip.setT(shown_frame_index+1);
                    overlay_deflections_ip.updateAndRepaintWindow();
                }                    
                //Overlay points = grid_plotter.getPoints(gxy[0], gxy[1]);
                boolean[] sub_cu_flag = grid_reader.get_CU_flag(candidate, cu_flag, shown_frame_index);
                Overlay arrows = grid_plotter.getPointsArrows(gxy[0], gxy[1], xy[0], xy[1], sub_cu_flag, zoom_in_image_grid, zoom_in_deflection);
                Overlay labels = null;
                boolean showlabel = checkbox_show_labels_grid.isSelected();
                boolean showdxy = checkbox_show_labels_grid.isIndeterminate();
                if(showlabel || showdxy){                        
                    labels = grid_reader.getLabels(candidate, shown_frame_index, zoom_in_image_grid, showdxy);
                }
                
                grid_plotter.plot_arrow_overlay(arrows,overlay_deflections_ip,labels);                                
            }
        }
    }
    
    private void showLabelsArrows(int[] candidate, int shown_frame_index){        
        if(shown_frame_index<0 || shown_frame_index>=grid_reader.frames) return;
        if(overlay_deflections_ip==null || !overlay_deflections_ip.isVisible()) return;
        double[][] xy = grid_reader.getXY(candidate, shown_frame_index);
        if(xy!=null){
            double[][] gxy = grid_reader.getGrids(candidate, shown_frame_index);
            if(gxy!=null){
                if(shown_frame_index<overlay_deflections_ip.getNSlices()){
                    overlay_deflections_ip.setSlice(shown_frame_index+1);
                    overlay_deflections_ip.updateAndRepaintWindow();
                }
                else if(shown_frame_index<overlay_deflections_ip.getNFrames()){
                    overlay_deflections_ip.setT(shown_frame_index+1);
                    overlay_deflections_ip.updateAndRepaintWindow();
                }                    
                //Overlay points = grid_plotter.getPoints(gxy[0], gxy[1]);
                Overlay arrows = grid_plotter.getPointsArrows(gxy[0], gxy[1], xy[0], xy[1], zoom_in_image_grid, zoom_in_deflection);
                Overlay labels = null;
                boolean showlabel = checkbox_show_labels_grid.isSelected();
                boolean showdxy = checkbox_show_labels_grid.isIndeterminate();
                if(showlabel || showdxy){                        
                    labels = grid_reader.getLabels(candidate, shown_frame_index, zoom_in_image_grid, showdxy);
                }
                
                grid_plotter.plot_arrow_overlay(arrows,overlay_deflections_ip,labels);                
            }
        }
    }
    
    public static double[][] getSubSet(int[] candidates, double[]x, double[] y, boolean[] active){
        if(candidates==null) return null;
        int n = candidates.length;
        int na = 0;        
        for(int jj=0; jj<n; jj++) {
            int c=candidates[jj];
            if(active[c]) na++;                
        }
        if(na<1) return null;
        
        double[][] XY = new double[2][na];
        na = 0;
        for(int jj=0; jj<n; jj++) {
            int c=candidates[jj];
            if(active[c]) {
                XY[0][na] = x[c];
                XY[1][na] = y[c];
                na++;
            }
        }
        
        return XY;
    }    
        
    private void updateDeflectionMap(double[] xc, double[] yc, double[] dx, double[] dy, boolean[] active_d){
            if(!load_suc_grid) return;
            int shown_frame_index = int_converter.fromString(textfield_frame_no_grid.getText()) - 1;
            if(shown_frame_index>=grid_reader.frames) return;
            if(overlay_deflections_ip!=null && overlay_deflections_ip.isVisible()){
                Roi roi = overlay_deflections_ip.getRoi();
                if(roi==null || !roi.isArea()){
                    roi = new Roi(0,0,overlay_deflections_ip.getWidth(), overlay_deflections_ip.getHeight());
                }
                
                candidate_grid = grid_reader.getSelections(roi, zoom_in_image_grid);   
                if(this.checkbox_remove_large.isSelected()){
                    String text = this.textfield_large_threshold.getText();
                    if(text!=null && !text.isEmpty()){
                        double threshold = double_converter.fromString(text);
                        if(threshold>0) candidate_grid = grid_reader.remove_large_deflection(candidate_grid, dx, dy, active_d, threshold);
                    }
                }
                showLabelsArrows(xc, yc, dx, dy,active_d, candidate_grid, shown_frame_index);            
            }
    }
    
    private void showLabelsArrows(double[] xc, double[] yc, double[] dx, double[] dy, boolean[] active_d, int[] candidate, int shown_frame_index){        
        if(shown_frame_index<0 || shown_frame_index>=grid_reader.frames) return;
        double[][] xy = getSubSet(candidate, xc, yc, active_d);
        if(xy!=null){
            double[][] dxy = getSubSet(candidate, dx, dy, active_d);
            double[][] gxy = GridDataLoader.getGrids(xy[0], xy[1], dxy[0], dxy[1]);
            if(gxy!=null){
                if(shown_frame_index<overlay_deflections_ip.getNSlices()){
                    overlay_deflections_ip.setSlice(shown_frame_index+1);
                    overlay_deflections_ip.updateAndRepaintWindow();
                }
                else if(shown_frame_index<overlay_deflections_ip.getNFrames()){
                    overlay_deflections_ip.setT(shown_frame_index+1);
                    overlay_deflections_ip.updateAndRepaintWindow();
                }                    
                //Overlay points = grid_plotter.getPoints(gxy[0], gxy[1]);
                Overlay arrows = grid_plotter.getPointsArrows(gxy[0], gxy[1], xy[0], xy[1], zoom_in_image_grid, zoom_in_deflection);
                Overlay labels = null;
                boolean showlabel = checkbox_show_labels_grid.isSelected();
                boolean showdxy = checkbox_show_labels_grid.isIndeterminate();
                if(showlabel || showdxy){                        
                    labels = grid_reader.getLabels(candidate, xc, yc, dx, dy, zoom_in_image_grid, showdxy);
                }
                else{
                
                }
                //Overlay labels = showlabel? grid_reader.getLabels(shown_frame_index) : null; 
                grid_plotter.plot_arrow_overlay(arrows,overlay_deflections_ip,labels);                
                //overlay_deflections_ip.setOverlay(arrows);
                //overlay_deflections_ip.updateAndDraw();
            }
        }
    }
    
    @FXML private void handleButtonShowLineGrid(ActionEvent event){
        SwingUtilities.invokeLater(new Runnable() {
            @Override public void run() {
                SetGridDialogResult dlg_rst = showSetGridDialog();
                if(dlg_rst.suc){
                    updateShowLineGrid(dlg_rst.s,dlg_rst.a, dlg_rst.b);
                    //updateLabelFileInfo();
                }
            }
        });        
    }
    
    
    private void updateShowLineGrid(double spacing, double oblique, double grid_angle){
        if(!load_suc_grid) return;
        
        boolean shownin = checkbox_shownin_grid.isSelected();
        if(shownin && WindowManager.getImageCount()>0){
            ImagePlus ip = WindowManager.getCurrentImage();//IJ.getImage();
            if(ip!=null){
                    overlay_deflections_ip = ip;
                    need_recreate_map_grid = false;
            }
        }
        else{
            overlay_deflections_ip = created_deflection_map_ip;
        }
        
        if(overlay_deflections_ip!=null){            
                Overlay lines = null;
                int shown_frame_index = int_converter.fromString(textfield_frame_no_grid.getText()) - 1;
                double[][] gxy = grid_reader.getGrids(shown_frame_index);
//                if(Double.isNaN(grid_reader.grid_oblique))
//                    lines = GetLineGrid(gxy[0],gxy[1],grid_reader.spacing,zoom_in_image_grid);
//                else{
//                    lines = GetLineGrid(gxy[0],gxy[1],grid_reader.spacing,grid_reader.grid_oblique, grid_reader.grid_angle, zoom_in_image_grid);
//                }
                if(!Double.isNaN(grid_reader.grid_oblique)){
                    //lines = GetLineGrid(gxy[0],gxy[1],spacing,oblique,grid_angle,zoom_in_image_grid);
                    GridCreation grid_creator = new GridCreation(gxy[0],gxy[1]);
                    grid_creator.label_grid(spacing,oblique, grid_angle);
                    int[] label_x = grid_creator.getIX();
                    int[] label_y = grid_creator.getIY();
                    boolean[] is_labeled = grid_creator.IsLabeled();
                    lines = grid_plotter.plotGrid(gxy[0],gxy[1],label_x,label_y,is_labeled, zoom_in_image_grid);   
                    show2DIndexTable(label_x, label_y, is_labeled, shown_frame_index);                    
                }
                if(lines!=null) overlay_deflections_ip.setOverlay(lines);            
        }
    }
    
    private void show2DIndexTable(int[] index_x, int[] index_y, boolean[] is_labeled, int frame){
        ResultsTable rt = new ResultsTable();        
        rt.setNaNEmptyCells(true);
        rt.showRowNumbers(false);
        
        int n = candidate_grid==null ? grid_reader.cens : candidate_grid.length;
        
        //for(int frame=0; frame<grid_reader.frames; frame++){
            for(int i=0; i<n; i++){
                rt.incrementCounter();
                int c = candidate_grid==null ? i : candidate_grid[i];
                if(grid_reader.active[c][frame]){                
                    rt.addValue("Frame Index", frame+1);
                    rt.addValue("Pillar Index", c+1);
                    if(is_labeled[c]){
                        rt.addValue("IX", index_x[c]);
                        rt.addValue("IY", index_y[c]);                    
                    }
                    
                    double x = grid_reader.trackX[c][frame];
                    double y = grid_reader.trackY[c][frame];
                    double dx = grid_reader.DX[c][frame];
                    double dy = grid_reader.DY[c][frame];                    
                    double gx = x-dx;
                    double gy = y-dy;
                    
                    rt.addValue("X", x);
                    rt.addValue("Y", y);                    
                    rt.addValue("GX", gx);
                    rt.addValue("GY", gy);                    
                    rt.addValue("DX", dx);
                    rt.addValue("DY", dy);  
                    rt.addValue("DIS", Math.sqrt(dx*dx+dy*dy));  
                }
            }
        //}
        
        rt.show("Data(Index2D) for Frame=" + (frame+1));        
    }
    
    private void plot_current_pillar_deflections(int shown_pillar_index, boolean type){
        //int shown_pillar_index = int_converter.fromString(textfield_pillar_no_grid.getText()) - 1; 
        if(shown_pillar_index>=grid_reader.cens) return;
        double[] xaxis = grid_reader.getFrameSeries();        
        double[][] dxy = grid_reader.getPillarDXY(shown_pillar_index);
        String name = "pillar track#" + (shown_pillar_index+1);
        
        if(type){
            Plot ploter_ipillar = grid_plotter.get_plot_driftXY(xaxis,dxy[0],dxy[1],name);  
            ploter_ipillar.setLegend("deflection X Axis\ndeflection Y Axis\ndeflection",Plot.AUTO_POSITION|Plot.LEGEND_TRANSPARENT);
            ploter_ipillar.show();        
        }
        else{
            Plot plot_time = new Plot(name,"frame #", "deflection(nm)");				
            plot_time.addPoints(xaxis, dxy[2], Plot.LINE);			
            plot_time.draw();	
            plot_time.show();
        }
    }
    
    private void updateDeflectionPlot(){
        if(!load_suc_grid) return;
        int shown_pillar_index = int_converter.fromString(textfield_pillar_no_grid.getText()) - 1;        
        if(shown_pillar_index>=grid_reader.cens) return;
        int shown_frame_index = int_converter.fromString(textfield_frame_no_grid.getText()) - 1;
        if(shown_frame_index>=grid_reader.frames) return;
        
        if(overlay_deflections_ip!=null && overlay_deflections_ip.isVisible()){
            double[] xy = grid_reader.getXY(shown_pillar_index);
            if(xy!=null){
                PointRoi roi = new PointRoi(xy[0]*zoom_in_image_grid, xy[1]*zoom_in_image_grid);
                overlay_deflections_ip.setRoi(roi);
                overlay_deflections_ip.updateAndDraw();
                int[] candidate = {shown_pillar_index};                
                showLabelsArrows(candidate, shown_frame_index);
            }
        }
    }
    
    private Overlay GetLineGrid(double[] x, double[] y, double seperation, double grid_oblique, double grid_angle, double zoom){
        GridCreation grid_creator = new GridCreation(x,y);
        grid_creator.label_grid(seperation, grid_oblique, grid_angle);
        Overlay grid = grid_plotter.plotGrid(x, y, grid_creator.getIX(), grid_creator.getIY(), grid_creator.IsLabeled(), zoom);
        //test line grid without pinchusion model:      
        /*
        double a = grid_oblique*Math.PI/180;
        Rectangle rect = new Rectangle(0,0,grid_creator.getMapWidth(),grid_creator.getMapHeight());
        double[] xy0 = grid_creator.getFisrtlabeledXY();
        double x0=0;
        double y0=0;
        if(xy0!=null){ x0=xy0[0]; y0=xy0[1];}
        IJ.log("x0=" + x0 + " y0=" + y0);
        Overlay grid = grid_plotter.plotGrid(x0, y0, seperation*Math.cos(a), -seperation*Math.sin(a), rect, zoom);
        */
        return grid;
        //ip.setOverlay(grid);
    }
    
    private Overlay GetLineGrid(double[] x, double[] y, double seperation, double zoom){
        GridCreation grid_creator = new GridCreation(x,y);
        double[] sa = grid_creator.estimate_square_grid_properties(seperation, seperation*0.2);
        if(sa==null) return null;        
        grid_creator.label_grid(sa[0], sa[1], 90);
        Overlay grid = grid_plotter.plotGrid(x, y, grid_creator.getIX(), grid_creator.getIY(), grid_creator.IsLabeled(), zoom);
        return grid;
        //ip.setOverlay(grid);
    }
    
    @FXML private void handleButtonDeflectionMap(ActionEvent event){
        boolean shownin = checkbox_shownin_grid.isSelected();
        if(shownin && WindowManager.getImageCount()>0){
            ImagePlus ip = WindowManager.getCurrentImage();//IJ.getImage();
            if(ip!=null){
                    overlay_deflections_ip = ip;
                    need_recreate_map_grid = false;
            }
        }
        else{
            overlay_deflections_ip = created_deflection_map_ip;
        }
        
        zoom_in_image_grid = int_converter.fromString(textfield_zoomin_image_grid.getText());
        if(overlay_deflections_ip==null || !overlay_deflections_ip.isVisible() || need_recreate_map_grid){
            File input_file = new File(grid_fname);	
            String title = input_file.getName();
            overlay_deflections_ip = created_deflection_map_ip = grid_reader.create_empty_map(zoom_in_image_grid, 5, title);
            overlay_deflections_ip.show();
            need_recreate_map_grid = false;
        }        
        
        int shown_frame_index = int_converter.fromString(textfield_frame_no_grid.getText()) - 1;
        if(shown_frame_index>=grid_reader.frames) return;
        if(overlay_deflections_ip!=null && overlay_deflections_ip.isVisible()){
            Roi roi = overlay_deflections_ip.getRoi();
            if(roi==null || !roi.isArea()){
                roi = new Roi(0,0,overlay_deflections_ip.getWidth(), overlay_deflections_ip.getHeight());
            }            
            candidate_grid = grid_reader.getSelections(roi, zoom_in_image_grid);    
            if(this.checkbox_remove_large.isSelected()){
                String text = this.textfield_large_threshold.getText();
                if(text!=null && !text.isEmpty()){
                    double threshold = double_converter.fromString(text);
                    if(threshold>0) candidate_grid = grid_reader.remove_large_deflection(candidate_grid, shown_frame_index, threshold);
                }                
            }
            //showLabelsArrows(candidate_grid, shown_frame_index);
            draw_highlight_cu(shown_frame_index);
        }
        //updateDeflectionMap();
    }
    
    private void updateLabelFileInfo(){
        Platform.runLater(new Runnable() {
            @Override public void run() {
                label_grid_load_info.setText(String.format("Frames:%d Pillars:%d Pixel Size:%.2fnm Spacing:%.2f Diameter:%.2f Grid Oblique:%.2f\u00b0 Grid Angle:%.2f\u00b0",
                        grid_reader.frames,
                        grid_reader.cens,
                        grid_reader.lattice,
                        grid_reader.spacing,
                        grid_reader.diameter,
                        grid_reader.grid_oblique,
                        grid_reader.grid_angle)
                );
            }
        });
    }
    
    private boolean showReCreateGridDialog(boolean set_grid){
        double threshold_deflection = threshold_deflection_dxy;
        if(Double.isNaN(threshold_deflection_dxy)){
            String text = textfield_large_threshold.getText();
            if(!text.isEmpty()) threshold_deflection =  double_converter.fromString(text);
        }
        
        GenericDialogPlus gd = new GenericDialogPlus("Settings for Grid Recreation");
        //gd.setSize(300, 100);
        double s = grid_reader.spacing;
        double a = grid_reader.grid_oblique;
        double b = grid_reader.grid_angle;
        if(set_grid){
            //estimate the current separtion and oblique angle for the square grid.
            double[] sa = settingGrid();    
            if(sa!=null){
                int len = sa.length;
                if(len>0){
                    s = sa[0];
                    if(len>1){
                        a = sa[1];
                        if(len>2) b = sa[2];
                    }
                }
            }
                //gd.setSize(300, 100);
            gd.addNumericField("        spacing:", s, 2, 10, "pixel");
            gd.addNumericField("        grid_oblique:", a, 2, 10, "degree");		
            gd.addNumericField("        grid_angle:", b, 2, 10, "degree");	
        }
        
        gd.addNumericField("        deflection threshold:", threshold_deflection, 2, 10, "nm");
        gd.addCheckbox(" process all frames?", false);
        String new_grid_fname = grid_fname;
        if(grid_fname.endsWith(".dxy")) new_grid_fname = grid_fname.replaceAll(".dxy", ".grid");
        
        gd.addFileField("save to:", new_grid_fname, 50);
        gd.showDialog();
        if (gd.wasCanceled()) return false;
        
        if(set_grid){
            s = gd.getNextNumber();
            a = gd.getNextNumber();
            b = gd.getNextNumber();
            if(s<0){
                IJ.showMessage("the spacing must be larger than 0");
                return false;
            }            
        }
        
        threshold_deflection = gd.getNextNumber();
        if(threshold_deflection<=0){
            IJ.showMessage("the threshold must be larger than 0");
            return false;
        }
        boolean is_all_frames = gd.getNextBoolean();
        new_grid_fname = gd.getNextString();
        
        if(set_grid){            
            grid_reader.spacing = s;
            grid_reader.grid_oblique = a;
            grid_reader.grid_angle = b;
        }
        threshold_deflection_dxy = threshold_deflection;
        
        return true;
    }
    
    private SetGridDialogResult showSetGridDialog(){
        double s = grid_reader.spacing;
        double a = grid_reader.grid_oblique;
        double b = grid_reader.grid_angle;
        SetGridDialogResult result = new SetGridDialogResult(s,a,b,false);   
        
        //estimate the current separtion and oblique angle for the square grid.
        double[] sa = settingGrid();    
        if(sa!=null){
            int len = sa.length;
            if(len>0){
                s = sa[0];
                if(len>1){
                    a = sa[1];
                    if(len>2) b = sa[2];
                }
            }
        }
        
        GenericDialog gd = new GenericDialog("Settings for Grid Properties");
        //gd.setSize(300, 100);
        gd.addNumericField("        spacing:", s, 2, 10, "pixel");
        gd.addNumericField("        grid_oblique:", a, 2, 10, "degree");		
        gd.addNumericField("        grid_angle:", b, 2, 10, "degree");		
        gd.showDialog();
               
        if (gd.wasCanceled()) return result;
        
        s = gd.getNextNumber();
        a = gd.getNextNumber();
        b = gd.getNextNumber();
        result = new SetGridDialogResult(s,a,b,true);  
//        grid_reader.spacing = gd.getNextNumber();
//        grid_reader.grid_oblique = gd.getNextNumber();
//        grid_reader.grid_angle = gd.getNextNumber();
//        return true;

        return result;
        
    }
    
    private double[] getFirstPoint(PointRoi roi){
        if(roi==null) return null;
        FloatPolygon p = roi.getFloatPolygon();
        int n = p.npoints;
        if(n<1) return null;
        double minx = p.xpoints[0];
        double miny = p.ypoints[0];
        
        double[] xy = {minx, miny};
        if(n<2) return xy;
        
        double mind2 = minx*minx + miny*miny;                 
        for(int i=1; i<n; i++){

            double dx = p.xpoints[i];
            double dy = p.ypoints[i];
            double d2 = dx*dx+dy*dy;
            if(d2<mind2){                    
                mind2 = d2;
                minx = dx;
                miny = dy;
            }
        }
        
        xy[0] = minx;
        xy[1] = miny;
        return xy;
    }
    
    private boolean showPrevewGridDialog(){
        int n = WindowManager.getImageCount();
        if(n<1){
            IJ_showMessage("there is no image opened!");
            return false;
        }
        
        ImagePlus ip = IJ.getImage();
        if(ip==null || !ip.isVisible()){
            IJ_showMessage("there is no image opened to draw on!");
            return false;
        }
        
        String text = textfiled_spacing.getText();
        double s = text.isEmpty() ? Double.NaN : double_converter.fromString(text);
        text = textfiled_grid_oblique.getText();
        double a = text.isEmpty() ? Double.NaN : double_converter.fromString(text);
        text = textfiled_grid_angle.getText();
        double b = text.isEmpty() ? Double.NaN : double_converter.fromString(text);
        double x0 = 0;
        double y0 = 0;
        
        boolean show = true;
        Roi roi = ip.getRoi();
        if(PointRoi.class.isInstance(roi)){
            double[] xy0 = getFirstPoint((PointRoi)roi);
            if(xy0!=null){
                x0 = xy0[0];
                y0 = xy0[1];
                show = false;
            }
        }
        
        if(show){
            GenericDialog gd = new GenericDialog("Settings for Grid Properties");
            //gd.setSize(300, 100);
            gd.addNumericField("        spacing:",  s, 2, 10, "");
            gd.addNumericField("        grid_oblique:",  a, 2, 10, "");	
            gd.addNumericField("        grid_angle:",  b, 2, 10, "");	
            gd.addNumericField("        x0:",  x0, 2, 10, "");
            gd.addNumericField("        y0:",  y0, 2, 10, "");	

            gd.showDialog();
            if (gd.wasCanceled()) return false;

            s = gd.getNextNumber();        
            a = gd.getNextNumber();
            b = gd.getNextNumber();
            x0 = gd.getNextNumber();
            y0 = gd.getNextNumber();
        
            if(s<0){
                IJ_showMessage("the spacing must be larger than zero!");
                return false;
            }
        }
        
        GridAnalysisPlotter.plotGrid(x0, y0, s, a, b, ip);
//        double aa = a*Math.PI/180;
//        double dx = s*Math.cos(aa);
//        double dy = -s*Math.sin(aa);
//
//        GridAnalysisPlotter.plotGrid(x0, y0, dx, dy, ip);
        
        return true;
    }
    
    private double[] settingGrid(){
        if(!load_suc_grid) return null;
        int shown_frame_index = int_converter.fromString(textfield_frame_no_grid.getText()) - 1;
        if(shown_frame_index>=grid_reader.frames) return null;
        double[][] gxy = grid_reader.getGrids(shown_frame_index);
        double s = GridCreation.estimate_seperation(gxy[0], gxy[1]);  
        if(Double.isNaN(s)) s = grid_reader.spacing;
        
        double[] sa = {s};
        if(grid_reader.grid_angle==90)
        {
            sa = GridCreation.estimate_square_grid_properties(gxy[0], gxy[1], s, s*0.2);
        }
        else if(grid_reader.IX!=null){
            int [][] ixy = grid_reader.getGridsIXY(shown_frame_index);
            sa = GridCreation.esimate_grid_properties(gxy[0], gxy[1], ixy[0], ixy[1], s, s*0.2);
        }
        return sa;        
    }
    
    private double estimate_current_grid_spacing(){
        if(!load_suc_grid) return Double.NaN;
        int shown_frame_index = int_converter.fromString(textfield_frame_no_grid.getText()) - 1;
        if(shown_frame_index>=grid_reader.frames) return Double.NaN;
        double[][] gxy = grid_reader.getGrids(shown_frame_index);
        double s = GridCreation.estimate_seperation(gxy[0], gxy[1]);  
        return s;        
    }
   
    @FXML private void handleButtonLoadFileGrid(ActionEvent event){            
        
        String fpath = textfile_grid_fname.getText();
        List<ContractionUnit> cu = Hunt_CU_PlugIn.load(fpath);
        if(cu!=null) {
            cu_list_grid = cu;
            if(load_suc_grid){
                ContractionUnit.show_table_list(cu_list_grid);
                int f = (int)this.slider_frame_no_grid.getValue()-1;
                draw_highlight_cu(f);
            }
            return;
        }
        
        load_suc_grid =  grid_reader.data_loader(fpath) && grid_reader.frames>0;                
        accordion_grid.setDisable(!load_suc_grid);  
        accordion_grid.setExpanded(load_suc_grid);
        if(load_suc_grid){
          
            button_recreate_grid.disableProperty().set(grid_reader.getFileVersion()==pillar_tracking_FD.fileversion2);
            updateLabelFileInfo();
            
            cu_hunter_grid = new Hunt_CU_PlugIn();
            cu_list_grid = null;            
            contraction_unit_list = null;
            overlay_deflections_ip = null;            
            //label_grid_load_info.setText("Frames:" + grid_reader.frames + "  Pillars:" + grid_reader.cens + "  Pixel Size(nm):" + grid_reader.lattice + "  Spacing:" + grid_reader.spacing + "  Diameter:" + grid_reader.diameter);
            slider_frame_no_grid.maxProperty().set(grid_reader.frames);
            slider_pillar_no_grid.maxProperty().set(grid_reader.cens);
            //File input_file = new File(fpath);		            
            grid_fname = fpath;//input_file.getName();
            need_recreate_map_grid = true;
        }        
    }
    
    @FXML private void handleButtonGetDataFrameDrift(ActionEvent event){
        int frame = int_converter.fromString(textfield_frame_no_drift.getText()) - 1;
        if(frame>=reader.nframes) return;
        
        ResultsTable rt = new ResultsTable();        
        rt.setNaNEmptyCells(true);
        rt.showRowNumbers(false);
        
        int n = selected_pillars_indexs==null ? reader.npillars : selected_pillars_indexs.length;
        boolean having_dxy = (reader.DX!=null);
        for(int i=0; i<n; i++){
            rt.incrementCounter();
            int c = selected_pillars_indexs==null ? i : selected_pillars_indexs[i];
            {                
                rt.addValue("Pillar Index", c+1);
                rt.addValue("Raw X", reader.raw_tracksX[c][frame]);
                rt.addValue("Raw Y", reader.raw_tracksY[c][frame]);
                
                if(reader.correted_data){
                    rt.addValue("corrected X", reader.correct_tracksX[c][frame]);
                    rt.addValue("corrected Y", reader.correct_tracksY[c][frame]);
                }
                
                if(having_dxy){
                    double dx = reader.DX[c][frame];
                    double dy = reader.DY[c][frame];
                    rt.addValue("DX", dx);
                    rt.addValue("DY", dy);  
                    rt.addValue("DIS", Math.sqrt(dx*dx+dy*dy));  
                }
                else if(reader.correted_data){
                    double[] xy = reader.getDeflections(c,frame,0);  
                    rt.addValue("DX",  xy[0]);
                    rt.addValue("DY",  xy[1]);  
                    rt.addValue("DIS", xy[2]);                                
                }
            }
        }
        rt.show("Data for Frame=" + (frame+1));
        
    }
    
    @FXML private void handleButtonGetALLDataFramesDrift(ActionEvent event){
        
        ResultsTable rt = new ResultsTable();        
        rt.setNaNEmptyCells(true);
        rt.showRowNumbers(false);
        
        int n = selected_pillars_indexs==null ? reader.npillars : selected_pillars_indexs.length;
        boolean having_dxy = (reader.DX!=null);
        
        for(int frame=0; frame<reader.nframes; frame++){
            for(int i=0; i<n; i++){
                rt.incrementCounter();
                int c = selected_pillars_indexs==null ? i : selected_pillars_indexs[i];
                {                
                    rt.addValue("Frame Index", frame+1);
                    rt.addValue("Pillar Index", c+1);
                    rt.addValue("Raw X", reader.raw_tracksX[c][frame]);
                    rt.addValue("Raw Y", reader.raw_tracksY[c][frame]);

                    if(reader.correted_data){
                        rt.addValue("corrected X", reader.correct_tracksX[c][frame]);
                        rt.addValue("corrected Y", reader.correct_tracksY[c][frame]);
                    }

                    if(having_dxy){
                        double dx = reader.DX[c][frame];
                        double dy = reader.DY[c][frame];
                        rt.addValue("DX", dx);
                        rt.addValue("DY", dy);  
                        rt.addValue("DIS", Math.sqrt(dx*dx+dy*dy));  
                    }
                    else if(reader.correted_data){
                        double[] xy = reader.getDeflections(c,frame,0);  
                        rt.addValue("DX",  xy[0]);
                        rt.addValue("DY",  xy[1]);  
                        rt.addValue("DIS", xy[2]);                                
                    }
                }
            }
        }
        
        rt.show("Data for ALL Frames");        
    }
    
    @FXML private void handleButtonGetDataFrame(ActionEvent event){
        int frame = int_converter.fromString(textfield_frame_no_grid.getText()) - 1;
        if(frame>=grid_reader.frames) return;
        //double[][] xy = grid_reader.getFrameXY(shown_frame_index);
        ResultsTable rt = new ResultsTable();        
        rt.setNaNEmptyCells(true);
        rt.showRowNumbers(false);

        int n = candidate_grid==null ? grid_reader.cens : candidate_grid.length;
        for(int i=0; i<n; i++){
            rt.incrementCounter();
            int c = candidate_grid==null ? i : candidate_grid[i];
            if(grid_reader.active[c][frame]){                
                rt.addValue("Pillar Index", c+1);
                rt.addValue("X", grid_reader.trackX[c][frame]);
                rt.addValue("Y", grid_reader.trackY[c][frame]);
                double dx = grid_reader.DX[c][frame];
                double dy = grid_reader.DY[c][frame];
                rt.addValue("DX", dx);
                rt.addValue("DY", dy);  
                rt.addValue("DIS", Math.sqrt(dx*dx+dy*dy));  
            }
        }
        rt.show("Data for Frame=" + (frame+1));
        
    }
    
    @FXML private void handleButtonGetDataALLFrames(ActionEvent event){

        ResultsTable rt = new ResultsTable();        
        rt.setNaNEmptyCells(true);
        rt.showRowNumbers(false);
        
        int n = candidate_grid==null ? grid_reader.cens : candidate_grid.length;
        
        for(int frame=0; frame<grid_reader.frames; frame++){
            for(int i=0; i<n; i++){
                rt.incrementCounter();
                int c = candidate_grid==null ? i : candidate_grid[i];
                if(grid_reader.active[c][frame]){                
                    rt.addValue("Frame Index", frame+1);
                    rt.addValue("Pillar Index", c+1);
                    rt.addValue("X", grid_reader.trackX[c][frame]);
                    rt.addValue("Y", grid_reader.trackY[c][frame]);
                    double dx = grid_reader.DX[c][frame];
                    double dy = grid_reader.DY[c][frame];
                    rt.addValue("DX", dx);
                    rt.addValue("DY", dy);  
                    rt.addValue("DIS", Math.sqrt(dx*dx+dy*dy));  
                }
            }
        }
        
        rt.show("Data for ALL Frames");        
    }
    
    @FXML private void handleMenuPlotPillarDXYDrift(ActionEvent event){
        //int c = int_converter.fromString(textfield_pillar_no.getText()) - 1;
        plot_pillar_deflection(true);
    }
    
    @FXML private void handleMenuPlotPillarDeflectionDrift(ActionEvent event){
        //int c = int_converter.fromString(textfield_pillar_no.getText()) - 1;
        plot_pillar_deflection(false);
    }
    
    
    @FXML private void handleMenuPlotPillarDXY(ActionEvent event){
        int c = int_converter.fromString(textfield_pillar_no_grid.getText()) - 1;
        plot_current_pillar_deflections(c, true);
    }
    
    @FXML private void handleMenuPlotPillarDeflection(ActionEvent event){
        int c = int_converter.fromString(textfield_pillar_no_grid.getText()) - 1;
        plot_current_pillar_deflections(c, false);
    }
    
    @FXML private void handleButtonGetDataPillar(ActionEvent event){
        int c = int_converter.fromString(textfield_pillar_no_grid.getText()) - 1;
        if(c>=grid_reader.cens) return;
        //double[][] xy = grid_reader.getFrameXY(shown_frame_index);
        ResultsTable rt = new ResultsTable();        
        rt.setNaNEmptyCells(true);
        rt.showRowNumbers(false);
        //rt.incrementCounter();
        for(int f=0; f<grid_reader.frames; f++){
            rt.incrementCounter();
            if(grid_reader.active[c][f]){                
                rt.addValue("Frame Index", f+1);
                rt.addValue("X", grid_reader.trackX[c][f]);
                rt.addValue("Y", grid_reader.trackY[c][f]);
                double dx = grid_reader.DX[c][f];
                double dy = grid_reader.DY[c][f];
                rt.addValue("DX", dx);
                rt.addValue("DY", dy);  
                rt.addValue("DIS", Math.sqrt(dx*dx+dy*dy));               
            }
        }
        rt.show("Data for Pillar=" + (c+1));
    }
    
    @FXML private void handleButtonGetDataPillarDrift(ActionEvent event){
        int c = int_converter.fromString(textfield_pillar_no_drift.getText()) - 1;
        if(c>=reader.npillars) return;
        //double[][] xy = grid_reader.getFrameXY(shown_frame_index);
        ResultsTable rt = new ResultsTable();        
        rt.setNaNEmptyCells(true);
        rt.showRowNumbers(false);
        //rt.incrementCounter();
        boolean having_dxy = reader.DX!=null;
        for(int f=0; f<reader.nframes; f++){
            rt.incrementCounter();
            rt.addValue("Frame Index", f+1);
            rt.addValue("Raw X", reader.raw_tracksX[c][f]);
            rt.addValue("Raw Y", reader.raw_tracksY[c][f]);
            if(reader.correted_data){                
                double x = reader.correct_tracksX[c][f];
                double y = reader.correct_tracksY[c][f];
                rt.addValue("X", x);
                rt.addValue("Y", y);  
            }
            
            if(having_dxy){
                double dx = reader.DX[c][f];
                double dy = reader.DY[c][f];
                rt.addValue("DX", dx);
                rt.addValue("DY", dy);  
                rt.addValue("DIS", Math.sqrt(dx*dx+dy*dy)); 
            }
            else if(reader.correted_data){
                double[] xy = reader.getDeflections(c,f,0);  
                rt.addValue("DX",  xy[0]);
                rt.addValue("DY",  xy[1]);  
                rt.addValue("DIS", xy[2]);  
                //rt.addValue("DIS", Math.sqrt(dx*dx+dy*dy));               
            }
        }
        rt.show("Data(Drfit Analysis)for Pillar=" + (c+1));
    }
    
    private GridDeflection[] recreate_grid_frame(int f, double deflection_threshold){
        if(f<0 || f>=grid_reader.frames) return null;
        boolean has_labeled = (grid_reader.IX != null);  
        //int f = current_frame_index;
        double lattice = grid_reader.lattice;
        double threshold = deflection_threshold/lattice;
        double[][] xy = grid_reader.getFrameXY(f);                    
        GridCreation grid_creator = null;
        if(has_labeled){
            int[][] ixy = grid_reader.getFrameIXY(f);
            double[][] dxy = grid_reader.getFrameDXY(f);
            boolean[] dxy_active = grid_reader.getFrameActive(f);
            boolean[] label_flag = grid_reader.getFrameActive(f);
            grid_creator = new GridCreation(xy[0], xy[1], ixy[0], ixy[1],dxy[0],dxy[1], dxy_active, label_flag);
            //grid_creator.label_grid(spacing, oblique, grid_angle);
            double spacing = grid_reader.spacing;
            if(spacing<0) spacing = 2;
            grid_creator.create_grid_points(spacing/2.0, threshold);
        }
        else{                        
            grid_creator = new GridCreation(xy[0], xy[1]);
            double spacing = grid_reader.spacing;                                                
            grid_creator.label_grid(spacing, grid_reader.grid_oblique, grid_reader.grid_angle);                        
            grid_creator.create_grid_points(spacing/2.0, threshold);
        }                    
        double[] dx = grid_creator.getDX();
        double[] dy = grid_creator.getDY();
        boolean[] active = grid_creator.getActive();
        int npillars = active.length;
        GridDeflection[] defletions = new GridDeflection[npillars];
        double[] xc = xy[0];
        double[] yc = xy[1];
        //updateDeflectionMap(xy[0],xy[1],dx,dy,active);        
        for(int i=0; i<npillars; i++){
            defletions[i] = new GridDeflection(xc[i],yc[i],dx[i],dy[i],active[i]);
            dx[i] *= lattice;
            dy[i] *= lattice;
        }
        grid_reader.updateFrameDXY(f,dx, dy,active);
        return defletions;
    }
    
    private void recreate_grid(){
        boolean has_labeled = (grid_reader.IX != null);   
        double threshold_deflection = threshold_deflection_dxy;
        if(Double.isNaN(threshold_deflection_dxy)){
            String text = textfield_large_threshold.getText();
            if(!text.isEmpty()) threshold_deflection =  double_converter.fromString(text);
        }
        GridProperty gp = new GridProperty(grid_reader.spacing, grid_reader.grid_oblique, grid_reader.grid_angle);
        boolean set_grid = !has_labeled;
        if(set_grid){
            //estimate the current separtion and oblique angle for the square grid.            
            double[] sa = settingGrid();    
            if(sa!=null){
                int len = sa.length;
                if(len>0){
                    gp.spacing = sa[0];
                    if(len>1){
                        gp.oblique = sa[1];
                        if(len>2) gp.grid_angle = sa[2];
                    }
                }
            }
        }
        ReCreate_Grid_Plugin grid_recreator = new ReCreate_Grid_Plugin();
        grid_recreator.threshold = threshold_deflection;
        grid_recreator.grid_property = gp;
        grid_recreator.grid_fname = grid_fname;        
        if(grid_recreator.showDialog(set_grid)){
            if(set_grid){
                grid_reader.spacing = grid_recreator.grid_property.spacing;
                grid_reader.grid_oblique = grid_recreator.grid_property.oblique;
                grid_reader.grid_angle = grid_recreator.grid_property.grid_angle;
            }
            threshold_deflection_dxy = grid_recreator.threshold;
            
            if(!grid_recreator.process_all_frames){
                int frame = int_converter.fromString(textfield_frame_no_grid.getText()) - 1;
                if(frame<0 || frame>=grid_reader.frames) return;
                current_frame_index = frame;   
                recreate_grid_frame(current_frame_index, threshold_deflection_dxy);
            }
            else{
                int frames = grid_reader.frames;
                for(int f=0; f<frames; f++){
                    recreate_grid_frame(f, threshold_deflection_dxy);
                    IJ.showProgress(f, frames);
                }         
                //save the data to file
                String fname = grid_recreator.grid_fname;    
                GridFileHeader header = new GridFileHeader(grid_reader.lattice, grid_reader.diameter, 
                        grid_reader.spacing, grid_reader.grid_oblique, grid_reader.grid_angle, 
                        threshold_deflection_dxy);
                Grid_Creation_Plugin.writeDXY_version2(fname,header,
                        grid_reader.trackX,grid_reader.trackY,
                        grid_reader.CX,grid_reader.CY,
                        grid_reader.DX,grid_reader.DY,
                        grid_reader.IX,grid_reader.IY
                );
                grid_fname = fname;
                this.textfile_grid_fname.setText(grid_fname);
            }
            
            updateDeflectionMap();
        }
    }
    
    @FXML private void handleButtonReCreateGrid(ActionEvent event){        
        Thread worker = new Thread() {
                public void run() {		            	            
                         try {
                             recreate_grid();
                         } catch (Exception ex) { }
                }
        };

        worker.start(); // So we don't hold up the dispatch thread.
    }
    
    @FXML private void handleButtonCreateTiffMovie(ActionEvent event){                
        if(overlay_deflections_ip==null || !overlay_deflections_ip.isVisible()) return;
        Thread worker = new Thread() {
            public void run() {		            	            
                 try {
                   // Something that takes a long time . . . in real life,                   
                   int nframes = grid_reader.frames;
                   GenericDialog gd = new GenericDialog("Create Movie");   
                   gd.addNumericField(" start frame:", 1, 0, 10, ">0");
                   gd.addNumericField(" end frame:", nframes, 0, 10, "<"+(nframes+1));
                   gd.showDialog();
                   if(gd.wasCanceled()) return;
                   int start_frame = (int)gd.getNextNumber();
                   int end_frame = (int)gd.getNextNumber();
                   if(start_frame<1) start_frame = 1;
                   if(end_frame>nframes) end_frame = nframes;
                   
                   if(end_frame<start_frame){
                       IJ.showMessage("end frame must be larger than start frame");
                       return;
                   }
                   
                   int deflection_image_id = overlay_deflections_ip.getID();
                   IJ.selectWindow(deflection_image_id);                   
                   //ScreenGrabber sg = new ScreenGrabber();
                   //ImagePlus ip = sg.captureImage();                   
                   ImagePlus ip = MyCaptureImage.captureImage(overlay_deflections_ip, 0);
                   ImageStack stack = new ImageStack(ip.getWidth(), ip.getHeight());
                   int num_frames = end_frame - start_frame + 1;
                   for(int i=start_frame; i<=end_frame; i++){
                       //textfield_frame_no_grid.setText(""+i);                       
                       slider_frame_no_grid.adjustValue(i);      
                       if(i>0 && i<=nframes) updateDeflectionMap();
                       ip = MyCaptureImage.captureImage(overlay_deflections_ip, 10); 
                       stack.addSlice("f="+i, ip.getProcessor());
                       IJ.showProgress(i-start_frame, num_frames);
                   }
                   ip = new ImagePlus("Movie_"+ip.getTitle(), stack);
                   ip.show();                   
                 } catch (Exception ex) {
                 }
               }
             };
        worker.start(); // So we don't hold up the dispatch thread.        
    }
    
    private void draw_frame_tracks(int frame){
        if(frame>0 && frame<=reader.nframes){
            int frame_window = int_converter.fromString(textfield_frame_window.getText());
            if(frame_window>0){
               int end_f = frame-1;                    
               if(end_f<overlay_tracks_ip.getNSlices()){
                   overlay_tracks_ip.setSlice(end_f+1);
                   overlay_tracks_ip.updateAndRepaintWindow();
               }
               else if(end_f<overlay_tracks_ip.getNFrames()){
                   overlay_tracks_ip.setT(end_f+1);
                   overlay_tracks_ip.updateAndRepaintWindow();
               }

               int start_f = end_f - frame_window;                   
               draw_list_pillars(selected_pillars_indexs,start_f,end_f);
               overlay_tracks_ip.updateAndRepaintWindow();
               //IJ.wait(10);
           }
        }
    }
    
    @FXML private void handleButtonCreateTiffMovieDrift(ActionEvent event){                
        if(overlay_tracks_ip==null || !overlay_tracks_ip.isVisible()) return;
        Thread worker = new Thread() {
            public void run() {		            	            
                 try {
                   // Something that takes a long time . . . in real life,                   
                   int nframes = reader.nframes;
                   GenericDialog gd = new GenericDialog("Create Movie");   
                   gd.addNumericField(" start frame:", 1, 0, 10, ">0");
                   gd.addNumericField(" end frame:", nframes, 0, 10, "<"+(nframes+1));
                   gd.showDialog();
                   if(gd.wasCanceled()) return;
                   int start_frame = (int)gd.getNextNumber();
                   int end_frame = (int)gd.getNextNumber();
                   if(start_frame<1) start_frame = 1;
                   if(end_frame>nframes) end_frame = nframes;
                   
                   if(end_frame<start_frame){
                       IJ.showMessage("end frame must be larger than start frame");
                       return;
                   }
                   
                   int deflection_image_id = overlay_tracks_ip.getID();
                   IJ.selectWindow(deflection_image_id);                   
                   //ScreenGrabber sg = new ScreenGrabber();
                   ImagePlus ip = MyCaptureImage.captureImage(overlay_tracks_ip, 0);                   
                   ImageStack stack = new ImageStack(ip.getWidth(), ip.getHeight());
                   int num_frames = end_frame - start_frame + 1;
                   for(int i=start_frame; i<=end_frame; i++){                       
                       //textfield_frame_no_drift.setText(""+i);
                       draw_frame_tracks(i);
                       ip = MyCaptureImage.captureImage(overlay_tracks_ip, 100); 
                       stack.addSlice("f="+i, ip.getProcessor());
                       IJ.showProgress(i-start_frame, num_frames);
                   }
                   ip = new ImagePlus("Movie_"+ip.getTitle(), stack);
                   ip.show();                   
                 } catch (Exception ex) {
                 }
               }
             };
        worker.start(); // So we don't hold up the dispatch thread.        
    }
 
    // </editor-fold>
    
    private void load_settings(String filePath){
        DriftDataLoader drift_loader = new DriftDataLoader();
        if(drift_loader.load_settings(filePath)){
            FileHeaderDrift header = drift_loader.file_header;                                
            textfiled_pixel_size.setText(Double.toString(header.lattice));
            textfiled_spacing.setText(Double.toString(header.spacing));
            textfiled_grid_oblique.setText(Double.toString(header.oblique));
            checkbox_square_grid.setSelected(header.grid_angle==90);
            textfiled_grid_angle.setText(Double.toString(header.grid_angle));
            textfiled_gauss_sigma.setText(Double.toString(header.sigma_PSF));
            textfiled_diameter.setText(Double.toString(header.diameter));            
            checkbox_dark_pillar.setSelected(header.dark_pillars);
            localization_algorithm = header.use_metric_CG ? pillar_tracking.localization_algorithm_CG : pillar_tracking.localization_algorithm_Levmar;
            choicebox_localization_algorithm.getSelectionModel().select(localization_algorithm);
            textfiled_searchwindow.setText(Integer.toString(header.kernel_w));
            textfiled_constrain_radius.setText(Integer.toString(header.box_constrian_R));
            textfiled_max_drift.setText(Integer.toString((int)header.catch_radius));
            checkbox_apply_enhancement.setSelected(false);
            checkbox_apply_pillar_recontruction.setSelected(false);
            meun_apply_rank_filter.setSelected(header.apply_mean_rank_filter);
            if(header.use_enhancer){
                image_enhancer_plugin = drift_loader.enhancer;
                if(image_enhancer_plugin.isGPUOn()){
                    if(IJ.isWindows()){
                        boolean suc = image_enhancer_plugin.gpu_lib.prepare();
                        if(!suc) image_enhancer_plugin.setGPUOn(false);
                    }
                    else image_enhancer_plugin.setGPUOn(false);
                }
            }
            else image_enhancer_plugin = null;
            
            if(header.use_fft){
                if(fft==null) fft = new HighPassFFT();
                fft.setOffCenterRadius(drift_loader.mask_radius);
                fft.step_radius_off_center = drift_loader.start_radius;
                fft.end_radius_off_center = drift_loader.end_radius;                
                int[][] points = drift_loader.fft_points;
                fft.setOffCenterPoints(points[0], points[1]);
            }                
            else fft = null;
        }
        else{
            GridDataLoader grid_loader = new GridDataLoader();
            if(grid_loader.load_settings(filePath)){
                textfiled_pixel_size.setText(Double.toString(grid_loader.lattice));
                textfiled_spacing.setText(Double.toString(grid_loader.spacing));
                textfiled_grid_oblique.setText(Double.toString(grid_loader.grid_oblique));
                checkbox_square_grid.setSelected(grid_loader.grid_angle==90);
                textfiled_grid_angle.setText(Double.toString(grid_loader.grid_angle));
                //textfiled_gauss_sigma.setText(Double.toString(grid_loader.sigma_PSF));
                textfiled_diameter.setText(Double.toString(grid_loader.diameter));
                checkbox_square_grid.setSelected(grid_loader.grid_angle==90);
                //checkbox_dark_pillar.setSelected(grid_loader.dark_pillars);
            }
        }
    }
    
    @Override
    public void initialize(URL location, ResourceBundle resources) {  
        
        textfiled_spacing.setTextFormatter(new TextFormatter(double_converter));
        textfiled_grid_oblique.setTextFormatter(new TextFormatter(double_converter));
        textfiled_diameter.setTextFormatter(new TextFormatter(double_converter));
        textfiled_gauss_sigma.setTextFormatter(new TextFormatter(double_converter));
        textfiled_pixel_size.setTextFormatter(new TextFormatter(double_converter));
        textfiled_unit_nm.setTextFormatter(new TextFormatter(double_converter));
        textfiled_unit_pixel.setTextFormatter(new TextFormatter(double_converter));
        
        textfiled_pixel_size.textProperty().addListener((observable, oldValue, newValue) -> {                       
            //if(newValue.isEmpty()) return;
           try{
                settings.pixel_size = newValue.isEmpty()? Double.NaN : double_converter.fromString(newValue);

                if(this.toggle_converter_forward.isSelected()){
                    String text = textfiled_unit_nm.getText();
                    if(!text.isEmpty()) {
                     double nm = double_converter.fromString(text);
                     textfiled_unit_pixel.setText(String.format("%.3f", settings.convertNanoToPixel(nm)));
                    }
                }
                else if(this.toggle_converter_back.isSelected()){
                    String text = textfiled_unit_pixel.getText();
                    if(!text.isEmpty()) {
                     double pixel = double_converter.fromString(text);
                     textfiled_unit_nm.setText(String.format("%.3f", settings.convertPixelToNano(pixel)));
                    }
                }
           }
           catch(Exception ex){}
        });
        
        textfiled_unit_nm.textProperty().addListener((observable, oldValue, newValue) -> {                       
           if(newValue.isEmpty()) return;         
           if(this.toggle_converter_forward.isSelected()){
               try{ 
                    double nm = double_converter.fromString(newValue);
                    textfiled_unit_pixel.setText(String.format("%.3f", settings.convertNanoToPixel(nm)));
               }
               catch(Exception ex){}
           }           
        });
        
        textfiled_unit_pixel.textProperty().addListener((observable, oldValue, newValue) -> {                       
           if(newValue.isEmpty()) return;         
           if(this.toggle_converter_back.isSelected()){              
               try{
                    double pixel = double_converter.fromString(newValue);
                    textfiled_unit_nm.setText(String.format("%.3f", settings.convertPixelToNano(pixel)));
               }
               catch(Exception ex){}
           }
        });        
        
        //tracking
        textfiled_searchwindow.setTextFormatter(new TextFormatter(int_converter));
        textfiled_constrain_radius.setTextFormatter(new TextFormatter(int_converter));
        textfiled_max_drift.setTextFormatter(new TextFormatter(int_converter));
        textfiled_num_threads.setTextFormatter(new TextFormatter(int_converter));
        textfiled_num_threads.setText(int_converter.toString(ThreadUtil.getNbCpus()));
        this.checkmenu_show_fft.visibleProperty().bind(this.checkbox_apply_pillar_recontruction.selectedProperty());
        this.checkmenu_show_recontruction.visibleProperty().bind(this.checkbox_apply_pillar_recontruction.selectedProperty());
        
        this.checkmenu_drift_correction.visibleProperty().bind(this.checkbox_apply_pillar_recontruction.selectedProperty().not());
        this.checkmenu_create_grid_poly_fit.visibleProperty().bind(this.checkbox_apply_pillar_recontruction.selectedProperty().not());
        
        this.combo_images.setOnMouseClicked(new EventHandler<MouseEvent>() {
            @Override
            public void handle(MouseEvent mouseEvent) {
                if(mouseEvent.getButton().equals(MouseButton.PRIMARY)){
                    if(mouseEvent.getClickCount() > 0){
                        updateComboImages();
                    }
                }
            }
        });
        
        
        
        this.combo_images.getSelectionModel().selectedIndexProperty().addListener((observable, oldValue, newValue) -> {   
            int n = WindowManager.getImageCount();
            if(n<1) return;
            
            if(list_image_ids==null) return;            
            int index = newValue.intValue();
            if(index<0 || index>=list_image_ids.size()) return;                

            int id = list_image_ids.get(index);
            ImagePlus ip = WindowManager.getImage(id);
            if(ip==null || !ip.isVisible()) return;   
            FileInfo fi = ip.getOriginalFileInfo();
            if(fi != null && fi.directory != null && fi.fileName!=null){
                String fname = ip.getShortTitle();
                if(fname==null) fname = fi.fileName;
                String info = fi.directory + fname;
                //String fileNameWithOutExt = FilenameUtils.removeExtension(info);
                if(localization_algorithm == pillar_tracking.localization_algorithm_CG)
                    this.textfiled_output_fname.setText(info+".cent");            
                else this.textfiled_output_fname.setText(info+".bin");            
            }
            else this.textfiled_output_fname.setText(""); 
        });
        
        this.choicebox_localization_algorithm.setItems(FXCollections.observableArrayList("Center of Mass (CM)", "Gaussian Fit (GF)"));
        this.choicebox_localization_algorithm.getSelectionModel().select(this.localization_algorithm);
        this.choicebox_localization_algorithm.getSelectionModel().selectedIndexProperty().addListener(new 
            ChangeListener<Number>(){
                public void changed(ObservableValue ov, Number value, Number new_value){
                    int selected = new_value.intValue();
                    if(selected>=0){
                        localization_algorithm = selected;
                        String ouput_fname = textfiled_output_fname.getText();
                        if(ouput_fname!=null && !ouput_fname.isEmpty()){
                            if(ouput_fname.endsWith(".bin") && localization_algorithm == pillar_tracking.localization_algorithm_CG){
                                textfiled_output_fname.setText(ouput_fname.replaceAll(".bin", ".cent"));
                            }
                            else if(ouput_fname.endsWith(".cent") && localization_algorithm != pillar_tracking.localization_algorithm_CG){
                                textfiled_output_fname.setText(ouput_fname.replaceAll(".cent", ".bin"));
                            }
                        }                        
                    }
            }            
        });
        
        checkbox_pixel_unit.selectedProperty().addListener(new ChangeListener<Boolean>() {
           @Override
           public void changed(ObservableValue<? extends Boolean> observable, Boolean oldValue, Boolean newValue) {
               if(newValue){
                    if(WindowManager.getImageCount()<1){
                        checkbox_pixel_unit.setSelected(false);
                        pixel_unit_change(false);
                    }
                    else pixel_unit_change(true);
               }
               else{
                   pixel_unit_change(false);
               }
           }
       });
        
       checkbox_apply_enhancement.selectedProperty().addListener(new ChangeListener<Boolean>() {
           @Override
           public void changed(ObservableValue<? extends Boolean> observable, Boolean oldValue, Boolean newValue) {
               if(newValue){
                   if(image_enhancer_plugin==null) setKernelforEnhancement(true);                   
               }
           }
       });
       
       this.checkbox_apply_pillar_recontruction.selectedProperty().addListener(new ChangeListener<Boolean>() {
           @Override
           public void changed(ObservableValue<? extends Boolean> observable, Boolean oldValue, Boolean newValue) {
               if(newValue){
                   if(fft==null || !fft.hasOffCenterPoints()){
                       IJ_showMessage("The mask has been not set yet!");
                       checkbox_apply_pillar_recontruction.selectedProperty().set(false);
                   }
                   else setMaskFilterRadius();
               
               }
           }
       });
       
       checkbox_show_current_kernel.selectedProperty().addListener(new ChangeListener<Boolean>() {
           @Override
           public void changed(ObservableValue<? extends Boolean> observable, Boolean oldValue, Boolean newValue) {
               if(newValue){
                   if(image_enhancer_plugin!=null){        
                        ImagePlus current_psf_ip = image_enhancer_plugin.getPSF();                    
                        current_psf_ip.show();                                                     
                   }
                   else{
                       checkbox_show_current_kernel.selectedProperty().set(false);
                   }
               }
           }
       });
       
               
        checkbox_shownin.selectedProperty().addListener(new ChangeListener<Boolean>() {
           @Override
           public void changed(ObservableValue<? extends Boolean> observable, Boolean oldValue, Boolean newValue) {
               if(newValue){
                    if(WindowManager.getImageCount()<1){
                        checkbox_shownin.setSelected(false);
                    }
               }
               else{
                   //need_recreate_map_drift = true;
               }
           }
       });
        
        checkbox_show_track_head.selectedProperty().addListener(new ChangeListener<Boolean>() {
           @Override
           public void changed(ObservableValue<? extends Boolean> observable, Boolean oldValue, Boolean newValue) {
               DriftAnalysisPlotter.SHOW_TRACK_HEADER = newValue;
               handleButtonTracksMap(null);
           }
       });
         checkbox_show_track_head.setSelected(DriftAnalysisPlotter.SHOW_TRACK_HEADER);
        
       checkbox_show_labels.selectedProperty().addListener(new ChangeListener<Boolean>() {
           @Override
           public void changed(ObservableValue<? extends Boolean> observable, Boolean oldValue, Boolean newValue) {
               updateCheckboxShowLabels();
           }
       });
       
       checkbox_still_pillars.selectedProperty().addListener(new ChangeListener<Boolean>() {
           @Override
           public void changed(ObservableValue<? extends Boolean> observable, Boolean oldValue, Boolean newValue) {
               handleButtonTracksMap(null);
           }
       });
       
        label_drift_load_info.setText("");
        label_grid_load_info.setText("");
       
       slider_zoomin_images.disableProperty().bind(checkbox_flatten.selectedProperty().not());
       slider_zoomin_tracks.valueProperty().addListener(new ChangeListener<Number>() {
           @Override
           public void changed(ObservableValue<? extends Number> observable, Number oldValue, Number newValue) {
               int setvalue = newValue.intValue();
               slider_zoomin_tracks.setValue(setvalue);
               if(Math.abs(setvalue-oldValue.intValue())>=1.0){
                   zoom_in_tracks = setvalue;
                   textfield_zoomin_tracks.setText(int_converter.toString(setvalue));
                   updateSliderZoominTracks();
               }
           }
       });
       
       slider_zoomin_images.valueProperty().addListener(new ChangeListener<Number>() {
           @Override
           public void changed(ObservableValue<? extends Number> observable, Number oldValue, Number newValue) {
               int setvalue = newValue.intValue();
               slider_zoomin_images.setValue(setvalue);
               if(Math.abs(setvalue-oldValue.intValue())>=1.0){
                   zoom_in_image = setvalue;
                   textfield_zoomin_images.setText(int_converter.toString(setvalue));
               }
           }
       });
       
       slider_pillar_no_drift.valueProperty().addListener(new ChangeListener<Number>() {
           @Override
           public void changed(ObservableValue<? extends Number> observable, Number oldValue, Number newValue) {
               int setvalue = newValue.intValue();
               slider_pillar_no_drift.setValue(setvalue);
               if(Math.abs(setvalue-oldValue.intValue())>=1.0){                   
                   textfield_pillar_no_drift.setText(int_converter.toString(setvalue));
               }
           }
       });
       
       slider_frame_no_drift.valueProperty().addListener(new ChangeListener<Number>() {
           @Override
           public void changed(ObservableValue<? extends Number> observable, Number oldValue, Number newValue) {
               int setvalue = newValue.intValue();
               slider_frame_no_drift.setValue(setvalue);
               if(Math.abs(setvalue-oldValue.intValue())>=1.0){                   
                   textfield_frame_no_drift.setText(int_converter.toString(setvalue));                   
               }
           }
       });
       
       textfield_frame_window.setTextFormatter(new TextFormatter(int_converter));
       textfield_frame_window.setText("10");
       textfield_frame_no_drift.setTextFormatter(new TextFormatter(int_converter));
       textfield_frame_no_drift.textProperty().addListener((observable, oldValue, newValue) -> {                       
           if(newValue==null || newValue.isEmpty()) return;
           try{
                int num = int_converter.fromString(newValue);
                int slide = (int)slider_frame_no_drift.getValue();
                int max = (int)slider_frame_no_drift.getMax();
                if(num!=slide){
                    slider_frame_no_drift.setValue(num);                    
                }
                if(num>0 && num<=max){
                     int frame_window = int_converter.fromString(textfield_frame_window.getText());
                     if(frame_window>0){
                        int end_f = num-1;                    
                        if(end_f<overlay_tracks_ip.getNSlices()){
                            overlay_tracks_ip.setSlice(end_f+1);
                            overlay_tracks_ip.updateAndRepaintWindow();
                        }
                        else if(end_f<overlay_tracks_ip.getNFrames()){
                            overlay_tracks_ip.setT(end_f+1);
                            overlay_tracks_ip.updateAndRepaintWindow();
                        }

                        int start_f = end_f - frame_window;                   
                        draw_list_pillars(selected_pillars_indexs,start_f,end_f);
                    }
                }
                else{
                    if(num<=0) textfield_frame_no_drift.setText(int_converter.toString(1));
                    else if(num>max) textfield_frame_no_drift.setText(int_converter.toString(max));
                }
           }
           catch(Exception ex){}
        });
       
       textfield_frame_no_drift.setText("1");
       
       textfield_pillar_no_drift.setTextFormatter(new TextFormatter(int_converter));
       textfield_pillar_no_drift.textProperty().addListener((observable, oldValue, newValue) -> {                       
           if(newValue==null || newValue.isEmpty()) return;
           try{
                int num = int_converter.fromString(newValue);
                int slide = (int)slider_pillar_no_drift.getValue();
                int max = (int)slider_pillar_no_drift.getMax();
                if(num!=slide) slider_pillar_no_drift.setValue(num);  
                if(num>0 && num<=reader.npillars){
                    selected_pillars_indexs = new int[]{num-1};
                    updateDriftPlotWindow();
                    updateSliderZoominTracks();
                }
                else{
                    if(num<=0) textfield_pillar_no_drift.setText(int_converter.toString(1));
                    else if(num>max) textfield_pillar_no_drift.setText(int_converter.toString(max));
                }
           }
           catch(Exception ex){}
        });
       textfield_pillar_no_drift.setText("1");
       
       //grid
        checkbox_shownin_grid.selectedProperty().addListener(new ChangeListener<Boolean>() {
           @Override
           public void changed(ObservableValue<? extends Boolean> observable, Boolean oldValue, Boolean newValue) {
               if(newValue){
                    if(WindowManager.getImageCount()<1){
                        checkbox_shownin_grid.setSelected(false);
                    }
                    
               }
               else{
                   //need_recreate_map_drift = true;
               }
           }
       });
       
       slider_zoomin_deflection.valueProperty().addListener(new ChangeListener<Number>() {
           @Override
           public void changed(ObservableValue<? extends Number> observable, Number oldValue, Number newValue) {
               int setvalue = newValue.intValue();
               slider_zoomin_deflection.setValue(setvalue);
               if(Math.abs(setvalue-oldValue.intValue())>=1.0){
                   zoom_in_deflection = setvalue;
                   textfield_zoomin_deflection.setText(int_converter.toString(setvalue));
                   updateDeflectionMap();
               }
           }
       });
       
       slider_zoomin_image_grid.valueProperty().addListener(new ChangeListener<Number>() {
           @Override
           public void changed(ObservableValue<? extends Number> observable, Number oldValue, Number newValue) {
               int setvalue = newValue.intValue();
               slider_zoomin_image_grid.setValue(setvalue);
               if(Math.abs(setvalue-oldValue.intValue())>=1.0){
                   //zoom_in_image_grid = setvalue;
                   textfield_zoomin_image_grid.setText(int_converter.toString(setvalue));
                   need_recreate_map_grid = true;
               }
           }
       });
       
       slider_frame_no_grid.valueProperty().addListener(new ChangeListener<Number>() {
           @Override
           public void changed(ObservableValue<? extends Number> observable, Number oldValue, Number newValue) {
               int setvalue = newValue.intValue();
               slider_frame_no_grid.setValue(setvalue);
               if(Math.abs(setvalue-oldValue.intValue())>=1.0){                   
                   textfield_frame_no_grid.setText(int_converter.toString(setvalue));
               }
           }
       });
       
       slider_pillar_no_grid.valueProperty().addListener(new ChangeListener<Number>() {
           @Override
           public void changed(ObservableValue<? extends Number> observable, Number oldValue, Number newValue) {
               int setvalue = newValue.intValue();
               slider_pillar_no_grid.setValue(setvalue);
               if(Math.abs(setvalue-oldValue.intValue())>=1.0){                   
                   textfield_pillar_no_grid.setText(int_converter.toString(setvalue));
//                   updateDeflectionPlot();
               }
           }
       });
       
       checkbox_show_labels_grid.selectedProperty().addListener(new ChangeListener<Boolean>() {
           @Override
           public void changed(ObservableValue<? extends Boolean> observable, Boolean oldValue, Boolean newValue) {
               updateDeflectionMap();
           }
       });
       
       checkbox_show_labels_grid.indeterminateProperty().addListener(new ChangeListener<Boolean>() {
           @Override
           public void changed(ObservableValue<? extends Boolean> observable, Boolean oldValue, Boolean newValue) {
               updateDeflectionMap();
           }
       });
       
       
       checkbox_remove_large.selectedProperty().addListener(new ChangeListener<Boolean>() {
           @Override
           public void changed(ObservableValue<? extends Boolean> observable, Boolean oldValue, Boolean newValue) {
               updateDeflectionMap();               
           }
       });      
       
       accordion_grid.setDisable(true);
       accordion_drift.setDisable(true);
       accordion_grid.setExpanded(false);
       accordion_drift.setExpanded(false);
       
       textfield_large_threshold.setTextFormatter(new TextFormatter(double_converter));
       textfield_large_threshold.disableProperty().bind(checkbox_remove_large.selectedProperty().not());
       textfield_large_threshold.textProperty().addListener((observable, oldValue, newValue) -> {                       
           if(newValue!=null && !newValue.isEmpty()){
               try{ 
                    double threshold = double_converter.fromString(newValue);
                    if(threshold>0) updateDeflectionMap();
               }
               catch(Exception ex){}
           }
        });
       
       textfield_frame_no_grid.setTextFormatter(new TextFormatter(int_converter));
       textfield_frame_no_grid.textProperty().addListener((observable, oldValue, newValue) -> {                       
           if(newValue==null || newValue.isEmpty()) return;
           try{ 
                int num = int_converter.fromString(newValue);
                int slide = (int)slider_frame_no_grid.getValue();
                int max = (int)slider_frame_no_grid.getMax();
                if(num!=slide) slider_frame_no_grid.setValue(num);  
                if(num>0 && num<=max) updateDeflectionMap();
                else{
                    if(num<1) textfield_frame_no_grid.setText(int_converter.toString(1));
                    else if(num>max) textfield_frame_no_grid.setText(int_converter.toString(max));
                }
           }
           catch(Exception ex){}
        });
       textfield_frame_no_grid.setText("1");
       
       textfield_pillar_no_grid.setTextFormatter(new TextFormatter(int_converter));
       textfield_pillar_no_grid.textProperty().addListener((observable, oldValue, newValue) -> {                       
           if(newValue==null || newValue.isEmpty()) return;
           try{
                int num = int_converter.fromString(newValue);
                //slider_pillar_no_grid.setValue(num);
                int slide = (int)slider_pillar_no_grid.getValue();
                int max = (int)slider_pillar_no_grid.getMax();
                if(num!=slide) slider_pillar_no_grid.setValue(num);
                if(num>0 && num<=max) updateDeflectionPlot();
                else{
                    if(num<1) textfield_pillar_no_grid.setText(int_converter.toString(1));
                    else if(num>max) textfield_pillar_no_grid.setText(int_converter.toString(max));
                }
           }
           catch(Exception ex){}
        });
       textfield_pillar_no_grid.setText("1");       
       
       textfield_zoomin_image_grid.setText(int_converter.toString((int)zoom_in_image_grid));
       textfield_zoomin_deflection.setText(int_converter.toString((int)zoom_in_deflection));       
       
       choice_crosshair.setItems(FXCollections.observableArrayList(cross_size));
       choice_arrows.setItems(FXCollections.observableArrayList(arrow_size));
       choice_labels.setItems(FXCollections.observableArrayList(label_size));
       choice_crosshair.getSelectionModel().selectedIndexProperty().addListener((observable, oldValue, newValue) -> {                       
           grid_plotter.cross_size = cross_size.get(newValue.intValue());
           updateDeflectionMap();
        });
       
       choice_arrows.getSelectionModel().selectedIndexProperty().addListener((observable, oldValue, newValue) -> {                       
           grid_plotter.arrow_head_size = arrow_size.get(newValue.intValue());
           updateDeflectionMap();
        });
       
       choice_labels.getSelectionModel().selectedIndexProperty().addListener((observable, oldValue, newValue) -> {                       
           grid_reader.label_size = label_size.get(newValue.intValue());
           updateDeflectionMap();
        });       
       choice_crosshair.getSelectionModel().select(1);
       choice_arrows.getSelectionModel().select(2);
       choice_labels.getSelectionModel().select(2);
       
       color_crosshair.valueProperty().addListener((observable, oldValue, newValue) -> {                       
           grid_plotter.cross_color = new Color((float)newValue.getRed(),
                                             (float) newValue.getGreen(),
                                             (float) newValue.getBlue(),
                                             (float) newValue.getOpacity());
           updateDeflectionMap();
        });
       
       color_arrows.valueProperty().addListener((observable, oldValue, newValue) -> {                       
           grid_plotter.arrow_color = new Color((float)newValue.getRed(),
                                             (float) newValue.getGreen(),
                                             (float) newValue.getBlue(),
                                             (float) newValue.getOpacity());
           updateDeflectionMap();
        });
       
       color_labels.valueProperty().addListener((observable, oldValue, newValue) -> {                       
           grid_reader.label_color = new Color((float)newValue.getRed(),
                                             (float) newValue.getGreen(),
                                             (float) newValue.getBlue(),
                                             (float) newValue.getOpacity());
           updateDeflectionMap();
        });
       
       color_crosshair.setValue(javafx.scene.paint.Color.WHITE);
       color_arrows.setValue(javafx.scene.paint.Color.CYAN);
       color_labels.setValue(javafx.scene.paint.Color.YELLOW);
       
       //about
       Image image = new Image(getClass().getResourceAsStream("logo-home.png"));
       link_mbi.setText("http://mbi.nus.edu.sg/");
       link_mbi.setContentDisplay(ContentDisplay.GRAPHIC_ONLY);
       link_mbi.setGraphic(new ImageView (image));
       link_mbi.setOnAction((ActionEvent e)->{
            if (Desktop.isDesktopSupported()) {
                    try {
		        URI uri = new URI(link_mbi.getText());
		        Desktop.getDesktop().browse(uri);
		    } catch (Exception ex) {	       
	    	}
            }     
       });       
       link_mbi.setVisible(true);       
       link_copyright.setText("Source Code on Github under GPL3.0 License");
       link_copyright.setOnAction((ActionEvent e)->{
            if (Desktop.isDesktopSupported()) {
                    try {
		        URI uri = new URI("https://github.com/scottreen/PillarTracker");
		        Desktop.getDesktop().browse(uri);
		    } catch (Exception ex) { }
            }    
       });       
       link_copyright.setVisible(true);
       //this.pane_properties.getChildren().add(link_copyright);
       
      link_version.setVisible(true);
      link_version.setText("Version:" + Pillar_Tracker_Plugin.VERSION);
      link_version.setOnAction((ActionEvent e)->{
           if (Desktop.isDesktopSupported()) {
			try {
		        URI uri = new URI("https://imagej.net/PillarTracker");
		        Desktop.getDesktop().browse(uri);
		    } catch (Exception ex) { }
            }    
       });        
      
      link_update_gpu_lib.setVisible(true);
      link_update_gpu_lib.setOnAction((ActionEvent e)->{
           SwingUtilities.invokeLater(new Runnable() {
                @Override public void run() {                    
                   PillarTrackerGPULibrary lib = new PillarTrackerGPULibrary();
                   lib.updateLibrary();
             }
            });
       });       
      
       this.textfile_drift_fname.setText(drift_fname);
       this.textfile_grid_fname.setText(grid_fname);
       
       javax.swing.filechooser.FileFilter filter = new FileNameExtensionFilter("Binary File|*.bin", "bin");
       drift_fileChooser.setDialogTitle("Open Drift File");       
       drift_fileChooser.setFileFilter(filter);
       javax.swing.filechooser.FileFilter cu_filter = new FileNameExtensionFilter("Contractile Unit File|*.hcu", "hcu");
       drift_fileChooser.addChoosableFileFilter(cu_filter);
       
       fileChooser_track_fname.setDialogTitle("Save Pillar Tracking File(Include Drift information)");       
       fileChooser_track_fname.setFileFilter(filter);
       
       javax.swing.filechooser.FileFilter filter_grid = new FileNameExtensionFilter("Binary File|*.bin, *.dxy", "bin", "dxy");       
       grid_fileChooser.setDialogTitle("Open Grid File");       
       grid_fileChooser.setFileFilter(filter_grid);
       grid_fileChooser.addChoosableFileFilter(cu_filter);

       pane_properties.setOnDragOver(new EventHandler<DragEvent>() {
            @Override
            public void handle(DragEvent event) {
                Dragboard db = event.getDragboard();
                if (db.hasFiles()) {
                    event.acceptTransferModes(TransferMode.LINK);  
                } else {
                    event.consume();
                }
            }
        });
       
       pane_properties.setOnDragDropped(new EventHandler<DragEvent>() {
            @Override 
            public void handle(DragEvent event) {
                Dragboard db = event.getDragboard();                
                boolean success = false;
                if (db.hasFiles()) {
                    //try
                    {
                        List<File> list = db.getFiles();
                        if(list.size()==1){
                            success = true;
                            Path p = list.get(0).toPath();
                            String filePath = p.toString(); 
                            load_settings(filePath);
                        }
                    }
                    //catch(Exception ex){success = true;}
                }
                event.setDropCompleted(success);
                event.consume();
            }
        });
       
       tab_drift.setOnDragOver(new EventHandler<DragEvent>() {
            @Override
            public void handle(DragEvent event) {
                Dragboard db = event.getDragboard();
                if (db.hasFiles()) {
                    event.acceptTransferModes(TransferMode.LINK);  
                } else {
                    event.consume();
                }
            }
        });
       
       
       tab_drift.setOnDragDropped(new EventHandler<DragEvent>() {
            @Override 
            public void handle(DragEvent event) {
                Dragboard db = event.getDragboard();                
                boolean success = false;
                if (db.hasFiles()) {
                    //try
                    {
                        List<File> list = db.getFiles();
                        if(list.size()==1){
                            success = true;
                            Path p = list.get(0).toPath();
                            String filePath = p.toString(); //.getAbsolutePath();
                            textfile_drift_fname.setText(filePath);      
                        }
                    }
                    //catch(Exception ex){success = true;}
                }
                event.setDropCompleted(success);
                event.consume();
            }
        });
       
       tab_grid.setOnDragOver(new EventHandler<DragEvent>() {
            @Override
            public void handle(DragEvent event) {
                Dragboard db = event.getDragboard();
                if (db.hasFiles()) {
                    event.acceptTransferModes(TransferMode.LINK);  
                } else {
                    event.consume();
                }
            }
        });
       
       tab_grid.setOnDragDropped(new EventHandler<DragEvent>() {
            @Override 
            public void handle(DragEvent event) {
                Dragboard db = event.getDragboard();
                boolean success = false;
                if (db.hasFiles() ) {
                    List<File> list = db.getFiles();
                    if(list.size()==1){
                        success = true;
                        Path p = list.get(0).toPath();
                        String filePath = p.toString(); 

                        textfile_grid_fname.setText(filePath);    
                    }  
                }
                event.setDropCompleted(success);
                event.consume();
            }
        });
       
       // square grid??
       textfiled_grid_angle.disableProperty().bind(this.checkbox_square_grid.selectedProperty());
       
       checkbox_square_grid.selectedProperty().addListener(new ChangeListener<Boolean>() {
           @Override
           public void changed(ObservableValue<? extends Boolean> observable, Boolean oldValue, Boolean newValue) {
               String grid_angle_str = newValue ? "90" : "" + settings.grid_angle;
               textfiled_grid_angle.textProperty().set(grid_angle_str);               
           }
       });
       
       this.checkbox_square_grid.selectedProperty().set(true);
       this.button_hunt_cu.disableProperty().bind(this.checkbox_raw.disableProperty());
    }
   
}
