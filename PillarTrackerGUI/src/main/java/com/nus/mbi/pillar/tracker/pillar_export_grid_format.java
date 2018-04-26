package com.nus.mbi.pillar.tracker;

//Xu Xiaochun @ Mechnobiology Institute, Singapore
//contact:mbixxc@nus.edu.sg
//update history:
//2016-10-11:

import ij.*;
import ij.plugin.*;
import java.io.*;
import java.nio.*;
import fiji.util.gui.*; 		//using the GenericDialogPlus;
import ij.io.OpenDialog;

public class pillar_export_grid_format implements PlugIn{
	
//input file->E:\Yangbo\2015-06-25\75fps_cell2_dish2_large-output\\75fps_cell2_dish2_large_X11-pillar-tracks-big.bin
String tracks_file = "E:\\Yangbolix\\2015-06-25\\75fps_cell2_dish2_large-output\\75fps_cell2_dish2_large_X11-pillar-tracks-big.bin";

public final int checksum=444666; // version 1 double precision floats

// image and pillar params
public double pixel_size = 40.4; // pixel size in nm
public double diameter=500.;
public double spacing=1000.;

int npillars = 0;
int nframes = 0;
double[][] raw_tracksX = null;
double[][] raw_tracksY = null;
double[][] correct_tracksX = null;
double[][] correct_tracksY = null;
double[]   driftX = null;
double[]   driftY = null;
boolean correted_data = false;
boolean load_data_suc = false;
boolean debug = true;
boolean little_endian = false; //Only used for the old binary format

public void run(String arg) {
		IJ.freeMemory();
		boolean checkonce = arg.contains("checkonce");
		boolean consistent;	
		do {
			consistent=true;
			if (!showDialog()) return;			
			File f = new File(tracks_file);
			if (!f.isFile()) {IJ.showMessage("the file doesn't exsit"); consistent=false;}
			else load_data_suc = load_tracks();			
		} while (!consistent && !checkonce);
		
		if(!consistent || !load_data_suc) return;
			
		IJ.log("input file->" + tracks_file);
		if(npillars<1) {IJ.log("the file has no pillar"); return;}
		if(nframes<1) {IJ.log("the file has no frame"); return;}
		else if(nframes<2) {IJ.log("the file has only one frame"); return;}		

		String output_fname = tracks_file + ".grid";
		IJ.log("output file->" + output_fname);
		save_grid_format(output_fname);		
}

	public void save_grid_format(String output_fname){
		int frames = nframes;
        int cens = npillars;
        
        try {                    
                FileOutputStream pout;
                DataOutputStream printout;
                pout = new FileOutputStream(output_fname);
                printout = new DataOutputStream(pout);

                // write credentials into the text file
                printout.writeInt(checksum);
                printout.writeDouble(pixel_size);
                printout.writeDouble(diameter);
                printout.writeDouble(spacing);
                printout.writeInt(frames);
                printout.writeInt(cens);
                byte[] active_lin=new byte[cens*frames];
                int j=0;
                for(int y=0; y<cens; y++) for(int x=0; x<frames; x++) {                       
                        if(Double.isNaN(raw_tracksX[y][x]) || Double.isNaN(raw_tracksY[y][x]))
                        	active_lin[j]=0; 
                        else active_lin[j]=1;
                        
                        j++;
                }
                printout.write(active_lin, 0, cens*frames);
                for(int y=0; y<cens; y++) {
                        String s="";
                        for(int x=0; x<frames; x++) {                                
                                if(active_lin[y*frames+x]>0) {
                                        s=s+"("+raw_tracksX[y][x]+","+raw_tracksY[y][x]+")	";
                                        printout.writeDouble(raw_tracksX[y][x]);
                                        printout.writeDouble(raw_tracksY[y][x]);
                                } else {
                                        s=s+"	";
                                }                                
                        }
                        IJ.log(s);
                        //printout.println(sf);
                }
                printout.close();
        }
        catch(FileNotFoundException fe)
         {
                IJ.log("FileNotFoundException : " + fe);
         }
        catch(IOException ioe)
        {
                IJ.log("IOException : " + ioe);
        }		
	}


	public boolean load_tracks(String filename, double pixel_size)
	{		
		int [] info=new int[2];
		double[][] matrix_raw = null;
		double[][] matrix_corrected = null;
		double[][] matrix_driftXY = null;
		try{
			FileInputStream fis = new FileInputStream(filename);			
			info = readfileinfo(fis,little_endian);		
			matrix_raw = readfile2matrix(fis, info[0], info[1], little_endian);
			matrix_corrected = readfile2matrix(fis, info[0], info[1], little_endian);
			matrix_driftXY = readfile2matrix(fis, 2, info[1], little_endian);
			fis.close();		
		}
		catch(Exception e){
			matrix_corrected = null;
			IJ.log("read file failed->" + e.getMessage());
		}
		
		if(matrix_raw==null || info[0]%2 != 0) return false;
		correted_data = (matrix_corrected!=null);

		tracks_file = filename;
		this.pixel_size = pixel_size;
		
		npillars = info[0]/2;
		nframes = info[1];
		raw_tracksX = new double[npillars][nframes];
		raw_tracksY = new double[npillars][nframes];		
		correct_tracksX = new double[npillars][nframes];
		correct_tracksY = new double[npillars][nframes];
		driftX = new double[nframes];
		driftY = new double[nframes];
		
		for(int i=0; i<nframes; i++){
			for(int j=0; j<npillars; j++){
				raw_tracksX[j][i] = pixel_size*matrix_raw[2*j][i];
				raw_tracksY[j][i] = pixel_size*matrix_raw[2*j+1][i];				
				if(correted_data){
					correct_tracksX[j][i] = pixel_size*matrix_corrected[2*j][i];
					correct_tracksY[j][i] = pixel_size*matrix_corrected[2*j+1][i];			
				}
			}

			if(correted_data){
				driftX[i] = pixel_size*matrix_driftXY[0][i];
				driftY[i] = pixel_size*matrix_driftXY[1][i];
			}
		}
		
		IJ.log("check the last values");
		IJ.log("	" + raw_tracksX[npillars-1][nframes-1] + "	 " + raw_tracksY[npillars-1][nframes-1]);
		if(correted_data){
			IJ.log("	" + (correct_tracksX[npillars-1][nframes-1] - raw_tracksX[npillars-1][nframes-1]) + "	 " + (correct_tracksY[npillars-1][nframes-1]-raw_tracksY[npillars-1][nframes-1]));
			IJ.log("	" + driftX[nframes-1] + "	 " + driftY[nframes-1]);
		}	
		
		return true;
	}
	
	public boolean load_tracks()
	{
		return load_tracks(tracks_file, pixel_size);		
	}
	
	public int[] readfileinfo(FileInputStream fis, boolean use_little_endian) throws Exception
	{	
		//double [][] matrix = null;
		//FileInputStream fis = new FileInputStream(tracks_file);
		int[]info = new int[2];
		byte[] header = new byte[8];
		int bytes = fis.read(header);
		if(bytes<=0) return null;
		//  create a byte buffer and wrap the array
		ByteBuffer bb = ByteBuffer.wrap(header);
		//  if the file uses little endian as apposed to network
		//  (big endian, Java's native) format,
		//  then set the byte order of the ByteBuffer
		if(use_little_endian) bb.order(ByteOrder.LITTLE_ENDIAN);
	
		for (int r = 0; r < info.length; r++){
			info[r] = bb.getInt();
			IJ.log("	" + info[r]);
		}
	
		return info;
	}
	
	public double[][] readfile2matrix(FileInputStream fis, int ncols, int nrows, boolean use_little_endian) throws Exception
	{
		//int ncols = info[0];
		//int nrows = info[1];		
		byte[] data = new byte[ncols*nrows*8];		
		int bytes = fis.read(data);	
		IJ.log("reading bytes="+bytes);
		if(bytes<=0) return null;
		ByteBuffer bb = ByteBuffer.wrap(data);	
		if(use_little_endian) bb.order(ByteOrder.LITTLE_ENDIAN);
		
		double[][] matrix = new double[ncols][nrows];
		for(int i=0; i<nrows; i++){
			for(int j=0; j<ncols; j++){							
					matrix[j][i] = bb.getDouble();
					//IJ.log("	 " + matrix[j][i]);			
			}
		}
	
		//fis.close();
		//IJ.log("loading matrix file done");
		return matrix;
	}


	boolean showDialog(){		
		
		String lastname = OpenDialog.getDefaultDirectory();//OpenDialog.getLastDirectory(); //
		if(lastname!=null) {
			String lastfilename = OpenDialog.getLastName();
			if(lastfilename!=null && lastfilename.toLowerCase().endsWith(".bin")) tracks_file = lastname + lastfilename;//"pillar-tracks.bin";//
			else tracks_file = lastname + "pillar-tracks.bin";//
		}
		int str_len = tracks_file.length();
		if(str_len<50) str_len = 50;
		
		GenericDialogPlus gd = new GenericDialogPlus("Pillar Tracks Exporter");
		gd.addFileField("Binary File Path:", tracks_file, 50);		
		gd.addNumericField("Pixel Size(nm):", pixel_size, 2);
		gd.addMessage("pillar settings");
		gd.addNumericField("pillar diameter in nm:", diameter, 0);
		gd.addNumericField("pillar spacing in nm:", spacing, 0);
		
		gd.showDialog();
		if (gd.wasCanceled()) return false;

		tracks_file = gd.getNextString();		
		pixel_size = gd.getNextNumber();
		diameter=(double)gd.getNextNumber();
		spacing=(double)gd.getNextNumber();
		return true;
	}
}