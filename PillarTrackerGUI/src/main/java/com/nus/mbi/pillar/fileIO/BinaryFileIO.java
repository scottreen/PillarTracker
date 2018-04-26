package com.nus.mbi.pillar.fileIO;

import java.io.DataInputStream;
import java.io.FileInputStream;
import java.nio.ByteBuffer;

/**
 *
 * @author xiaochun
 */
public class BinaryFileIO {
    public static int readInt(FileInputStream fis) throws Exception
    {	            
            int size_int = Integer.BYTES;
            byte[] header = new byte[size_int];
            int bytes = fis.read(header);
            if(bytes<=0) return Integer.MIN_VALUE;
            //  create a byte buffer and wrap the array
            ByteBuffer bb = ByteBuffer.wrap(header); 
            return bb.getInt();
    }
    
    public static byte[] readBytes(FileInputStream fis, int size) throws Exception
    {	    
            byte[] header = new byte[size];
            int bytes = fis.read(header);            
            if(bytes<=0) return null;           
            return header;
    }
    
    public static double readDouble(FileInputStream fis) throws Exception
    {	            
            int size_double = Double.BYTES;
            byte[] header = new byte[size_double];
            int bytes = fis.read(header);
            if(bytes<=0) return Double.NaN;
            //  create a byte buffer and wrap the array
            ByteBuffer bb = ByteBuffer.wrap(header); 
            return bb.getDouble();
    }
    
    public static boolean readBoolean(FileInputStream fis) throws Exception
    {	            
            DataInputStream din = new DataInputStream(fis);
            return din.readBoolean();
    }
    
    public static double[][] readfile2matrix(FileInputStream fis, int ncols, int nrows) throws Exception
    {
            int size_double = Double.BYTES;		
            byte[] data = new byte[ncols*nrows*size_double];		
            int bytes = fis.read(data);	
            //IJ.log("reading bytes="+bytes);
            if(bytes<=0) return null;
            ByteBuffer bb = ByteBuffer.wrap(data);           

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
    
}
