package com.nus.mbi.pillar.solver;

import com.sun.jna.Library;

/**
 *
 * @author xiaochun
 */
public interface GPUCrossCorrelationLibrary extends Library{
        
        void MinimumErrorAMPGPU(
		float[] imagef,
		float[] amp_map,		
		int img_width,
		int img_height, 
		float[] psf, 
		int[] psf_dx,
		int[] psf_dy, 
		int psf_size		
		);
        
        void CrossCorrelationGPU(
		float[] imagef,
		float[] amp_map,		
		int img_width,
		int img_height, 
		float[] psf, 
		int[] psf_dx,
		int[] psf_dy, 
		int psf_size
		);
        
        
	float checkGPU();	
}
