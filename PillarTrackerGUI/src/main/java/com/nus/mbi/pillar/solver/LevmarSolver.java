package com.nus.mbi.pillar.solver;

import com.nus.mbi.pillar.tracker.LocalizationWindow;

/**
 * rewrite the lev-mar algorithm as following cite, under GPL license
 * @misc{lourakis04LM,
    author={M.I.A. Lourakis},
    title={levmar: Levenberg-Marquardt nonlinear least squares algorithms in {C}/{C}++},
    howpublished={[web page] \verb+http://www.ics.forth.gr/~lourakis/levmar/+},
    year={Jul. 2004},
    note={[Accessed on 31 Jan. 2005.]},
    }   
 * @author xiaochun
 */
public class LevmarSolver {     
    public final static double EPSILON   =    1E-12;
    public final static double ONE_THIRD  =   0.3333333334; /* 1.0/3.0 */
    public final static int    LM_OPTS_SZ   = 	 4; /* max(4, 5) */
    public final static int    LM_ERROR    =     -1;
    public final static double LM_INIT_MU    =	 1E-03;
    public final static double LM_STOP_THRESH = 1E-17;    
    public final static double DBL_MAX = Double.MAX_VALUE;
    public final static double DBL_MIN = Double.MIN_VALUE;

    public double[] opts = new double[LM_OPTS_SZ];
    
    public LevmarSolver(){ 
	/* optimization control parameters; passing to levmar NULL instead of opts reverts to defaults */
	opts[0]=1E-03; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;	
    }
    
    /* evaluate the least square function */
    double fevl(double[] p, double[] I, int n, double[] winx, double[] winy, double sigma2)
    {
        int i;  

        double evl = 0;
        for(i=0; i<n; ++i){		
            double dx = winx[i] - p[2];
            double dy = winy[i] - p[3];
            double df = I[i] - (p[0]*Math.exp(-(dx*dx + dy*dy)/sigma2) + p[1]);
            evl += (df*df);
        }
        return evl;
    }

    /* Jacobian of expfunc() */
    void fjac(double[] p, double[] I, double[] jacTjac, double[] jacTe, int n, double[] winx, double[] winy, double sigma2)
    {   
       int m = 4;
       int i, j, k;
       double[] g = {0, 0, 0, 0};
       for(i=m*m; i-->0; ) jacTjac[i]=0.0;
       for(i=m; i-->0; ) jacTe[i]=0.0;
       /* fill Jacobian row by row 
        * Since J^T J is symmetric, its computation can be sped up by computing
        * only its upper triangular part and copying it to the lower part
        */
       for(i=0; i<n; ++i){
            double dx = winx[i] - p[2];
            double dy = winy[i] - p[3];
            double e = Math.exp(-(dx*dx + dy*dy)/sigma2);
            double df = I[i] - (p[0]*e + p[1]);

            double t = 2*e*p[0]/sigma2;
            g[0] = e;
            g[1] = 1;
            g[2] = t*dx;
            g[3] = t*dy;	

            for(j=0; j<m; j++){			
                    jacTe[j] += (df*g[j]);
                    for(k=j+1; k-->0;) /* j<=i computes lower triangular part only */	
                            jacTjac[j*m+k] += (g[j]*g[k]);			
            }		
        }

        for(i=m; i-->0; ) /* copy to upper part */
                for(j=i+1; j<m; ++j)
                        jacTjac[i*m+j]=jacTjac[j*m+i];   

    }   
    
    void updatePSF(double[] p, double[] psf, int n, double[] winx, double[] winy, double sigma2)
    {
        int i;          
        for(i=0; i<n; ++i){		
            double dx = winx[i] - p[2];
            double dy = winy[i] - p[3];
            psf[i] = Math.exp(-(dx*dx + dy*dy)/sigma2);
        }
    }
    
    /* evaluate the least square function */
    double fevl(double[] p, double[] I, double[] psf, int n, double[] winx, double[] winy, double sigma2)
    {
        int i;  

        double evl = 0;
        for(i=0; i<n; ++i){		
            double e = psf[i]; 
            double df = I[i] - (p[0]*e + p[1]);
            evl += (df*df);
        }
        return evl;
    }

    /* compute Jacobian of least square function */
    void fjac(double[] p, double[] I, double[] psf, double[] jacTjac, double[] jacTe, int n, double[] winx, double[] winy, double sigma2)
    {   
       int m = 4;
       int i, j, k;
       double[] g = {0, 0, 0, 0};
       for(i=m*m; i-->0; ) jacTjac[i]=0.0;
       for(i=m; i-->0; ) jacTe[i]=0.0;
       /* 
        * Since Jacobian matrix is symmetric, 
        * copye the lower triangular part to the upper part
        */
       for(i=0; i<n; ++i){
            double dx = winx[i] - p[2];
            double dy = winy[i] - p[3];
            double e = psf[i];//Math.exp(-(dx*dx + dy*dy)/sigma2);
            double df = I[i] - (p[0]*e + p[1]);

            double t = 2*e*p[0]/sigma2;
            g[0] = e;
            g[1] = 1;
            g[2] = t*dx;
            g[3] = t*dy;	

            for(j=0; j<m; j++){			
                    jacTe[j] += (df*g[j]);
                    for(k=j+1; k-->0;) /* only lower triangular part is computed */	
                            jacTjac[j*m+k] += (g[j]*g[k]);			
            }		
        }

        for(i=m; i-->0; ) /* copy to upper part */
                for(j=i+1; j<m; ++j)
                        jacTjac[i*m+j]=jacTjac[j*m+i];   

    }   
    
    int AX_EQ_B_INVERT4(double[] m, double[] b, double[] x)
    {
        double[] inv = new double[16];    
        inv[0] = m[5]  * m[10] * m[15] - 
                 m[5]  * m[11] * m[14] - 
                 m[9]  * m[6]  * m[15] + 
                 m[9]  * m[7]  * m[14] +
                 m[13] * m[6]  * m[11] - 
                 m[13] * m[7]  * m[10];

        inv[4] = -m[4]  * m[10] * m[15] + 
                  m[4]  * m[11] * m[14] + 
                  m[8]  * m[6]  * m[15] - 
                  m[8]  * m[7]  * m[14] - 
                  m[12] * m[6]  * m[11] + 
                  m[12] * m[7]  * m[10];

        inv[8] = m[4]  * m[9] * m[15] - 
                 m[4]  * m[11] * m[13] - 
                 m[8]  * m[5] * m[15] + 
                 m[8]  * m[7] * m[13] + 
                 m[12] * m[5] * m[11] - 
                 m[12] * m[7] * m[9];

        inv[12] = -m[4]  * m[9] * m[14] + 
                   m[4]  * m[10] * m[13] +
                   m[8]  * m[5] * m[14] - 
                   m[8]  * m[6] * m[13] - 
                   m[12] * m[5] * m[10] + 
                   m[12] * m[6] * m[9];

        inv[1] = -m[1]  * m[10] * m[15] + 
                  m[1]  * m[11] * m[14] + 
                  m[9]  * m[2] * m[15] - 
                  m[9]  * m[3] * m[14] - 
                  m[13] * m[2] * m[11] + 
                  m[13] * m[3] * m[10];

        inv[5] = m[0]  * m[10] * m[15] - 
                 m[0]  * m[11] * m[14] - 
                 m[8]  * m[2] * m[15] + 
                 m[8]  * m[3] * m[14] + 
                 m[12] * m[2] * m[11] - 
                 m[12] * m[3] * m[10];

        inv[9] = -m[0]  * m[9] * m[15] + 
                  m[0]  * m[11] * m[13] + 
                  m[8]  * m[1] * m[15] - 
                  m[8]  * m[3] * m[13] - 
                  m[12] * m[1] * m[11] + 
                  m[12] * m[3] * m[9];

        inv[13] = m[0]  * m[9] * m[14] - 
                  m[0]  * m[10] * m[13] - 
                  m[8]  * m[1] * m[14] + 
                  m[8]  * m[2] * m[13] + 
                  m[12] * m[1] * m[10] - 
                  m[12] * m[2] * m[9];

        inv[2] = m[1]  * m[6] * m[15] - 
                 m[1]  * m[7] * m[14] - 
                 m[5]  * m[2] * m[15] + 
                 m[5]  * m[3] * m[14] + 
                 m[13] * m[2] * m[7] - 
                 m[13] * m[3] * m[6];

        inv[6] = -m[0]  * m[6] * m[15] + 
                  m[0]  * m[7] * m[14] + 
                  m[4]  * m[2] * m[15] - 
                  m[4]  * m[3] * m[14] - 
                  m[12] * m[2] * m[7] + 
                  m[12] * m[3] * m[6];

        inv[10] = m[0]  * m[5] * m[15] - 
                  m[0]  * m[7] * m[13] - 
                  m[4]  * m[1] * m[15] + 
                  m[4]  * m[3] * m[13] + 
                  m[12] * m[1] * m[7] - 
                  m[12] * m[3] * m[5];

        inv[14] = -m[0]  * m[5] * m[14] + 
                   m[0]  * m[6] * m[13] + 
                   m[4]  * m[1] * m[14] - 
                   m[4]  * m[2] * m[13] - 
                   m[12] * m[1] * m[6] + 
                   m[12] * m[2] * m[5];

        inv[3] = -m[1] * m[6] * m[11] + 
                  m[1] * m[7] * m[10] + 
                  m[5] * m[2] * m[11] - 
                  m[5] * m[3] * m[10] - 
                  m[9] * m[2] * m[7] + 
                  m[9] * m[3] * m[6];

        inv[7] = m[0] * m[6] * m[11] - 
                 m[0] * m[7] * m[10] - 
                 m[4] * m[2] * m[11] + 
                 m[4] * m[3] * m[10] + 
                 m[8] * m[2] * m[7] - 
                 m[8] * m[3] * m[6];

        inv[11] = -m[0] * m[5] * m[11] + 
                   m[0] * m[7] * m[9] + 
                   m[4] * m[1] * m[11] - 
                   m[4] * m[3] * m[9] - 
                   m[8] * m[1] * m[7] + 
                   m[8] * m[3] * m[5];

        inv[15] = m[0] * m[5] * m[10] - 
                  m[0] * m[6] * m[9] - 
                  m[4] * m[1] * m[10] + 
                  m[4] * m[2] * m[9] + 
                  m[8] * m[1] * m[6] - 
                  m[8] * m[2] * m[5];

        double det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

        if (det == 0) return 0;

        det = 1.0 / det;
        for (int i = 0; i<4; i++){
                x[i] = 0;
                int ioff = i*4; 
                for(int j=0; j<4; j++){
                        x[i] += (inv[ioff+j]*b[j]);
                }
                x[i] =  x[i]*det;
        }

        return 1;
    }

    /*
    common function for least square levmar algorithm.
    */
    int levmar_double(
      double[] p,         /* I/O: initial parameter estimates. On output has the estimated solution */
      double[] x,          /* I: measurement vector. NULL implies a zero vector */         
      int n,              /* I: measurement vector dimension */
      int itmax,          /* I: maximum number of iterations */
      double[] opts,    /* I: minim. options [\mu, \epsilon1, \epsilon2, \epsilon3]. Respectively the scale factor for initial \mu,
                           * stopping thresholds for ||J^T e||_inf, ||Dp||_2 and ||e||_2. Set to NULL for defaults to be used
                           */      
      double[] winx, 
      double[] winy,
      double sigma2)     
    {
        final int m = 4; /* I: parameter vector dimension (i.e. #unknowns) */
        int i, k;
        int issolved;
        /* temp work arrays */
        double[] jacTe = new double[m];      /* J^T e_i mx1 */
        double[] jacTjac = new double[m*m];    /* mxm */
        double[] Dp = new double[m];         /* mx1 */
        double[] diag_jacTjac = new double[m];   /* diagonal of J^T J, mx1 */
        double[] pDp = new double[m];        /* p + Dp, mx1 */
        double[] psf = new double[n];
        double[] newpsf = new double[n];
        double mu,  /* damping constant */
               tmp; /* mainly used in matrix & vector multiplications */
        double p_eL2, jacTe_inf, pDp_eL2; /* ||e(p)||_2, ||J^T e||_inf, ||e(p+Dp)||_2 */
        double p_L2, Dp_L2=DBL_MAX, dF, dL;
        double tau, eps1, eps2, eps2_sq, eps3;
        double init_p_eL2;
        int nu=2, nu2, stop=0, nfev, njev=0, nlss=0;
        //final int nm=n*m;

        mu=jacTe_inf=0.0; /* -Wall */

        if(n<m){
          //fprintf(stderr, LCAT(LEVMAR_DER, "(): cannot solve a problem with fewer measurements [%d] than unknowns [%d]\n"), n, m);
          return LM_ERROR;
        }

        if(opts!=null){
                tau=opts[0];
                eps1=opts[1];
                eps2=opts[2];
                eps2_sq=opts[2]*opts[2];
                eps3=opts[3];
        }
        else{ // use default values
                tau=LM_INIT_MU;
                eps1=LM_STOP_THRESH;
                eps2=LM_STOP_THRESH;
                eps2_sq=LM_STOP_THRESH*LM_STOP_THRESH;
                eps3=LM_STOP_THRESH;
        }

        /* compute e=x - f(p) and its L2 norm */  
        updatePSF(p, psf, n, winx, winy, sigma2);
        p_eL2 = fevl(p, x, psf, n, winx, winy, sigma2); nfev=1; 
        
        init_p_eL2=p_eL2;
        if(Double.isInfinite(p_eL2)) stop=7;

        for(k=0; k<itmax && stop==0; ++k){
          /* Note that p and e have been updated at a previous iteration */
          if(p_eL2<=eps3){ /* error is small */
            stop=6;
            break;
          }

          /* J^T J, J^T e */    
          fjac(p, x,  psf, jacTjac, jacTe, n, winx, winy, sigma2); ++njev;

          /* Compute ||J^T e||_inf and ||p||^2 */
          for(i=0, p_L2=jacTe_inf=0.0; i<m; ++i){
            if(jacTe_inf < (tmp=Math.abs(jacTe[i]))) jacTe_inf=tmp;

            diag_jacTjac[i]=jacTjac[i*m+i]; /* save diagonal entries so that augmentation can be later canceled */
            p_L2+=p[i]*p[i];
          }
          //p_L2=sqrt(p_L2);

          /* check for convergence */
          if((jacTe_inf <= eps1)){
            Dp_L2=0.0; /* no increment for p in this case */
            stop=1;
            break;
          }

          /* compute initial damping factor */
          if(k==0){
            for(i=0, tmp=DBL_MIN; i<m; ++i)
              if(diag_jacTjac[i]>tmp) tmp=diag_jacTjac[i]; /* find max diagonal element */
            mu=tau*tmp;
          }

          /* determine increment using adaptive damping */
          while(true){
            /* augment normal equations */
            for(i=0; i<m; ++i) jacTjac[i*m+i]+=mu;

            issolved=AX_EQ_B_INVERT4(jacTjac, jacTe, Dp); ++nlss;
            if(issolved>0){
                /* compute p's new estimate and ||Dp||^2 */
                for(i=0, Dp_L2=0.0; i<m; ++i){
                   pDp[i]=p[i] + (tmp=Dp[i]);
                   Dp_L2+=tmp*tmp;
                }
                
                if(Dp_L2<=eps2_sq*p_L2){ /* relative change in p is small, stop */                 
                   stop=2;
                   break;
                }
                
                if(Dp_L2>=(p_L2+eps2)/(EPSILON*EPSILON)) { /* almost singular */                
                  stop=4;
                  break;
                }
                
                /* compute ||e(pDp)||_2 */   
                updatePSF(pDp, newpsf, n, winx, winy, sigma2);
                pDp_eL2 = fevl(pDp, x,  newpsf, n, winx, winy, sigma2); ++nfev; /* evaluate function at p + Dp */
                if(Double.isInfinite(pDp_eL2)){ /* sum of squares is not finite, most probably due to a user error.
                                          * This check makes sure that the inner loop does not run indefinitely.   
                                          * Thanks to Steve Danauskas for reporting such cases
                                          */
                  stop=7;
                  break;
                }

                for(i=0, dL=0.0; i<m; ++i)
                  dL+=Dp[i]*(mu*Dp[i]+jacTe[i]);

                dF=p_eL2-pDp_eL2;

                if(dL>0.0 && dF>0.0){ /* reduction in error, increment is accepted */
                  tmp=2.0*dF/dL-1.0;
                  tmp=1.0-tmp*tmp*tmp;
                  mu=mu*( (tmp>=ONE_THIRD)? tmp : ONE_THIRD );
                  nu=2;            
                  for(i=0 ; i<m; ++i) p[i]=pDp[i];/* update p's estimate */
                  for(i=0 ; i<n; ++i) psf[i]=newpsf[i];/* update psf's estimate */
                  p_eL2=pDp_eL2; /* update ||e||_2 */
                  break;
                }
            }

            /* if this point is reached, either the linear system could not be solved or
             * the error did not reduce; in any case, the increment must be rejected
             */
            mu*=nu;
            nu2=nu<<1; // 2*nu;
            if(nu2<=nu){ /* nu has wrapped around (overflown). Thanks to Frank Jordan for spotting this case */
              stop=5;
              break;
            }
            nu=nu2;

            for(i=0; i<m; ++i) /* restore diagonal J^T J entries */
              jacTjac[i*m+i]=diag_jacTjac[i];
          } /* inner loop */
        }

        if(k>=itmax) stop=3;

        for(i=0; i<m; ++i) /* restore diagonal J^T J entries */
          jacTjac[i*m+i]=diag_jacTjac[i];

        return (stop!=4 && stop!=7)?  k : LM_ERROR;
    }
    
    /*direct solver for a and b*/
    double direct_solver(double[] win, double[] win_x, double[] win_y, int win_size, double sigma2, double[] v)
    {	
        double x0 = v[2];
        double y0 = v[3];
        //IJ.log("x0=" + x0 + "y0=" + y0 + "\n");

        double amp = 0; 
        double bkg = 0;
        double sigmaY = 0.0;
        double sigmaY2 = 0.0;		
        double sigmaXY = 0.0;		
        double sigmaX = 0;
        double sigmaX2 = 0;


        double min_bkg = win[0];
        for(int i=0; i<win_size; i++)
        {			
            if(min_bkg>win[i]) min_bkg = win[i];
            sigmaX += win[i];
            sigmaX2 += win[i]*win[i];			

            double dx = win_x[i]-x0;
            double dy = win_y[i]-y0;
            double psfi = Math.exp(-(dx*dx+dy*dy)/sigma2);			
            sigmaXY += psfi*win[i];
            sigmaY +=  psfi;
            sigmaY2 += psfi*psfi;		
        }	

        //return 0;		
        double err = sigmaX2;
        double det = sigmaY2*win_size-sigmaY*sigmaY;	
        if(det>1e-20 || det<-1e-20){
                amp = (sigmaXY*win_size-sigmaX*sigmaY)/det;
                bkg = (sigmaX - amp*sigmaY)/win_size;
                err = amp*amp*sigmaY2 + win_size*bkg*bkg + sigmaX2 + 2*amp*bkg*sigmaY - 2*amp*sigmaXY - 2*bkg*sigmaX;
        }		
        
        v[0] = amp;
        v[1] = bkg;

        double m10 = 0;
        double m01 = 0;
        double m00 = 0;
//        if(bkg<min_bkg && bkg>=0) min_bkg = bkg;
        if(bkg<min_bkg) min_bkg = bkg;
        
        for(int i=0; i<win_size; i++)
        {			
                double w = win[i] - min_bkg;
                m10 += win_x[i]*w;
                m01 += win_y[i]*w;
                m00 += w;
        }

        v[2] = m10/m00;
        v[3] = m01/m00;
        //err = err/2;
        return err;
    }
    
    
    static double direct_solver(double[] win, double[] win_x, double[] win_y, int win_size, double sigma2, double[] v, boolean dark_object)
    {	
        double x0 = v[2];
        double y0 = v[3];
        //IJ.log("x0=" + x0 + "y0=" + y0 + "\n");

        double amp = 0; 
        double bkg = 0;
        double sigmaY = 0.0;
        double sigmaY2 = 0.0;		
        double sigmaXY = 0.0;		
        double sigmaX = 0;
        double sigmaX2 = 0;


        double min_bkg = win[0];
        double max_bkg = win[0];
        for(int i=0; i<win_size; i++)
        {			
            double wini = win[i];
            if(min_bkg>wini) min_bkg = wini;
            if(max_bkg<wini) max_bkg = wini;
            sigmaX += wini;
            sigmaX2 += wini*wini;			

            double dx = win_x[i]-x0;
            double dy = win_y[i]-y0;
            double psfi = Math.exp(-(dx*dx+dy*dy)/sigma2);			
            sigmaXY += psfi*wini;
            sigmaY +=  psfi;
            sigmaY2 += psfi*psfi;		
        }	

        //return 0;		
        double err = sigmaX2;
        double det = sigmaY2*win_size-sigmaY*sigmaY;	
        if(det>1e-20 || det<-1e-20){
                amp = (sigmaXY*win_size-sigmaX*sigmaY)/det;
                bkg = (sigmaX - amp*sigmaY)/win_size;
                err = amp*amp*sigmaY2 + win_size*bkg*bkg + sigmaX2 + 2*amp*bkg*sigmaY - 2*amp*sigmaXY - 2*bkg*sigmaX;
        }		
        
        v[0] = amp;
        v[1] = bkg;

        double m10 = 0;
        double m01 = 0;
        double m00 = 0;
//        if(bkg<min_bkg && bkg>=0) min_bkg = bkg;
        if(bkg<min_bkg) min_bkg = bkg;
        if(bkg>max_bkg) max_bkg = bkg;
        
        for(int i=0; i<win_size; i++)
        {			
                double w = dark_object ? (max_bkg - win[i]) : (win[i] - min_bkg);
                m10 += win_x[i]*w;
                m01 += win_y[i]*w;
                m00 += w;
        }

        v[2] = m10/m00;
        v[3] = m01/m00;
        //err = err/2;
        return err;
    }
    
    //load window, dx0, dy0
    static int load_window(double[] pixels, int width, int height, int xc, int yc, int window_r, double[] win, double[] winx, double[] winy)
    {
        //IJ.log("x0=" + x0 + "y0=" + y0 + "\n");
        int k=0;
        for(int x=-window_r; x<=window_r; x++){
            for(int y=-window_r; y<=window_r; y++){
                int i=y+yc;
                int j=x+xc;
                if(i>=0 && i<height && j>=0 && j<width){					
                        winx[k] = x;
                        winy[k] = y;
                        double t = pixels[i*width+j];					
                        win[k] = t;				

                        k++;
                }
            }
        }
        
        return k;
    }
    
    public boolean estimate_gaussian_fit(int xc, int yc, double[] pixels, int width, int height, int kernel_width, double gauss_sigma, double box_constrian_limit, boolean dark_object, double[] para)
    {						
        int window_r = kernel_width/2;
        int window_w = window_r*2 + 1;
        int window_s = window_w*window_w;
        double[] win = new double[window_s];
        double[] winx= new double[window_s];
        double[] winy= new double[window_s];            
        int win_size = load_window(pixels,width,height,xc,yc,window_r, win, winx, winy);            
        double sigma2 = gauss_sigma*gauss_sigma*2;

        double[] par_seed = {para[0], para[1], 0, 0};
        direct_solver(win, winx, winy, win_size, sigma2, par_seed, dark_object);
        //printf("seed par in %d:\n	a=%f	b=%f	x0=%f	y0=%f\n",iter, par_seed[0],par_seed[1], par_seed[2],par_seed[3]);
        double[] new_para = {par_seed[0], par_seed[1],par_seed[2], par_seed[3]}; /* starting value */		
        //double[] new_para = {par_seed[0], par_seed[1],0, 0}; /* starting value */		
        /* invoke the optimization function */
        int ret=levmar_double(new_para, win, win_size, 100, opts, winx, winy, sigma2); // with analytic Jacobian		

        boolean suc = false;
        if(ret>=0){           
//            if(new_para[1]>=0 && new_para[2] != par_seed[2] && new_para[3] != par_seed[3]){ //check the direct solver?
//            if(new_para[2] != par_seed[2] && new_para[3] != par_seed[3]){ //check the direct solver?
//                    suc = true;                    
//            }
            suc=true;
        }
        
        for(int i=0; i<4; i++) new_para[i] = suc? new_para[i] : par_seed[i]; 
        
        //check the constrains and dark pillars
        if(suc){
            suc = false;
            if((dark_object==false && new_para[0]>0) || (dark_object && new_para[0]<0)){
                if(Math.abs(new_para[2])<box_constrian_limit && Math.abs(new_para[3])<box_constrian_limit){
                    for(int i=0; i<4; i++) para[i] = new_para[i]; 
                    suc = true;
                }
            }                
        }
                   
        return suc;		
    }
    
    public boolean estimate_gaussian_fit(LocalizationWindow window, double gauss_sigma, double box_constrian_limit, boolean dark_object, double[] para)
    {						
        double[] win = window.win;
        double[] winx= window.winx;
        double[] winy= window.winy;
        int win_size = window.win_size;
        double sigma2 = gauss_sigma*gauss_sigma*2;

        double[] par_seed = {para[0], para[1], 0, 0};
        direct_solver(win, winx, winy, win_size, sigma2, par_seed, dark_object);
        //printf("seed par in %d:\n	a=%f	b=%f	x0=%f	y0=%f\n",iter, par_seed[0],par_seed[1], par_seed[2],par_seed[3]);
        double[] new_para = {par_seed[0], par_seed[1],par_seed[2], par_seed[3]}; /* starting value */		
        //double[] new_para = {par_seed[0], par_seed[1],0, 0}; /* starting value */		
        /* invoke the optimization function */
        int ret=levmar_double(new_para, win, win_size, 100, opts, winx, winy, sigma2); // with analytic Jacobian		

        boolean suc = false;
        if(ret>=0){           
//            if(new_para[1]>=0 && new_para[2] != par_seed[2] && new_para[3] != par_seed[3]){ //check the direct solver?
//            if(new_para[2] != par_seed[2] && new_para[3] != par_seed[3]){ //check the direct solver?
//                    suc = true;                    
//            }
            suc=true;
        }
        
        for(int i=0; i<4; i++) new_para[i] = suc? new_para[i] : par_seed[i]; 
        
        //check the constrains and dark pillars
        if(suc){
            suc = false;
            if((dark_object==false && new_para[0]>0) || (dark_object && new_para[0]<0)){
                if(Math.abs(new_para[2])<box_constrian_limit && Math.abs(new_para[3])<box_constrian_limit){
                    for(int i=0; i<4; i++) para[i] = new_para[i]; 
                    suc = true;
                }
            }                
        }
                   
        return suc;		
    }
}

