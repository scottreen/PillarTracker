package com.nus.mbi.pillar.solver;

/**
 *
 * @author xiaochun
 */
public class Matrix_Util {
       /*
	   Recursive definition of determinate using expansion by minors.
	*/
	public static double Determinant(double[][] a,int n)
	{
	   int i,j,j1,j2;
	   double det = 0;  
	
	   if (n < 1) { /* Error */
	
	   } else if (n == 1) { /* Shouldn't get used */
	      det = a[0][0];
	   } else if (n == 2) {
	      det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
	   } else {
	      det = 0;
	      for (j1=0;j1<n;j1++) {
		double[][] m = new double[n-1][n-1];
	         for (i=1;i<n;i++) {
	            j2 = 0;
	            for (j=0;j<n;j++) {
	               if (j == j1)
	                  continue;
	               m[i-1][j2] = a[i][j];
	               j2++;
	            }
	         }
	         det += Math.pow(-1.0,j1+2.0) * a[0][j1] * Determinant(m,n-1);
	      }
	   }
	   return(det);
	}
	
	/*
	   Find the cofactor matrix of a square matrix
	*/
	public static void CoFactor(double[][] a,int n,double[][] b)
	{
	   int i,j,ii,jj,i1,j1;
	   double det;
	   double[][] c = new double[n-1][n-1];
	
	   for (j=0;j<n;j++) {
	      for (i=0;i<n;i++) {
	
	         /* Form the adjoint a_ij */
	         i1 = 0;
	         for (ii=0;ii<n;ii++) {
	            if (ii == i)
	               continue;
	            j1 = 0;
	            for (jj=0;jj<n;jj++) {
	               if (jj == j)
	                  continue;
	               c[i1][j1] = a[ii][jj];
	               j1++;
	            }
	            i1++;
	         }
	
	         /* Calculate the determinate */
	         det = Determinant(c,n-1);
	
	         /* Fill in the elements of the cofactor */
	         b[i][j] = Math.pow(-1.0,i+j+2.0) * det;
	      }
	   }
	}
	
	/*
	   Transpose of a square matrix, do it in place
	*/
	public static void Transpose(double[][] a,int n)
	{
	   int i,j;
	   double tmp;
	
	   for (i=1;i<n;i++) {
	      for (j=0;j<i;j++) {
	         tmp = a[i][j];
	         a[i][j] = a[j][i];
	         a[j][i] = tmp;
	      }
	   }
	}
        
        public static double[][] transpose(double[][] a){
                int m = a.length;
                int n = a[0].length;	
                int i,j;	   
                double[][] b= new double[n][m]; 
                for(i=0;i<n;i++) for(j=0;j<m;j++) b[i][j] = a[j][i];
                return b;
        }
	
	public static void InverseMatrix(double[][] a,double[][] b, int n)
	{
		if(n==1) b[0][0] = 1/a[0][0];
		else{
			double[][] c = new double[n][n];	
			for (int i=0;i<n;i++) for(int j=0; j<n; j++) c[i][j]=a[i][j];	
			double det = Determinant(c,n);		
			double[][] d = new double[n][n];				
			CoFactor(c,n,d);
			Transpose(d,n);
			for(int i=0; i<n; i++) for(int j=0; j<n; j++) b[i][j] = d[i][j]/det;
		}
	}

	public static double[][] MatrixTans(double[][]a){
	   int m = a.length;
           int n = a[0].length;	
	   int i,j;	   
	   double[][] b= new double[n][m]; 
	   for(i=0;i<n;i++) for(j=0;j<m;j++) b[i][j] = a[j][i];
	   return b;
	}

        public static double[][] MatrixMul(double[][] A, double[][] B)
        {
            int Am = A.length;
            int An = A[0].length;

            int Bm = B.length;
            int Bn = B[0].length;

            int Cm=Am;
            int Cn=Bn;
            if (An != Bm)
            {
                Cm = Cn = 0;
                //throw new Exception("Matrix Mulitiply error!");
                return null;
            }
            
            double[][] C = new double[Cm][Cn];

            for (int i = 0; i < Cm; i++)
            {
                for (int j = 0; j < Cn; j++)
                {
                    C[i][j] = 0;
                    for (int k = 0; k < An; k++)
                    {
                        C[i][j] += (A[i][k]*B[k][j]);
                    }
                }
            }
            
            return C;
        }

        public static double[][] MatrixSub(double[][] A, double[][] B)
        {
            int Am = A.length;
            int An = A[0].length;

            int Bm = B.length;
            int Bn = B[0].length;

            if (An != Bn && Am != Bm) return null;
            
            double[][] C = new double[Am][An];

            for (int i = 0; i < Am; i++)
            {
                for (int j = 0; j < An; j++)
                {
                    C[i][j] = A[i][j]-B[i][j];
                }
            }

            return C;
        }

        public static double[][] MatrixAdd(double[][] A, double[][] B)
        {
            int Am = A.length;
            int An = A[0].length;

            int Bm = B.length;
            int Bn = B[0].length;

            if (An != Bn && Am != Bm) return null;
            
            double[][] C = new double[Am][An];

            for (int i = 0; i < Am; i++)
            {
                for (int j = 0; j < An; j++)
                {
                    C[i][j] = A[i][j]+B[i][j];
                }
            }

            return C;
        }

        public static double[] MatrixMul(double[][] A, double[] b)
        {
            int Am = A.length;
            int An = A[0].length;

            int Bn = b.length;

            int Cm = Am;
            if (An != Bn)
            {
                Cm = 0;
                return null;
            }

            double[] c = new double[Cm];

            for (int i = 0; i < Cm; i++)
            {
                c[i] = 0;
                for (int k = 0; k < An; k++)
                {
                    c[i] += (A[i][k] * b[k]);
                }
            }

            return c;
        }

	//matrix A transpose multiply with vector b
        public static double[] MatrixMulTrans(double[][] A, double[] b)
        {
            int Am = A.length;
            int An = A[0].length;

            int Bn = b.length;

            int Cn = An;
            if (Am != Bn)
            {
                Cn = 0;
                return null;
            }

            double[] c = new double[Cn];

            for (int i = 0; i < Cn; i++)
            {
                c[i] = 0;
                for (int k = 0; k < Am; k++)
                {
                    c[i] += (A[k][i] * b[k]);
                }
            }

            return c;
        }

        public static double[] VectorAdd(double[] a, double[] b)
        {
            int an = a.length;
            int bn = b.length;
            if(an!= bn) return null;
            double[] c = new double[an];
            for(int i=0; i<an; i++) c[i] = a[i] + b[i];
            return c;
        }

        public static double VectorMul(double[] a, double[] b)
        {
            int an = a.length;
            int bn = b.length;
            if (an != bn) return 0;
            double c = 0;
            for (int i = 0; i < an; i++) c += (a[i] + b[i]);
            return c;
        }

	public static double getAbs(double[] vector, int size)
	{
		double sum = 0;
		for(int i=0; i<size; i++)
		{
			double t = vector[i];
			sum+=(t*t);
		}
	
		return sum;
	}
	
	public static double getMax(double[][] M, int n)
	{
		double m = M[0][0];
		for(int i=0; i<n; i++)
		{
			//for(int j=0; j<n; j++)
			{
				double t = M[i][i];
				if(m<t) m = t;
			}
		}
		return m;
	}
	
	public static double[][] IdentityMatrix(int n){
		double[][] I = new double[n][n];	
		for(int i=0; i<n; i++)	for(int j=0; j<n; j++)	I[i][j]= (i==j)?1:0;
		return I;
	}

}
