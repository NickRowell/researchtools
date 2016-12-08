/**
 * Copyright (c) 2011 University of Dundee
 *
 * Name:
 * QR.java
 *
 * Purpose:
 * This class simply provides an easier interface to JLAPACK routine, which
 * itself interfaces to LAPACK routines.
 *
 * It is for calculating the QR decomposition with column pivoting. This
 * method is useful because it works on Jama matrices.
 * 
 * Language:
 * Java
 *
 * Author:
 * Nicholas Rowell (nickrowell@computing.dundee.ac.uk)
 *
 */

package numeric.linearalgebra;

// Must include lapack as a library to get java API to basic fortran routines,
// then lapack_simple to get the simplified Java API. Must also include
// f2jutil to get intW class.
//import org.netlib.lapack.DGEQP3;
//import org.netlib.util.intW;

import Jama.*;

public class QR {


    /**
     * This method interfaces to fortran linear algebra package LAPACK using
     * JLAPACK routines. It performs housekeeping operations so that API
     * to LAPACK is in terms of Jama.Matrix objects and not arrays.
     */
    public static Matrix[] qr(Matrix A){

        int m = A.getRowDimension();               // row dimension
        int n = A.getColumnDimension();            // column dimension
        int k = Math.min(m,n);
        // Hack up data required by lapack API
        int[] PIVOT   = new int[n];                // Pivot matrix in sparse vector form
        double[] TAU  = new double[k];             // The scalar factors of the elementary reflectors, whatever the hell that means.
//        double[] WORK = new double[3*n+1];         // Work vector, I think fortran uses elements of this to store temporary values. Or something.
//        intW STATUS   = new intW(-1);              // Status integer returned by fortran routine

        // Get a double array representing matrix
        double[][] a = A.getArrayCopy();

        // Pass to JLAPACK function, which in turn calls fortan routine
//        DGEQP3.DGEQP3(m, n, a, PIVOT, TAU, WORK, WORK.length, STATUS);

        // a, PIVOT and TAU are now overwritten with P,Q,R data

        // Build output matrices:
        Matrix P = new Matrix(new double[n][n]);
        Matrix R = new Matrix(new double[m][n]);
        Matrix Q = new Matrix(new double[m][m]);

        // Set elements of P from pivot vector
        for(int p = 0; p<n;p++)
            P.set(PIVOT[p]-1, p, 1);

        // Set elements of R
        // loop over rows of a
        for(int r=0; r<m; r++)
            // for each consecutive row, read off one less element
            for(int c=r; c<n; c++){

                R.set(r, c, a[r][c]);

                // Reset elements of A so elementary reflectors
                // can be read straight off array.
                if(r==c) a[r][c] = 1;
                else     a[r][c] = 0;

            }

        // Set elements of Q
        //
        // Turn A back into a Matrix so that extracting vectors v is easier
        Matrix A_OUT = new Matrix(a);
        //
        // Store elementary reflectors
        Matrix[] H = new Matrix[k];
        for(int h=0; h<H.length; h++){

            // Read reflector vector v off array A
            Matrix v = A_OUT.getMatrix(0,m-1,h,h);

            H[h] = Matrix.identity(m,m).minus(v.times(v.transpose()).times(TAU[h]));

        }

        // Q = H(0)*H(1)*....*H(k)
        // Take copy of first elementary reflector
        Q = H[0].copy();
        // Multiply all others together
        for(int h=1; h<H.length; h++)
            Q = Q.times(H[h]);

        return new Matrix[]{P,Q,R};

    }

}
