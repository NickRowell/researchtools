package misc;

import Jama.*;

// Must include lapack as a library to get java API to basic fortran routines,
// then lapack_simple to get the simplified Java API. Must also include
// f2jutil to get intW class.
//import org.netlib.lapack.DGEQP3;
//import org.netlib.util.intW;

/**
 * Can probably delete this class. Not sure whether a test function is
 * necessary for the QR class. Note that this is not a working test function
 * anyway.
 *
 * @author nickrowell
 */
public class QRTester{


    public static void main(String[] args){

        // Test LAPACK QR decomposition algorithm

        // input matrix - rows have varying linear dependency
        double[][] d = new double[][]{
                                      {1,0,0},
                                      {2,0,0},
                                      {0,0,4},
                                      {0,3,0},
                                      {0,1,0},
                                      {0,5,0},
                                      {1,5,0},
                                      {0,0,1},
                                      {0,2,0},
                                      {0,0,2},
                                      {0,0,2},
        };

        // Make a matrix
        Matrix A = new Matrix(d);

        // Now take transpose so that columns have varying linear dependency
        A = A.transpose();

        int m = A.getRowDimension();                                // row dimension
        int n = A.getColumnDimension();                                 // column dimension
        // Hack up data required by lapack API
        int[] PIVOT = new int[n];                  // Pivot matrix in sparse vector form
        double[] TAU = new double[Math.min(m,n)];  // The scalar factors of the elementary reflectors, whatever the hell that means.
        double[] WORK = new double[3*n+1];         // Work vector, I think fortran uses elements of this to store temporary values. Or something.
//        intW STATUS = new intW(-1);                // Status integer returned by fortran routine


        // Get a double array representing matrix
        double[][] a = A.getArrayCopy();

        // Pass to LAPACK function
//        DGEQP3.DGEQP3(m, n, a, PIVOT, TAU, WORK, WORK.length, STATUS);

        // Print status flag
        //System.out.println("STATUS = "+STATUS.val);
        //System.out.println("Optimal work vector size = "+WORK[0]+", used "+(LDA+1));

        // Now hack array A into Q,R and check against MATLAB algorithm

        // Make output matrices:
        Matrix P = new Matrix(new double[n][n]);
        Matrix R = new Matrix(new double[m][n]);
        Matrix Q = new Matrix(new double[m][m]);
        

        // Set elements of P
        for(int p = 0; p<n;p++)
            P.set(PIVOT[p]-1, p, 1);
        
        // Set elements of R
        // loop over rows
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
        // Turn A back into a Matrix so that getting vectors v is easier
        Matrix A_OUT = new Matrix(a);
        //
        // Store elementary reflectors
        Matrix[] H = new Matrix[Math.min(m,n)];
        for(int h=0; h<H.length; h++){

            // Read reflector vector v off array A
            Matrix v = A_OUT.getMatrix(0,m-1,h,h);

            H[h] = Matrix.identity(m,m).minus(v.times(v.transpose()).times(TAU[h]));

        }

        // Take copy of first elementary reflector
        Q = H[0].copy();
        // Multiply all others together
        for(int h=1; h<H.length; h++)
            Q = Q.times(H[h]);
        



        System.out.println("Input matrix A = "); A.print(2,0);


        System.out.println("Pivot matrix P = "); P.print(2, 0);
        

        //System.out.println("Matrix Q = "); Q.print(5, 5);


        //System.out.println("Matrix R = "); R.print(5,5);


        System.out.println("A * P = "); A.times(P).print(2, 0);


        //System.out.println("QR = "); Q.times(R).print(5, 5);
        

        // Exit
        System.exit(0);








    }



}

