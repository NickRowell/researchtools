/**
 * Copyright (c) 2011 University of Dundee
 *
 * Name:
 * Ellipse.java
 *
 * Purpose:
 * Class represents ellipses. These are defined in terms of the conic 
 * conic parameters and more human-readable semi-axes and centre etc.
 * Class also provides methods to fit ellipses to sets of data points by
 * various techniques.
 * 
 * Language:
 * Java
 *
 * Author:
 * Nicholas Rowell (nickrowell@computing.dundee.ac.uk)
 *
 */

package numeric.fitting;

import java.text.DecimalFormat;
import java.util.ArrayList;
import Jama.*;


public class Ellipse {

    // Conic parameters of ellipse a*x**2 + b*x*y + c*y**2 + d*x + e*y + f = 0
    public double a;
    public double b;
    public double c;
    public double d;
    public double e;
    public double f;

    // More human-readable ellipse parameters

    /** x coordinate of centre */
    public double x0=0;
    /** y coordinate of centre */
    public double y0=0;
    /** major semi-axis */
    public double r_a=0;
    /** minor semi-axis */
    public double r_b=0;
    /** Anti-clockwise angle of rotation from x axis */
    public double theta=0;    

    // Other...

    /** Condition number of matrix used in ellipse fitting */
    public double cond;
    /** Tilt angle. This is the angle between the line normal to the image
     * plane and the line normal to the ellipse, assuming the ellipse is the
     * result of the orthographic projection of a tilted circle.
     */
    public double tilt=0;
    /**
     * ID tag. Useful during RANSAC to identify individual ellipses after
     * list has been shuffled.
     */
    public int tag=-1;


    public Ellipse(){};

    public Ellipse(double[] p){

        a = p[0];
        b = p[1];
        c = p[2];
        d = p[3];
        e = p[4];
        f = p[5];

        //+++ Ellipse specified in terms of a general conic +++//
        //
        // a*x**2 + b*x*y + c*y**2 + d*x + e*y + f
        //

        // Transform coordinates to remove xy cross term
        double phi = 0.5*Math.atan2(b,(c-a));

        // Intermediate variables
        double A = a*Math.cos(phi)*Math.cos(phi)-b*Math.sin(phi)*Math.cos(phi)
                        + c*Math.sin(phi)*Math.sin(phi);
        double B = a*Math.sin(phi)*Math.sin(phi)+b*Math.sin(phi)*Math.cos(phi)
                        + c*Math.cos(phi)*Math.cos(phi);
        double C = d*Math.cos(phi) - e*Math.sin(phi);
        double D = d*Math.sin(phi) + e*Math.cos(phi);
        double F = f;

        //+++ Conic in terms of transformed coefficients +++//
        //
        // A*x**2 + B*y**2 + C*x + D*y + F = 0
        //
        // Get semi-axis lengths and centre in rotated frame

        r_a = Math.max(Math.sqrt((C*C/(4*A) + D*D/(4*B) - F)/A),
                            Math.sqrt((C*C/(4*A) + D*D/(4*B) - F)/B));
        r_b = Math.min(Math.sqrt((C*C/(4*A) + D*D/(4*B) - F)/A),
                            Math.sqrt((C*C/(4*A) + D*D/(4*B) - F)/B));

        tilt = Math.toDegrees(Math.acos(r_b/r_a));

        double x = -C/(2.0*A);     // Centre in rotated frame
        double y = -D/(2.0*B);

        // Transform centre back to original frame
        x0 = x*Math.cos(-phi)  - y*Math.sin(-phi);
        y0 = x*Math.sin(-phi)  + y*Math.cos(-phi);

        //+++ Now set theta - angle between major axis and x direction   +++//
        //+++ Formulae from http://mathworld.wolfram.com/Ellipse.html    +++//
        //+++ except for the absolute values that i found were necessary +++//
        if(Math.abs(b)<1E-9 && Math.abs(a)<Math.abs(c))
            theta=0;
        if(Math.abs(b)<1E-9 && Math.abs(a)>Math.abs(c))
            theta=0.5*Math.PI;
        if(Math.abs(b)>1E-9 && Math.abs(a)<Math.abs(c))
            theta=0.5*Math.atan(b/(a-c));
        if(Math.abs(b)>1E-9 && Math.abs(a)>Math.abs(c))
            theta=0.5*Math.PI + 0.5*Math.atan(b/(a-c));

    }


    public double[] getCentre(){ return new double[]{x0,y0};}


    @Override
    /** Return matlab command for plotting ellipse */
    public String toString(){
        return ("ellipse("+r_a+","+r_b+","+theta+","+x0+","+y0+")");
    }

    /** Return gnuplot command for plotting ellipse */
    public String toStringG(int obj){
        DecimalFormat xpx = new DecimalFormat("0.0");
        return ("set obj "+obj+" ellipse center "+xpx.format(x0)+", "+xpx.format(y0)+"  size "+xpx.format(2*r_a)+", "+xpx.format(2*r_b)+"  angle "+xpx.format(Math.toDegrees(theta))+" behind fs empty bo -1 lw 1 # tilt = "+xpx.format(tilt) + ", cond = "+cond);
    }


    /**
     * There now follows the various fitting methods that implement the
     * ellipse fitting algorithm.
     */



    /**
     * Implementation of "Numerically Stable Direct Least Squares Fitting of Ellipses"
     * This method will fit an arbitrary ellipse to a set of points.
     * @param track
     * @return
     */
    public static Ellipse fitEllipse(ArrayList<double[]> track){

        //+++ Read coordinate pairs into two arrays +++//
        double[] x = new double[track.size()];
        double[] y = new double[track.size()];

        for(int c=0; c<track.size(); c++){
            double[] xy = track.get(c);
            x[c] = xy[1];
            y[c] = xy[2];

        }

        // Calculate mean and range of each coordinate
//        double x_m = mean(x);
//        double y_m = mean(y);
//        double x_s = range(x)/2.0;
//        double y_s = range(y)/2.0;

        // Now normalise positions
//        for(int c=0; c<track.size(); c++){
//            x[c] = (x[c] - x_m)/x_s;
//            y[c] = (y[c] - y_m)/y_s;
//        }


        //+++ Split design matrix into quadratic and linear parts +++//

        //+++ Build design matrix +++//
        double[][] d1 = new double[x.length][3];
        double[][] d2 = new double[x.length][3];

        for(int i=0; i<x.length; i++){

            d1[i][0] = x[i]*x[i];
            d1[i][1] = x[i]*y[i];
            d1[i][2] = y[i]*y[i];

            d2[i][0] = x[i];
            d2[i][1] = y[i];
            d2[i][2] = 1.0;
        }

        Matrix D1 = new Matrix(d1);
        Matrix D2 = new Matrix(d2);

        //+++ Scatter matrix is split into four sub-matrices +++//
        Matrix S1 = D1.transpose().times(D1);
        Matrix S2 = D1.transpose().times(D2);
        Matrix S3 = D2.transpose().times(D2);

        //+++ Build constraint matrix +++//
        Matrix C1     = new Matrix(new double[][]{{0, 0, 2},
                                                  {0,-1, 0},
                                                  {2, 0, 0}});
        // Inverse is used - specify this directly
        Matrix C1_inv = new Matrix(new double[][]{{0,   0, 0.5},
                                                  {0,  -1,   0},
                                                  {0.5, 0,   0}});

        // Intermediate step...
        Matrix T = S3.inverse().times(S2.transpose()).times(-1);

        // Start to build reduced scatter matrix M...
        Matrix M = S1.plus(S2.times(T));

        // Premultiply by C1^{-1}
        M = C1_inv.times(M);

        // Eigenvalue decomposition of M:
        EigenvalueDecomposition evd = new EigenvalueDecomposition(M);

        // Minimal non-negative eigenvalue is solution. However, for almost
        // singular matrix M, numerical instability can lead to wrong answers.
        // This trick gets around that.
        double[] values = new double[3];

        Matrix evecs = evd.getV();

        for(int evec = 0; evec<3; evec++){

            // Build column vector a1 from eigenvector evec
            Matrix a1 = evecs.getMatrix(0, 2, evec, evec);

            // Get aT*C*a
            values[evec] = a1.transpose().times(C1.times(a1)).get(0, 0);

        }

        // Look for single positive value of aT*C*a
        int index = -1;
        int N=0;
        for(int e=0; e<3; e++){
            //+++ Detect positive eigenvalues +++//
            if (values[e] > 0){
                index = e;
                N++;
            }
        }

        if(N!=1){
            System.err.println("Ellipse solution failed! "+N+
                                " positive eigenvalues.");
            return new Ellipse();
        }

        // Now get parameter vector a1
        Matrix a1 = evecs.getMatrix(0, 2, index, index);

        // Get remaining parameters a2
        Matrix a2 = T.times(a1);

        // Build array from these matrices
        double[] a = new double[]{a1.get(0, 0), a1.get(1, 0), a1.get(2, 0),
                                  a2.get(0, 0), a2.get(1, 0), a2.get(2, 0)};

        // Now un-normalise parameters
//        double A = a[0]/(x_s*x_s);
//        double B = a[1]/(x_s*y_s);
//        double C = a[2]/(y_s*y_s);
//        double D = -2*a[0]*x_m/(x_s*x_s) - a[1]*y_m/(x_s*y_s) + a[3]/x_s;
//        double E = -2*a[2]*y_m/(y_s*y_s) - a[1]*x_m/(x_s*y_s) + a[4]/y_s;
//        double F = a[0]*x_m*x_m/(x_s*x_s) + a[1]*x_m*y_m/(x_s*y_s) + a[2]*y_m*y_m/(y_s*y_s) - a[3]*x_m/x_s - a[4]*y_m/y_s + a[5];

        //Ellipse out = new Ellipse(new double[]{A,B,C,D,E,F});
        Ellipse out = new Ellipse(a);

        out.cond = M.cond();

        return out;
    }



    /**
     * This method follows a similar approach to fitEllipse(), except that
     * the ellipse form used in this method is restricted so that it's
     * centre lies on the y-axis and it has zero rotation of the semi axes.
     * @param track
     * @return
     */
    public static Ellipse fitConstrainedEllipse(ArrayList<double[]> track){

        //+++ Read coordinate pairs into two arrays +++//
        double[] x = new double[track.size()];
        double[] y = new double[track.size()];

        for(int c=0; c<track.size(); c++){
            double[] xy = track.get(c);
            x[c] = xy[1];
            y[c] = xy[2];
        }

//        // Calculate mean and range of each coordinate
//        double x_m = mean(x);
//        double y_m = mean(y);
//        double x_s = range(x)/2.0;
//        double y_s = range(y)/2.0;
//
//        // Now normalise positions
//        for(int c=0; c<track.size(); c++){
//            x[c] = (x[c] - x_m)/x_s;
//            y[c] = (y[c] - y_m)/y_s;
//        }

        //+++ Build design matrix +++//
        double[][] d1 = new double[x.length][2];
        double[][] d2 = new double[x.length][2];

        for(int i=0; i<x.length; i++){

            d1[i][0] = x[i]*x[i];
            d1[i][1] = y[i]*y[i];

            d2[i][0] = y[i];
            d2[i][1] = 1.0;
        }

        Matrix D1 = new Matrix(d1);
        Matrix D2 = new Matrix(d2);

        //+++ Scatter matrix is split into four sub-matrices +++//
        Matrix S1 = D1.transpose().times(D1);
        Matrix S2 = D1.transpose().times(D2);
        Matrix S3 = D2.transpose().times(D2);

        //+++ 2x2 constraint matrix +++//
        double[][] c = new double[][]{{ 0, -2},
                                      { -2, 0}};

        //+++ Build constraint matrix +++//
        Matrix C1     = new Matrix(new double[][]{{ 0, 2},{ 2, 0}});
        // Inverse is used - specify this directly
        Matrix C1_inv = new Matrix(new double[][]{{ 0, 0.5},{0.5, 0}});

        // Intermediate step...
        Matrix T = S3.inverse().times(S2.transpose()).times(-1);

        // Start to build reduced scatter matrix M...
        Matrix M = S1.plus(S2.times(T));

        // Premultiply by C1^{-1}
        M = C1_inv.times(M);

        // Eigenvalue decomposition of M:
        EigenvalueDecomposition evd = new EigenvalueDecomposition(M);

        // Minimal non-negative eigenvalue is solution. However, for almost
        // singular matrix M, numerical instability can lead to wrong answers.
        // This trick gets around that.
        double[] values = new double[2];

        Matrix evecs = evd.getV();

        for(int evec = 0; evec<2; evec++){

            // Build column vector a1 from eigenvector evec
            Matrix a1 = evecs.getMatrix(0, 1, evec, evec);

            // Calculate aT*C*a
            values[evec] = a1.transpose().times(C1.times(a1)).get(0, 0);

        }

        // Look for single positive value of aT*C*a
        int index = -1;
        int N=0;
        for(int e=0; e<2; e++){
            //+++ Detect positive values +++//
            if (values[e] > 0){
                index = e;
                N++;
            }
        }

        if(N!=1){
            System.err.println("Ellipse solution failed! "+N+
                               " positive eigenvalues.");
            return new Ellipse();
        }

        // Now get parameter vector a1
        Matrix a1 = evecs.getMatrix(0, 1, index, index);

        // Get remaining parameters a2
        Matrix a2 = T.times(a1);

        // Build array from these matrices
        double[] a = new double[]{a1.get(0, 0), 0.0,          a1.get(1, 0),
                                  0.0,          a2.get(0, 0), a2.get(1, 0)};

        // Now un-normalise parameters
//        double A = a[0]/(x_s*x_s);
//        double B = a[1]/(x_s*y_s);
//        double C = a[2]/(y_s*y_s);
//        double D = -2*a[0]*x_m/(x_s*x_s) - a[1]*y_m/(x_s*y_s) + a[3]/x_s;
//        double E = -2*a[2]*y_m/(y_s*y_s) - a[1]*x_m/(x_s*y_s) + a[4]/y_s;
//        double F = a[0]*x_m*x_m/(x_s*x_s) + a[1]*x_m*y_m/(x_s*y_s) + a[2]*y_m*y_m/(y_s*y_s) - a[3]*x_m/x_s - a[4]*y_m/y_s + a[5];

        //Ellipse out = new Ellipse(new double[]{A,B,C,D,E,F});
        Ellipse out = new Ellipse(a);

        out.cond = M.cond();

        return out;

    }



    /**
     * Fit an ellipse to a set of points, with the constraint that the centre 
     * lies along the line y*sin(theta) + x*cos(theta) = r.
     *
     * @param track
     * @param theta Measured from positive x axis. Ranges [-PI:+PI]
     * @param r
     * @return
     */
    public static Ellipse fitEllipse(ArrayList<double[]> track, 
                                     double theta,
                                     double r){

        //+++ Origin of transformed frame is defined by the point +++//
        //+++ of closest approach between line and origin         +++//
        double x_0 = r*Math.cos(theta);
        double y_0 = r*Math.sin(theta);

        //+++ Now transform points as they are read into arrays.  +++//
        //+++ Note that reverse transformation is used because we +++//
        //+++ want to shift the coordinate axes not the points.   +++//
        ArrayList<double[]> track_prime = new ArrayList<double[]>();

        double x[] = new double[track.size()];
        double y[] = new double[track.size()];  // Transformed coordinates

        for(int c=0; c<track.size(); c++){

            double[] xy = track.get(c);

            // Rotation and translation
            x[c] =  (xy[1]-x_0)*Math.cos(theta) + (xy[2]-y_0)*Math.sin(theta);
            y[c] = -(xy[1]-x_0)*Math.sin(theta) + (xy[2]-y_0)*Math.cos(theta);

            track_prime.add(new double[]{xy[0],x[c],y[c]});
        }

        // Calculate mean and range of each coordinate
//        double x_m = mean(x);
//        double y_m = mean(y);
//        double x_s = range(x)/2.0;
//        double y_s = range(y)/2.0;
//
//        // Now normalise positions
//        for(int c=0; c<track.size(); c++){
//
//            double[] fxy = track_prime.get(c);
//
//            fxy[1] = (fxy[1] - x_m)/x_s;
//            fxy[2] = (fxy[2] - y_m)/y_s;
//
//            track_prime.set(c, fxy);
//        }


        Ellipse e_prime = fitConstrainedEllipse(track_prime);

        //+++ Ellipse in rotated frame is given by +++//
        //
        // A*x**2 + C*y**2 + E*y + F = 0
        //
        // i.e. B=D=0 in this frame.
        //
        // Recover full conic parameters in original frame
        double ct = Math.cos(theta);
        double st = Math.sin(theta);
        double A = e_prime.a;
        double B = e_prime.c;
        double C = e_prime.e;
        double D = e_prime.f;

        double[] a = new double[6];

        a[0] = A*ct*ct + B*st*st;
        a[1] = 2*A*ct*st - 2*B*ct*st;
        a[2] = A*st*st + B*ct*ct;
        a[3] =  -2*A*x_0*ct*ct - 2*A*y_0*st*ct - 2*B*x_0*st*st
                    + 2*B*y_0*st*ct - C*st;
        a[4] =  -2*A*x_0*ct*st - 2*A*y_0*st*st
                    + 2*x_0*B*st*ct - 2*y_0*B*ct*ct + C*ct;
        a[5] = A*x_0*x_0*ct*ct + 2*A*x_0*y_0*ct*st + A*y_0*y_0*st*st
                    + B*x_0*x_0*st*st - 2*x_0*y_0*B*st*ct + B*y_0*y_0*ct*ct
                    + C*x_0*st - C*y_0*ct + D;

        // un-normalise parameters
//        A = a[0]/(x_s*x_s);
//        B = a[1]/(x_s*y_s);
//        C = a[2]/(y_s*y_s);
//        D = -2*a[0]*x_m/(x_s*x_s) - a[1]*y_m/(x_s*y_s) + a[3]/x_s;
//        double E = -2*a[2]*y_m/(y_s*y_s) - a[1]*x_m/(x_s*y_s) + a[4]/y_s;
//        double F = a[0]*x_m*x_m/(x_s*x_s) + a[1]*x_m*y_m/(x_s*y_s) + a[2]*y_m*y_m/(y_s*y_s) - a[3]*x_m/x_s - a[4]*y_m/y_s + a[5];

        //+++ Create ellipse with transformed parameters +++//
        //Ellipse out = new Ellipse(new double[]{A,B,C,D,E,F});
        Ellipse out = new Ellipse(a);

        out.cond = e_prime.cond;

        return out;
    }


    /**
     * Get the mean of an array of values.
     * @param x
     * @return
     */
    public static double mean(double[] x){

        double sum = 0;

        for(int i=0; i<x.length; i++)
            sum += x[i];

        return sum/(double)x.length;
    }


    /**
     * Get the range of values in an array.
     * @param x
     * @return
     */
    public static double range(double[] x){

        double MIN = x[0], MAX = x[0];

        for(int i=0; i<x.length; i++){

            if(x[i] < MIN) MIN = x[i];
            if(x[i] > MAX) MAX = x[i];

        }

        return (MAX - MIN);

    }



}
