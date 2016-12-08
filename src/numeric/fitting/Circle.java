/**
 * Copyright (c) 2011 University of Dundee
 *
 * Name:
 * Circle.java
 *
 * Purpose:
 * This class represents geometric circle objects. It also provides methods to
 * fit circles to sets of 2D coordinates. These include direct least squares
 * fitting based on minimization of algebraic distance, least squares fitting
 * with RANSAC outlier rejection, and non-linear fitting using Levenberg
 * Marquardt algorithm to minimize the sum of the geometric distance of points
 * from the circle.
 * 
 * Language:
 * Java
 *
 * Author:
 * Nicholas Rowell (nickrowell@computing.dundee.ac.uk)
 *
 */

package numeric.fitting;

import java.util.ArrayList;
import java.util.List;
import java.util.Collections;
import java.awt.image.BufferedImage;

import Jama.*;
import numeric.minimisation.nllsq.algo.LevenbergMarquardt;

public class Circle extends LevenbergMarquardt {

    /**
     * Size of finite steps used in Levenberg-Marquardt algorithm.
     */
    private static final double dr=1.0, dx0=1.0, dy0=1.0;
    
    /**
     * Threshold for zero-finding in floating point arithmetic.
     */
    private static final double EPSILON = 1e-9;
    
    /**
     * Radius.
     */
    double r;
    
    /**
     * x coordinate of centre.
     */
    double x0;
    
    /**
     * y coordinate of centre.
     */
    double y0;

    /**
     * Parameter a of the conic section parameterization:
     * 
     * a(x*x + y*y) + dx + ey + f = 0
     */
    double a;
    /**
     * Parameter d of the conic section parameterization:
     * 
     * a(x*x + y*y) + dx + ey + f = 0
     */
    double d;
    
    /**
     * Parameter e of the conic section parameterization:
     * 
     * a(x*x + y*y) + dx + ey + f = 0
     */
    double e;
    
    /**
     * Parameter f of the conic section parameterization:
     * 
     * a(x*x + y*y) + dx + ey + f = 0
     */
    double f;
    
    /**
     * List of points to fit circle to.
     */
    List<double[]> data = new ArrayList<double[]>();

    /**
     * List of inliers found by RANSAC. Subset of data.
     */
    List<double[]> inliers = new ArrayList<double[]>();

    /**
     * List of outliers found by RANSAC. Subset of data.
     */
    List<double[]> outliers = new ArrayList<double[]>();

    /**
     * Default constructor.
     */
    public Circle() {
    	this(1.0, 0.0, 0.0);
    }

    /**
     * Constructor based on conventional circle parameters. The conic section
     * parameters are derived from these.
     * 
     * @param R
     * 	The radius.
     * @param X0
     * 	The X coordinate of the centre.
     * @param Y0
     * 	The Y coordinate of the centre.
     */
    public Circle(double R, double X0, double Y0) {

        r = R;
        x0 = X0;
        y0 = Y0;

        // Set conic section parameters. These are not unique; here a is
        // arbitrarily chosen to be 1.
        a = 1.0;
        d = -2*x0;
        e = -2*y0;
        f = x0*x0 + y0*y0 - r*r;
    }

    /**
     * Constructor based on conic section parameters. The conventional
     * parameters are derived from these.
     * 
     * @param a
     * 	Parameter a of the conic section parameterization.
     * @param d
     * 	Parameter d of the conic section parameterization.
     * @param e
     * 	Parameter e of the conic section parameterization.
     * @param f
     * 	Parameter f of the conic section parameterization.
     */
    public Circle(double a, double d, double e, double f){

        this.a = a;
        this.d = d;
        this.e = e;
        this.f = f;

        x0 = -this.d/(2*this.a);
        y0 = -this.e/(2*this.a);
        r  = Math.sqrt((this.d*this.d + this.e*this.e)/(4*this.a*this.a) - this.f/this.a);
    }

    /**
     * Get the radius.
     * @return
     * 	The radius.
     */
    public double getR() {
    	return r;
    }
    
    /**
     * Get the X coordinate of the centre.
     * @return
     * 	The X coordinate of the centre.
     */
    public double getX0() {
    	return x0;
    }

    /**
     * Get the Y coordinate of the centre.
     * @return
     * 	The Y coordinate of the centre.
     */
    public double getY0() {
    	return y0;
    }

    /**
     * Direct least squares fitting of a {@link Circle} to a set of 2D points. All
     * points are added to the {@link Circle}'s internal data and inlier lists.
     * The least squares solution is that which minimises the algebraic form of
     * the cost function (rather than the geometric form).
     * 
     * @param points
     * 	The {@link List} of (x,y) points to fit the {@link Circle} to.
     * @return
     * 	The best fitting circle (algebraic, not geometric).
     */
    public static Circle fitCircle(List<double[]> points){

        // Use this to make a Circle
        Circle out = new Circle();

        // Copy points from argument to internal data and inliers lists.
        for(double[] point: points){
            out.data.add(point);
            out.inliers.add(point);
        }

        // Now use points to calculate parameters of Circle. If false is
        // returned by call to calculateCircle, then no parameters could be
        // determined and null should be returned from this method.

        return out.calculateCircle() ? out : null;
    }

    /**
     * Least squares fitting of a circle to a set of datapoints, with RANSAC
     * outlier rejection to separate data points into inliers and outliers.
     * This is useful when there are catastrophic outliers included in dataset,
     * points not drawn from the known distribution of the data.
     *
     * @param points
     * 	List of all data points
     * @param consensusThreshold
     * 	Threshold required for Circle acceptance
     * @param outlierThreshold
     * 	Distance in pixels for a point to be an outlier.
     * @param setSize
     * 	Number of samples selected for initial model.
     * @return
     * 	The best fitting {@link Circle}.
     */
    public static Circle fitCircleWithRANSAC(List<double[]> points, double consensusThreshold,
    		double outlierThreshold, int setSize, int maxLoop) {

    	// Sanity checks
        if(points.size()<setSize || setSize<3) {
            throw new RuntimeException("Insufficient points for Circle!");
        }
        if(consensusThreshold < 0 || consensusThreshold > 1) {
            throw new RuntimeException("Consensus threshold below zero or " +
                                       "greater than one!");
        }

        // Main RANSAC loop

        int N = 0;

        while(N++ < maxLoop) {

            // Shuffle data points list to randomize order
            Collections.shuffle(points);

            // Get Circle fitted to initial random set picked from start of
            // points List.
            Circle ransac = fitCircle(points.subList(0, setSize));

            // If Circle is null, skip rest of loop. This can occur if two or
            // more of the initial set members are too close and circle is not
            // sufficiently constrained.
            if(ransac==null) {
            	continue;
            }

            // Now loop over remaining data points and divide into inliers
            // and outliers.
            for(double[] data : points.subList(setSize, points.size())) {

                // Add all points to internal data list
                ransac.data.add(data);

                // Then also add each to either inliers or outliers depending
                // on geometric distance from the circumference. Note that the
                // geometric distance is negative if point lies inside circle.
                if(Math.abs(ransac.getGeometricDistance(data)) < outlierThreshold) {
                    ransac.inliers.add(data);
                }
                else {
                    ransac.outliers.add(data);
                }

            }

            // Now update model using all inliers to constrain parameters
            ransac.calculateCircle();

            // Now check consensus level found by the current Circle
            double consensus = (double)ransac.inliers.size()/(double)ransac.data.size();

            // Consensus threshold exceeded - model is a good enough fit
            if(consensus > consensusThreshold) {
                return ransac;
            }

        }

        // No sufficiently large consensus found in allowed number of loops.
        return null;
    }

    /**
     * Use current inliers to get the least squares solution for the parameters
     * of the {@link Circle}. The algorithm minimizes the algebraic distance between
     * the points and the {@link Circle}, i.e. finds (a,d,e,f) which minimizes the
     * squared sum of a(x_i*x_i + y_i*y_i) + dx_i + ey_i + f = 0 over all i.
     *
     * @return boolean
     * 	Specifies whether inliers set could be used to calculate {@link Circle}
     * parameters (true) or not (false).
     */
    public boolean calculateCircle() {

        if(inliers.size()<3) {
            throw new RuntimeException("Need 3 or more points for Circle!");
        }

        // Build matrix X as 2D double array initially
        double[][] x = new double[inliers.size()][4];

        int N=0;

        for(double[] xy : inliers){
            x[N][0] = xy[0]*xy[0] + xy[1]*xy[1];
            x[N][1] = xy[0];
            x[N][2] = xy[1];
            x[N][3] = 1.0;
            N++;
        }

        // Convert to Matrix. Vector of algebraic distances of each point from 
        // circle is given by X*(a,d,e,f)^T. The solution is found by squaring
        // this and setting derivative equal to zero, i.e.
        //
        // X^T * X * (a,d,e,f)^T = 0
        //
        Matrix X = new Matrix(x);

        // Make design matrix X^T * X
        Matrix D = X.transpose().times(X);
        
        // Get SVD
        SingularValueDecomposition svd = D.svd();
                
        // Get right singular vector corresponding to lowest singular value
        Matrix rsv = svd.getV().getMatrix(new int[]{0,1,2,3}, new int[]{3});

        // Check that second lowest singular value is not zero. This would
        // indicate that no solution has been found, i.e. multiple points lie
        // at the same coordinates.
        if(svd.getS().get(2, 2)<EPSILON) {
            // Failure!
            return false;
        }

        // Set new Circle parameters
        a = rsv.get(0,0);
        d = rsv.get(1,0);
        e = rsv.get(2,0);
        f = rsv.get(3,0);

        x0 = -d/(2*a);
        y0 = -e/(2*a);
        r  = Math.sqrt((d*d + e*e)/(4*a*a) - f/a);

        // Success!
        return true;
    }

    /**
     * Get geometric distance of point (X[0],X[1]) to the circumference of this
     * Circle. Value is positive if the point lies outside the Circle, and
     * negative if it lies inside.
     */
    public double getGeometricDistance(double[] X) {

        // distance of point from circle centre
        double l = Math.sqrt((X[0]-x0)*(X[0]-x0) + (X[1]-y0)*(X[1]-y0));

        // Subtract radius of circle to get geometric distance
        l -= r;

        return l;
    }

    /**
     * Draw a {@link Circle} onto a {@link BufferedImage}.
     * @param image
     * 	The {@link BufferedImage}.
     * @param rgb
     * 	The RGB colour to draw the circle.
     */
    public void drawCircle(BufferedImage image, int rgb) {

        // Angular step between points on circumference that lie one
        // pixel apart
        double ang = 2*Math.asin(0.5/r);

        // Loop round circumference of circle:
        for (double theta = 0; theta < 2 * Math.PI; theta += ang) {

            double i = x0 + r * Math.sin(theta);
            double j = y0 + r * Math.cos(theta);

            try {
            	image.setRGB((int)i, (int)j, rgb);
            }
            catch (ArrayIndexOutOfBoundsException aioobe) {
            	// Ignored
            }

        }
        return;
    }

    /**
     * Get the coordinate of the point on the {@link Circle} that
     * lies closest to the given (x,y) point.
     * @param xy
     * 	The (x,y) point.
     * @return
     * 	The coordinate of the point on the {@link Circle} that
     * lies closest to the given (x,y) point.
     */
    private double[] getClosestPointOnCircle(double[] xy) {
        
        // Angle towards point XY from circle centre
        double theta = Math.atan2(xy[1] - y0, xy[0] - x0);
        
        // Closest point on circumference is point lying distance r along 
        // this direction
        return new double[]{x0+r*Math.cos(theta), y0+r*Math.sin(theta)};   
    }

    
    
    
    
    /** Levenberg-Marquardt algorithm implemented methods */

    /** Required by LMA */
    public int getDataN() {
    	return 2*inliers.size();
    }
    
    /** Required by LMA */
    public int getParametersN() {
    	return 3;
    }
    
    /** Required by LMA */
    public Matrix getModel() {

        Matrix model = new Matrix(getDataN(),1);
        
        // For each inlying data point, get the closest point on
        // the circumference
        for(int p=0; p<inliers.size(); p++) {
    
            // Get closest point on circumference to inlying point p.
            double[] xy = getClosestPointOnCircle(inliers.get(p));
            
            model.set(2*p, 0, xy[0]);    // x coordinate of model point
            model.set(2*p+1, 0, xy[1]);  // y coordinate of model point
        }

        return model;
    }

    public Matrix getData() {

        Matrix DATA = new Matrix(getDataN(),1);
        
        // For each inlying data point, get the closest point on
        // the circumference
        
        for(int p=0; p<inliers.size(); p++){
    
            // Get coordinates of this datapoint
            double[] xy = inliers.get(p);
            
            DATA.set(2*p, 0, xy[0]);    // x coordinate of model point
            DATA.set(2*p+1, 0, xy[1]);  // y coordinate of model point
            
        }

        return DATA;
    }
    
    
    public boolean updateData(Matrix delta) {
        
        for(int p=0; p<inliers.size(); p++){
    
            // Get coordinates of this datapoint
            double[] xy = inliers.get(p);
            
            xy[0] += delta.get(2*p, 0);     // update x coordinate
            xy[1] += delta.get(2*p+1, 0);   // update y coordinate
            
            // Store updated point back to list
            inliers.set(p, xy);
            
        }

        return true;
    }
    
    // Get covariance matrix for data points. This would be set by the user
    // when the data points are generated, or something like that.
    public Matrix getCovariance(){
    
        Matrix COVAR = new Matrix(getDataN(),getDataN());
        
        // For each inlying data point, get the closest point on
        // the circumference
        
        for(int p=0; p<inliers.size(); p++){
    
            // Diagonal matrix
            COVAR.set(2*p, 2*p, 1.0);
            COVAR.set(2*p+1, 2*p+1, 1.0);
            
        }

        return COVAR;   
        
    }
    
    

    /** Required by LMA
     *
     * Parameters are organised in the adjustment matrix delta in the following
     * manner:
     *
     * 0-1: Coordinates of circle centre
     *   2: Radius
     *
     */
//    public boolean updateParameters(Matrix delta){
//
//        x0 += delta.get(0, 0);
//        y0 += delta.get(1, 0);
//        r  += delta.get(2, 0);
//
//        a = 1.0;
//        d = -2*x0;
//        e = -2*y0;
//        f = x0*x0 + y0*y0 - r*r;
//        
//        return true;
//    }
    
    /**
     * Override setParameters() method.
     *  
     */
    @Override
    public void setParameters(Matrix delta){
        
        x0 = delta.get(0, 0);
        y0 = delta.get(1, 0);
        r  = delta.get(2, 0);

        a = 1.0;
        d = -2*x0;
        e = -2*y0;
        f = x0*x0 + y0*y0 - r*r;
        
        return;
    }
    
    
    @Override
    public Matrix getParameters(){
        return new Matrix(new double[][]{{x0},{y0},{r}});
    }

    @Override
    public Matrix getParametersSteps(){
        return new Matrix(new double[][]{{dx0},{dy0},{dr}});
    } 
    
    
    /**
     * Required by LMA.
     */
    public void printParameters(){

        System.out.println("LMA: "+toString());

    }

    @Override
    public String toString(){
        return "Circle: (x0,y0) =  ("+x0+","+y0+"), r =  "+r+
                ", "+data.size()+" points and "+inliers.size()+" inliers.";
    }

}