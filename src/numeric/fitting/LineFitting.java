/**
 * Copyright (c) 2011 University of Dundee
 *
 * Name:
 * LineFitting.java
 *
 * Purpose:
 * This class provides a method used to perform total least squares fitting
 * of a straight line to a bunch of data points.
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

public class LineFitting {

    /**
     * Total least squares fitting of straight line to points.
     *
     * @param points
     * @return
     */
    public static double[] tls(double[][] points){

        double x_M = 0, y_M = 0, xy_M = 0, x2_M = 0, y2_M = 0;

        for(int p=0; p<points.length; p++)
        {
            x_M  += points[p][0];
            y_M  += points[p][1];
            xy_M += points[p][0]*points[p][1];
            x2_M += points[p][0]*points[p][0];
            y2_M += points[p][1]*points[p][1];
        }

        x_M  /= points.length;
        y_M  /= points.length;
        xy_M /= points.length;
        x2_M /= points.length;
        y2_M /= points.length;

        double two_T = Math.atan2(2*(x_M*y_M - xy_M),  (x_M*x_M - x2_M) - (y_M*y_M - y2_M));
        double theta = two_T/2.0;
        double r     = x_M * Math.cos(theta) + y_M * Math.sin(theta);

        // Return parameters of normal representation of line
        return new double[]{r,theta};

    }

    /**
     * Total least squares fitting of straight line to points, with an
     * ArrayList interface so data can be dynamically expanded within
     * RANSAC loop without having to
     *
     * @param points
     * @return
     */
    public static double[] tls(ArrayList<double[]> points){

        double x_M = 0, y_M = 0, xy_M = 0, x2_M = 0, y2_M = 0;

        for(double[] p: points)
        {
            x_M  += p[0];
            y_M  += p[1];
            xy_M += p[0]*p[1];
            x2_M += p[0]*p[0];
            y2_M += p[1]*p[1];
        }

        x_M  /= points.size();
        y_M  /= points.size();
        xy_M /= points.size();
        x2_M /= points.size();
        y2_M /= points.size();

        double two_T = Math.atan2(2*(x_M*y_M - xy_M),  (x_M*x_M - x2_M) - (y_M*y_M - y2_M));
        double theta = two_T/2.0;
        double r     = x_M * Math.cos(theta) + y_M * Math.sin(theta);

        // Return parameters of normal representation of line
        return new double[]{r,theta};

    }


    /**
     * Weighted total least squares fitting of straight line to points.
     *
     * @param points
     * @return
     */
    public static double[] wtls(double[][] points, double[] weights)
    {

        double w = 0, xw_M = 0, yw_M = 0, xyw_M = 0, x2w_M = 0, y2w_M = 0;

        for(int p=0; p<points.length; p++){
            w     += 1.0 / weights[p];
            xw_M  += points[p][0] / weights[p];
            yw_M  += points[p][1] / weights[p];
            xyw_M += points[p][0]*points[p][1] / weights[p];
            x2w_M += points[p][0]*points[p][0] / weights[p];
            y2w_M += points[p][1]*points[p][1] / weights[p];
        }

        w     /= points.length;
        xw_M  /= points.length;
        yw_M  /= points.length;
        xyw_M /= points.length;
        x2w_M /= points.length;
        y2w_M /= points.length;

        double two_T = Math.atan2(2*(xw_M*yw_M/w - xyw_M),  (xw_M*xw_M/w - x2w_M) - (yw_M*yw_M/w - y2w_M));
        double theta = two_T/2.0;
        double r     = xw_M * Math.cos(theta)/w + yw_M * Math.sin(theta)/w;

        // Return parameters of normal representation of line
        return new double[]{r,theta};

    }

    /**
     * Weighted total least squares fitting of straight line to points.
     *
     * @param points
     * @return
     */
    public static double[] wtls(ArrayList<double[]> points)
    {

        double w = 0, xw_M = 0, yw_M = 0, xyw_M = 0, x2w_M = 0, y2w_M = 0;

        for(double[] p: points)
        {
            w     +=  1.0/p[2];
            xw_M  += p[0]/p[2];
            yw_M  += p[1]/p[2];
            xyw_M += p[0]*p[1]/p[2];
            x2w_M += p[0]*p[0]/p[2];
            y2w_M += p[1]*p[1]/p[2];
        }

        w     /= points.size();
        xw_M  /= points.size();
        yw_M  /= points.size();
        xyw_M /= points.size();
        x2w_M /= points.size();
        y2w_M /= points.size();

        double two_T = Math.atan2(2*(xw_M*yw_M/w - xyw_M),  (xw_M*xw_M/w - x2w_M) - (yw_M*yw_M/w - y2w_M));
        double theta = two_T/2.0;
        double r     = xw_M * Math.cos(theta)/w + yw_M * Math.sin(theta)/w;

        // Return parameters of normal representation of line
        return new double[]{r,theta};

    }

}