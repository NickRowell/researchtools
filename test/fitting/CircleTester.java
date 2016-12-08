package fitting;

import java.util.ArrayList;
import java.util.Random;

import fitting.Circle;

/**
 *
 * @author nickrowell
 */
public class CircleTester {


    public static void main(String[] args){

        ArrayList<double[]> points = new ArrayList<double[]>();

        // Randomly add points on the circumference of a Circle to list
        int N=6;

        Random noise = new Random();

        // True circle parameters
        double x0 = 1, y0 = 0, r = 10;

        // STD of Gaussian noise added to each coordinate
        double std = 0.1;

        for(int n=0; n<N; n++){

            // Select random angle, measured clockwise from y axis
            double angle = Math.random()*2.0*Math.PI;

            double[] point = new double[]{x0 + r * Math.sin(angle),
                                          y0 + r * Math.cos(angle)};

            // Now add noise to each coordinate
            point[0] += noise.nextGaussian()*std;
            point[1] += noise.nextGaussian()*std;

            points.add(point);
            
        }

        // Fit a circle to all datapoints
        Circle test1 = Circle.fitCircle(points);

        if(test1!=null)
            System.out.println(test1.toString());
       
        
        // Fit a circle to data using RANSAC to identify outliers
        Circle test2 = Circle.fitCircleWithRANSAC(points, 0.5, 1, 5, 5000);

        if(test2!=null)
            System.out.println(test2.toString());

        // Now use LM algorithm to improve fit for each circle
        if(test1!=null) test1.fit(1000, false);
        if(test2!=null) test2.fit(1000, false);

        System.out.println("Improved fits:");
        if(test1!=null) System.out.println(test1.toString());
        if(test2!=null) System.out.println(test2.toString());
        
        System.out.println("\nCovariance matrices on parameters by covariance propagation:");
        if(test1!=null) test1.getFourthOrderCovariance().print(5, 5);
        if(test2!=null) test2.getFourthOrderCovariance().print(5, 5);
        
        System.out.println("Covariance matrices on parameters by Gnuplot method:");
        if(test1!=null) test1.getParameterCovariance().print(5, 5);
        if(test2!=null) test2.getParameterCovariance().print(5, 5);
    }

}
