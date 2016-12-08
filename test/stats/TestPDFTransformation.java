package stats;

import Jama.Matrix;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import data.Histogram2D;

/**
 * This class tests the transformation law for multi-variate probability
 * density functions.
 * 
 * @author nickrowell
 */
public class TestPDFTransformation{
    
    static Histogram2D pXY = new Histogram2D(-10,10,1,
                                             -10,10,1,
                                             true);
    
    static Histogram2D pUV = new Histogram2D(-10,10,1,
                                             -10,10,1,
                                             true);
    
    public static void main(String[] args) throws IOException{
    
        // Draw lots of random points from P(X,Y)
        for(int i=0; i<1000000; i++){
            
            double[] XY = drawXY();
            
            // Add to (X,Y) histogram
            pXY.add(XY[0], XY[1]);
            // Transform to (U,V) and add to (U,V) histogram
            pUV.add(getU(XY[0],XY[1]), getV(XY[0],XY[1]));
        
        }
        
        // Output file for results:
        BufferedWriter out = new BufferedWriter(new FileWriter(new File("results/tests/transformPDF_in")));
        
        // Print input pdf P(X,Y)
        out.write(pXY.print(true));
        out.flush();
        
        // Switch writer to transformed PDF file
        out = new BufferedWriter(new FileWriter(new File("results/tests/transformPDF_out")));
        
        // Print output pdf P(U,V)
        out.write(pUV.print(true));
        out.flush();
        
    }
    
    /**
     * Draw values of random variables X and Y.
     * 
     * X and Y are distributed as bivariate Gaussian.
     * 
     * @return 
     */
    private static double[] drawXY(){
    
        Matrix cov  = new Matrix(new double[][]{{4,0},{0,1}});
        Matrix mean = new Matrix(new double[][]{{0},{0}});
        
        Matrix rand = Statistics.drawRandomVector(cov, mean);
    
        return new double[]{rand.get(0,0), rand.get(1,0)};
    }
    
    // For simple rotation of the input PDF, this is the angle
    static double angle = Math.toRadians(90);
    
    /**
     * U = U(X,Y)
     */
    private static double getU(double X, double Y){
        return Math.cos(angle) * X + Math.sin(angle) * Y;
    }
    
    /**
     * V = V(X,Y)
     */
    private static double getV(double X, double Y){
        return -Math.sin(angle) * X + Math.cos(angle) * Y;
    }
    
    
    
}