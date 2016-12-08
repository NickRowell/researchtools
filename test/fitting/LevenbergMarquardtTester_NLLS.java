package fitting;

import Jama.*;

import java.util.List;
import java.util.ArrayList;

import minimisation.nllsq.LevenbergMarquardt;

/**
 * Fit function to datapoints. The results can then be compared to the output
 * of the Gnuplot 'fit' algorithm to verify the behaviour of this 
 * implementation against Gnuplot.
 * @author nickrowell
 */
public class LevenbergMarquardtTester_NLLS extends LevenbergMarquardt{
    
    // Initial parameters of quadratic that will be fitted to data.
    // f(x) = A*x*x + B*x + C
    public double A=1;
    public double B=1;
    public double C=1;
    
    // Size of finite steps to be used with each parameter
    double dA=0.1, dB=0.1, dC=0.1;
    
    // Data points for the fit: List of (x,y) points. Only y values are
    // considered in the fit. The x values are used to compute the model.
    public List<double[]> data = new ArrayList<double[]>();
    
    
    public LevenbergMarquardtTester_NLLS(List<double[]> dat){
        
        data = dat;
    
    }
    
    
    public static void main(String[] args){
            
        List<double[]> input = new ArrayList<double[]>();
        
        input.add(new double[]{0, 0.0});
        input.add(new double[]{1, 1.0});
        input.add(new double[]{2, 4.0});
        input.add(new double[]{3,9.0});
        input.add(new double[]{4,16.0});

        LevenbergMarquardtTester_NLLS test = new LevenbergMarquardtTester_NLLS(input);
                        
        // Print out Gnuplot script that runs a fitting algorithm using same
        // data and input parameters:
        System.out.println("# Gnuplot datafile 'tmp':");        
         for (double[] ds : input) 
            System.out.println(ds[0] + " " + ds[1]);       
        System.out.println("# Gnuplot script:");
        System.out.println("A = "+test.A);
        System.out.println("B = "+test.B);
        System.out.println("C = "+test.C);
        System.out.println("f(x) = A*x**2 + B*x + C");
        System.out.println("plot f(x), 'tmp' ");
        System.out.println("fit f(x) 'tmp' via A,B,C \n\n\n");    
        

        // Perform Levenberg-Marquardt fitting algorithm on data & function
        test.fit(1000, true);
        
        System.out.println("Fitted parameters:");
        
        test.printParameters();
        
        System.out.println("Asymptotic standard error:");
        test.getAsymptoticStandardError().print(5,5);
        
        System.out.println("Parameter correlation:");
        test.getParameterCorrelation().print(5, 5);
        
        System.out.println("Parameter covariance:");
        test.getParameterCovariance().print(5, 5);       
        
        System.out.println("Parameter covariance by first order propagation:");
        test.getFourthOrderCovariance().print(5,5);
        
    }
    
    // Print parameters on successful update
    @Override
    public void updateCallback(){
        System.out.println("LMA: Successful update: (A,B,C) = "
                           + "("+A+","+B+","+C+")");
    }
    
    
    
    // Get predicted y coordinates given x coordinate and current parameters
    public Matrix getModel(){
        
        Matrix MOD = new Matrix(getDataN(),1);
        
        for(int p=0; p<data.size(); p++)
            MOD.set(p, 0, A*data.get(p)[0]*data.get(p)[0] +
                          B*data.get(p)[0] + C);

        return MOD;
    }
    
    // Get observed y coordinates
    public Matrix getData(){
        
        Matrix DATA = new Matrix(getDataN(),1);
        
        for(int p=0; p<data.size(); p++)
            DATA.set(p, 0, data.get(p)[1]);

        return DATA;
    }
    
    // Update data points temporarily
    public boolean updateData(Matrix delta){
    
        for(int p=0; p<data.size(); p++){
            data.get(p)[1] += delta.get(p,0);
        }
        
        return true;   
        
    }
    
    
    public Matrix getCovariance(){
        Matrix cov = new Matrix(new double[][]{{0.5,  0,  0,  0,  0},
                                               {  0,0.5,  0,  0,  0},
                                               {  0,  0,0.5,  0,  0},
                                               {  0,  0,  0,0.5,  0},
                                               {  0,  0,  0,  0,0.5}});
        return cov;
    
    }
    
    
    
    public Matrix getParameters(){return new Matrix(new double[][]{{A},{B},{C}});}
    
    public Matrix getParametersSteps(){
        return new Matrix(new double[][]{{dA},{dB},{dC}});
    }
    
    
    public void setParameters(Matrix param){
        A = param.get(0,0);
        B = param.get(1,0);
        C = param.get(2,0);
    }        
    
    public int getParametersN(){ return 3;}
    public int getDataN(){ return 5;}
    public void printParameters(){ System.out.println("A = "+A+"\nB = "+B+"\nC = "+C);}
    
}
