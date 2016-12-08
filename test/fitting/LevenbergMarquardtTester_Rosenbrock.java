package fitting;

import minimisation.nllsq.LevenbergMarquardt;
import Jama.*;


/**
 * This class tests the performance of the LevenbergMarquardt class on
 * the Rosenbrock function.
 *
 * There are two parameters (x & y) and a single datapoint (f(x,y) - Rosenbrock
 * function to minimize).
 *
 * @author nickrowell
 */
public class LevenbergMarquardtTester_Rosenbrock extends LevenbergMarquardt{

    // Parameters
    private double x;
    private double y;

    // Finite step sizes in each parameter
    private double dx=1E-4, dy=1E-4;
    
    // Data
    private double z = 0.0;
    
    public LevenbergMarquardtTester_Rosenbrock(double X, double Y){
        x = X;
        y = Y;
    }

    @Override
    public String toString(){
        return x+","+y;
    }



    public static void main(String[] args){

        // Initialise parameters to some values
        LevenbergMarquardtTester_Rosenbrock test = new LevenbergMarquardtTester_Rosenbrock(-2,-0.5);

        test.fit(100000,true);

        test.printParameters();
        
        System.out.println("Function should converge to (1,1)");

        System.out.println("Parameter covariance by fourth order propagation:");
        test.getFourthOrderCovariance().print(10, 10);
        
    }

    /** 
     * Single data point is zero - the minimum of the Rosenbrock function. The
     * LM algorithm will adjust the parameters (x & y) until it finds the
     * point in the Rosenbrock function where the difference is smallest.
     */
    public Matrix getData(){
    
        return new Matrix(new double[][]{{z}});
    
    }
    
    /**
     * Update the data point to test parameter sensitivity.
     */
    public boolean updateData(Matrix delta){
                
        z += delta.get(0,0);
        return true;
    }
    

    /** Model is Rosenbrock function */
    public Matrix getModel(){

        Matrix R = new Matrix(1,1);
        
        // Define Rosenbrock function
        R.set(0,0, (1-x)*(1-x) + 100*(y-x*x)*(y-x*x));

        return R;

    }    
    
    /** 
     * Inverse covariance weighting. Rosenbrock function is very wide and flat
     * so use small variance on depth of minimum to provide realistic imitation
     * of behaviour on real data.
     */
    public Matrix getCovariance(){
        
        return new Matrix(new double[][]{{0.00001}});
    
    }
    


    /** Implement updateParameter() function */
//    public boolean updateParameters(Matrix delta){
//
//        x += delta.get(0,0);
//        y += delta.get(1,0);
//
//        return true;
//
//    }

    /** Implement getParameters() method */
    public Matrix getParameters(){
        return new Matrix(new double[][]{{x},{y}});
    }
    
    /** Implement getParametersSteps() method */
    public Matrix getParametersSteps(){
        return new Matrix(new double[][]{{dx},{dy}});
    }    
    
    
    public void setParameters(Matrix param){
        x = param.get(0,0);
        y = param.get(1,0);
    }
     
    /** Get number of parameters */
    public int getParametersN(){ return 2;}

    /** Get number of data points */
    public int getDataN(){ return 1;}

    /** Implement printParameters() method */
    public void printParameters(){
        System.out.println("Fitted parameters = "+this.toString());
    }


}

