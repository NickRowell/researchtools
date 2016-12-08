package numeric.minimisation.gui.infra;

import Jama.Matrix;

/**
 * Parent class for CostFunction types. Defines main interface.
 * @author nickrowell
 */
public abstract class CostFunction
{
    public static enum Type
    {
        Rosenbrock(new Rosenbrock()),
        Parabola(new Parabola());
    
        private CostFunction costfunction;
        
        Type(CostFunction pcostfcn){ costfunction = pcostfcn;}
        
        public CostFunction getCostFunction(){ return costfunction;}
        
        @Override
        public String toString(){ return costfunction.getNameString();}
    
    }
    
    // Cost function is the sum-of-squared residuals
    public double getCostFunction(double x, double y)
    {
        
        Matrix s_minus_s0 = getModelVector(x,y).minus(getDataVector());

        // f(x,y)
        Matrix f = s_minus_s0.transpose().times(s_minus_s0);
        
        return 0.5 * f.get(0, 0);
    }
    
    
    
    
    
    // Observation vector is 'fixed': we know the observation s(x) but not
    // the parameters x that produce it. That's what we're trying to deduce
    // by fitting the model to the data.
    public abstract Matrix getDataVector();
    
    // Model vector is calculated using the current parameters values.
    public abstract Matrix getModelVector(double x, double y);
    
    // Jacobian of model vector wrt the parameters
    public abstract Matrix getJacobian(double x, double y);
    
    
    
    public abstract Matrix getJacobianAtSolution();
    
    
    // Stuff related to GNUplot graphing:
    
    /** Suggest a suitable starting point. */
    public abstract double getSuggestedStartX();
    public abstract double getSuggestedStartY();
    
    /** Get GNUplot script that writes a table of the CostFunction values. */
    @Override
    public abstract String toString();
    
    // Specify some values to tweak GNUplot graphing:
    public abstract double[] getPlotRangeX();
    public abstract double[] getPlotRangeY();
    public abstract String getNameString();
}
