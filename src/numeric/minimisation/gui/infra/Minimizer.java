package numeric.minimisation.gui.infra;


/**
 * Parent class representing all minimization algorithms.
 * 
 * @author nickrowell
 */
public abstract class Minimizer
{
    /** Each Minimizer has a corresponding Type. */
    public static enum Type
    {
        GaussNewton(new GaussNewton()),
        SteepestDescent(new SteepestDescent()),
        MeanJacobianPseudoInverses(new MeanJacobianPseudoInverses()),
        PseudoInverseMeanJacobian(new PseudoInverseMeanJacobian());
    
        private Minimizer minimizer;
        
        Type(Minimizer pminimizer){ minimizer = pminimizer;}
        
        public Minimizer getMinimizer(){ return minimizer;}
        
        @Override
        public String toString(){ return minimizer.getNameString();}
    
    }
    
    /**
     * Main minimization function - get a single step for the parameter.
     * @param costFcn   Cost function to be minimized
     * @param x         Current X coordinate
     * @param y         Current Y coordinate
     * @return 			The step to apply to the current coordinates
     */
    public abstract double[] getStep(CostFunction costFcn, double x, double y);
    
    /**
     * Get an appropriate name for this minimizer.
     * @return
     * 		String containing a human-readable name for the minimizer
     */
    public abstract String getNameString();
    
}
