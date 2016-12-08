package numeric.minimisation.gui.infra;

import Jama.Matrix;

/**
 *
 * @author nickrowell
 */
public class SteepestDescent
extends Minimizer
{

    private double gain = 0.001;
    
    @Override
    public double[] getStep(CostFunction costFcn, double x, double y)
    {
        
        Matrix J = costFcn.getJacobian(x, y);
        
        Matrix s_minus_s0 = costFcn.getModelVector(x,y).minus(costFcn.getDataVector());
        
        // Steepest descent step...
        Matrix step = J.transpose().times(s_minus_s0);
        
        double dx = -gain*step.get(0, 0);
        double dy = -gain*step.get(1, 0);
                
        return new double[]{dx,dy};
    }

    @Override
    public String getNameString()
    {
        return "Steepest Descent";
    }
    
}
