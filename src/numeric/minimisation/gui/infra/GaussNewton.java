package numeric.minimisation.gui.infra;

import Jama.Matrix;

/**
 *
 * @author nickrowell
 */
public class GaussNewton
extends Minimizer
{

    private double gain = 0.01;
    
    @Override
    public double[] getStep(CostFunction costFcn, double x, double y)
    {
        Matrix J = costFcn.getJacobian(x, y);
        
        Matrix s_minus_s0 = costFcn.getModelVector(x,y).minus(costFcn.getDataVector());
        
        Matrix step = J.inverse().times(s_minus_s0);
        
        double dx = -gain*step.get(0, 0);
        double dy = -gain*step.get(1, 0);
                
        return new double[]{dx,dy};
    }

    @Override
    public String getNameString()
    {
        return "Gauss-Newton";
    }
    
}
