package numeric.minimisation.gui.infra;

import Jama.Matrix;

/**
 *
 * @author nickrowell
 */
public class PseudoInverseMeanJacobian
extends Minimizer
{

    private double gain = 1.0;
    
    @Override
    public double[] getStep(CostFunction costFcn, double x, double y)
    {
        // Jacobian at current position
        Matrix J2 = costFcn.getJacobian(x, y);
        
        // Jacobian of cost function at solution
        Matrix J1 = costFcn.getJacobianAtSolution();
        
        Matrix s_minus_s0 = costFcn.getModelVector(x,y).minus(costFcn.getDataVector());
        
        // Mean of Jacobian pseudoinverses...
        Matrix step = (J1.plus(J2).times(0.5)).inverse().times(s_minus_s0);
        
        double dx = -gain*step.get(0, 0);
        double dy = -gain*step.get(1, 0);
                
        return new double[]{dx,dy};
    }

    @Override
    public String getNameString()
    {
        return "Pseudoinverse of mean Jacobian";
    }
    
}
