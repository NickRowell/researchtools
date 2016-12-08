package numeric.minimisation.gui.infra;

import Jama.Matrix;

/**
 *
 * @author nickrowell
 */
public class Rosenbrock 
extends CostFunction
{
    
    // Location of solution (minimum cost function)
    public double x0 = 1;
    public double y0 = 1;
    
    @Override
    /**
     * The data vector is what is obtained when the parameters solution is
     * processed through the model.
     */
    public Matrix getDataVector()
    {
        return getModelVector(x0,y0);
    }
    
    @Override
    public Matrix getModelVector(double x, double y)
    {
        return new Matrix(new double[][]{{Math.sqrt(2)*(1-x)},
                                         {Math.sqrt(2)*10*(y-x*x)}});
    }
    
    @Override
    /**
     * 
     * Jacobian has the form 
     * 
     *                      | ds_x/dx ds_x/dy |
     *                  J = |                 |
     *                      | ds_y/dx ds_y/dy |
     * 
     */
    public Matrix getJacobian(double x, double y)
    {
        return new Matrix(new double[][]{{-Math.sqrt(2),      0},
                                         {-Math.sqrt(2)*20*x, Math.sqrt(2)*10}});
    }
    
    @Override
    public Matrix getJacobianAtSolution()
    {
        return getJacobian(x0,y0);
    }
    
    
    
    
    @Override
    public String toString()
    {
        return "f(x,y) = (1-x)**2 + 100*((y-x**2)**2)";
    }
    
    @Override
    public double[] getPlotRangeX(){ return new double[]{-2,2};}
    @Override
    public double[] getPlotRangeY(){ return new double[]{-1,3};}
    
    @Override
    public double getSuggestedStartX() { return -0.9;}
    @Override
    public double getSuggestedStartY() { return -0.7;}
    
    @Override
    public String getNameString()
    {
        return "Rosenbrock";
    }
    
    
}
