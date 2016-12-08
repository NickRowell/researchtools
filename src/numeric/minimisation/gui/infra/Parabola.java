package numeric.minimisation.gui.infra;

import Jama.Matrix;

/**
 *
 * @author nickrowell
 */
public class Parabola 
extends CostFunction
{
    // Location of solution (minimum cost function)
    public double x0 = 0;
    public double y0 = 0;
    
    @Override
    public Matrix getDataVector()
    {
        return getModelVector(x0,y0);
    }
    
    @Override
    public Matrix getModelVector(double x, double y)
    {
        return new Matrix(new double[][]{{x}, {y}});
    }
    
    @Override
    public Matrix getJacobian(double x, double y)
    {
        return new Matrix(new double[][]{{1,0},
                                         {0,1}});
    }
    
    @Override
    public Matrix getJacobianAtSolution()
    {
        return getJacobian(x0,y0);
    }
    
    
    
    
    // Format cost function equation to a string for insertion into gnuplot script
    @Override
    public String toString()
    {
        return "f(x,y) = x**2 + y**2";
    }
    
    @Override
    public double[] getPlotRangeX(){ return new double[]{-1,1};}
    @Override
    public double[] getPlotRangeY(){ return new double[]{-1,1};}
    
    @Override
    public double getSuggestedStartX() { return 1.0;}
    @Override
    public double getSuggestedStartY() { return 1.0;}
    
    @Override
    public String getNameString()
    {
        return "Parabola";
    }
    
    
}
