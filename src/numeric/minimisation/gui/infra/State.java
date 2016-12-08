/**
 * 
 * Name:
 *  State
 * 
 * Purpose:
 *  Contains current execution state of program.
 * 
 * Language:
 * Java
 *
 * Author:
 * Nicholas Rowell
 * 
 */
package numeric.minimisation.gui.infra;


import java.io.File;
import java.util.LinkedList;
import java.util.List;

public class State 
{
    
    /**
     * Working directory for temporary files.
     */
    public File outputDirectory = new File("/home/nrowell/Temp");
    
    /** 
     * Current Minimization algorithm.
     */
    public Minimizer minimizer = Minimizer.Type.GaussNewton.getMinimizer();
    
    /**
     * Objective/cost function.
     */
    public CostFunction costFunction = CostFunction.Type.Rosenbrock.getCostFunction();
    
    /**
     * Current parameter values.
     */
    public double x = costFunction.getSuggestedStartX();
    public double y = costFunction.getSuggestedStartY();
    
    /**
     * List of previous steps.
     */
    public List<double[]> steps = new LinkedList<>();
    
    /**
     * File containing tabulated cost function.
     */
    public File costFcnTable = new File(outputDirectory,"costfcntab.dat");
    
    /**
     * File containing contours of cost function.
     */
    public File costFcnConts = new File(outputDirectory,"costfcnconts.dat");

}
