package numeric.minimisation.gui.exec;

import infra.gui.IPanel;
import infra.io.Gnuplot;
import infra.os.OSChecker;

import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import javax.swing.*;

import numeric.minimisation.gui.dev.Version;
import numeric.minimisation.gui.infra.CostFunction;
import numeric.minimisation.gui.infra.Minimizer;
import numeric.minimisation.gui.infra.State;


/**
 * 
 * 
 * @author nrowell
 */
public class Main 
extends JFrame
{
    
    /**
	 * The serial version UID.
	 */
	private static final long serialVersionUID = 1L;

	/** Handle to main State object. */
    private State state;
    
    /** Starting coordinate text entry fields. */
    JTextField startX, startY;
    
    /**
     * IPanel to display plot of cost function contours and current/previous solution.
     */
    IPanel plotPanel;
    
    
    /** Default constructor. */
    public Main()
    {
        state = new State();
        initComponents();
    }
    
    /** Creates new form SyntheticSFRForm. */
    public Main(State pstate) 
    {
        state = pstate;
        initComponents();
    }

    /** Set up GUI components. */
    private void initComponents()
    {
       
        /** Right hand panel with plot of cost function contours and parameter steps. */
        plotPanel  = new IPanel();
        // Initialise graph
        writeCostFcnContours();
        plotGraph();
        
        /** Panel containing algorithm selection and minimisation controls. */
        JPanel entryPanel = new JPanel(new GridLayout(4,1));
        
        JLabel minimizerLabel = new JLabel("Algorithm:");
        
        // Configure minimiser type combo box
        final JComboBox<Minimizer.Type> minimizerComboBox = new JComboBox<>();
        minimizerComboBox.setSelectedItem(state.minimizer);
        minimizerComboBox.setModel(new DefaultComboBoxModel<Minimizer.Type>(Minimizer.Type.values()));
        minimizerComboBox.addActionListener(new java.awt.event.ActionListener() {
            @Override
            public void actionPerformed(java.awt.event.ActionEvent evt) 
            {
                state.minimizer = ((Minimizer.Type)minimizerComboBox.getSelectedItem()).getMinimizer();
            }
        });
        
        
        JLabel costFcnLabel = new JLabel("Function:");
        
        // Configure cost function type combo box
        final JComboBox<CostFunction.Type> costFcnComboBox = new JComboBox<>();
        costFcnComboBox.setSelectedItem(state.costFunction);
        costFcnComboBox.setModel(new DefaultComboBoxModel<CostFunction.Type>(CostFunction.Type.values()));
        costFcnComboBox.addActionListener(new java.awt.event.ActionListener() {
            @Override
            public void actionPerformed(java.awt.event.ActionEvent evt) 
            {
                // Update CostFunction in main State object
                state.costFunction = ((CostFunction.Type)costFcnComboBox.getSelectedItem()).getCostFunction();
                
                // Write CostFunction's suggested start point to the text boxes
                startX.setText(Double.toString(state.costFunction.getSuggestedStartX()));
                startY.setText(Double.toString(state.costFunction.getSuggestedStartY()));
                writeCostFcnContours();
                plotGraph();
            }
        });
        
        JPanel algorithmPanel = new JPanel(new GridLayout(2,2));
        
        
        algorithmPanel.add(minimizerLabel);
        algorithmPanel.add(minimizerComboBox);
        
        algorithmPanel.add(costFcnLabel);
        algorithmPanel.add(costFcnComboBox);
        
        
        JButton resetButton = new JButton();
        resetButton.setText("Reset");
        resetButton.addActionListener(new ActionListener() 
        {
            @Override
            public void actionPerformed(ActionEvent evt)
            {
                // Parse values from text fields:
                StartPointVerifier startPoint = new StartPointVerifier(startX.getText(), startY.getText());
                
                if(startPoint.valid)
                {
                    state.x = startPoint.x;
                    state.y = startPoint.y;
                    state.steps.clear();
                    plotGraph();
                }
                else
                {
                    // Present JOptionPane with error message
                    errorMessage(startPoint.message.toString());
                }
                
                
                
            }
        });
        
        
        
        JPanel startingPointPanel = new JPanel(new FlowLayout());
        
        startX  = new JTextField(Double.toString(state.x));
        startY  = new JTextField(Double.toString(state.y));
        
        startingPointPanel.add(new JLabel("Starting point (x,y) = ("));
        startingPointPanel.add(startX);
        startingPointPanel.add(new JLabel(","));
        startingPointPanel.add(startY);
        startingPointPanel.add(new JLabel(")"));
        startingPointPanel.add(resetButton);
        
        
        JButton stepButton = new JButton();
        stepButton.setText("One step");
        stepButton.addActionListener(new ActionListener() 
        {
            @Override
            public void actionPerformed(ActionEvent evt)
            {
                // Perform one optimization step
                
                // Push current parameters onto stack
                state.steps.add(new double[]{state.x, state.y});
                
                // Perform one optimization step
                double[] step = state.minimizer.getStep(state.costFunction, state.x, state.y);
                
                // Update parameters
                state.x += step[0];
                state.y += step[1];
                
                // Replot graph with new step
                plotGraph();
                
            }
        });
        
        JButton allStepsButton = new JButton();
        allStepsButton.setText("1000 steps");
        allStepsButton.addActionListener(new ActionListener()
        {
            @Override
            public void actionPerformed(ActionEvent evt)
            {
                // Perform one optimization step
                
                
                for(int i=0; i<1000; i++)
                {
                    
                    // Push current parameters onto stack
                    state.steps.add(new double[]{state.x, state.y});

                    // Perform one optimization step
                    double[] step = state.minimizer.getStep(state.costFunction, state.x, state.y);

                    // Update parameters
                    state.x += step[0];
                    state.y += step[1];

                    
                }
                // Replot graph with new step
                plotGraph();
                
            }
        });
        
        entryPanel.add(algorithmPanel);
        entryPanel.add(startingPointPanel);
        entryPanel.add(stepButton);
        entryPanel.add(allStepsButton);
        
        
        
        // Add to main panel
        setLayout(new FlowLayout());
        add(entryPanel);
        add(plotPanel);
        
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setTitle("Nonlinear minimization algorithms "+Version.string());
       
        pack();
        // Set initial position of frame
        setLocationRelativeTo(null);
        setVisible(true);
    }
    
    
    
    /**
     * Write the cost function data table and contours file that will be used
     * to plot the graph.
     */
    private void writeCostFcnContours()
    {
        StringBuilder script = new StringBuilder();
        
        // Get values from cost function for configuring GNUplot script
        String fxy = state.costFunction.toString();
        String xmin = Double.toString(state.costFunction.getPlotRangeX()[0]);
        String xmax = Double.toString(state.costFunction.getPlotRangeX()[1]);
        String ymin = Double.toString(state.costFunction.getPlotRangeY()[0]);
        String ymax = Double.toString(state.costFunction.getPlotRangeY()[1]);
        
        // Write table of cost function values
        script.append("reset").append(OSChecker.newline);
        script.append(fxy).append(OSChecker.newline);
        script.append("set xrange [").append(xmin).append(":").append(xmax).append("]").append(OSChecker.newline);
        script.append("set yrange [").append(ymin).append(":").append(ymax).append("]").append(OSChecker.newline);
        script.append("set isosample 250, 250").append(OSChecker.newline);
        script.append("set table '").append(state.costFcnTable.getAbsolutePath()).append("'").append(OSChecker.newline);
        script.append("splot f(x,y)").append(OSChecker.newline);
        script.append("unset table").append(OSChecker.newline);
        
        // Write cost function contours
        script.append("set contour base").append(OSChecker.newline);
        script.append("set cntrparam level incremental -3, 0.5, 3").append(OSChecker.newline);
        script.append("unset surface").append(OSChecker.newline);
        script.append("set table '").append(state.costFcnConts.getAbsolutePath()).append("'").append(OSChecker.newline);
        script.append("splot f(x,y)").append(OSChecker.newline);
        script.append("unset table").append(OSChecker.newline);
        
        // Write script to file
        File scriptf = writeGnuplotScript(script.toString());
        
        // Execute Gnuplot script
        Gnuplot.executeScript(scriptf);
    }
    
    
    /** Plot a graph of current cost function and parameter steps. */
    private void plotGraph()
    {
        // Build a Gnuplot script
        StringBuilder script = new StringBuilder();
        
        script.append("reset").append(OSChecker.newline);
        script.append("set terminal pngcairo enhanced color size 410,350").append(OSChecker.newline);
        
        // CostFunction should suggest a range here, also taking into account
        // steps taken by parameters.
        String xmin = Double.toString(state.costFunction.getPlotRangeX()[0]);
        String xmax = Double.toString(state.costFunction.getPlotRangeX()[1]);
        String ymin = Double.toString(state.costFunction.getPlotRangeY()[0]);
        String ymax = Double.toString(state.costFunction.getPlotRangeY()[1]);
        
        script.append("set xrange [").append(xmin).append(":").append(xmax).append("]").append(OSChecker.newline);
        script.append("set yrange [").append(ymin).append(":").append(ymax).append("]").append(OSChecker.newline);
        
        script.append("unset key").append(OSChecker.newline);
        script.append("set palette rgbformulae 33,13,10").append(OSChecker.newline);
        script.append("p '").append(state.costFcnTable.getAbsolutePath()).append("' with image,\\").append(OSChecker.newline);
        script.append("'").append(state.costFcnConts.getAbsolutePath()).append("' w l lt -1 lw 1.5,\\").append(OSChecker.newline);
        
        // Plot parameters steps:
        script.append("'-' w l lt -1 lc rgb 'red' lw 1.5,\\").append(OSChecker.newline);
        
        // Plot current parameters:
        script.append("'-' w p pt 5 ps 1 lc rgb 'green' lw 1.5").append(OSChecker.newline);
        
        // Inline all parameter steps:
        for(double[] step : state.steps)
        {
            script.append(Double.toString(step[0])).append("\t");
            script.append(Double.toString(step[1])).append(OSChecker.newline);
        }
        script.append(Double.toString(state.x)).append("\t");
        script.append(Double.toString(state.y)).append(OSChecker.newline);
        script.append("e").append(OSChecker.newline);
        
        // Inline current parameter
        script.append(Double.toString(state.x)).append("\t");
        script.append(Double.toString(state.y)).append(OSChecker.newline);
        script.append("e").append(OSChecker.newline);
        
        // Write script to file
        File scriptf = writeGnuplotScript(script.toString());
        
        // Execute script file
        BufferedImage img = Gnuplot.executeScript(scriptf);
        
        // Get image from Gnuplot
        plotPanel.setImage(img);
        // Reset size of panel with new image
        plotPanel.setPreferredSize();
        // Repaint etc.
        plotPanel.repaint();
        
        repaint();
        
    }
    
    
    private File writeGnuplotScript(String script)
    {
        // Write script to temp file
        File scriptf = new File(state.outputDirectory,"script.p");
        scriptf.deleteOnExit();
        
        try
        {
            BufferedWriter out = new BufferedWriter(new FileWriter(scriptf));
            out.write(script);
            out.close();
        }
        catch(IOException ioe)
        { 
            System.err.println("IOException Gnuplot.executeScript(): "+ioe.getMessage());
        }
        return scriptf;
        
    }
    
    
    /**
     * Shows a message in a pop-up JOptionPane with an error theme.
     * 
     * @param message Message to display
     */
    private void errorMessage(String message)
    {
        // Present JOptionPane with error message
        JOptionPane.showMessageDialog(this, message, "Input Error", JOptionPane.ERROR_MESSAGE);
    }
    
    /** Test function to check appearance of form. */
    public static void main(String[] args){
        
        // Create and display the form
        java.awt.EventQueue.invokeLater(
                new Runnable() 
                    {
                        @Override
                        public void run() 
                        {
                            new Main(new State());
                        }
                    });
    
    }
    
}

class StartPointVerifier 
{
    
    /** Validity flag for entries. */
    public boolean valid = true;
    /** Contains any error message produced by checking form entries. */
    public String message = "";
    
    // Parsed and checked starting point coordinates
    public double x, y;
    
    public StartPointVerifier(String sx, String sy)
    {
            
        StringBuilder msg_builder = new StringBuilder();
        
        // Attempt to parse x coordinate
        try
        {
            x = Double.parseDouble(sx);
        }
        catch(NumberFormatException nfe)
        {
            msg_builder.append("Could not parse x as a double: ").append(sx).append("\n");
            valid = false;
        }
        try
        {
            y = Double.parseDouble(sy);
        }
        catch(NumberFormatException nfe)
        {
            msg_builder.append("Could not parse y as a double: ").append(sy).append("\n");
            valid = false;
        }
        
        message = msg_builder.toString();
    }

}