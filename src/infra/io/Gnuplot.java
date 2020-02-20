/**
 * 
 * Name:
 *  Gnuplot.java
 * 
 * Purpose:
 *  Utility class used to handle Gnuplot processes.
 * 
 * Language:
 * Java
 *
 * Author:
 * Nicholas Rowell
 * 
 */
package infra.io;

import java.awt.BorderLayout;
import java.awt.image.BufferedImage;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.swing.JFrame;

import infra.gui.IPanel;
import infra.os.OSChecker;
import util.GuiUtil;

/**
 * Class provides a simple, clean interface to GNUplot and methods to execute scripts, return plots
 * as BufferedImages and other utilities.
 * 
 * @author nrowell
 * @version $Id$
 */
public class Gnuplot 
{
	/**
	 * The Logger.
	 */
	private static final Logger log = Logger.getLogger( Gnuplot.class.getName() );
    
	/**
	 * Spawns a Gnuplot process, executes the given File as a Gnuplot script and attempts
	 * to read the output as a BufferedImage.
	 * @param script
	 * 		The Gnuplot script file to execute
	 * @return
	 * 		The BufferedImage produced as a result of using Gnuplot to execute the script
	 */
    public static BufferedImage executeScript(File script)
    {
        try
        {
            // Run GNUplot from Java, with script file just written.
            final Process proc = Runtime.getRuntime().exec(getGnuplotCommand(script));

            // Capture error stream from Gnuplot process
            StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "GNUPLOT");
            
            // Capture image from Gnuplot output
            ImageGobbler imageGobbler = new ImageGobbler(proc.getInputStream());
            
            // Start errorgobbler
            errorGobbler.start();
            imageGobbler.start();
            
            // Wait for GNUplot process to complete.
            proc.waitFor();
            
            // wait for error (messages) thread output to finish
            errorGobbler.join();
            
            // wait for output (image related) thread to finish
            imageGobbler.join();
            
            // Get image from imageGobbler
            return imageGobbler.img;
            
        }
        catch(IOException | InterruptedException e) {
        	log.log(Level.SEVERE, "Experienced "+e+" when executing Gnuplot script!", e);
        	return null;
        }
    }
    
    /**
	 * Spawns a Gnuplot process, executes the given File as a Gnuplot script and does not
	 * attempt to read any output (e.g. an image) from the process. This is useful for
	 * scripts that don't produce a plot to the output stream, i.e. where the output has
	 * been redirected to file or cases where there is no output.
	 * 
	 * @param script
	 * 		The Gnuplot script file to execute
	 * @return
	 * 		The exit code returned by GNUplot process
	 */
    public static int executeScriptNoOutput(File script)
    {
        try
        {
            // Run GNUplot from Java, with script file just written.
            final Process proc = Runtime.getRuntime().exec(getGnuplotCommand(script));

            // Capture error stream from Gnuplot process
            StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "GNUPLOT");
            
            // Start errorgobbler
            errorGobbler.start();
            
            // Wait for GNUplot process to complete.
            int exitCode = proc.waitFor();
            
            // wait for error (messages) thread output to finish
            errorGobbler.join();
            
            return exitCode;
        }
        catch(IOException | InterruptedException e) {
        	log.log(Level.SEVERE, "Experienced "+e+" when executing Gnuplot script!", e);
        	// Return exit code indicating failure
        	return 1;
        }
    }
    
    
    /**
     * Writes the script to a temporary file, executes it, deletes the file and returns the
     * output as a BufferedImage.
     * 
     * @param script
     * 		String containing the Gnuplot script to execute
     * @return
     * 		BufferedImage containing the output of the script
     * @throws IOException 
     */
    public static BufferedImage executeScript(String script) throws IOException
    {
    	File scriptFile = File.createTempFile("gnuplot", null);
    	BufferedWriter scriptOut = new BufferedWriter(new FileWriter(scriptFile));
    	scriptOut.write(script);
    	scriptOut.close();
    	BufferedImage plot = Gnuplot.executeScript(scriptFile);
    	scriptFile.delete();
    	return plot;
    }
    
    /**
     * Creates a JFrame and displays the given BufferedImage.
     * 
     * @param image
     * 	BufferedImage to display.
     */
    public static void displayImage(BufferedImage image) {
    	
    	final IPanel ipanel = new IPanel(image);
		
		// Create and display the form
        java.awt.EventQueue.invokeLater(
                new Runnable() 
                    {
                        @Override
                        public void run() 
                        {
                            final JFrame tester = new JFrame();
                            tester.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
                            tester.setLayout(new BorderLayout());
                            tester.add(ipanel, BorderLayout.CENTER);
                            
                            // Add right-click menu allowing displayed image to be saved to file
                            GuiUtil.addRightClickMenuSaveJPanelImage(tester, ipanel);
                            
                            tester.pack();
                            tester.setVisible(true);
                        }
                    });
    }
    
    /**
     * Get a simple Gnuplot script that produces no output. This is useful for testing if
     * GNUplot works and can be found on the system path.
     * 
     * @return
     * 	A simple GNUplot script that produces no output.
     */
    private static String getTestScript()
    {
        return  "# Gnuplot script used simply to test if Gnuplot is working "+
                OSChecker.newline + "exit";
    }
    
    /**
     * Tests if Gnuplot is available by attempting to execute a simple script.
     * 
     * @param tempDir
     * 	Working directory to which a script file will be temporarily written.
     * @return
     * 	Boolean indicating if GNUplot executed successfully.
     */
    public static boolean isGnuplotWorking(File tempDir)
    {
    	
    	// Write test script to temporary directory:
    	File script = null;
		try {
			script = File.createTempFile("gnuplot", null, tempDir);
			BufferedWriter out = new BufferedWriter(new FileWriter(script));
	        out.write(getTestScript());
	        out.close();
		} catch (IOException e) {
			log.log(Level.SEVERE, "Experienced IOException when writing temporary script file!", e);
            return false;
		}
        
    	// Execute test script and get exit code
    	int exitCode = executeScriptNoOutput(script);
    	
    	// Cleanup
    	script.delete();
    	
    	if(exitCode==0)
        {
            // Gnuplot found & works
        	log.log(Level.INFO, "Gnuplot works: exit code "+exitCode);
        	return true;
        }
        else
        {
            log.log(Level.SEVERE, "Gnuplot failed! Exit code "+exitCode);
            return false;
        }
    	
    }
    
    /**
     * Assembles a GNUplot command as an array of Strings ready to be executed using
     * {@link Runtime#ggetRuntime()#exec(String[])}.
     * 
     * @param script
     * 	The File containing the script to execute.
     * @return
     * 	Array of Strings containing the commands to launch GNUplot on the script using
     * the Runtime.
     */
    private static String[] getGnuplotCommand(File script) {
    	
    	// On Windows platforms, must surround path to script with double
        // quotes in case it contains spaces, i.e. "Program Files\gnuplot".
        // Linux doesn't seem to have a problem handling these cases when 
        // the command is invoked from Java using 
        // Runtime.getRuntime().exec(...), but bizarrely it doesn't like 
        // the quotes, so for linux platforms we dont use them.
        String scriptArg = script.getAbsolutePath();

        switch (OSChecker.getOS())
        {
            case WINDOWS:
                scriptArg = "\"" + scriptArg + "\"";
			default:
				break;
        }

        String[] command = new String[]{"gnuplot", scriptArg};
        
        return command;
    }
}
