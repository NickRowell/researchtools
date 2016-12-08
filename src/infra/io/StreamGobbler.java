/**
 * 
 * Name:
 *  StreamGobbler.java
 * 
 * Purpose:
 *  Class used to handle IO stream from Gnuplot.
 * 
 * Language:
 * Java
 *
 * Author:
 * Nicholas Rowell
 * 
 */
package infra.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.logging.Level;
import java.util.logging.Logger;

public class StreamGobbler 
extends Thread
{
	/**
	 * The Logger.
	 */
	private static final Logger log = Logger.getLogger( StreamGobbler.class.getName() );
	
	/**
	 * The InputStream to gobble.
	 */
    InputStream is;
    
    /**
     * A string appended to any output in order to identify the source.
     */
    String type;

    public StreamGobbler(InputStream is, String type) 
    {
        this.is = is;
        this.type = type;
    }

    @Override
    public void run() 
    {
        try 
        {
            InputStreamReader isr = new InputStreamReader(is);

            BufferedReader br = new BufferedReader(isr);
            String line = null;
            while ((line = br.readLine()) != null) 
            {
                log.log(Level.INFO, type + "> " + line);
            }

        } 
        catch (IOException ioe) 
        {
            ioe.printStackTrace();
        }
    }
}
