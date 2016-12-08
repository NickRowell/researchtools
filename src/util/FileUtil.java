package util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Utilities related to file handling.
 *
 *
 * @author nrowell
 * @version $Id$
 */
public class FileUtil {
	
	/**
	 * The logger.
	 */
	protected static final Logger logger = Logger.getLogger(FileUtil.class.getName());
	
	/**
	 * Reads the text file and returns the final N lines in sequence.
	 * @param file
	 * 	The {@link File} containing text to read
	 * @param N
	 * 	The number of lines to read from the end of the file
	 * @return
	 * 	The final N lines of text
	 */
	public static String[] getTail(File file, int N) {
		
		// Implement a simple ring buffer
		String[] buffer = new String[N];
		int writePos=0;
		
		try (BufferedReader in = new BufferedReader(new FileReader(file))) {
            while(in.ready()) {
            	
            	buffer[writePos] = in.readLine();
            	// Increment the write position
            	writePos++;
            	// Wrap around
            	writePos = writePos%N;
            }
        }
        catch(IOException e) {
        	logger.log(Level.SEVERE, "Exception reading file "+file, e);
			return new String[0];
        }
		
		// Read the contents out from the buffer in the order they were written to it
		String[] out = new String[N];
		for(int i=0; i<N; i++) {
			out[i] = buffer[(writePos+i)%N];
		}
		return out;
	}
	
}