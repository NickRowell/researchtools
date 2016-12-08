package util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.LinkedList;
import java.util.List;

/**
 * Class provides utilities for parsing data from Strings and Files.
 *
 * @author nrowell
 * @version $Id$
 */
public class ParseUtil {
	
	/**
	 * Delimiter representing any amount of whitespace.
	 */
	public static final String whitespaceDelim = "\\s+";

	/**
	 * A {@link NumberFormat} used to format large numbers for display as strings.
	 */
	public static final NumberFormat expNumFormat = new DecimalFormat("0.######E0");
	
	/**
	 * Attempt to parse a double value from the String.
	 * @param variable
	 * 	The String to parse
	 * @param name
	 * 	The name of the variable (for inclusion in error message)
	 * @return
	 * 	The parsed value of the variable
	 */
	public static double parseNumber(String variable, String name) {
    	try {
            return Double.parseDouble(variable);
        }
        catch(NumberFormatException nfe) {
        	throw new NumberFormatException("Could not parse "+name+" from: "+variable);
        }
	}
    
    public static double parseAndCheckGreaterThan(String variable, String name, double x) {
    	double v = parseNumber(variable, name);
        if(v<=x) {
        	throw new IllegalArgumentException(name + " must be greater than "+x+"!");
        }
    	return v;
    }
    
    public static double parseAndCheckGreaterThanOrEqualTo(String variable, String name, double x) {
    	double v = parseNumber(variable, name);
        if(v<x) {
        	throw new IllegalArgumentException(name + " must be greater than or equal to "+x+"!");
        }
    	return v;
    }
    
    public static double parseAndCheckLessThan(String variable, String name, double x) {
    	double v = parseNumber(variable, name);
        if(v>=x) {
        	throw new IllegalArgumentException(name + " must be less than "+x+"!");
        }
    	return v;
    }
    
    public static double parseAndCheckWithinRangeInclusive(String variable, String name, double min, double max) {
    	double v = parseNumber(variable, name);
        if(v < min || v > max) {
        	throw new IllegalArgumentException(name + " must be within range ["+min+":"+max+"]!");
        }
    	return v;
    }
	
    /**
     * Reads the contents of the text file into a List of Strings.
     * @param resourceLocation
     * 	The location of the file as a resource.
     * @return
     * 	List<String> containing each line from the file.
     */
    public static List<String> parseResource(String resourceLocation) {
    	
    	// Open reader on file
        InputStream is = (new ParseUtil()).getClass().getClassLoader().getResourceAsStream(resourceLocation);
        
        // Read entire WDLF data file into a List of Strings
        List<String> lines = new LinkedList<>();
        String line;
        
        try(BufferedReader in = new BufferedReader(new InputStreamReader(is))) {
			while((line=in.readLine()) != null)
			{
				lines.add(line);
			}
		} catch (IOException e) {
			System.out.println("Error reading "+resourceLocation);
			e.printStackTrace();
		}
    	
        return lines;
    }
    
    /**
	 * Parse a file containing tabulated values.
	 * 
	 * @param in
	 * 	A {@link BufferedReader} open on the file to parse.
	 * @param comments
	 * 	A {@link List<String>} containing comment characters/strings. Lines that start with, or whose first column contains,
	 * an entry in this list is ignored as a comment line.
	 * @return
	 * 	A 2D array containing the tabulated data parsed from the file. The leading index is the column number, the
	 * trailing index is the row number; so the array data[0] contains the column 0 entries for each row in sequence.
	 * @throws IOException
	 * 	If there's a problem reading the file.
	 */
	public static double[][] parseFile(BufferedReader in, String delimiter, List<String> comments) throws IOException {
		
		String record;
		
		List<double[]> dataLists = new LinkedList<>();

        // Loop over table and read all rows into lists
		int lineNum = 0;
        while ((record = in.readLine()) != null) {
        	
        	lineNum++;
        	
        	// Trim off leading & trailing whitespace
        	// TODO: remove leading/trailing delimiters in general
        	record = record.trim();
        	
        	// Avoid completely blank lines
        	if (record.length()==0) {
                continue;
        	}
        	
            // Avoid any commented out lines
        	String firstChar = record.substring(0, 1);
            if (comments.contains(firstChar)) {
                continue;
            }
            
            // Split the String into columns
            String[] parts = record.split(delimiter);
            
            // Check if first column contains a comment flag
            if(comments.contains(parts[0])) {
            	continue;
            }
            
            // Parse the values and store in array
            double[] values = new double[parts.length];
            
            for(int i=0; i<parts.length; i++) {
            	values[i] = parseNumber(parts[i], "Line "+lineNum+", Column "+i);
            }
            dataLists.add(values);
        }
        
        // Convert lists to 2D array of doubles
        double[][] data = new double[dataLists.get(0).length][dataLists.size()];
        
        // Sanity check that all lines have the same number of entries
        int n=dataLists.get(0).length;
        
        for(int i=0; i<dataLists.size(); i++) {
        	
        	double[] dataList = dataLists.get(i);
        	
        	if(dataList.length != n) {
        		throw new RuntimeException("Different number of entries in consecutive records! "+
        									"Expected "+n+", found "+dataList.length);
        	}
        	
        	n = dataList.length;
        	
        	for(int j=0; j<n; j++) {
        		data[j][i] = dataList[j];
        	}
        }
        return data;
	}
	
}