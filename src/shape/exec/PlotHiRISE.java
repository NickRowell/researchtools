package shape.exec;

import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Logger;

import javax.imageio.ImageIO;

import shape.Equirectangular;

/**
 * This class provides an application that processes HiRISE RDR images and writes transformed
 * areocentric image data to file.
 * 
 * @author nrowell
 */
public class PlotHiRISE {

	/**
	 * The {@link Logger}.
	 */
	private static final Logger logger = Logger.getLogger(PlotHiRISE.class.getName());
	
	/**
	 * Main application entry point.
	 * 
	 * @param args
	 * 	The command line arguments (ignored).
	 * @throws IOException
	 * 	If there's a problem loading the image and label data
	 */
	public static void main(String[] args) throws IOException {
		
		// HiRISE JPEG image file
		File imageFile = new File("/home/nrowell/Temp/HiRISE/AEB_000001_0000_RED.NOMAP.browse.jpg");
		
		// Associated HiRISE label file containing calibration coefficients etc
		File labelFile = new File("/home/nrowell/Temp/HiRISE/AEB_000001_0000_RED.LBL.txt");
		
		// File to store the transformed image data ready for plotting
		File output = new File("/home/nrowell/Temp/HiRISE/AEB_000001_0000_RED.txt");
		
		// Load the transformation from the label file
		Equirectangular eq = getHiriseProjection(labelFile);
		
		logger.info("Loaded HiRISE Equirectangular projection coefficients: " + eq.toString());
		
		// Load the image
		BufferedImage im = ImageIO.read(imageFile);
		
		logger.info("Loaded HiRISE image " + imageFile.getName());
		
		// Write the projected image data to file
		BufferedWriter out = new BufferedWriter(new FileWriter(output));
		
		// Write header
		out.write("# HiRISE image " + imageFile.getName() + "\n");
		out.write("# LINE_PROJECTION_OFFSET   = " + eq.LINE_PROJECTION_OFFSET + " [pixels]\n");
		out.write("# SAMPLE_PROJECTION_OFFSET = " + eq.SAMPLE_PROJECTION_OFFSET + " [pixels]\n");
		out.write("# MAP_SCALE                = " + eq.MAP_SCALE + " [metres/pixel]\n");
		out.write("# CENTER_LONGITUDE         = " + eq.CENTER_LONGITUDE + " [radians]\n");
		out.write("# CENTER_LATITUDE          = " + eq.CENTER_LATITUDE + " [radians]\n");
		out.write("# A_AXIS_RADIUS            = " + eq.A_AXIS_RADIUS + " [km]\n");
		out.write("# C_AXIS_RADIUS            = " + eq.C_AXIS_RADIUS + " [km]\n");
		out.write("# Columns:\n");
		out.write("# 1) Areocentric longitude [deg]\n");
		out.write("# 2) Areocentric latitude [deg]\n");
		out.write("# 3) Pixel value [-]\n");
		
		// Loop over the image pixels; the sample coordinate spans the image width and the line
		// coordinate spans the image height.
		for(int line = 0; line < im.getHeight(); line++) {
			for(int sample = 0; sample < im.getWidth(); sample++) {
				
				// Pixel value in ARGB format
				int pix = im.getRGB(sample, line);
				
				// Get red channel in 0:255
				int r = (pix & 0x00FF0000) >> 16; 
			
				// Corresponding longitude and latitude
				double[] lonLat = eq.getLongLat(sample, line);
				
				out.write(Math.toDegrees(lonLat[0]) + "\t" + Math.toDegrees(lonLat[1]) + "\t" + r + "\n");
			}
		}
		
		out.close();
	}
	
	/**
	 * Reads a HiRISE label file associated with an image, and creates an {@link Equirectangular}
	 * instance configured with the appropriate calibration coefficients.
	 * 
	 * @param labelFile
	 * 	The {@link File} containing HiRISE calibration coefficients.
	 * @return
	 * 	An {@link Equirectangular} instance configured with the calibration coefficients contained
	 * in the given label file.
	 * @throws IOException 
	 */
	public static Equirectangular getHiriseProjection(File labelFile) throws IOException {
		
		BufferedReader in = new BufferedReader(new FileReader(labelFile));
		
		String line = null;
		
		// Map containing the calibration coefficients to be read from the file
		Map<String, Double> coefficients = new HashMap<>();
		coefficients.put("LINE_PROJECTION_OFFSET", Double.NaN);
		coefficients.put("SAMPLE_PROJECTION_OFFSET", Double.NaN);
		coefficients.put("MAP_SCALE", Double.NaN);
		coefficients.put("CENTER_LONGITUDE", Double.NaN);
		coefficients.put("CENTER_LATITUDE", Double.NaN);
		coefficients.put("A_AXIS_RADIUS", Double.NaN);
		coefficients.put("C_AXIS_RADIUS", Double.NaN);
		
		while((line = in.readLine()) != null) {
			
			// Trim leading & trailing whitespace
			line = line.trim();
			
			// Ignore comments
			if(line.startsWith("/*")) {
				continue;
			}
			
			// Check if any of our coefficients are declared on this line
			for(String coeff : coefficients.keySet()) {
				if(line.contains(coeff)) {
					// Parse the parameter from the line
					
					// Split the string on whitespace
					String[] parts = line.split("\\s+");
					
					// Assume line follows format COEFFICIENT_NAME = COEFFICIENT_VALUE UNITS
					double value = Double.parseDouble(parts[2]);
					
					// Store the parsed value
					coefficients.put(coeff, value);
				}
			}
		}
		
		in.close();
		
		// Check we found all the coefficients
		if(coefficients.containsValue(Double.NaN)) {
			
			// Missing keys
			Set<String> missingKeys = new HashSet<>();
			
			for(Entry<String,Double> entry : coefficients.entrySet()) {
				if(entry.getValue() == Double.NaN) {
					missingKeys.add(entry.getKey());
				}
			}
			
			throw new RuntimeException("Couldn't locate calibration coefficients: " + missingKeys);
		}
		
		// Create the Equirectangular instance
		Equirectangular eq = new Equirectangular(coefficients.get("LINE_PROJECTION_OFFSET"),
				coefficients.get("SAMPLE_PROJECTION_OFFSET"),  coefficients.get("MAP_SCALE"),
				coefficients.get("CENTER_LONGITUDE"), coefficients.get("CENTER_LATITUDE"),
				coefficients.get("A_AXIS_RADIUS"), coefficients.get("C_AXIS_RADIUS"));
		
		return eq;
	}
	
}