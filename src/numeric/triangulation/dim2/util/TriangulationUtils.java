package numeric.triangulation.dim2.util;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import numeric.geom.dim2.Polygon2d;
import numeric.geom.dim2.Triangle2d;

/**
 * Utilities associated with 2D triangulations.
 *
 * @author nrowell
 * @version $Id$
 */
public class TriangulationUtils {
	
	/**
	 * Assigns colours to a tessellation of polygons such that no two
	 * neighbouring polygons share the same colour.
	 * 
	 * @param neighboursMap
	 * 	Mapping of the {Polygon2d}s to be coloured to the set of neighbouring {Polygon2d}s.
	 * @param colours
	 * 	The set of colours to assign (requires at least 4), encoded as 24bit integers.
	 * @return
	 * 	A {@link Map} of {@link Polygon2d} to the corresponding colour.
	 */
	public static Map<Polygon2d, Integer> colourTessellation(Map<Polygon2d, Set<Polygon2d>> neighboursMap, int[] colours) {
		
		// Sanity check
		if(colours.length < 4) {
			throw new RuntimeException("Too few colours ("+colours.length+"), require at least 3!");
		}
		
		// Now assign colours randomly to each polygon
		Map<Polygon2d, Integer> colourMap = new HashMap<>();
		
		// Random instance used to select random colours
		Random rnd = new Random();
		   
		for(Polygon2d poly : neighboursMap.keySet()) {
			
			// Set of colours not used by the neighbours of this polygon
			List<Integer> availableColours = new LinkedList<>();
			
			// Add all the colours initially
			for(int colour : colours) {
				availableColours.add(colour);
			}
			
			// Loop over the polygon neighbours and eliminate colours already used
			for(Polygon2d neighbour : neighboursMap.get(poly)) {
				// Found coloured neighbour
				if(colourMap.containsKey(neighbour)) {
					Integer neighbourColour = colourMap.get(neighbour);
					availableColours.remove(neighbourColour);
				}
			}
			
			// Randomly select an available colour for this polygon
			int i = rnd.nextInt(availableColours.size());
			int colour = availableColours.get(i);
			
			colourMap.put(poly, colour);
		}
		
		return colourMap;
	}
	
	/**
	 * Assigns colours to triangles in a triangulation such that no two
	 * neighbouring triangles share the same colour.
	 * 
	 * @param neighboursMap
	 * 	Mapping of the {Triangle2d}s to be coloured to the set of neighbouring {Triangle2d}s.
	 * @param colours
	 * 	The set of colours to assign (requires at least 4), encoded as 24bit integers.
	 * @return
	 * 	A {@link Map} of {@link Triangle2d} to the corresponding colour.
	 */
	public static Map<Triangle2d, Integer> colourTriangulation(Map<Triangle2d, Set<Triangle2d>> neighboursMap, int[] colours) {
		
		// Sanity check
		if(colours.length < 4) {
			throw new RuntimeException("Too few colours ("+colours.length+"), require at least 3!");
		}
		
		// Now assign colours randomly to each triangle
		Map<Triangle2d, Integer> colourMap = new HashMap<>();
		
		// Random instance used to select random colours
		Random rnd = new Random();
		   
		for(Triangle2d tri : neighboursMap.keySet()) {
			
			// Set of colours not used by the neighbours of this triangle
			List<Integer> availableColours = new LinkedList<>();
			
			// Add all the colours initially
			for(int colour : colours) {
				availableColours.add(colour);
			}
			
			// Loop over the Triangle neighbours and eliminate colours already used
			for(Triangle2d neighbour : neighboursMap.get(tri)) {
				// Found coloured neighbour
				if(colourMap.containsKey(neighbour)) {
					Integer neighbourColour = colourMap.get(neighbour);
					availableColours.remove(neighbourColour);
				}
			}
			
			// Randomly select an available colour for this triangle
			int i = rnd.nextInt(availableColours.size());
			int colour = availableColours.get(i);
			
			colourMap.put(tri, colour);
		}
		
		return colourMap;
	}
	
}