/**
 * Copyright (c) 2011 University of Dundee
 *
 * Name:
 * ContourDetection.java
 *
 * Purpose:
 * Contour detection algorithm processes 2D array of floating point values and
 * a threshold value, and finds contour lines where the threshold is crossed.
 * 
 * Coordinates of contours vary from (0,0) at the top left corner of array to
 * (N+1,M+1) at bottom right corner of array, for N*M arrays.
 * 
 * To do: linear interpolation of crossing points rather than just
 *        midpoint of neighbouring cells.
 * 
 * Language:
 * Java
 *
 * Author:
 * Nicholas Rowell (nickrowell@computing.dundee.ac.uk)
 *
 */
package numeric.data;


import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import numeric.data.ContourDetection.DIRECTION;

public class ContourDetection {
    
    // Represent possible search orders
    static enum ROTATION{ CLOCKWISE, ANTICLOCKWISE};
    
    // Possible states
    static enum STATE{ABOVE_THRESHOLD,
                      BELOW_THRESHOLD};
    
    // There are eight directions. The order of these in the list sets initial
    // direction of search once a new edge has been started.
    static enum DIRECTION{ W(0,-1), 
                           SW(1,-1),
                           S(1,0), 
                           SE(1,1), 
                           E(0,1), 
                           NE(-1,1), 
                           N(-1,0),
                           NW(-1,-1);
    
                           private final int i_off;
                           private final int j_off;
                           
                           DIRECTION(int di, int dj){
                               i_off = di;
                               j_off = dj;
                           }
                           
                           DIRECTION getOpposite(){
                               return DIRECTION.values()[(ordinal()+4)%8];
                           }
                           
                           boolean isHorizontal(){
                               return ordinal()==0 || ordinal()==4;
                           }
                           boolean isVertical(){
                               return ordinal()==2 || ordinal()==6;
                           }
                           
                           boolean isLateral(){
                               return isVertical() || isHorizontal();
                           }
                           
                           boolean isDiagonal(){
                               return !isLateral();
                           }
                           
                           boolean isEast(){return ordinal()==4;}
                           boolean isWest(){return ordinal()==0;}
                           
                        };

    
    // Record current search state, i.e. whether we are above or below the
    // contour threshold.
    private STATE state;
    
    // Contour level
    private double threshold;
    
    // Data structure representing vertical boundaries between cells in grid.
    // When one is found to be part of a contour, it is indicated in this
    // array so that it is not re-found later on.
    boolean[][] BOUNDARY;
    
    // Copy of input array padded with a one-pixel border of below-threshold
    // values, so that contours always form closed loops and we know the
    // initial state of the search algorithm.
    double[][] padded;
    
    // List of contours
    List<Contour> contours = new ArrayList<Contour>();
    
    // In order to transform from array index coordinates to physically 
    // meaningful coordinates, we record physical coordinate of top left of
    // cell [0][0], and the step size between cells.
    double i0, j0, di, dj;
    
    // First index iterates over rows.
    // Second index iterates over columns.
    public ContourDetection(double[][] array, double _threshold){
        
        // Give default values for coordinate transformation.
        this(array, _threshold, 0, 0, 1, 1);
        
    }
    
    // First index iterates over rows.
    // Second index iterates over columns.
    public ContourDetection(double[][] array, double _threshold, double _i0,
                            double _j0, double _di, double _dj){
        
        i0 = _i0;
        j0 = _j0;
        di = _di;
        dj = _dj;
        
        threshold = _threshold;
        
        // Pad boundary with below-threshold values. This is
        // necessary to ensure that all contours form closed loops.
        padded = new double[array.length+2][array[0].length+2];
        
        for(int i=0; i<padded.length; i++){
            for(int j=0; j<padded[i].length; j++){
                
                // Are we in a border pixel?
                if(i==0 || i == padded.length-1 || 
                   j==0 || j == padded[i].length-1)
                    padded[i][j] = threshold - 1.0;
                else
                    padded[i][j] = array[i-1][j-1];
            }
        }
        
        // Data structure to represent boundaries between neighbouring cells.
        // Values initialised to false.
        BOUNDARY  = new boolean[array.length][array[0].length+1];
        
        // Initial state of contour search. Because we pad with below-threshold
        // values, we know that the initial state is SEARCHING_BELOW_THRESHOLD.
        state = STATE.BELOW_THRESHOLD;
        
        // Loop over all pixels, looking for contour boundaries
        for(int i=0; i<padded.length; i++){
            for(int j=0; j<padded[i].length; j++){     
                
                // Have we crossed a contour?
                if((padded[i][j] > threshold && state == STATE.BELOW_THRESHOLD) ||
                   (padded[i][j] < threshold && state == STATE.ABOVE_THRESHOLD)){
                    
                    // We ALWAYS cross from the -j direction. Contour lies between
                    // cells [i][j-1] and [i][j] in the padded array.
                    // This corresponds to cell [i-1][j-2] and [i-1][j-1] in the
                    // original array.
                    
                    // Has this boundary been previously logged as a contour?
                    if(BOUNDARY[i-1][j-1]){
                        // Reset state and take no further action.
                        switchState();
                    }
                    
                    else{
                        
                        // Log this point as a contour, then proceed with
                        // contour growing.
                        BOUNDARY[i-1][j-1] = true;
                        switchState();
                        
                        // Initialise new contour and log starting point & search direction
                        Contour contour = new Contour(new int[]{i,j}, DIRECTION.W);
                        
                        // Log first contour crossing at midpoint of cells
                        
                        // Reverse direction is west... general form for
                        // midpoint from current location and search direction:
                        contour.addPoint(j-0.5+((double)DIRECTION.W.j_off)/2.0,
                                         i-0.5+((double)DIRECTION.W.i_off)/2.0);                    
                        
                        
                        if(state==STATE.ABOVE_THRESHOLD)
                            // Grow contour starting search in SW direction.
                            growContour(contour, DIRECTION.SW, i, j);
                        else if(state==STATE.BELOW_THRESHOLD)
                            // Grow contour starting search in S direction.
                            growContour(contour, DIRECTION.S, i, j);
                                    
                        // Add contour to list
                        contours.add(contour);
                        
                    }
                    
                }
                
            }
            
        }
    
    }
    
    private void switchState(){
        if(state == STATE.ABOVE_THRESHOLD) state = STATE.BELOW_THRESHOLD;
        else state = STATE.ABOVE_THRESHOLD;
    }
    
    public String printContours(){
        
        StringBuilder conts = new StringBuilder();
        for(Contour contour : contours)
            conts.append(contour.toString(i0, di, j0, dj)).append("\n\n");
        return conts.toString();
    }
     
    /**
     * 
     * @param contour       List to add edge points to.
     * @param edges         edge map
     * @param isDone        array storing 'done' status for each pixel
     * @param initSearch    Initial search direction.
     * @param i             Coordinates of current point
     * @param j 
     */
    private void growContour(Contour contour,
                             DIRECTION initSearch,
                             int i, int j){
        
        // When above the threshold (on inside of contour), we do 8-neighbour
        // search, incrementing direction by one each loop.
        int step = 1;
        
        // When below the threshold, we do 4-neighbour search.
        if(state==STATE.BELOW_THRESHOLD) step = 2;
        
        
        // Search for edge sections surrounding this point. 
        // Incrementing/decrementing 'search' sets anticlockwise/clockwise
        // direction of traversal.
        for(int search = initSearch.ordinal(); ; search+=step){
            
            // get search direction from integer
            DIRECTION d = DIRECTION.values()[search%8];
            
            // Coordinates of the pixel in this direction
            int x = i+d.i_off;
            int y = j+d.j_off;
            
            // If this is the original start point, then we are done...
            if(contour.isStart(new int[]{i,j}, d)){
                // Final point to connect up line to form closed loop.
                contour.addPoint(j-0.5+((double)d.j_off)/2.0,
                                 i-0.5+((double)d.i_off)/2.0);
                return;
            }
            
            // Contour is crossed in this direction
            if((padded[x][y] > threshold && state == STATE.BELOW_THRESHOLD) ||
               (padded[x][y] < threshold && state == STATE.ABOVE_THRESHOLD)){
            
                // Contour crossed
                
                // Is this a horizontal search direction? If so, mark the 
                // boundary as having been checked.
                if(d.isWest()) BOUNDARY[i-1][j-1] = true;
                if(d.isEast()) BOUNDARY[i-1][j]   = true;
                
                // Don't mark diagonal points
                if(d.isLateral())
                    // Log midpoint midpoint and search on...
                    contour.addPoint(j-0.5+((double) d.j_off)/2.0,
                                     i-0.5+((double) d.i_off)/2.0);
            }
            
            // Contour is not crossed in this direction
            else if((padded[x][y] > threshold && state == STATE.ABOVE_THRESHOLD) ||
                    (padded[x][y] < threshold && state == STATE.BELOW_THRESHOLD)){
            
                // Not a contour crossing - jump to this cell and continue search.
                
                // Set initial search direction from next cell, in terms of
                // rotational offset from transit direction.
                DIRECTION opposite = d.getOpposite();
                
                DIRECTION initialSearch;
                
                if(opposite.isLateral())
                    // If horizontal/vertical direction, advance search direction by two.
                    initialSearch = DIRECTION.values()[(opposite.ordinal() + 2) % 8];
                else
                    // If diagonal direction, advance search direction by one.
                    initialSearch = DIRECTION.values()[(opposite.ordinal() + 1) % 8];
                
                growContour(contour, initialSearch, x, y);
                return;
            }
            
        }
        
    }
    
}


class Contour{
    
    // Indices of starting element in array.
    int[] start_point = new int[2];
    DIRECTION start_dir;
    
    
    List<double[]> points;

    public Contour(int[] _start_point, DIRECTION _start_dir){
        
        start_point[0] = _start_point[0];
        start_point[1] = _start_point[1];
        
        start_dir = _start_dir;
        
        points = new LinkedList<double[]>();
    
    }

    public void addPoint(double x, double y){ points.add(new double[]{x,y}); }

    
    public boolean isStart(int[] point, DIRECTION d){
        return point[0]==start_point[0] && point[1]==start_point[1] && d.equals(start_dir);
    }
    
    
    @Override
    public String toString(){
    
        StringBuilder out = new StringBuilder();
        
        for(double[] point : points)
            out.append(point[0]).append("\t").append(point[1]).append("\n");
    
        return out.toString();
    }
    
    /**
     * Print contour coordinates, with units transformed from array index
     * coordinates to user defined scale.
     * @param i0    i coordinate of top left cell
     * @param di    step size in i
     * @param j0    j coordinate of top left cell
     * @param dj    step size in j
     * @return 
     */
    public String toString(double i0, double di, double j0, double dj){
    
        StringBuilder out = new StringBuilder();
        
        for(double[] point : points)
            out.append(point[0]*di + i0).append("\t").append(point[1]*dj + j0).append("\n");
    
        return out.toString();
    }    
    
    
    
}