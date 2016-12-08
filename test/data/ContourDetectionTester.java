/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package data;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import data.ContourDetection;

/**
 *
 * @author nickrowell
 */
public class ContourDetectionTester {
    
    // Output file for test data
    static BufferedWriter out;
    static{
        try{
           out = new BufferedWriter(new FileWriter(new File("/home/nickrowell/tmp3")));
        }
        catch(IOException ignored){}
    }  
   
    
    
    public static void main(String[] args) throws IOException{
        
        // This is the input array for which edge detection is to be performed
        //
        // In this layout, second index iterates over columns.
        // First index iterates over rows.
        //
        // West is left (column index decremented)
        // South is down (row index incremented)
        // 
        //
        // [17*17]
        //
        // First index iterates over rows (up-down, or Y)
        // Second index iterates over columns (left-right, or X)
        // 
        // origin is at top left
        //
//        double[][] region = new double[][]{{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},   // [0][0:16]
//                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},   // [1][0:16]
//                                           {0,0,0,1,0,1,1,1,1,1,1,1,1,1,1,0,0},   //   ...
//                                           {0,0,1,0,1,0,1,1,1,1,1,1,1,1,1,0,0},
//                                           {0,0,0,1,0,1,1,1,1,0,1,1,1,1,1,0,0},
//                                           {0,0,1,1,1,1,1,1,1,1,0,1,1,1,1,0,0},
//                                           {0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0},   // [6][0:16]
//                                           {0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0},   // [7][0:16]
//                                           {0,0,1,1,1,1,1,0,0,0,0,0,0,1,1,0,0},   // [8][0:16]
//                                           {1,1,1,1,1,1,0,0,0,1,1,1,1,1,1,0,0},
//                                           {0,0,1,1,1,1,0,0,0,1,1,1,1,1,1,0,0},
//                                           {0,0,1,1,1,1,1,0,0,1,1,1,1,1,0,0,0},
//                                           {1,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0},
//                                           {0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0},
//                                           {0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0},
//                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
//                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
        
        
        double[][] region = new double[][]{{0,0,0,0,1,1,1,0},
                                           {0,1,0,0,1,0,1,1},
                                           {0,0,1,0,1,1,0,1},
                                           {0,0,0,0,1,1,1,1}};
         
         
        ContourDetection contourDetection = new ContourDetection(region,0.5,-2,-3,5,5);
            
        out.write(contourDetection.printContours()); 
        out.close();
    
    }
   
    
}
