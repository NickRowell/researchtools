/**
 * Copyright (c) 2012 University of Dundee
 *
 * Name: TestFourier2D.java
 *
 * Purpose: 
 * This class is designed for testing Fourier transform algorithms on images
 * using the JTransforms library.
 * 
 * Language: Java
 *
 * Author: Nicholas Rowell (nickrowell@computing.dundee.ac.uk)
 *
 */
package transforms;

import images.Rendering;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.imageio.ImageIO;
import javax.swing.ImageIcon;
import javax.swing.JLabel;

import transforms.Direction;
import transforms.FFT_2D;
import misc.Complex;


public class TestFourier2D extends javax.swing.JFrame {

    // Images used in testing fourier transform
    BufferedImage untransformed = ImageIO.read(new File("test/transforms/untransformed.png"));
    BufferedImage magnitude     = ImageIO.read(new File("test/transforms/magnitude.gif"));
    BufferedImage logmagnitude  = ImageIO.read(new File("test/transforms/logmagnitude.gif"));
    BufferedImage phase         = ImageIO.read(new File("test/transforms/phase.gif"));
    
    /**
     * Creates new form TestFourier2D
     */
    public TestFourier2D() throws IOException {
        initComponents();
        
        // Load raw images in JLabels for adding to Swing panels in GUI
        JLabel untrans = new JLabel(new ImageIcon(untransformed));
        JLabel mag     = new JLabel(new ImageIcon(magnitude));
        JLabel logmag  = new JLabel(new ImageIcon(logmagnitude));
        JLabel pha     = new JLabel(new ImageIcon(phase));

        // Place in Frame
        add(untrans);
        add(mag);
        add(logmag);
        add(pha);
        
    }

    /**
     * 
     */
    @SuppressWarnings("unchecked")
    private void initComponents() {

        transformButton = new javax.swing.JButton();

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        getContentPane().setLayout(new java.awt.FlowLayout());

        transformButton.setText("transform >>");
        transformButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                transformButtonActionPerformed(evt);
            }
        });
        getContentPane().add(transformButton);

        pack();
    }

    private void transformButtonActionPerformed(java.awt.event.ActionEvent evt) {
        
        // Read raw image into Complex array
        int[][] fint = Rendering.imToArray(untransformed);
        
        int w = fint.length;
        int h = fint[0].length;
        
        // Read to a 1D row-packed real-only double array
        double[] X = new double[w*h];
        for(int i=0; i<w; i++)
        {
        	for(int j=0; j<h; j++)
        	{
        		X[i*h + j] = fint[i][j];
        	}
        }
        
        // Take Fourier transform
        Complex[] x = FFT_2D.transformReal(X, w, h, Direction.FORWARD);
        
        // Build magnitude and log-magnitude arrays
        int[][] mag    = new int[w][h];
        int[][] logmag = new int[w][h];
        
        double max_mag = Integer.MIN_VALUE;
        double max_logmag = Integer.MIN_VALUE;
        
        for(int i=0; i<w; i++)
        {
        	for(int j=0; j<h; j++)
        	{
        		mag[i][j]    = (int)Math.rint(x[i*h + j].abs());
        		logmag[i][j] = (int)Math.rint(Math.log(x[i*h + j].abs()));
        		
        		if(mag[i][j] > max_mag) max_mag = mag[i][j];
        		if(logmag[i][j] > max_logmag) max_logmag = logmag[i][j];
        	}
        }
        
        // Normalise arrays to 0:255 range
        for(int i=0; i<w; i++)
        {
        	for(int j=0; j<h; j++)
        	{
        		mag[i][j]    *= (255.0/max_mag);
        		logmag[i][j] *= (255.0/max_logmag);
        		
        	}
        }
        
        
        // Get visualisations of transform results
        JLabel myMagnitude = new JLabel(new ImageIcon(Rendering.arrayToIm(mag)));
        JLabel myLogMagnitude = new JLabel(new ImageIcon(Rendering.arrayToIm(logmag)));
        JLabel myPhase = new JLabel(new ImageIcon(Rendering.arrayToIm(logmag)));
        
        // Add to frame
        add(myMagnitude);
        add(myLogMagnitude);
        add(myPhase);
        
        repaint();
        
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        /*
         * Set the Nimbus look and feel
         */
        /*
         * If Nimbus (introduced in Java SE 6) is not available, stay with the
         * default look and feel. For details see
         * http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException ex) {
            java.util.logging.Logger.getLogger(TestFourier2D.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(TestFourier2D.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(TestFourier2D.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(TestFourier2D.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }

        /*
         * Create and display the form
         */
        java.awt.EventQueue.invokeLater(new Runnable() {

            public void run() {
                try {
                    new TestFourier2D().setVisible(true);
                } catch (IOException ex) {
                    Logger.getLogger(TestFourier2D.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        });
    }
    private javax.swing.JButton transformButton;
}
