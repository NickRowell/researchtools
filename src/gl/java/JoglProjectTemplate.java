package gl.java;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Toolkit;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.GLAutoDrawable;
import com.jogamp.opengl.GLCapabilities;
import com.jogamp.opengl.GLEventListener;
import com.jogamp.opengl.GLProfile;
import com.jogamp.opengl.awt.GLJPanel;
import com.jogamp.opengl.util.Animator;

import Jama.Matrix;
import gl.java.util.GLUtils;
import gl.java.util.PerspectiveCamera;
import numeric.geom.dim3.Quaternion;
import numeric.geom.dim3.Triangle3d;
import numeric.geom.dim3.Vector3d;

/**
 * TODO:
 * 
 * 
 * 1) Use shaders to do the rendering
 * 1.5) Include normals for vertices, light source?
 * 2) Define camera position and orientation & model orientation as Vector3d and Matrix fields.
 *    User interaction updates these. They are converted to the GL versions whenever drawing.
 * 
 * 5) Set the frame size according to the pinhole camera image size; don't allow resizes.
 * 
 * Tidy up.
 * Make the interface nice and clean as possible.
 * 
 *
 * @author nrowell
 * @version $Id$
 */
public class JoglProjectTemplate extends GLJPanel implements GLEventListener, MouseListener, MouseWheelListener, MouseMotionListener {
	
	Animator animator;
	  
	/**
	 * The serial version UID.
	 */
	private static final long serialVersionUID = 5765608097799901659L;
	
	// Elements to include:
	
	// Vertex & Fragment shaders
	
	/**
	 * Pinhole camera model; encapsulates the camera intrinsic matrix.
	 */
	PerspectiveCamera cam = new PerspectiveCamera(10, 10, 512, 512);
	
	/**
	 * Camera orientation represented as a 3x3 orthonormal matrix.
	 * The standard pinhole camera model is used, where the camera
	 * boresight points along +Z axis, Y is 'down' and X is 'right' in the image.
	 * 
	 * *** Matrix rotates vectors from the external frame to the camera frame
	 * 
	 */
	private Matrix r_ext_cam;
	
	/**
	 * Camera position (relative to external frame origin, expressed in external frame).
	 * Represented as 3x1 column matrix.
	 * 
	 * *** position of camera origin relative to external frame origin, expressed in external frame
	 * 
	 */
	private Matrix t_ext_cam;
	
	/**
	 * Model orientation (position is fixed at the origin).
	 * 
	 * *** Matrix rotates vectors from the external frame to the model frame
	 * 
	 */
    private Matrix r_ext_model;
    
    /**
     * Model position (relative to external frame origin, expressed in external frame).
     * Represented as 3x1 column matrix.
     * 
     * *** position of model origin relative to external frame origin, expressed in external frame
     * 
     */
    private Matrix t_ext_model;
	
	/**
	 * List of {@link Triangle<Vector3d>}s to draw. The vertices are expressed in the model frame.
	 */
	List<Triangle3d<Vector3d>> triangles = new LinkedList<>();
	
	/**
	 * Main constructor for the {@link JoglProjectTemplate}.
	 */
	public JoglProjectTemplate(GLCapabilities glcaps) {
		
		super(glcaps);

		triangles.addAll(GLUtils.getCube(1.0));
		
        // Initialise model orientation (identity)
        r_ext_model = Matrix.identity(3, 3);
        
        // Initialise model position (at external frame origin)
        t_ext_model = new Matrix(new double[][]{{0}, {0}, {0}});
        
        // Initialise camera orientation
        r_ext_cam = new Matrix(new double[][]{{1, 0, 0},
            								  {0, 0,-1},
            								  {0, 1, 0}});
        
        // Initialise camera position
        t_ext_cam = new Matrix(new double[][]{{0}, {-30}, {0}});
        
        addGLEventListener(this);      // Capture window resize etc and display stuff
        addMouseListener(this);        // Capture mouse clicks (click to drag)
        addMouseMotionListener(this);  // Capture mouse drags (rotate the model)
        addMouseWheelListener(this);	// Capture middle wheel scrolls (zoom in & out)
        
        setPreferredSize(new Dimension(cam.width, cam.height));
        
        animator = new Animator(this);
        animator.start();
	}
	
	
	@Override
	public void display(GLAutoDrawable drawable) {

        /**
         * OpenGL view matrix is constant (camera fixed position/orientation).
         * 
         * OpenGL has camera boresight pointing along -Z axis, Y axis 'up' and
         * X axis 'right'.
         * 
         * Basis vectors of camera frame expressed in model frame can be read off 
         * the columns, so:
         * 
         * Column 1: X axis of camera frame points along +X axis of model frame
         * Column 2: Y axis of camera frame points along -Z axis of model frame
         * Column 3: Z axis of camera frame points along +Y axis of model frame
         * 
         * Position of model origin in camera frame can be read from first three
         * elements of fourth row:
         * 
         * Model origin at (0,0,-range) in camera frame
         * 
         * Alternatively, the rows give the basis vectors of the model frame 
         * expressed in the camera frame, so
         * 
         * Row 1: X axis of model frame points along +X axis of camera frame
         * Row 2: Y axis of model frame points along +Z axis of camera frame
         * Row 3: Z axis of model frame points along -Y axis of camera frame
         * 
         */
//         glVMatrix = new Matrix(new double[][]{{1, 0, 0, 0},
//                                               {0, 0, 1, 0},
//                                               {0,-1, 0, 0},
//                                               {0, 0,-range, 1}});
         
         /** 
          * Corresponding pinhole camera version.
          * 
          * Axes are oriented differently:
          * Camera boresight points along +Z axis, Y is 'down' and X is 'right'
          * 
          * It is also transposed wrt OpenGL convention because matrix operations
          * are applied in a different order.
          * 
          */
//         camVMatrix = new Matrix(new double[][]{{1, 0, 0, 0},
//                                                {0, 0, 1, 0},
//                                                {0,-1, 0, range},
//                                                {0, 0, 0, 1}});
        
        // Convert the camera position and orientation to the convention required by GL.
		
		// Matrix that converts between the GL and pinhole camera conventions for the camera axes
		Matrix gl2Cam = new Matrix(new double[][]{{1, 0, 0}, {0, -1, 0}, {0, 0, -1}});
		
		// We need to construct the matrix that transforms positions from the model frame to
		// the camera frame
		
		// Position of model origin in camera frame
		Matrix t_cam_model = r_ext_cam.times(t_ext_model.minus(t_ext_cam));
		
		// Orientation of model in camera frame; this is the matrix that rotates vectors from
		// the model frame to the camera frame.
		Matrix r_mod_cam = r_ext_cam.times(r_ext_model.transpose());
		
		// Now, compute the minimum and maximum distance to any of the vertices. These values
		// will be used to configure the near and far plane positions.
		double near =  Double.MAX_VALUE;
		double far  = -Double.MAX_VALUE;
		// Loop over all the Triangles
		for(Triangle3d<Vector3d> tri : triangles) {
			
			// Loop over all the vertices
			for(Vector3d vert : new Vector3d[]{tri.v0, tri.v1, tri.v2}) {
				
				// Create a Matrix from this vertex
				Matrix r0 = new Matrix(new double[][]{{vert.x}, {vert.y}, {vert.z}});
				
				Matrix r0_cam = r_mod_cam.times(r0).plus(t_cam_model);
				
				// Get the distance to this point along the Z axis
				double rz = r0_cam.get(2, 0);
				
				if(rz > 0) {
					// Only interested in minimum distance to points in front of camera
					near = Math.min(rz, near);
				}
				far = Math.max(rz, far);
			}
		}
		
		// Transform these to the GL convention for the camera frame axes
		t_cam_model = gl2Cam.times(t_cam_model);
		r_mod_cam = gl2Cam.times(r_mod_cam);
        
        // Build the GL modelview matrix
        Matrix glMVMatrix = Matrix.identity(4, 4);
        glMVMatrix.setMatrix(new int[]{0,1,2}, new int[]{0,1,2}, r_mod_cam);
        glMVMatrix.setMatrix(new int[]{0,1,2}, new int[]{3}, t_cam_model);
        // GL likes matrices transposed relative to common sense
        glMVMatrix = glMVMatrix.transpose();
        
        
        GL2 gl = (GL2) drawable.getGL();
	        
        // Set colour that screen is cleared to
        gl.glClearColor(0.5f, 0.69f, 1.0f, 1.0f);
        
        // Clear the screen to the underlying colour
        gl.glClear( GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT);
        
        // Load model view matrix
        gl.glMatrixMode(GL2.GL_MODELVIEW);
        gl.glLoadIdentity();
        gl.glLoadMatrixf(GLUtils.toFloatArray(glMVMatrix), 0);

        // Load perspective projection matrix
        gl.glMatrixMode(GL2.GL_PROJECTION);
        gl.glLoadIdentity();
        gl.glLoadMatrixf(cam.getGlProjectionMatrix(near, far), 0);
        
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        
        gl.glPolygonMode(GL.GL_FRONT_AND_BACK, GL2.GL_LINE);
        gl.glBegin(GL.GL_TRIANGLES);
        for(Triangle3d<Vector3d> triangle : triangles)
        {
        	// Set colour
        	gl.glColor3f((float)Math.random(), 0, 1);
        	gl.glVertex3f((float)triangle.v0.x,(float)triangle.v0.y,(float)triangle.v0.z);
        	gl.glVertex3f((float)triangle.v1.x,(float)triangle.v1.y,(float)triangle.v1.z);
        	gl.glVertex3f((float)triangle.v2.x,(float)triangle.v2.y,(float)triangle.v2.z);
        }
        gl.glEnd();
        gl.glPolygonMode(GL.GL_FRONT_AND_BACK, GL2.GL_FILL);
        
        gl.glFlush();
	}

	@Override
	public void dispose(GLAutoDrawable arg0) {
		// TODO Auto-generated method stub
		animator.stop();
	}

	@Override
	public void init(GLAutoDrawable drawable) {
		
        GL2 gl = (GL2) drawable.getGL();

        // Load initial projection matrix
        gl.glMatrixMode(GL2.GL_PROJECTION);
        gl.glLoadIdentity();
        gl.glLoadMatrixf(cam.getGlProjectionMatrix(1.0, 10.0), 0);
        gl.glClearDepth(1.0);
        gl.glEnable(GL.GL_DEPTH_TEST);
        gl.glDepthFunc(GL.GL_LEQUAL);
        
        
	}

	@Override
	public void reshape(GLAutoDrawable arg0, int arg1, int arg2, int arg3, int arg4) {
		// TODO Auto-generated method stub
		
	}
	
    // Mouse support.
    
    int prevMouseX = 0;
    int prevMouseY = 0;
    boolean mouseLButtonDown = false;
    boolean mouseCButtonDown = false;
    boolean mouseRButtonDown = false;
    
    @Override
    public void mouseClicked(MouseEvent e) 
    {
    }
    
    @Override
    public void mousePressed(MouseEvent e)
    {
        
        prevMouseX = e.getX();
        prevMouseY = e.getY();
        
        if ((e.getModifiers() & MouseEvent.BUTTON1_MASK) != 0)
        {
            mouseLButtonDown = true;
        }
        if ((e.getModifiers() & MouseEvent.BUTTON2_MASK) != 0)
        {
            mouseCButtonDown = true;
        }
        if ((e.getModifiers() & MouseEvent.BUTTON3_MASK) != 0)
        {
            mouseRButtonDown = true;
        }
        
    }
    
    @Override
    public void mouseReleased(MouseEvent e) 
    {
        
        if ((e.getModifiers() & MouseEvent.BUTTON1_MASK) != 0) 
        {
            mouseLButtonDown = false;
        }
        if ((e.getModifiers() & MouseEvent.BUTTON2_MASK) != 0) 
        {
            mouseCButtonDown = false;
        }
        if ((e.getModifiers() & MouseEvent.BUTTON3_MASK) != 0) 
        {
            mouseRButtonDown = false;
        }
    }
    
    @Override
    public void mouseDragged(MouseEvent e) 
    {
        int x = e.getX();
        int y = e.getY(); 
            
        // Drags with the left button rotate the asteroid/database
        if (mouseLButtonDown) 
        {
            // Mouse motion vector lies in XZ plane (camera positioned on -z axis)
            Vector3d mouse = new Vector3d((float)(prevMouseX - x), 0, (float)(y - prevMouseY ));
            
            // Check that we have recorded motion:
            if(x==prevMouseX && y==prevMouseY)
            {
                // Attempting to calculate rotation will cause an error
            }
            else
            {
                // Magnitude of motion - one pixel = one degree.
                double angle = Math.toRadians(mouse.norm());

                // Rotation axis is perpendicular to motion of mouse and +y axis
                Vector3d axis = mouse.normalise().cross(new Vector3d(0,1,0));

                // Rotation quaternion
                Quaternion rot = new Quaternion(axis, angle);
                r_ext_model = r_ext_model.times(rot.toMatrix());
            }
        }
        
        // Drags with the middle button
        if (mouseCButtonDown) 
        {
            // Zoom according to displacement along Y axis
            // Smaller motions = zoom value closer to 1.0, which equals smaller
            // change to range.
            float zoom = 1.0f - ((float)(y - prevMouseY))/(float)cam.height;
            
            // How does zoom affect range?
            if(Math.signum(zoom) > 0.0) {
            	t_ext_cam.timesEquals(Math.abs(zoom));
            }
            if(Math.signum(zoom) < 0.0) {
            	t_ext_cam.timesEquals(1.0/Math.abs(zoom));
            }
        }
        // Drags with the right button zoom in/out
        if (mouseRButtonDown) 
        {
            
        }
        
        prevMouseX = x;
        prevMouseY = y;
        
    }

	@Override
	public void mouseWheelMoved(MouseWheelEvent e) {
		if (e.getWheelRotation() < 0) {
			// Mouse wheel moved up
			t_ext_cam.timesEquals(1.1);
		} else {
			// Mouse wheel moved down
			t_ext_cam.timesEquals(1.0/1.1);
		}
	}

	@Override
	public void mouseEntered(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}


	@Override
	public void mouseMoved(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}
	


	@Override
	public void mouseExited(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}

	/**
	 * Main application entry point.
	 * @param args
	 * 	The args - ignored.
	 */
	public static void main(String[] args) {

		// Highly recommended when running JOGL applications on Linux.
		GLProfile.initSingleton();
		
		final JFrame frame = new JFrame("JOGL Project Template");

	    // Get the default GLProfile - which represents the profile (version of OpenGL) best suited to
	    // the running platform.
	    GLProfile prof = GLProfile.getDefault();
	    final GLCapabilities glcaps = new GLCapabilities(prof);
	      
		SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
            	frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
                frame.setLayout(new BorderLayout());
                frame.add(new JoglProjectTemplate(glcaps), BorderLayout.CENTER);
                frame.setSize(512, 512);
                frame.pack();
                Dimension screenSize =
           	         Toolkit.getDefaultToolkit().getScreenSize();
           	      Dimension frameSize = frame.getSize();
           	      if (frameSize.width > screenSize.width ) {
           	         frameSize.width = screenSize.width;
           	      }
           	      if (frameSize.height > screenSize.height) {
           	         frameSize.height = screenSize.height;
           	      }
           	      frame.setLocation ((screenSize.width - frameSize.width ) >> 1, (screenSize.height - frameSize.height) >> 1);
                  frame.setVisible(true);
            }
        });
	}


}