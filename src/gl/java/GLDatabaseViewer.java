package GUI;

import Clustering.Cluster;
import Jama.Matrix;
import Triangulation.IntersectionTest;
import Triangulation.Triangle;
import Utils.*;
import gl.java.util.GLUtils;
import images.ImageUtil;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.FloatBuffer;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import javax.media.opengl.GL;
import javax.media.opengl.GLAutoDrawable;
import javax.media.opengl.GLEventListener;
import javax.media.opengl.GLJPanel;
import javax.media.opengl.glu.GLU;
import javax.swing.BorderFactory;
import javax.swing.border.BevelBorder;
import uk.ac.dundee.spacetech.pangu.ClientLibrary.ClientConnection;
import uk.ac.dundee.spacetech.pangu.ClientLibrary.ConnectionFactory;
import uk.ac.dundee.spacetech.pangu.ClientLibrary.ValidPoint;
import uk.ac.dundee.spacetech.pangu.ClientLibrary.Vector2D;

/**
 * Class provides a visualisation of point cloud and clusters located within
 * it. We draw clusters by projecting their covariance matrices into the image
 * plane and drawing a confidence ellipse.
 * 
 * We translate and project the point cloud using inbuilt OpenGL methods.
 * The clusters are rendered as ellipses by projecting the covariance ellipsoids
 * into the image plane and drawing a line of constant probability density.
 * This projection is done manually, which means we need to store the 
 * view and projection matrices in both opengl and standard pinhole camera
 * formats, which are rotated wrt each other and have different conventions
 * for the camera boresight and image plane.
 * 
 * +++ Concurrent access of clusters & triangles lists +++
 * This class maintains it's own internal lists of Cluster and Triangle objects
 * that it uses to render frames at high speed. From time to time as the rest 
 * of the software runs, ProcPointCloudPanel will update these internal lists
 * with new output from the database generation algorithm.
 * 
 * We need to ensure that the present class is not busy iterating over one of
 * the lists at the same time that ProcPointCloudPanel attempts to clear and
 * refill the list.
 * 
 * The way to do this is to use Collections.synchronizedList() to make the 
 * lists stored in this class synchronised for atomic operations (.add(), 
 * .clear() etc), which ensures that atomic operations on the lists performed 
 * by the threads running inside this class don't interleave with atomic
 * operations performed by other threads on the same lists.
 * It is still vital that the user manually synchronise on the
 * lists when iterating over them however - hence the synchonized blocks 
 * surrounding the iterations in the display() method below.
 * 
 * 
 * 
 * 
 * @author nickrowell
 */
public class GLDatabaseViewer
extends GLJPanel
implements GLEventListener, MouseListener, MouseMotionListener
{
    /** Reference to main State object. */
    private final State state;
    
    /** FloatBuffer used to transfer (fixed) point cloud vertex data to GPU. */
    private final FloatBuffer pointCloudVertices;
    private final FloatBuffer pointCloudColours;
    private final int N_VERTS;
    /** Number of bytes in a java float. */
    private static final int FLOAT_SIZE_BYTES = 4;
    
    // Local copies of the database objects, updated each time the database is
    // recalculated. These objects are used to render the database. This avoids
    // concurrency problems that could arise if we tried to render the main
    // database objects while they were being recalculated.
    final List<Cluster> clusters;
    final List<Triangle2d> triangles;
    
    /** GL Utilities instance. */
    private final GLU glu = new GLU();
    // Texture handle (used to draw PANGU image into frame)
    int mTextureID;    
    
    /** View matrix in OpenGL standard frame. */
    private Matrix glVMatrix;
    /** Corresponding pinhole camera version. */
    private Matrix camVMatrix;
    /** Model matrix. */
    private Matrix mMatrix;
     /** Near and far plane distances. */
    private double near;
    private double far;
     /** Camera projection matrix (for performing manual projections). */
    CameraProperties camera;   
    
    /** 
     * Handle to triangle currently selected in viewer by right mouse click.
     * When we select the triangle, we check that it is part of only one
     * tetrahedron (this is a condition for triangles that are part of the
     * external hull of the triangulation), so when it comes to deleting the
     * tetrahedron we can access it via SELECTED.tetras.get(0)
     */
    public Triangle2d selected_triangle;
    
    // List of all selected clusters
    public List<Cluster> selected_clusters = new LinkedList<Cluster>();
    
    // Configure various database visualisation options
    public boolean DRAW_INTERNAL_TRIANGLES = false;
    public boolean DRAW_OBSCURED_CLUSTERS  = true;
    public boolean DRAW_OBSCURED_MESH      = true;
    
    /** PANGU server socket connection. */
    ClientConnection pangu;
    /** Byte array contains the current PANGU image. */
    byte[] raw;
    /** Boolean flag indicates whether to draw an image of the asteroid. */
    public boolean show_asteroid = false;
    
    
    
    /** Boolean flag used to indicate that next frame should be written to disk. */
    public boolean save_image = false;
    /** Integer used to increment file names. */
    int filename_index = 0;
    
    // Size of (square) visualisation - can differ from PANGU server.
    private final int imsize;
    
    public GLDatabaseViewer(State astate, int pimsize)
    {
        
        super();
        
        state  = astate;
        imsize = pimsize;
        
        clusters = Collections.synchronizedList(new LinkedList<Cluster>());
        triangles = Collections.synchronizedList(new LinkedList<Triangle2d>());
        
        // Use same FOV as PANGU server but a different image size
        camera = new CameraProperties(state.fov, state.fov, imsize, imsize);
        
        // Configure visible panel
        setPreferredSize(new Dimension(imsize, imsize));
        setMinimumSize(new Dimension(imsize, imsize));
        setLayout(new BorderLayout());
        setBorder(BorderFactory.createBevelBorder(BevelBorder.LOWERED));
        setVisible(true);
        
        // Copy point cloud to local storage
        int[] N = {0};
        pointCloudVertices = state.pointCloud.syncToFloatBuffer(N);
        pointCloudVertices.position(0);
        N_VERTS = N[0];
        
        // Make a colour buffer that specifies the colour of each vertex
        float[] colours_arr = new float[N_VERTS*3];
        for(int i=0; i<colours_arr.length; i+=3)
        {
            colours_arr[i+0] = 0;
            colours_arr[i+1] = 0;
            colours_arr[i+2] = 0;
        }
        pointCloudColours = ByteBuffer.allocateDirect(colours_arr.length * FLOAT_SIZE_BYTES).order(ByteOrder.nativeOrder()).asFloatBuffer();
        pointCloudColours.put(colours_arr);
        pointCloudColours.position(0);
        
        setViewMatrixAndClipPlanes(state.range);
        
        // Model matrix initialised to identity.
        mMatrix = new Matrix(new double[][]{{ 1, 0, 0, 0},
                                            { 0, 1, 0, 0},
                                            { 0, 0, 1, 0},
                                            { 0, 0, 0, 1}});
        
        // Open socket connection to PANGU server.
        try
        {
            pangu = ConnectionFactory.makeConnection("localhost", 10363);
            updateAsteroidImage();
        }
        catch(IOException ioe)
        {
            System.err.println("Couldn't connect to PANGU server on localhost "
                    + "port 10363");
            System.exit(1);
        }
        
        addGLEventListener(this);      // Capture window resize etc and display stuff
        addMouseListener(this);        // Capture mouse clicks
        addMouseMotionListener(this);  // Capture mouse drags
        
    }
    
    public void cleanUp()
    {
        try
        {
            pangu.stop();
        } catch (IOException ex)
        {
            // Ignored
        }
    }
    
    /**
     * There are two view matrices representing the OpenGL and pinhole camera
     * conventions for the same viewpoint. We stay at a fixed orientation and
     * line of sight, but can zoom in & out along the camera boresight 
     * according to the range passed into this method.
     * @param range 
     */
    private void setViewMatrixAndClipPlanes(double range)
    {
        
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
        glVMatrix = new Matrix(new double[][]{{1, 0, 0, 0},
                                              {0, 0, 1, 0},
                                              {0,-1, 0, 0},
                                              {0, 0,-range, 1}});
        
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
        camVMatrix = new Matrix(new double[][]{{1, 0, 0, 0},
                                               {0, 0, 1, 0},
                                               {0,-1, 0, range},
                                               {0, 0, 0, 1}});
        
        // Should be appropriate at all ranges
        near = range/8.0;
        // Ad-hoc value for far plane distance
        far  = 2*range + 1000;
        
    }
    
    /** Take local copies of database objects from State variables. */
    public final void updateVisualisationObjects()
    {
        clusters.clear();
        clusters.addAll(state.database.clusters);
        triangles.clear();
        triangles.addAll(state.database.del.tris);
        selected_triangle = null;
    }

    @Override
    public void init(GLAutoDrawable drawable)
    {
        
        GL gl = drawable.getGL();
        
        // Load initial projection matrix
        gl.glMatrixMode(GL.GL_PROJECTION);
        gl.glLoadIdentity();
        glu.gluPerspective(camera.fovy, camera.aspect, near, far);
        
        // New stuff for DatabasePruningViewer
        // Create texture
        int[] textures = new int[1];

        // Get a suitable name (an integer) for a texture object.
        gl.glGenTextures(1, textures, 0);

        // Get texture name out of array (integer other than zero).
        mTextureID = textures[0];
        
        // Create new texture object: specify target (2D texture) and previously
        // unused texture name (mTextureID)
        gl.glBindTexture(GL.GL_TEXTURE_2D, mTextureID);
        gl.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_LINEAR);
        gl.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_NEAREST);
        gl.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_S, GL.GL_REPEAT);
        gl.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_T, GL.GL_REPEAT);
        
    }

    public void dispose(GLAutoDrawable glad)
    {
    }

    @Override
    public void display(GLAutoDrawable drawable)
    {
        GL gl = drawable.getGL();
        
        // Set colour that screen is cleared to
        gl.glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        // Clear the screen to the underlying colour
        gl.glClear( GL.GL_DEPTH_BUFFER_BIT | GL.GL_COLOR_BUFFER_BIT);
        
        // Chain transformations (model rotation plus view roto-translation)
        Matrix glMVMatrix = mMatrix.times(glVMatrix);
        // Pinhole camera version. Note operations applied in a different order,
        // so mMatrix is transposed.
        Matrix camMVMatrix = camVMatrix.times(mMatrix.transpose());
        

        if(show_asteroid)
        {
            drawAsteroidImage(gl);
            drawDatabase(gl, camMVMatrix, glMVMatrix);
        }
        else
        {
            drawPointCloud(gl, glMVMatrix);
            drawDatabase(gl, camMVMatrix, glMVMatrix);
        }
        
        gl.glFlush();
        
        if(save_image)
        {
            // Allocate float buffer to store pixels
            FloatBuffer fb = FloatBuffer.wrap(new float[imsize*imsize*3]);
            
            gl.glPixelStorei(GL.GL_PACK_ALIGNMENT, 1);
            // Copy current frame buffer pixels in RGB float format
            gl.glReadPixels(0, 0, imsize, imsize, GL.GL_RGB, GL.GL_FLOAT, fb);
            
            // Convert to 24-bit integer pixel array. Must reverse rows because
            // OpenGL has origin at lower left of image, whereas the upper left
            // is assumed normally for images.
            int[] pixels = new int[imsize * imsize];
            for (int r=0; r < imsize; r++)
                for (int c=0; c < imsize; c++)
                {
                    // Index in opengl image
                    int i = r*imsize + c;
                    
                    // Index for pixel in vertically-flipped image
                    int j = (imsize-1-r)*imsize + c;
                    
                    // Need to map floats from [0:1] to [0:255] range and cast to ints
                    pixels[j] = ((int)(fb.get(i*3+0)*255.0f) << 16) +
                                ((int)(fb.get(i*3+1)*255.0f) <<  8) +
                                ((int)(fb.get(i*3+2)*255.0f));
                    
                }
            // Use Java BufferedImage to write image to disk
            try
            {
                BufferedImage image = new BufferedImage(imsize,imsize,BufferedImage.TYPE_INT_RGB);
                image.setRGB(0, 0, imsize, imsize, pixels, 0, imsize);
                String filename = String.format("image_%05d.png", filename_index);
                javax.imageio.ImageIO.write(image, "png", new File(state.output,filename));
            }
            catch (IOException ioe)
            {
                System.out.println("Save image exception: " +ioe.getMessage());
            }
            
            filename_index++;
            save_image = false;
        }
        
    }
    
    /**
     * GL drawing commands for the point cloud.
     * @param gl         GL instance
     * @param glMVMatrix GL model view matrix
     */
    public void drawPointCloud(GL gl, Matrix glMVMatrix)
    {
        gl.glMatrixMode(GL.GL_MODELVIEW);
        gl.glLoadIdentity();
        gl.glLoadMatrixf(GLUtils.toFloatArray(glMVMatrix), 0);
        
        // Load initial perspective projection matrix
        gl.glMatrixMode(GL.GL_PROJECTION);
        gl.glLoadIdentity();
        glu.gluPerspective(camera.fovy, camera.aspect, near, far);
        
        
        
        
        // Draw point cloud from vertex array
        gl.glEnableClientState(GL.GL_VERTEX_ARRAY);
        gl.glEnableClientState(GL.GL_COLOR_ARRAY);
        gl.glPointParameterf(GL.GL_POINT_SIZE, 25);
        gl.glEnable(GL.GL_POINT_SMOOTH);
        gl.glVertexPointer(3, GL.GL_FLOAT, 3*FLOAT_SIZE_BYTES, pointCloudVertices);
        gl.glColorPointer(3, GL.GL_FLOAT, 3*FLOAT_SIZE_BYTES, pointCloudColours);
        gl.glDrawArrays(GL.GL_POINTS, 0, N_VERTS);
        gl.glDisable(GL.GL_POINT_SMOOTH);
        gl.glDisableClientState(GL.GL_COLOR_ARRAY);
        gl.glDisableClientState(GL.GL_VERTEX_ARRAY);
    }
    
    public void drawAsteroidImage(GL gl)
    {
        
        // Orthogonal projection to draw 2D texture
        gl.glMatrixMode(GL.GL_PROJECTION);
        gl.glLoadIdentity();
        glu.gluOrtho2D(0.0, state.imsize, 0.0, state.imsize);

        gl.glMatrixMode(GL.GL_MODELVIEW);
        gl.glLoadIdentity();
        
        // Number of bytes taken up by PPM header
        int header = ImageUtil.getPPMHeaderSize(raw);
        // Allocate byte buffer of the correct size
        ByteBuffer imageBuf = ByteBuffer.allocate(raw.length-header);
        // Efficient transfer of data
        imageBuf.put(raw, header, raw.length-header);
        imageBuf.rewind();
        
        gl.glBindTexture(GL.GL_TEXTURE_2D, mTextureID);
        // Load texture from buffer
        gl.glTexImage2D(GL.GL_TEXTURE_2D, 0, GL.GL_RGB, state.imsize, state.imsize,
                0, GL.GL_RGB, GL.GL_UNSIGNED_BYTE, imageBuf);
        
        // Draw a textured quad
        gl.glEnable(GL.GL_TEXTURE_2D);
        gl.glTexEnvf(GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_MODULATE);
        gl.glBindTexture(GL.GL_TEXTURE_2D, mTextureID);
        gl.glPolygonMode(GL.GL_FRONT_AND_BACK, GL.GL_FILL);
        gl.glBegin(GL.GL_QUADS);
            gl.glColor3f(1,1,1);
            gl.glTexCoord2f(0.0f, 1.0f); gl.glVertex2f(0f, 0f);
            gl.glTexCoord2f(0.0f, 0.0f); gl.glVertex2f(0f, state.imsize);
            gl.glTexCoord2f(1.0f, 0.0f); gl.glVertex2f(state.imsize, state.imsize);
            gl.glTexCoord2f(1.0f, 1.0f); gl.glVertex2f(state.imsize, 0f);
        gl.glEnd();
        gl.glDisable(GL.GL_TEXTURE_2D);
        
    }
    
    // Note that to avoid the cluster & triangle lists being reset mid-way through
    // the display() method, we must synchronise on the lists for the entire
    // duration of the method.
    public void drawDatabase(GL gl, Matrix camMVMatrix, Matrix glMVMatrix)
    {
        // Perspective projection of database
        gl.glMatrixMode(GL.GL_MODELVIEW);
        gl.glLoadIdentity();
        gl.glLoadMatrixf(GLUtils.toFloatArray(glMVMatrix), 0);

        // Load initial perspective projection matrix
        gl.glMatrixMode(GL.GL_PROJECTION);
        gl.glLoadIdentity();
        glu.gluPerspective(camera.fovy, camera.aspect, near, far);
        
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        
        // Manually synchronize on clusters list until we are finished sending
        // all geometry to the rendering pipeline. This avoids the cluster list
        // being reset mid way through, which would reset the camera frame
        // positions to null and cause later operations to fail.
        // 
        // I don't think it's necessary to also synchronize on triangles list,
        // because the two are always refilled together.
        synchronized(clusters)
        { 
            // Get pure rotation (used to rotate cluster covariance matrices).
            // Note that OpenGL transformation matrices are transposed wrt the
            // conventional format.
            Matrix camMVRot = camMVMatrix.getMatrix(0, 2, 0, 2);
        
            // Now transform all Clusters to the camera frame
            for(Cluster cluster : clusters)
            {
                // Homogenous position vector of cluster in camera frame. Fourth
                // component of this is always 1.0
                Matrix pos = camMVMatrix.times(cluster.getHomogenousPosition());
                cluster.X_cam = new Vector3d(pos.get(0,0),pos.get(1,0),pos.get(2,0));
                cluster.cov_X_cam = camMVRot.times(cluster.cov_X_w).times(camMVRot.transpose());
                // Check if cluster is obscured from the current viewpoint
                cluster.isObscured = clusterObscured(cluster, triangles);
            }
        
            // Internal triangles drawn in wireframe
            gl.glPolygonMode(GL.GL_FRONT_AND_BACK, GL.GL_LINE);
            
            gl.glBegin(GL.GL_TRIANGLES);
            // Draw internal triangles first, so external triangles overwrite
            // the internal ones that share the same edges.
            if(DRAW_INTERNAL_TRIANGLES)
            {
                for(Triangle2d triangle : triangles)
                {
                    if(!triangle.isExternal())
                    {
                        // Set colour for internal triangles
                        gl.glColor3f(0.5f,0.5f,1);
                        gl.glVertex3f((float)triangle.v0.x,(float)triangle.v0.y,(float)triangle.v0.z);
                        gl.glVertex3f((float)triangle.v1.x,(float)triangle.v1.y,(float)triangle.v1.z);
                        gl.glVertex3f((float)triangle.v2.x,(float)triangle.v2.y,(float)triangle.v2.z);
                    }
                }
            }
            gl.glEnd();
            
            // External triangles drawn first in wireframe
            gl.glBegin(GL.GL_TRIANGLES);
            // Now (always) draw external triangles
            for(Triangle2d triangle : triangles)
            {
                // Only draw external triangles in the case where we are drawing
                // ALL surface, or where the triangle is unobscured.
                if(triangle.isExternal() && 
                   ((DRAW_OBSCURED_MESH || DRAW_INTERNAL_TRIANGLES) ||
                    (!triangleObscured(triangle))))
                {
                    // Set colour
                    gl.glColor3f(1, 0, 1);
                    gl.glVertex3f((float)triangle.v0.x,(float)triangle.v0.y,(float)triangle.v0.z);
                    gl.glVertex3f((float)triangle.v1.x,(float)triangle.v1.y,(float)triangle.v1.z);
                    gl.glVertex3f((float)triangle.v2.x,(float)triangle.v2.y,(float)triangle.v2.z);
                }
            }
            gl.glEnd();
            
            // External triangles drawn as transparent polygons
            gl.glPolygonMode(GL.GL_FRONT_AND_BACK, GL.GL_FILL);
            gl.glEnable(GL.GL_BLEND);
            gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);
            
            gl.glBegin(GL.GL_TRIANGLES);
            // Render external triangles
            for(Triangle2d triangle : triangles)
            {
                // Only draw external triangles in the case where we are drawing
                // ALL surface, or where the triangle is unobscured.
                if(triangle.isExternal() && 
                   ((DRAW_OBSCURED_MESH || DRAW_INTERNAL_TRIANGLES) ||
                    (!triangleObscured(triangle))))
                {
                    // Set colour
                    gl.glColor4f(1, 0, 1, 0.2f);
                    gl.glVertex3f((float)triangle.v0.x,(float)triangle.v0.y,(float)triangle.v0.z);
                    gl.glVertex3f((float)triangle.v1.x,(float)triangle.v1.y,(float)triangle.v1.z);
                    gl.glVertex3f((float)triangle.v2.x,(float)triangle.v2.y,(float)triangle.v2.z);
                }
            }
            gl.glEnd();
            gl.glDisable(GL.GL_BLEND);
            
            // Now draw currently selected triangle, if one exists.
            if(selected_triangle != null)
            {

                // Selected triangles drawn as transparent polygons
                gl.glPolygonMode(GL.GL_FRONT_AND_BACK, GL.GL_FILL);
                gl.glEnable(GL.GL_BLEND);
                gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);
                
                // Draw internal triangles first, so that external ones
                // overwrite internal ones
                gl.glBegin(GL.GL_TRIANGLES);
                for(Triangle2d triangle : selected_triangle.tetras.get(0).tris)
                {
                    if(!triangle.isExternal())
                    {
                        gl.glColor4f(1,0,0,0.2f);
                        gl.glVertex3f((float)triangle.v0.x,(float)triangle.v0.y,(float)triangle.v0.z);
                        gl.glVertex3f((float)triangle.v1.x,(float)triangle.v1.y,(float)triangle.v1.z);
                        gl.glVertex3f((float)triangle.v2.x,(float)triangle.v2.y,(float)triangle.v2.z);
                    }
                }
                
                for(Triangle2d triangle : selected_triangle.tetras.get(0).tris)
                {
                    if(triangle.isExternal())
                    {
                        gl.glColor4f(0,1,0,0.2f);
                        gl.glVertex3f((float)triangle.v0.x,(float)triangle.v0.y,(float)triangle.v0.z);
                        gl.glVertex3f((float)triangle.v1.x,(float)triangle.v1.y,(float)triangle.v1.z);
                        gl.glVertex3f((float)triangle.v2.x,(float)triangle.v2.y,(float)triangle.v2.z);
                    }
                }
                gl.glEnd();
                gl.glDisable(GL.GL_BLEND);
                
                // Now draw selected triangles in wireframe
                gl.glPolygonMode(GL.GL_FRONT_AND_BACK, GL.GL_LINE);
                
                // Draw internal triangles first, so that external ones
                // overwrite internal ones
                gl.glBegin(GL.GL_TRIANGLES);
                for(Triangle2d triangle : selected_triangle.tetras.get(0).tris)
                {
                    if(!triangle.isExternal())
                    {
                        // Get surface normal for triangle
                        //Vector3d norm = selected_triangle.tetras.get(0).getOutwardNormalForTriangle(triangle);
                        gl.glColor3f(1,0,0);
                        gl.glVertex3f((float)triangle.v0.x,(float)triangle.v0.y,(float)triangle.v0.z);
                        gl.glVertex3f((float)triangle.v1.x,(float)triangle.v1.y,(float)triangle.v1.z);
                        gl.glVertex3f((float)triangle.v2.x,(float)triangle.v2.y,(float)triangle.v2.z);
                    }
                }
                
                for(Triangle2d triangle : selected_triangle.tetras.get(0).tris)
                {
                    if(triangle.isExternal())
                    {
                        gl.glColor3f(0,1,0);
                        gl.glVertex3f((float)triangle.v0.x,(float)triangle.v0.y,(float)triangle.v0.z);
                        gl.glVertex3f((float)triangle.v1.x,(float)triangle.v1.y,(float)triangle.v1.z);
                        gl.glVertex3f((float)triangle.v2.x,(float)triangle.v2.y,(float)triangle.v2.z);
                    }
                }
                gl.glEnd();
            }
            
            // Switch to orthogonal projection to draw manually-projected clusters
            // into frame.
            gl.glMatrixMode(GL.GL_PROJECTION);
            gl.glLoadIdentity();
            glu.gluOrtho2D(0.0, imsize, 0.0, imsize);

            gl.glMatrixMode(GL.GL_MODELVIEW);
            gl.glLoadIdentity();
            
//            // Now draw cluster ellipses into image
//            for(Cluster cluster : clusters)
//            {
//                // Don't draw weak clusters
//                if(cluster.WEAK) continue;
//                
//                // Cluster selected by right click in cluster processing
//                boolean cluster_is_selected = selected_clusters.contains(cluster);
//                
//                // Cluster selected indirectly by being attached to a selected tetrahedron
//                if(selected_triangle != null)
//                {
//                    Cluster[] s = selected_triangle.tetras.get(0).verts;
//                    cluster_is_selected = (cluster == s[0] ||
//                                           cluster == s[1] ||
//                                           cluster == s[2] ||
//                                           cluster == s[3]);
//                }
//                
//                // Apply visibility test to clusters if surface mesh exists...
//                if(state.database.TRIANGULATION_EXISTS)
//                    // ... and we are culling obscured
//                    // objects AND NOT drawing internal triangles
//                    if(!DRAW_OBSCURED_CLUSTERS && !DRAW_INTERNAL_TRIANGLES)
//                        if(cluster.isObscured && !cluster_is_selected)
//                            continue;
//
//                // Project into image plane, with origin at top left.
//                float[] pos_im = camera.projectVector(cluster.X_cam);
//
//                // Optionally draw mean cluster position
//    //            gl.glBegin(GL2.GL_POINTS);
//    //            gl.glColor3f(1,0,0);
//    //            // OpenGL has origin at bottom left, so need to flip y coordinate
//    //            gl.glVertex2f((float)pos_im[0],height-(float)pos_im[1]);
//    //            gl.glEnd();
//
//                // Project cluster covariance into image plane to get covariance 
//                // matrix in pixels
//                Matrix cov_im = camera.projectMatrix(cluster.X_cam, cluster.cov_X_cam);
//
//                // Draw ellipse circumference as a series of vertices rendered
//                // as a GL_LINE_LOOP
//                float[] ellipse = Image.drawEllipse(pos_im, cov_im, state.confidence);
//                
//                // Mark selected clusters in some non-colour way
//                //
//                // 2014.08.11: Draw ALL clusters with thicker lines, in order to produce
//                // images suitable for putting in a paper.
////                if(cluster_is_selected)
////                {
////                    gl.glLineWidth(3.0f);
////                    gl.glEnable(GL.GL_LINE_SMOOTH);
////                    
////                }
//                gl.glLineWidth(3.0f);
//                gl.glEnable(GL.GL_LINE_SMOOTH);
//                
//                // Set cluster colour
//                if(state.database.TRIANGULATION_EXISTS)
//                {
//                    if(cluster.EXTERNAL)
//                    {   
//                        // If external clusters are connected to the surface
//                        // mesh by only one tetrahedron, we colour them green.
//                        if(cluster.tetras.size()==1)
//                            gl.glColor3f(0,1,0);
//                        // otherwise, regular external clusters are blue
//                        else gl.glColor3f(0,0,1);
//                    }  
//                    else
//                        // internal clusters red, and also those that are
//                        // separated from the mesh (e.g. by deletion of tetrahedra).
//                        gl.glColor3f(1,0,0);
//                }
//                else
//                {
//                    if(cluster_is_selected) gl.glColor3f(0,1,0);
//                    else                    gl.glColor3f(0,0,1);
//                }
//                
//                gl.glBegin(GL.GL_LINE_LOOP);
//                for(int i=0; i<ellipse.length; i+=2)
//                {
//
//                    // OpenGL has origin at bottom left, so need to flip y coordinate
//                    gl.glVertex2f((float)ellipse[i+0],imsize-(float)ellipse[i+1]);
//                }
//                gl.glEnd();
//                
//                // Turn off any special rendering options used for selected clusters
////                if(cluster_is_selected)
////                {
////                    gl.glLineWidth(1.0f);
////                    gl.glDisable(GL.GL_LINE_SMOOTH);
////                }
//                
//            }
            
        }
    }
    
    /**
     * Called in event that bit-depth of screen changes, or window is dragged
     * onto another screen in multiple-screen setups. Not required to 
     * implement this method.
     * @param drawable
     * @param modeChanged
     * @param deviceChanged 
     */
    @Override
    public void displayChanged(GLAutoDrawable drawable, boolean modeChanged, boolean deviceChanged)
    {
    }
    @Override
    public void reshape(GLAutoDrawable drawable, int x, int y, int new_width, int new_height) 
    {
    }
    
    private boolean triangleObscured(Triangle2d triangle)
    {
        // If triangle surface normal points away from the camera, then the
        // triangle is always obscured (lies on far side of asteroid)
        Vector3d camN = triangle.getNormalCamFrame();

        // Take dot product with camera-frame position vector
        // for any one of it's vertices.
        double d = camN.dot(triangle.v0.X_cam);

        // If d is positive, then surface normal points in the same direction
        // as vertex position vector, and we can't see the outward face of the
        // triangle. In other words, it's obscured.
        if (d>0) return true;
        
        // Triangle surface normal points towards camera - the triangle can 
        // still be obscured by intervening terrain however. To resolve this
        // case we need to check the visibility of each cluster that it is
        // attached to
        
        // Count how many of the attached landmarks are visible
        int N_CLUSTERS_VIS = 3;
        if(triangle.v0.isObscured) N_CLUSTERS_VIS--;
        if(triangle.v1.isObscured) N_CLUSTERS_VIS--;
        if(triangle.v2.isObscured) N_CLUSTERS_VIS--;
        
        // Can adjust this condition.
        return (N_CLUSTERS_VIS == 0);  // Zero visible clusters equals obscured triangle
        
    }
    
    /**
     * Main geometric test that determines if a given landmark is obscured by
     * any foreground triangles in the surface mesh.
     * @param cluster
     * @param triangles
     * @return 
     */
    private boolean clusterObscured(Cluster cluster, List<Triangle2d> triangles)
    {
        // Distance to feature along LOS
        double distance = cluster.X_cam.norm();
        // Unit vector towards feature, in camera frame
        Vector3d camPhat = cluster.X_cam.mult(1.0/distance);

        // Manually synchronize on triangles list before iterating
        synchronized(triangles)
        {
            // List of Triangles forming database...
            for(Triangle2d triangle : triangles)
            {
                // Only need to test triangles that form part of the external hull
                if(!triangle.isExternal()) continue;
                
                // If this triangle is attached to the present cluster, don't
                // check it for intersection
                if(cluster.tris.contains(triangle))
                {
                    continue;
                }

                // Otherwise, check if the LOS to the cluster intersects this
                // triangle.
                IntersectionTest intersection = new IntersectionTest(camPhat, triangle);

                // Note that this would indicate clusters as obscured in the case
                // where they and the triangle lie behind the camera, and the
                // cluster is on the near side of the triangle. Technically
                // unobscured, although the camera is not pointing in the right
                // direction to see it.
                if(intersection.occurs && intersection.distance < distance)
                    // Cluster is obscured by this triangle
                    return true;

            }
        }
        // Tried every triangle and found no intervening ones
        return false;
    
    }
    
    /**
     * Use line-of-sight intersection test to find the triangle lying along
     * the given line of sight direction.
     * @return
     */
    private Triangle2d lookupTriangle(Vector3d line_of_sight, List<Triangle2d> triangles)
    {

        // We will intersect two or more triangles if the LOS vector pierces
        // the database. Need to record the minimum distance in order to find
        // the triangle on the nearside.
        double min_distance = 2*state.range;
        Triangle2d target = null;
        
        // Manually synchronize on triangles list before iterating
        synchronized(triangles)
        {
            // List of Triangles forming database...
            for(Triangle2d triangle : triangles)
            {
                    // Only consider external triangles
                    if(!triangle.isExternal()) continue;
                    
                    // Check if the LOS intersects the given triangle
                    IntersectionTest intersection = new IntersectionTest(line_of_sight, triangle);

                    if(intersection.occurs && intersection.distance < min_distance)
                    {
                        min_distance = intersection.distance;
                        target = triangle;
                    }
            }
        }
        return target;
    }
    
    private void lookupCluster(int x, int y)
    {
        // Covariance-weighted distance to nearest cluster to click point
        double min_distance2 = 2 * state.confidence * state.confidence;
        Cluster target = null;
        
        // Distance threshold (squared) for a positive match. This ensures
        // that we only select a cluster when we click inside it's ellipse.
        double dist_threshold2 = state.confidence * state.confidence;
        
        // Manually synchronize on clusters list before iterating
        synchronized(clusters)
        {
            // List of clusters
            for(Cluster cluster : clusters)
            {
                
                // Cluster centre in image
                float[] pos_im = camera.projectVector(cluster.X_cam);
                
                // Vector between cluster centre and click point
                Matrix r = new Matrix(new double[][]{{pos_im[0]-x},{pos_im[1]-y}});
                
                // Cluster position covariance in image
                Matrix cov_im = camera.projectMatrix(cluster.X_cam, cluster.cov_X_cam);
                
                // Inverse covariance matrix:
                double inv_det = 1.0/cov_im.det();
                Matrix inv_cov = new Matrix(new double[][]{{cov_im.get(1, 1)*inv_det, -cov_im.get(0, 1)*inv_det},
                                                           {-cov_im.get(0, 1)*inv_det, cov_im.get(0, 0)*inv_det}});
                
                // Get covariance-weighted distance between this cluster and
                // the click point.
                double dist2 = (r.transpose().times(inv_cov.times(r))).get(0, 0);
                
                if(dist2 < dist_threshold2 && dist2 < min_distance2)
                {
                    min_distance2 = dist2;
                    target = cluster;
                }
            }
        }
        
        if(target != null)
        {
            // Clicking an already-selected cluster will remove it from the list
            if(selected_clusters.contains(target))
                selected_clusters.remove(target);
            else
                selected_clusters.add(target);
        }
        
    }
    
    /**
     * Sets the current camera and asteroid position & orientation on the PANGU
     * server to match the settings used to visualise the point cloud and
     * database. This ensures that the image obtained from PANGU matches the
     * current database view. This is called before requesting any image or
     * other measurements from the PANGU server.
     */
    private void setPANGUConfig()
    {
        
        // Extract rotation and translation from viewMatrix.
        Matrix viewRot = camVMatrix.getMatrix(0, 2, 0, 2);
        // Convert to Quaternion
        Quaternion att =  (new Quaternion(viewRot));
        // Extract translation
        Vector3d pos = new Vector3d(camVMatrix.get(0,3),
                                    camVMatrix.get(1,3),
                                    camVMatrix.get(2,3));
        // PANGU requires conjugate pose
        Pose conj = (new Pose(pos,att)).conjugate();
        
         // Get orientation of asteroid.
        Matrix modelRot = mMatrix.transpose().getMatrix(0, 2, 0, 2);
        Quaternion mAtt =  (new Quaternion(modelRot));
        // Translation of model to account for offset origin
        double x = mMatrix.get(3, 0);
        double y = mMatrix.get(3, 1);
        double z = mMatrix.get(3, 2);
        
        try
        {
            // Set camera position
            pangu.setViewpointByQuaternion(conj.position.x,
                                           conj.position.y,
                                           conj.position.z,
                                           conj.attitude.re,
                                           conj.attitude.getQ1(),
                                           conj.attitude.getQ2(),
                                           conj.attitude.getQ3());
            
            // Set position & rotation of asteroid
            pangu.setObjectPositionAttitude(0, x, y, z, mAtt.getRe(), 
                        mAtt.getQ1(), mAtt.getQ2(), mAtt.getQ3());
            
            
            // Set sun position
            pangu.setSunByDegrees(1e15,30.0,10.0);
            
        }
        catch(IOException ioe)
        {
            System.err.println("IOException in setPANGUConfig: "+ioe.getMessage());
        }      
        
    }
    
    /**
     * Configure PANGU for the current view, and get an image.
     */
    public final void updateAsteroidImage()
    {
        try
        {
            setPANGUConfig();
            raw = pangu.getImage();
        }
        catch(IOException ioe)
        {
            System.err.println("IOException in getPanguImage: "+ioe.getMessage());
        }
    }
    
    // Mouse support.
    
    int prevMouseX = 0;
    int prevMouseY = 0;
    boolean mouseLButtonDown = false;
    boolean mouseRButtonDown = false;
    
    @Override
    public void mouseClicked(MouseEvent e)
    {
        // Right mouse clicks select either triangles from the triangulation
        // solution or landmark clusters, depending on whether a triangulation
        // solution exists currently
        if ((e.getModifiers() & MouseEvent.BUTTON3_MASK) != 0)
        {
            
            if(state.database.TRIANGULATION_EXISTS)
            {
                Vector3d line_of_sight = camera.getUnitVector(e.getX(), e.getY());

                selected_triangle = lookupTriangle(line_of_sight, triangles);

                // Should always be linked to only one tetrahedron
                if(selected_triangle != null)
                    if(selected_triangle.triangles.size()!=1)
                        throw new RuntimeException("Selected triangle is connected"
                                + " to "+selected_triangle.triangles.size()+" tetrahedra!");
            }
            else
            {
                // Find cluster closest to click point (based on covariance-
                // weighted distance) and add to internal list of selected
                // clusters.
                lookupCluster(e.getX(), e.getY());
                
            }
            
        }
        else if((e.getModifiers() & MouseEvent.BUTTON1_MASK) != 0)
        {
            // Left mouse button clicked
            if(e.getClickCount() == 2)
            {
                // Double click - get intersection of LOS with PANGU model
                try
                {
                    // Update camera and model on PANGU server
                    setPANGUConfig();
                    
                    double i = (double)e.getX()/(double)camera.width;
                    double j = 1.0 - (double)e.getY()/(double)camera.height;
                    
                    ValidPoint point = pangu.lookupPoint(new Vector2D(i,j));
                    
                    if(point.valid)
                    {
                        // Get components of surface point
                        double x = point.point.i;
                        double y = point.point.j;
                        double z = point.point.k;
                        
                        // Now fix camera position & orientation so we orbit
                        // about this point rather than the origin.
                        mMatrix.set(3, 0, mMatrix.get(3, 0) - x);
                        mMatrix.set(3, 1, mMatrix.get(3, 1) - y);
                        mMatrix.set(3, 2, mMatrix.get(3, 2) - z);
                    }
                    else
                    {
                        // Reset position of model
                        mMatrix.set(3, 0, 0);
                        mMatrix.set(3, 1, 0);
                        mMatrix.set(3, 2, 0);
                    }
                    
                    // If we are in show_asteroid mode, also update asteroid image
                    if(show_asteroid) updateAsteroidImage();
                    
                }
                catch(IOException ioe)
                {
                    System.err.println("IOException thrown: "+ioe.getMessage());
                }
            }
        }
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
            Vector3d mouse = new Vector3d((float)(x - prevMouseX), 0, (float)(y - prevMouseY ));
            
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

                Matrix rotation = Matrix.identity(4, 4);

                rotation.setMatrix(new int[]{0,1,2}, new int[]{0,1,2}, rot.toMatrix());

                mMatrix = mMatrix.times(rotation);
                
                // If we are in show_asteroid mode, also update asteroid image
                if(show_asteroid) updateAsteroidImage();
            }
        }
        
        // Drags with the right button zoom in/out
        if (mouseRButtonDown) 
        {
            
            // Zoom according to displacement along Y axis
            // Smaller motions = zoom value closer to 1.0, which equals smaller
            // change to range.
            float zoom = 1.0f - ((float)(y - prevMouseY))/(float)imsize;
            
            // How does zoom affect range?
            if(Math.signum(zoom) > 0.0) state.range *= zoom;
            if(Math.signum(zoom) < 0.0) state.range /= zoom;
            
            setViewMatrixAndClipPlanes(state.range);
            
            // If we are in show_asteroid mode, also update asteroid image
            if(show_asteroid) updateAsteroidImage();
        }
        
        prevMouseX = x;
        prevMouseY = y;
        
    }
    
    @Override
    public void mouseMoved(MouseEvent e) 
    {
        // No actions when mouse simply moves without button pressed down
    }
    
    @Override
    public void mouseEntered(MouseEvent e) 
    {
        // Any actions when mouse enters frame?
    }
    
    @Override
    public void mouseExited(MouseEvent e) 
    {
        // Any actions when mouse leaves frame?
    }
    
    
}
