package gl.java.util;


import java.util.LinkedList;
import java.util.List;

import com.jogamp.opengl.GL2;

import Jama.Matrix;
import numeric.geom.dim3.Triangle3d;
import numeric.geom.dim3.Vector3d;

/**
 * Utilities class for functions used in GLJPanel classes to perform various
 * GL calls.
 * 
 * @author nickrowell
 */
public class GLUtils
{
    /**
     * Creates and returns a {@link List} of {@link Triangle<Vector3d>} objects that
     * define a cube centred on the origin with a width equal to the given value
     * @param width
     * 	The width of the cube (total length of each edge)
     * @return
     * 	A {@link List} of {@link Triangle<Vector3d>}.
     */
	public static List<Triangle3d<Vector3d>> getCube(double width) {
		
		double halfWidth = width/2.0;
		
		// Define the 8 corners of the cube
		Vector3d tlf = new Vector3d(-halfWidth,  halfWidth,  halfWidth);
		Vector3d trf = new Vector3d( halfWidth,  halfWidth,  halfWidth);
		Vector3d tlb = new Vector3d(-halfWidth, -halfWidth,  halfWidth);
		Vector3d trb = new Vector3d( halfWidth, -halfWidth,  halfWidth);
		Vector3d blf = new Vector3d(-halfWidth,  halfWidth, -halfWidth);
		Vector3d brf = new Vector3d( halfWidth,  halfWidth, -halfWidth);
		Vector3d blb = new Vector3d(-halfWidth, -halfWidth, -halfWidth);
		Vector3d brb = new Vector3d( halfWidth, -halfWidth, -halfWidth);
		
		// Build triangles
		List<Triangle3d<Vector3d>> tris = new LinkedList<>();
		
		// Top face
		tris.add(new Triangle3d<Vector3d>(tlf, trf, tlb));
		tris.add(new Triangle3d<Vector3d>(tlb, trf, trb));
		// Bottom face
		tris.add(new Triangle3d<Vector3d>(blf, brf, blb));
		tris.add(new Triangle3d<Vector3d>(blb, brf, brb));
		// Left face
		tris.add(new Triangle3d<Vector3d>(tlf, tlb, blb));
		tris.add(new Triangle3d<Vector3d>(blb, blf, tlf));
		// Right face
		tris.add(new Triangle3d<Vector3d>(trf, trb, brb));
		tris.add(new Triangle3d<Vector3d>(brb, brf, trf));
		// Front face
		tris.add(new Triangle3d<Vector3d>(tlf, trf, brf));
		tris.add(new Triangle3d<Vector3d>(brf, blf, tlf));
		// Back face face
		tris.add(new Triangle3d<Vector3d>(tlb, trb, brb));
		tris.add(new Triangle3d<Vector3d>(brb, blb, tlb));
		
		return tris;
	}
	
	/**
	 * Checks the current GL error status and throws a {@link RuntimeException} if an
	 * error is reported.
	 * @param op
	 * 	Prefix for the error message, if any.
	 * @param gl
	 * 	The {@link GL2} instance
	 */
    public static void checkGlError(String op, GL2 gl) {
        int error;
        while ((error = gl.glGetError()) != GL2.GL_NO_ERROR) {
            throw new RuntimeException(op + ": glError " + error);
        }
    }
    
    /**
     * 
     * @param shaderType
     * 	Shader type, e.g. GL.GL_VERTEX_SHADER, GL.GL_FRAGMENT_SHADER
     * @param source
     * 	Shader source as an array of Strings
     * @param gl
     * 	The {@link GL2} instance
     * @return
     * 	Non-zero integer by which the shader function can be referenced
     * 	
     */
    public static int loadShader(int shaderType, String[] source, GL2 gl) {
        
        // Make integer array of String lengths
        int[] length = new int[source.length];
        int index=0;
        for(String str: source) length[index++] = str.length();
        
        int shader = gl.glCreateShader(shaderType);
        
        if (shader != 0) {
            
            gl.glShaderSource(shader, source.length, source, length, 0);
            
            gl.glCompileShader(shader);
            
            int[] compiled = new int[1];
            
            gl.glGetShaderiv(shader, GL2.GL_COMPILE_STATUS, compiled, 0);
            
            if (compiled[0] == 0) {
                
                System.err.println("Could not compile shader "+shaderType);
                gl.glDeleteShader(shader);
                shader = 0;
            }
        }
        return shader;
    }  
    
    /**
     * 
     * @param vertexSource
     * 	Vertex shader source contained in an array of strings
     * @param fragmentSource
     * 	Fragment shader source contained in an array of strings
     * @param gl
     * 	The {@link GL2} instance
     * @return
     * 	Non-zero value by which the GL program can be referenced
     */
    public static int createProgram(String[] vertexSource, String[] fragmentSource, GL2 gl) {
        
        int vertexShader = loadShader(GL2.GL_VERTEX_SHADER, vertexSource, gl);
        if (vertexShader == 0) {
            System.err.println("loadShader(vertex) failed!");
            return 0;
        }

        int pixelShader = loadShader(GL2.GL_FRAGMENT_SHADER, fragmentSource, gl);
        if (pixelShader == 0) { 
            System.err.println("loadShader(fragment) failed!");          
            return 0;
        }      

        int program = gl.glCreateProgram();
        
        if (program != 0) {
            
            gl.glAttachShader(program, vertexShader);
            
            gl.glAttachShader(program, pixelShader);
            
            gl.glLinkProgram(program);
            
            int[] linkStatus = new int[1];
            
            gl.glGetProgramiv(program, GL2.GL_LINK_STATUS, linkStatus, 0);
            
            if (linkStatus[0] != GL2.GL_TRUE) {
                gl.glDeleteProgram(program);
                program = 0;
            }
        }
        return program;
    }
    
    /**
     * Inverse of {@link GLUtils#toFloatArray(Matrix)}.
     * @param in
     * 	The row-packed array of 16 floats
     * @return
     * 	A 4x4 {@link Matrix} in row-major order.
     */
    public static Matrix toJama4(float[] in){
    
        return new Matrix(new double[][]{{in[0],in[1],in[2],in[3]},
                                         {in[4],in[5],in[6],in[7]},
                                         {in[8],in[9],in[10],in[11]},
                                         {in[12],in[13],in[14],in[15]}});
    
    }
    
    /**
     * Converts the 4x4 Jama {@link Matrix} into a 1D row-packed array of floats.
     * @param in
     * 	The 4x4 {@link Matrix}
     * @return
     * 	The elements of the {@link Matrix} packed into a 1D row-packed array of floats.
     */
    public static float[] toFloatArray(Matrix in){
    
        double[] array = in.getRowPackedCopy();
        
        return new float[]{(float)array[0],(float)array[1],(float)array[2],(float)array[3],
                           (float)array[4],(float)array[5],(float)array[6],(float)array[7],
                           (float)array[8],(float)array[9],(float)array[10],(float)array[11],
                           (float)array[12],(float)array[13],(float)array[14],(float)array[15]};
    
    }
}