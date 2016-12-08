package numeric.triangulation.dim3.infra;

import numeric.geom.dim3.Sphere;
import numeric.geom.dim3.Vector3d;

/**
 * Class represents Tetrahedra, i.e. the basic unit in the 3D Delaunay triangulation.
 *
 * @author nrowell
 * @version $Id$
 */
public class Tetrahedron {
    
    /**
     * Array of (four) Triangles that define this Tetrahedron.
     */
    public Triangle3d[] tris = new Triangle3d[4];
    
    /**
     * Array of (four) Vertices that define this Tetrahedron.
     */
    public Vertex3d[] verts = new Vertex3d[4];
    
    /**
     * Circumscribing sphere for Tetrahedron.
     */
    public Sphere sphere;
    
    /**
     * Does this Tetrahedron connect to the external hull?
     */
    public boolean EXTERNAL;
    
    /**
     * Main constructor. Vertex v1 must lie opposite Triangle t1, and so on.
     * 
     * Also the winding order is important, though this is calculated by the
     * constructor. The vertices are arranged like this:
     * 
     * 
     * 
     *                          V3
     *                         /|\
     *                        / | \
     *                       /  |  \
     *                      /   |   \
     *                     /    |    \
     *                    /     |     \
     *                   /      |      \
     *                  /       |       \
     *                 /        |        \
     *                V1-----------------V2
     *                 \        |        /
     *                  \       |       /
     *                   \      |      /
     *                    \     |     /
     *                     \    |    /
     *                      \   |   /
     *                       \  |  /
     *                        \ | /
     *                         \|/
     *                          V0
     * 
     * @param t0    First Triangle.
     * @param t1    Second Triangle.
     * @param t2    Third Triangle.
     * @param t3    Fourth Triangle.
     * @param v0    Vertex lying opposite first triangle.
     * @param v1    Vertex lying opposite second triangle.
     * @param v2    Vertex lying opposite third triangle.
     * @param v3    Vertex lying opposite fourth triangle.
     */
    public Tetrahedron(Triangle3d t0, Triangle3d t1, Triangle3d t2, Triangle3d t3,
                       Vertex3d v0, Vertex3d v1, Vertex3d v2, Vertex3d v3, Sphere asphere) {
        tris[0] = t0;
        tris[1] = t1;
        tris[2] = t2;
        tris[3] = t3;
        
        verts[0] = v0;
        verts[1] = v1;
        verts[2] = v2;
        verts[3] = v3;
        
        sphere = asphere;
        
        // OK: now check the current winding order for vertices, and reverse
        // if necessary so that getOutwardNormalForTriangle(Triangle t) always
        // returns the outward pointing normal.
        
        // Check that surface normal for triangle 3 (opposite vertex 3) points
        // in the opposite direction to the vector v3-v2. If not, we need to 
        // reverse the order of vertices v1 & v2.
        if(getOutwardNormalForTriangle(t3).dot(v3.minus(v2)) > 0)
        {
            verts[1] = v2;
            verts[2] = v1;
            tris[1]  = t2;
            tris[2]  = t1;
        }
        
    }
    
    /**
     * Construct a Tetrahedron from a single Triangle that forms the base and a
     * single Vertex that forms the peak.
     */
    public Tetrahedron(Triangle3d tri, Vertex3d v)
    {
        // Check that this cluster isn't part of the triangle.
        if(tri.v0 == v || tri.v1 == v || tri.v2 == v)
            throw new RuntimeException("Vertex forms part of the Triangle!");
        
        Triangle3d t01v = new Triangle3d(tri.v0, tri.v1, v);
        Triangle3d t02v = new Triangle3d(tri.v0, tri.v2, v);
        Triangle3d t12v = new Triangle3d(tri.v1, tri.v2, v);
        
        tris[0] = tri;
        tris[1] = t01v;
        tris[2] = t02v;
        tris[3] = t12v;
        
        verts[0] = v;
        verts[1] = tri.v2;
        verts[2] = tri.v1;
        verts[3] = tri.v0;
        
        sphere = new Sphere(verts[0], verts[1], verts[2], verts[3]);
        
        if(getOutwardNormalForTriangle(tris[3]).dot(verts[3].minus(verts[2])) > 0)
        {
            verts[1] = tri.v1;
            verts[2] = tri.v2;
            tris[1]  = t02v;
            tris[2]  = t01v;
        }
        
    }
    
    /**
     * Does the Vector3d lie inside the circumscribing sphere for this
     * tetrahedron?
     */
    public boolean circumSphereContains(Vector3d vertex) {
        return sphere.centre.minus(vertex).norm2() < sphere.r2;
    }
    
    /**
     * Does the given point lie inside this tetrahedron?
     */
    public boolean contains(Vector3d v)
    {
        // Can immediately rule out points that lie outside of the
        // circumscribing sphere. This is an easy test.
        if(!circumSphereContains(v)) {
        	return false;
        }
        
        // Point lies within circumscribing sphere. Now apply more complex test
        // to see if the point lies within the tetrahedron itself.
        
        // Un-normalised outward facing normals for each triangle
        Vector3d N0 = (verts[2].minus(verts[1])).cross(verts[3].minus(verts[1]));
        Vector3d N1 = (verts[0].minus(verts[2])).cross(verts[3].minus(verts[2]));
        Vector3d N2 = (verts[0].minus(verts[3])).cross(verts[1].minus(verts[3]));
        Vector3d N3 = (verts[2].minus(verts[0])).cross(verts[1].minus(verts[0]));
        
        // Vectors between test point and each vertex in tetrahedron
        Vector3d vmv0 = v.minus(verts[0]);
        Vector3d vmv1 = v.minus(verts[1]);
        Vector3d vmv2 = v.minus(verts[2]);
        Vector3d vmv3 = v.minus(verts[3]);
        
        // Condition for point lying inside tetrahedron is that the dot products
        // of the normal vectors with these vectors are all negative.
        if(vmv0.dot(N1) > 0.0 || vmv1.dot(N2) > 0.0 || vmv2.dot(N3) > 0.0 || vmv3.dot(N0) > 0.0) {
        	return false;
        }
        
        return true;
    }
    
    /**
     * Gets the {@link Triangle} that lies opposite the given {@link Vertex3d} (which
     * must itself be part of this {@link Tetrahedron}), i.e. finds the single {@link Triangle}
     * that is part of this {@link Tetrahedron} and which doesn't have the given {@link Vertex3d}
     * as one it's vertices.
     * @param v
     * 	The {@link Vertex3d} (which must be part of this {@link Tetrahedron}).
     * @return
     * 	The {@link Triangle} that lies opposite the given {@link Vertex3d}.
     * @throws RuntimeException
     * 	If the {@link Vertex3d} does not form part of this {@link Tetrahedron}.
     */
    public Triangle3d getTriangleOppositeVertex(Vertex3d v) {
        if(v == verts[0]) {
        	return tris[0];
        }
        else if (v == verts[1]) {
        	return tris[1];
        }
        else if (v == verts[2]) {
        	return tris[2];
        }
        else if (v == verts[3]) {
        	return tris[3];
        }
        else {
            throw new RuntimeException("The given Vertex3d does not form part of this Tetrahedron!");
        }
    }
    
    /**
     * Get the outward-pointing surface normal for this {@link Triangle} (which
     * must itself be part of this {@link Tetrahedron}).
     * @param t
     * 	The {@link Triangle}.
     * @return
     * 	A {@link Vector3d} representing the outward-pointing surface normal for the
     * {@link Triangle}.
     * @throws RuntimeException
     * 	If the {@link Triangle} does not form part of this {@link Tetrahedron}.
     */
    public final Vector3d getOutwardNormalForTriangle(Triangle3d t) {
        if(t==tris[0]) {
            Vector3d N = (verts[2].minus(verts[1])).cross(verts[3].minus(verts[1]));
            return N.normalise();
        }
        else if(t==tris[1]) {
            Vector3d N = (verts[0].minus(verts[2])).cross(verts[3].minus(verts[2]));
            return N.normalise();
        }
        else if(t==tris[2]) {
            Vector3d N = (verts[0].minus(verts[3])).cross(verts[1].minus(verts[3]));
            return N.normalise();
        }
        else if(t==tris[3]) {
            Vector3d N = (verts[2].minus(verts[0])).cross(verts[1].minus(verts[0]));
            return N.normalise();
        }
        else {
            throw new RuntimeException("Triangle not found among this Tetrahedron's set!");
        }
    }
    
}