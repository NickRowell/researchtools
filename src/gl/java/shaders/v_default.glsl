precision mediump float;

// Vertex position & normal in model frame
attribute vec3 vertexPosition;
attribute vec3 vertexNormal;

// Calculate lighting and send vertex colour out
varying vec4 color_v2f;

// view and projection matrices are used to transform geometry to
// (Euclidean 3D) eye coordinates, then to clip space.
uniform mat4 mvMatrix;
uniform mat4 pMatrix;

// Sun position in model frame
uniform vec3 sun;

// Inverse transpose of the modelview matrix
uniform mat4 nMatrix;

void main(void)
{
	// Transform vertex position to eye coordinates
    vec4 V = mvMatrix * vec4(vertexPosition, 1.0);

    // Transform vertex position to clip space
    gl_Position = pMatrix * V;

    // Transform the normal to eye coordinates
    vec3 N = (nMatrix * vec4(vertexNormal,0)).xyz;

	// Vertex-Sun direction in model frame
    vec3 L = sun - vertexPosition;

	// Cosine of angle between vertex normal and light direction.
	// This determines the reflectance under the Lambert law.
    float cosA = max(dot(normalize(vertexNormal), normalize(L)),0.0);

    // Vertex shading based on Lambert reflectance law.
	color_v2f = vec4(cosA, cosA, cosA, 1.0);
}