/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package shape.vrml;

/**
 *
 * @author nickrowell
 */
// Represents points and rgb colours
class Triple<T extends Number>
{
    T x;
    T y;
    T z;

    public Triple(final T _x, final T _y, final T _z)
    {
        x = _x;
        y = _y;
        z = _z;
    }
    
}