package numeric.functions;

/**
 * Represents line segments in e.g. Linear interpolation objects.
 * @author nickrowell
 */
public class Line {
    
    // Standard parameters of line y = m*x + c
    public double m,c;

    // Range of this segment in the X axis
    public double min_x,max_x;
    
    // Range of this segment in the Y axis
    public double min_y,max_y;    
    
    /** Construct a Line that passes through (x1,y1) and (x2,y2). */
    public Line(double x1, double y1, double x2, double y2){
        
        // Sanity checking
        if(x2 <= x1) {
            throw new RuntimeException("x points ordered incorrectly or equal!");
        }
        
        m = (y2-y1)/(x2-x1);
        c = y2-((y2-y1)/(x2-x1))*x2;
        min_x = x1;
        max_x = x2;
        
        // Get limits on Y - gradient can be positive or negative so need
        // to check order of points before assigning min/max.
        min_y = Math.min(getY(min_x),getY(max_x));
        max_y = Math.max(getY(min_x),getY(max_x));
        
    }
    
    /** Alternative constructor based on line parameters and x range. */
    public Line(double M, double C, double[] x12){
        m = M;
        c = C;
        
        if(x12.length!=2)
            throw new RuntimeException("Incorrect number of x limits: "+x12.length);
        
        min_x = x12[0];
        max_x = x12[1];
        
        // Get limits on Y - gradient can be positive or negative so need
        // to check order of points before assigning min/max.
        min_y = Math.min(getY(min_x),getY(max_x));
        max_y = Math.max(getY(min_x),getY(max_x));               
    }
    
    /** Interpolate Y at given X */
    public final double getY(double X){ return m*X + c;}
    
    /** Interpolate X at given Y */
    public final double getX(double Y){ return (Y - c)/m;}
    
    /** Get gradient of line */
    public double getGradient(){ return m;}
    
    /** 
     * Check if this y value lies within range spanned by this line segment,
     * inclusive of end points.
     */
    public boolean containsY(double y){
        return y>=min_y && y<=max_y;
    }
    /** 
     * Check if this y value lies within range spanned by this line segment,
     * excluding points coincident with lower limit of range. This is useful
     * for linear interpolation where multiple Lines form a piecewise function.
     * This avoids multiple solutions when data lie exactly on the point
     * where two lines meet.
     */
    public boolean containsYRestricted(double y){
        return y>min_y && y<=max_y;
    }
    /** 
     * Check if this x value lies within range spanned by this line segment,
     * inclusive of end points.
     */
    public boolean containsX(double x){
        return x>=min_x && x<=max_x;
    }
    /** 
     * Check if this x value lies within range spanned by this line segment,
     * excluding points coincident with lower limit of range. This is useful
     * for linear interpolation where multiple Lines form a piecewise function.
     * This avoids multiple solutions when data lie exactly on the point
     * where two lines meet.
     */
    public boolean containsXRestricted(double x){
        return x>min_x && x<=max_x;
    }
    
    /**
     * Definite integral of this line segment wrt x between min_x and x.
     */
    public double integrateWrtX(double x){
        
        // if X is below range of this line segment, integral is zero
        if(x <= min_x) return 0.0;
        
        // Otherwise, set integration limits
        double a = min_x;
        double b = Math.min(max_x, x);  // Restrict upper limit to max_x
        
        return m * (b*b - a*a)/2.0 + c * (b - a);
        
    }
    
    /** Integral over total range. */
    public double integralWrtX(){
        return integrateWrtX(max_x);
    }
    
    
    @Override
    public String toString(){
        return "f(x) = "+m+" * x + "+c+", [ "+min_x+" : "+max_x+" ]";
    }

    
}
