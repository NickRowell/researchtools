package numeric.complex;

/**
 * Class to represent complex numbers.
 */
public class Complex {


    private double re;
    private double im;

    /** Basic constructor */
    public Complex(double r, double i){
        setRe(r);
        setIm(i);
    }

    /** Copy constructor */
    public Complex(Complex copyme){
        setRe(copyme.re);
        setIm(copyme.im);
    }

    final void setRe(double r){ re = r;}
    final void setIm(double i){ im = i;}

    /** Get real part of Complex number */
    public double getRe(){ return re;}
    
    /** Get imaginary part of Complex number */
    public double getIm(){ return im;}

    /** Basic arithmetic operations */

    public Complex add(Complex B){
        double r = this.re + B.re;
        double i = this.im + B.im;
        return new Complex(r,i);
    }

    public void addEquals(Complex B){
        this.re += B.re;
        this.im += B.im;
    }

    public Complex minus(Complex B){
        double r = this.re - B.re;
        double i = this.im - B.im;
        return new Complex(r,i);
    }

    public void minusEquals(Complex B){
        this.re -= B.re;
        this.im -= B.im;
    }


    public Complex times(Complex B){
        double r = this.re * B.re  - this.im * B.im;
        double i = this.im * B.re  + this.re * B.im;
        return new Complex(r,i);
    }

    public void timesEquals(Complex B){
        double r = this.re * B.re  - this.im * B.im;
        double i = this.im * B.re  + this.re * B.im;

        this.setRe(r);
        this.setIm(i);
    }

    public Complex times(double B){
        double r = this.re * B;
        double i = this.im * B;
        return new Complex(r,i);
    }    
    
    public void timesEquals(double B){
        double r = this.re * B;
        double i = this.im * B;
        
        this.setRe(r);
        this.setIm(i);        
    }     
    
    
    
    public Complex divide(Complex B){
        double r = (this.re * B.re  + this.im * B.im)/(B.re*B.re + B.im*B.im);
        double i = (this.im * B.re  - this.re * B.im)/(B.re*B.re + B.im*B.im);
        return new Complex(r,i);
    }

    public void divideEquals(Complex B){
        double r = (this.re * B.re  + this.im * B.im)/(B.re*B.re + B.im*B.im);
        double i = (this.im * B.re  - this.re * B.im)/(B.re*B.re + B.im*B.im);
        this.setRe(r);
        this.setIm(i);
    }

    /** Other operations */

    public Complex conjugate(){
        return new Complex(re,-im);
    }

    public Complex squared(){
        double r = this.re * this.re  - this.im * this.im;
        double i = this.im * this.re  + this.re * this.im;
        return new Complex(r,i);
    }

    public Complex sqrt(){
        double gamma = Math.sqrt((re + Math.sqrt(re*re + im*im))/2.0);
        double delta = Math.signum(im)*Math.sqrt((-re + Math.sqrt(re*re + im*im))/2.0);
        return new Complex(gamma,delta);
    }

    public double abs(){ return Math.sqrt(re*re + im*im);}

    public double arg(){ return Math.atan2(im,re);}

    /**
     * String representation of Complex, in form re +/- im.
     * @return 
     */
    @Override
    public String toString()
    {
        return String.format("%f %s %fi", re, im < 0 ? "-" : "+", Math.abs(im));
    }


    /**
     * Parse a Complex number from a String. Accepted format is
     *
     * (+/-/ ) a (+/-) b (anything following number e.g. i/*i/ * i)
     *
     * with any amount of white space between tokens.
     *
     */
    public static Complex parseComplex(String comp) throws NumberFormatException{

        char[] str = comp.toCharArray();

        // Set signs of real and imaginary parts. Default to positive.
        double SGN_RE = 1.0;
        double SGN_IM = 1.0;

        // Strings to hold real and imaginary values temporarily
        String re = "";
        String im = "";

        // Index of current character, stored in the first element of a
        // one-element array, so it can be manipulated by reference in methods.
        int[] p={0};

        // Eat characters until plus, minus, decimal point or numerals are found
        eatWhiteSpace(p, str);

        // Now check what character has been found...
        if((int)str[p[0]]=='+'){
            // Real part is positive
            SGN_RE = 1.0;
            // Advance pointer to next character
            p[0]++;
            // Eat characters until another plus, minus, decimal point or 
            // numeral is found
            eatWhiteSpace(p, str);
        }

        else if((int)str[p[0]]=='-'){
            // Real part is negative
            SGN_RE = -1.0;
            // Advance pointer to next character
            p[0]++;
            // Eat characters until another plus, minus, decimal point or
            // numeral is found
            eatWhiteSpace(p, str);
        }

        // Pointer should now be placed at first character of real value. If
        // this is not a numeral or decimal point, throw an exception.
        checkNumeral(p, str);

        // Read numerals and decimal points until no more are found, and store
        // the results in the String re for later parsing as a double.
        re = readNumeral(p, str);

        // Now eat characters until plus or minus is found. This is the
        // sign of the imaginary part.
        eatWhiteSpace(p, str);

        // Now check what character has been found...
        if((int)str[p[0]]=='+'){
            // Imaginary part is positive
            SGN_IM = 1.0;
            // Advance pointer to next character
            p[0]++;
            // Eat characters until another plus, minus, decimal point or 
            // numeral is found
            eatWhiteSpace(p, str);
        }

        else if((int)str[p[0]]=='-'){
            // Imaginary part is negative
            SGN_IM = -1.0;
            // Advance pointer to next character
            p[0]++;
            // Eat characters until another plus, minus, decimal point or
            // numeral is found
            eatWhiteSpace(p, str);
        }
        // Something other than a +/- is an error. Always require the sign
        // of the imaginary part at this point.
        else
            throw new NumberFormatException("Cannot parse Complex! Sign of" +
                                            " imaginary part missing!");
        
        // Pointer should now be placed at first character of imaginary value.
        // If this is not a numeral or decimal point, throw an exception.
        checkNumeral(p, str);

        // Read numerals and decimal points until no more are found, and store
        // the results in the String im for later parsing as a double.
        im = readNumeral(p, str);

        // Now parse real and imaginary values from Strings to doubles
        double real=0, imag=0;

        try{ real = Double.parseDouble(re);}
        catch(NumberFormatException nfe){
            // Catch the NumberFormatException then throw another one with
            // a more helpful message attached.
            throw new NumberFormatException("Cannot parse real value "+re);
        }
        try{ imag = Double.parseDouble(im);}
        catch(NumberFormatException nfe){
            // Catch the NumberFormatException then throw another one with
            // a more helpful message attached.
            throw new NumberFormatException("Cannot parse imaginary value "
                                            +im);
        }

        // Set signs of real and imaginary parts
        real *= SGN_RE;
        imag *= SGN_IM;

        // Now convert to a Complex and return
        return new Complex(real, imag);

    }

    private static void eatWhiteSpace(int[] p, char[] str){
         while(!((int)str[p[0]]=='+' || (int)str[p[0]]=='-' || (int)str[p[0]]=='.' ||
                ((int)str[p[0]] > 47 && (int)str[p[0]] < 58)))
            p[0]++;
    }

    private static void checkNumeral(int[] p, char[] str){
        if(!((int)str[p[0]]=='.' || ((int)str[p[0]]>47 && (int)str[p[0]]<58)))
            throw new NumberFormatException("Cannot parse Complex! Found "+
                                            "non-numeral "+str[p[0]]);
    }

    private static String readNumeral(int[] p, char[] str){

        String NUM = "";

        while((int)str[p[0]]=='.' || ((int)str[p[0]]>47 && (int)str[p[0]]<58))
        {
            NUM = NUM.concat(""+str[p[0]]);
            p[0]++;
        }

        return NUM;
    }

    /**
     * Format 2D array of doubles to Complex numbers with imaginary parts
     * set to zero.
     */
    public static Complex[][] toComplexArray(double[][] real){
        
        Complex[][] complex = new Complex[real.length][];

        for(int i=0; i< real.length; i++)
            complex[i] = toComplexArray(real[i]);
            
        return complex;
    }
    
    /**
     * Format 1D array of doubles to Complex numbers with imaginary parts
     * set to zero.
     */
    public static Complex[] toComplexArray(double[] real){
        
        Complex[] complex = new Complex[real.length];

        for(int i=0; i< real.length; i++)
            complex[i] = new Complex(real[i], 0.0);
        
        return complex;
    }   
    
    
    /**
     * Format 2D array of ints to Complex numbers with imaginary parts
     * set to zero.
     */ 
    public static Complex[][] toComplexArray(int[][] real){
        
        Complex[][] complex = new Complex[real.length][];

        for(int i=0; i< real.length; i++)
            complex[i] = toComplexArray(real[i]);
        
        return complex;
    }
    
    /**
     * Format 1D array of ints to Complex numbers with imaginary parts
     * set to zero.
     */     
    public static Complex[] toComplexArray(int[] real){
        
        Complex[] complex = new Complex[real.length];

        for(int i=0; i< real.length; i++)
            complex[i] = new Complex(real[i], 0.0);
            
        return complex;
    }
    
    
    
    
}
