package misc;

import java.io.*;

import misc.Complex;

/**
 * Test functions of the Complex class.
 * @author nickrowell
 */
public class ComplexTester {


    public static void main(String[] args) throws IOException{

        BufferedReader in = new BufferedReader(
                new InputStreamReader(System.in));


        while(true){

            System.out.print("\nEnter complex number: ");
            String comp_str = in.readLine();

            Complex comp = Complex.parseComplex(comp_str);

            System.out.println("Original string = "+comp_str+
                               "\nParsed complex = "+comp.toString());
        }

    }

}
