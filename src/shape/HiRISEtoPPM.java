package shape;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Scanner;

/**
 * Program written to parse a DTM produced by the HiRISE team and convert it to
 * a PANGU elevation model in ppm format that can be converted to pan using
 * dem2pan.
 * 
 * We extract only a (rectangular) sub-region of the terrain model, because the
 * full terrain models are typically huge and cannot be handled in full by PANGU
 * facilities. The sub-region is centred on a latitude and longitude specified
 * by the user.
 * 
 * @author nickrowell
 */
public class HiRISEtoPPM
{
    
    /** Handle to input IMG file. */
    File img_file;
    /** Final output PXN file. */
    File ppm_file;
    
    // Stats used to digitize height values:
    
    /** Maximum valid height. */
    float hmax;
    /** Minimum valid height. */
    float hmin;
    /** Height difference between consecutive levels in tga16. */
    float hres;
    
    /**
     * Values from IMG image header that are required to perform the conversion.
     */
    
    static enum PDS_SAMPLE_TYPE
    {
        PC_REAL(ByteOrder.LITTLE_ENDIAN),
        UNSIGNED_INTEGER(ByteOrder.LITTLE_ENDIAN),
        UNKNOWN(null);
        
        ByteOrder ENDIANNESS;
        
        PDS_SAMPLE_TYPE(ByteOrder pENDIANNESS){ ENDIANNESS = pENDIANNESS;}
    }
    
    PDS_SAMPLE_TYPE SAMPLE_TYPE = PDS_SAMPLE_TYPE.UNKNOWN;
    
    static PDS_SAMPLE_TYPE parseSampleType(String sample)
    {
        if(sample.equalsIgnoreCase("PC_REAL")) return PDS_SAMPLE_TYPE.PC_REAL;
        if(sample.equalsIgnoreCase("UNSIGNED_INTEGER")) return PDS_SAMPLE_TYPE.UNSIGNED_INTEGER;
        return PDS_SAMPLE_TYPE.UNKNOWN;
    }
    
    static enum PDS_PROJECTION_TYPE
    {
        EQUIRECTANGULAR,
        SINUSOIDAL,
        UNKNOWN;
    }
    
    PDS_PROJECTION_TYPE MAP_PROJECTION_TYPE = PDS_PROJECTION_TYPE.UNKNOWN;
    
    // Projection object that encapsulates projection algebra
    Projection projection;
    
    static PDS_PROJECTION_TYPE parseProjectionType(String proj)
    {
        if(proj.equalsIgnoreCase("\"EQUIRECTANGULAR\"")) return PDS_PROJECTION_TYPE.EQUIRECTANGULAR;
        if(proj.equalsIgnoreCase("EQUIRECTANGULAR")) return PDS_PROJECTION_TYPE.EQUIRECTANGULAR;
        if(proj.equalsIgnoreCase("\"SINUSOIDAL\"")) return PDS_PROJECTION_TYPE.SINUSOIDAL;
        if(proj.equalsIgnoreCase("SINUSOIDAL")) return PDS_PROJECTION_TYPE.SINUSOIDAL;
        return PDS_PROJECTION_TYPE.UNKNOWN;
    }
    
    // Enumerated types to represent integer valued parameters.
    enum IntPar
    {
        RECORD_BYTES(0),          // Size of IMG header [bytes]
        LINES(0),                 // Size of IMG data [lines*samples]
        LINE_SAMPLES(0),
        SAMPLE_BITS(0),           // Number of bits in each sample
        MISSING_DATA(0xFF7FFFFB), // Default value for missing data
        BANDS(1);                 // Number of channels per sample
        
        // Parameter value
        int value;
        // Indicates if value has been parsed from file (or is default)
        boolean set = false;
        
        IntPar(int def){ value = def;}
        
        void set(String sval, boolean verbose)
        {
            value = Integer.parseInt(sval);
            set   = true;
            if(verbose)
                System.out.println("Parsed "+this.toString()+" = "+value);
        }
        
        
    }
    
    
    // Enumerated types to represent double valued parameters.
    enum DoublePar
    {
        CENTER_LATITUDE(0),           // Centre of projection [degrees]
        CENTER_LONGITUDE(0),
        LINE_PROJECTION_OFFSET(0),    // Offset from samples origin [pixels]
        SAMPLE_PROJECTION_OFFSET(0),
        // Rotation of projection [degrees] - note we assume this is zero. The 
        // conversion will give dodgy results if not.
        MAP_PROJECTION_ROTATION(0),
        // Radius of reference ellipsoid model, usually set equal to the local
        // radius at the centre of projection to avoid confusion [km]
        A_AXIS_RADIUS(0),     // Equals equatorial radius [km]
        B_AXIS_RADIUS(0),
        C_AXIS_RADIUS(0),     // Equals polar radius [km]
        MAP_SCALE(0),         // Physical scale of map [m/pixel]
        MAP_RESOLUTION(0);    // Angular scale of pixels [pixel/deg]

        // Parameter value
        double value;
        // Indicates if value has been parsed from file (or is default)
        boolean set = false;
        
        DoublePar(int def){ value = def;}
        
        void set(String sval, boolean verbose)
        {
            value = Double.parseDouble(sval);
            set   = true;
            if(verbose)
                System.out.println("Parsed "+this.toString()+" = "+value);
        }
    }
    
    
    
    public static void main(String[] args) throws IOException
    {
        
        //File parent = new File("/home/nickrowell/Projects/PANGU_4/Documents/TN01/Planetary_Science/Available_Data/Shape_Models/Mars/HiRISE/Gale_Crater/mrohr_0001_0001/dtm/psp/orb_010500_010599/psp_010573_1755_psp_010639_1755");
        File parent = new File("/home/nickrowell/Projects/PANGU_4/Documents/TN01/Planetary_Science/Available_Data/Shape_Models/Mars/HRSC/mexhrs_2001/data/1927");
        


        // Name of input IMG file
        //File input_file = new File(parent, "dteec_010573_1755_010639_1755_u01.img");
        File input_file = new File(parent, "h1927_0000_nd4.img");
        
        
        // Name out output PXN file (will be written to same directory)
        //File ppm_file = new File(parent, "dteec_010573_1755_010639_1755_u01.ppm");
        File ppm_file = new File(parent, "h1927_0000_nd4.ppm");
        
        
        // Latitude and longitide of MSL landing site (extracted DEM will be
        // centred on this point).
        double lat = Math.toRadians(-4.5895);
        double lon = Math.toRadians(137.4417);
        
        HiRISEtoPPM imgToPXN = new HiRISEtoPPM(input_file, ppm_file, lat, lon, true);
    }
    
    
    
    public HiRISEtoPPM(File pimg_file, File pppm_file, double lat, double lon, boolean verbose) throws IOException
    {
        
        img_file  = pimg_file;
        ppm_file  = pppm_file;
        
        // Parse IMG model from file.
        try 
        {
            extractIMG(img_file, ppm_file, lat, lon, verbose);
        }
        catch (IOException ex) 
        {
            System.err.println("Error reading file "+img_file.toString()+
                               ":\n"+ex.getMessage());
            System.exit(1);
        }
        
        
    }
    
    private void extractIMG(File input, File output, double lat, double lon, boolean verbose) throws IOException
    {
        
        // Open reader on input file
        FileInputStream imgreader = new FileInputStream(input);
        
        readIMGHeader(imgreader, verbose);
        
        // Now we rationalise the parameters to decide if we can proceed to
        // interpret the data section.
        StringBuilder message = new StringBuilder();
        if(!rationaliseParameters(message))
        {
            throw new IOException(message.toString());
        }
        
        // Get indices in data array of sample closest to specfied point.
        int[] ls = projection.getLineSampleIndices(lon, lat);
        int s0 = ls[0];
        int l0 = ls[1];
        
        if(verbose)
            System.out.println("Centering extracted DTM on "+l0+", "+s0);
                
        // Set radius of region to extract (if it's odd we get a central pixel)
        int WIDTH  = 4097;
        int HEIGHT = 4097;
        
        // Temporary storage for each sample
        int nbytes = IntPar.SAMPLE_BITS.value/8;
        
        byte[] bsample = new byte[nbytes];
        
        if(verbose)
            System.out.println("First pass over the data: reading IMG and "
                    + "computing height stats...");
        
        // Record max & min heights (to allow later digitizing to 16-bit TGA)
        boolean initialised_heights = false;
        
        // Read all samples into a 2D array
        float[][] samples = new float[IntPar.LINES.value][IntPar.LINE_SAMPLES.value];
        
        // Read in the data in samples of 32 bits (4 bytes)
        for(int l=0; l<IntPar.LINES.value; l++)
        {
            for(int s=0; s<IntPar.LINE_SAMPLES.value; s++)
            {
                
                int nread = imgreader.read(bsample);      // Read one sample
                ByteBuffer bb = ByteBuffer.wrap(bsample); // Create a byte buffer and wrap the array
                
                bb.order(SAMPLE_TYPE.ENDIANNESS);         // Set endian-ness
                
                boolean isMissing = false;
                float   sample;
                
                switch(SAMPLE_TYPE)
                {
                        case PC_REAL:
                        {
                            // 32-bit samples: read 4 bytes into an integer
                            int data = bb.getInt();
                            isMissing = (data==IntPar.MISSING_DATA.value);
                            sample = Float.intBitsToFloat(data);
                            break;
                        }
                        case UNSIGNED_INTEGER:
                        {
                            // 8-bit samples: read one byte
                            byte data = bb.get();
                            isMissing = (data==IntPar.MISSING_DATA.value);
                            sample = (float)(data & 0xFF);
                            break;
                        }
                        default:
                            sample = 0;
                }
                
                
                
                if(!isMissing)                           // Skip invalid samples
                {
                    // Only use valid samples lying inside selected region to 
                    // determine height range, in order to maximise resolution
                    // of output PPM height.
                    if(Math.abs(l0-l) <= (HEIGHT-1)/2 && Math.abs(s0-s) <= (WIDTH-1)/2)
                    {
                        if(!initialised_heights)
                        {
                            hmax = sample; hmin = sample;           // Initialise values on
                            initialised_heights = true;             // first valid data
                        }
                        else
                        {
                            if(sample > hmax) hmax = sample;        // Update samples on all
                            if(sample < hmin) hmin = sample;        // remaining data
                        }
                    }
                    
                    samples[l][s] = sample;
                    
                }
                else
                {
                    samples[l][s] = Float.NaN;
                }
                
            }
        }
        
        // Scale height to restrict to range +/-32768 when digitised.
        float scale_factor = 32768/Math.max(Math.abs(hmax),Math.abs(hmin));
        
        System.out.println("hmax = "+hmax);
        System.out.println("hmin = "+hmin);
        
        // Metres-per-digital-level
        hres = 1.0f/scale_factor;
        
        
        imgreader.close();
        
        if(verbose)
            System.out.println("Second pass over data: digitizing samples and"
                    + " writing PPM file...");
        
        // Open writer for final PPM file
        FileOutputStream ppmwriter = new FileOutputStream(output);
        // Write header
        String ppmheader = "P6\n"+WIDTH+" "+HEIGHT+" 255\n";
        ppmwriter.write(ppmheader.getBytes());
        
        for(int l=l0-(HEIGHT-1)/2; l<=l0+(HEIGHT-1)/2; l++)
        {
            for(int s=s0-(WIDTH-1)/2; s<=s0+(WIDTH-1)/2; s++)
            {
                // Analogue height
                float sample;
                
                // Check if the desired sample falls in the image area
                if(l<0 || l>=IntPar.LINES.value || s<0 || s>=IntPar.LINE_SAMPLES.value)
                {
                    sample = hmin;
                }
                // Sample lies inside image area; final check for invalid data
                else
                {
                    sample = (Float.isNaN(samples[l][s])) ? hmin : samples[l][s];
                }
                
                // Digitize, rounding away from zero (see PANGU UM S.3.1.2.3)
                int tga16 = sample < 0 ?
                              (int)Math.floor(sample * scale_factor) :
                              (int)Math.ceil(sample * scale_factor);
                
                tga16 += 32768;
                
                // Write out 3 bytes per pixel (only 2 are used)
                byte[] bdata = new byte[3];
                bdata[0] = (byte) (tga16/256);
                bdata[1] = (byte) (tga16%256);
                bdata[2] = (byte) 0;
                
                ppmwriter.write(bdata);
                
            }
        }
        
        ppmwriter.close();
        
        // Print dem2pan command...
        System.out.println("/opt/pangu/source/pangu_testing/bin/linux_x86/dem2pan \\");
        System.out.println("-o dteec_010573_1755_010639_1755_u01.pan \\");
        System.out.println("-z "+hres+" -s "+DoublePar.MAP_SCALE.value+" \\");
        System.out.println("dteec_010573_1755_010639_1755_u01.ppm");

    }
    
    /**
     * Reads header section of an IMG file and returns the entire header text
     * as a string array (not yet...). On exit, the reader is positioned at the
     * start of the data section.
     * 
     * Note: could instead return header fields as a key/value map.
     * 
     * @param freader
     * @param verbose
     * @throws IOException 
     */
    void readIMGHeader(FileInputStream freader, boolean verbose) throws IOException
    {
        // Needs to have more elements than the longest ASCII line in header
        byte[] bytes = new byte[1024];
        int bytes_so_far = 0;
        
        // Flag used to indicate current comment status of text being parsed
        // (necessary to support multi-line comments)
        boolean COMMENT_ACTIVE = false;
        
        while(true)
        {
            // Read a single line from header as a byte array
            int nbytes = parseASCIILine(bytes, freader);
            
            // Must add one to account for carriage return, which is not
            // included in parsed ascii but does count towards number of bytes
            // in header.
            bytes_so_far += nbytes+1;  
            
            // Encode as a String using the default charset
            String line = new String(bytes, 0, nbytes);
            
//            if(verbose) System.out.println(line);
            
            // detect the end of the header text
            if(line.trim().equalsIgnoreCase("END")) break;
            
            // Use a Scanner to parse values from header
            Scanner scan = new Scanner(line);
            
            while(scan.hasNext())
            {
                // Non-empty line
                String token = scan.next();
                
                if(COMMENT_ACTIVE)
                {
                    // In COMMENT state, we only recognize COMMENT_END symbols
                    if(token.equalsIgnoreCase("*/"))
                    {
                        COMMENT_ACTIVE = false;
                    }
                    
                }
                else
                {
                    if(token.equalsIgnoreCase("/*"))
                    { 
                        COMMENT_ACTIVE = true;
                    }
                    if(token.equalsIgnoreCase("RECORD_BYTES"))
                    {
                        scan.next();
                        IntPar.RECORD_BYTES.set(scan.next(), verbose);
                    }
                    if(token.equalsIgnoreCase("BANDS"))
                    {
                        scan.next();
                        IntPar.BANDS.set(scan.next(), verbose);
                    }
                    if(token.equalsIgnoreCase("LINES"))
                    {
                        scan.next();
                        IntPar.LINES.set(scan.next(), verbose);
                    }
                    if(token.equalsIgnoreCase("LINE_SAMPLES"))
                    {
                        scan.next();
                        IntPar.LINE_SAMPLES.set(scan.next(), verbose);
                    }
                    if(token.equalsIgnoreCase("SAMPLE_BITS"))
                    {
                        scan.next();
                        IntPar.SAMPLE_BITS.set(scan.next(), verbose);
                        
                    }
                    
                    // Not figured out how to set this properly - for example:
                    //
                    // int HOLE = 0xFF7FFFFB;
                    //
                    // works, but
                    //
                    // int HOLE = Integer.parseInt("FF7FFFFB", 16);
                    //
                    // overflows.
                    //
                    if(token.equalsIgnoreCase("MISSING_CONSTANT"))
                    {
                        scan.next();
                        
                        // The MISSING_CONSTANT string has the format:
                        //
                        // 16#FF7FFFFB#
                        //
                        // where the first number gives the base for the number
                        // between the hash symbols, which is the actual value
                        // used to indicate missing data.
                        Scanner scan2 = new Scanner(scan.next());
                        
                        scan2.useDelimiter("#");
                        
                        String sbase = scan2.next();
                        String shole = scan2.next();
                        
                        int base = Integer.parseInt(sbase);
                        long hole = Long.parseLong(shole, base);

//                        LongPar.MISSING_DATA.set(hole, verbose);
                        
                    }
                    if(token.equalsIgnoreCase("A_AXIS_RADIUS"))
                    {
                        scan.next(); 
                        DoublePar.A_AXIS_RADIUS.set(scan.next(), verbose);
                        
                    }
                    if(token.equalsIgnoreCase("B_AXIS_RADIUS"))
                    {
                        scan.next(); 
                        DoublePar.B_AXIS_RADIUS.set(scan.next(), verbose);
                        
                    }
                    if(token.equalsIgnoreCase("C_AXIS_RADIUS"))
                    {
                        scan.next(); 
                        DoublePar.C_AXIS_RADIUS.set(scan.next(), verbose);
                        
                    }
                    if(token.equalsIgnoreCase("SAMPLE_TYPE"))
                    {
                        scan.next();
                        SAMPLE_TYPE = parseSampleType(scan.next());
                        if(verbose)
                            System.out.println("Parsed SAMPLE_TYPE = "+SAMPLE_TYPE.toString());
                    }
                    
                    if(token.equalsIgnoreCase("CENTER_LATITUDE"))
                    {
                        scan.next(); 
                        DoublePar.CENTER_LATITUDE.set(scan.next(), verbose);
                    }
                    if(token.equalsIgnoreCase("CENTER_LONGITUDE"))
                    {
                        scan.next(); 
                        DoublePar.CENTER_LONGITUDE.set(scan.next(), verbose);
                    }
                    if(token.equalsIgnoreCase("MAP_SCALE"))
                    {
                        scan.next(); 
                        DoublePar.MAP_SCALE.set(scan.next(), verbose);
                    }
                    if(token.equalsIgnoreCase("MAP_RESOLUTION"))
                    {
                        scan.next(); 
                        DoublePar.MAP_RESOLUTION.set(scan.next(), verbose);
                    }
                    if(token.equalsIgnoreCase("LINE_PROJECTION_OFFSET"))
                    {
                        scan.next(); 
                        DoublePar.LINE_PROJECTION_OFFSET.set(scan.next(), verbose);
                    }
                    if(token.equalsIgnoreCase("SAMPLE_PROJECTION_OFFSET"))
                    {
                        scan.next(); 
                        DoublePar.SAMPLE_PROJECTION_OFFSET.set(scan.next(), verbose);
                    }
                    if(token.equalsIgnoreCase("MAP_PROJECTION_TYPE"))
                    {
                        scan.next();
                        MAP_PROJECTION_TYPE = parseProjectionType(scan.next());
                        if(verbose)
                            System.out.println("Parsed MAP_PROJECTION_TYPE = "+MAP_PROJECTION_TYPE.toString());
                    }
                    if(token.equalsIgnoreCase("MAP_PROJECTION_ROTATION"))
                    {
                        scan.next(); 
                        DoublePar.MAP_PROJECTION_ROTATION.set(scan.next(), verbose);
                    }
                    
                }
                
            }
            
        }
        
        // Reached end of header text but not necessarily header section. The
        // variable bytes_so_far records the number of bytes read so far from
        // the header section. We now skip over any remaining header bytes to
        // reach the start of the data section.
        while(bytes_so_far < IntPar.RECORD_BYTES.value) { freader.read(); bytes_so_far++;}
    }
    
    
    /**
     * Reads single bytes from the input FileInputStream and stores them in the
     * given byte array, until an ascii '\n' is encountered. 
     * @param data
     * @param freader 
     */
    int parseASCIILine(byte[] data, FileInputStream freader) throws IOException
    {
        int bytes = 0;
        byte b;
        
        while((b = (byte)freader.read()) != '\n')
        {
            // Reached end of file
            if(b == -1) return bytes;
            
            data[bytes++] = b;
        }
        
        // Reached end of line
        return bytes;
    }
    
//    private void writePXN(File pxn_file, boolean verbose) throws IOException
//    {
//        
//        // Open writer on output file
//        BufferedWriter out = new BufferedWriter(new FileWriter(pxn_file));
//        
//        // Write PXN header
//        out.write("<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
//        out.write("<pangu_model ver=\"0\">\n"); 
//        
//        // Cease recursion and write out mesh for this node
//        out.write("<mesh>\n");
//        
////        new Vector3d(l, s, bb.getFloat());
//        
//        // Write out vertex positions
//        out.write("<positions>\n");
//        for(int l=0; l<LINE_LAST_PIXEL; l++)
//        {
//            for(int s=0; s<SAMPLE_LAST_PIXEL; s++)
//            {
//                // Here is where we will eventually invert the map projection.
//                Vector3d pixel = new Vector3d(l, s, samples[l][s]);
//                
//                out.write("<vec3> "+pixel.x+" "+pixel.y+" "+pixel.z+" </vec3>\n"); 
//            }
//        }
//        out.write("</positions>\n");
//        // Write out surface normals
//        out.write("<normals>\n");
//        for(int l=0; l<LINE_LAST_PIXEL; l++)
//        {
//            for(int s=0; s<SAMPLE_LAST_PIXEL; s++)
//            {
//                // Dummy normals along edges currently (can set these to normal
//                // of nearest triangular facet)
//                if(l==0 || l==(LINE_LAST_PIXEL-1))
//                {
//                    out.write("<vec3> 0 0 1 </vec3>\n");
//                    continue;
//                }
//                if(s==0 || s==(SAMPLE_LAST_PIXEL-1))
//                { 
//                    out.write("<vec3> 0 0 1 </vec3>\n");
//                    continue;
//                }
//                
//                
//                // Everywhere else: sample the surrounding vertices to estimate
//                // the normal at the present point.
//                //
//                // Use the 4 nearest neighbours. These vertices are laid out:
//                //
//                //           l-->
//                //
//                //   s        E
//                //   |    D   A   B
//                //   V        C
//                //
//                Vector3d A = new Vector3d(l, s, samples[l][s]);
//                Vector3d B = new Vector3d(l+1, s, samples[l+1][s]);
//                Vector3d C = new Vector3d(l, s+1, samples[l][s+1]);
//                Vector3d D = new Vector3d(l-1, s, samples[l-1][s]);
//                Vector3d E = new Vector3d(l, s-1, samples[l][s-1]);
//                
//                // Calculate normals of each facet
//                Vector3d ABC = Vector3d.getClockwiseSurfaceNormal(A, B, C);
//                Vector3d ACD = Vector3d.getClockwiseSurfaceNormal(A, C, D);
//                Vector3d ADE = Vector3d.getClockwiseSurfaceNormal(A, D, E);
//                Vector3d AEB = Vector3d.getClockwiseSurfaceNormal(A, E, B);
//                
//                // Take average
//                Vector3d norm = (ABC.add(ACD).add(ADE).add(AEB)).mult(0.25);
//                
//                out.write("<vec3> "+norm.x+" "+norm.y+" "+norm.z+" </vec3>\n");
//            }
//        }
//        out.write("</normals>\n");
//
//        // Optionally write a texcoord for each vertex based on projection
//        // of vertex onto corresponding face of ICQ model
////        if(TEXTURE)
////        {
////            out.write("<texcoords>\n");
////            for(int i=stride_i; i<=stride_i + Q; i++)
////            {
////                for(int j=stride_j; j<=stride_j + Q; j++)
////                {
////                    // Coordinate of this vertex in the ICQ face is (i,j)
////                    // i and j range from 0 to q-1: map this to range 0:1
////                    // for texture coordinate
////                    double u = (double)i / (double)(q-1);
////                    double v = (double)j / (double)(q-1);
////
////                    // Apply texture scale
////                    u *= tex_scale;
////                    v *= tex_scale;
////
////                    // Third texture coordinate always zero (only 2D
////                    // textures are supported).
////
////                    out.write("<vec3> "+u+" "+v+" 0 </vec3>\n"); 
////                }
////            }
////            out.write("</texcoords>\n");
////        }
//
//        // Write out triangular strips
//        out.write("<tristrips>\n");
//        for(int l=0; l<LINE_LAST_PIXEL; l++)
//        {
//            out.write("<tristrip>\t");
//            for(int s=0; s<=SAMPLE_LAST_PIXEL; s++)
//            {
//                out.write((l*(SAMPLE_LAST_PIXEL+1) + s) + " ");
//                out.write(((l+1)*(SAMPLE_LAST_PIXEL+1) + s) + " ");
//            }
//            out.write("\t</tristrip>\n");
//        }
//        out.write("</tristrips>\n");
//        
//        out.write("</mesh>\n");
//        
//        out.write("</pangu_model>\n");  // Closes pangu model
//        out.flush();
//        out.close();
//        
//    }
    
    
    /**
     * Checks the current parameter set to see if we have enough information to
     * interpret the data section.
     */
    boolean rationaliseParameters(StringBuilder message)
    {
        
        boolean valid = true;
        
        switch(SAMPLE_TYPE)
        {
            case PC_REAL:
            {
                // OK
                break;
            }
            case UNKNOWN:
            {
                message.append("\nUnknown SAMPLE_TYPE!"); 
                valid = false;
                break;
            }
        }
        
        switch(MAP_PROJECTION_TYPE)
        {
            case EQUIRECTANGULAR:
            {
                // Check that all necessary parameters have been set
                projection = new Equirectangular(DoublePar.LINE_PROJECTION_OFFSET.value,
                                                 DoublePar.SAMPLE_PROJECTION_OFFSET.value,
                                                 DoublePar.MAP_SCALE.value,
                                                 DoublePar.CENTER_LONGITUDE.value,
                                                 DoublePar.CENTER_LATITUDE.value,
                                                 DoublePar.A_AXIS_RADIUS.value,
                                                 DoublePar.C_AXIS_RADIUS.value);
                break;
            }
            case SINUSOIDAL:
            {
                // Check that all necessary parameters have been set
                projection = new Sinusoidal(DoublePar.LINE_PROJECTION_OFFSET.value,
                                            DoublePar.SAMPLE_PROJECTION_OFFSET.value,
                                            DoublePar.MAP_RESOLUTION.value,
                                            DoublePar.CENTER_LONGITUDE.value,
                                            DoublePar.CENTER_LATITUDE.value);
                break;
            }
            case UNKNOWN:
            {
                message.append("\nUnknown MAP_PROJECTION_TYPE!");
                valid = false;
                break;
            }
        }
        
        if(IntPar.SAMPLE_BITS.value%8!=0)
        {
            message.append("\nSAMPLE_BITS = "+IntPar.SAMPLE_BITS.value
                             + ": not a whole number of bytes!");
            valid = false;
        }
        
        return valid;
    }
    
    
    
}