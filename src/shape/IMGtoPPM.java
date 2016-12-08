package shape;

import Jama.Matrix;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

import numeric.dynamics.Pose;
import numeric.geom.dim3.Quaternion;
import numeric.geom.dim3.Vector3d;

/**
 * This program reads DTM files in IMG format, extracts a region and writes it
 * to a PPM file that can then be used by dem2pan to create a PANGU surface
 * model.
 * 
 * All parameters are entered manually rather than attempting to read them
 * automatically from the image header. This proved to be challenging due to
 * several factors, including that the units for the same quantity sometimes
 * differ between files (such as MAP_SCALE).
 * 
 * IMG files are typically huge areas and cannot be handled fully by PANGU (or
 * Java for that matter) due to the fact they can't be contained in memory. For
 * that reason, we only extract a sub-region.
 * 
 * To do:
 *  - figure out transformation to an established mars reference frame, so
 * that we can eventually use these models as a target for MSL landing simulations.
 *  - implement conversion to PXN file. This way, it is possible to correct for
 * the map projection distortion directly by calculating the 3D position of each
 * sample in the Mars centred reference frame (a DEM places all points on a 
 * regular grid).
 *  - implement more map projections (e.g. stereographic).
 * 
 * 
 * @author nickrowell
 */
public class IMGtoPPM
{
    
    // IO parameters
    
    /** Handle to input IMG file. */
    File img_file;
    /** Final output PXN file. */
    File ppm_file;
    
    // Data processing parameters
    
    /** Data type of samples {PC_REAL, UNSIGNED_INTEGER, ...}. */
    PDS_SAMPLE_TYPE SAMPLE_TYPE;
    /** Endian-ness for sample data type (usually LITTLE_ENDIAN). */
    ByteOrder ENDIANNESS;
    /** Map projection type {EQUIRECTANGULAR, SINUSOIDAL, ...}. */
    PDS_PROJECTION_TYPE MAP_PROJECTION_TYPE;
    /** Size of 'record': a basic unit of memory [bytes]. */
    int RECORD_BYTES;
    /** 
     * Indicates record number of start of image data [bytes]. 
     * Often (but not always! e.g. MOLA) labelled ^IMAGE. 
     */
    int IMAGE_POINTER;
    /** Size of IMG data [lines*samples]. Line number increases N/S. */
    int LINES;
    /** Size of IMG data [lines*samples]. Sample number increases E/W. */
    int LINE_SAMPLES;
    /** Number of bits in each sample (should be a multiple of 8). */
    int SAMPLE_BITS;
    /** Default value for missing data. */
    int MISSING_DATA;
    /** Offset applied to height data [metres]. */
    double OFFSET = 0;
    /** Scale factor applied to height data [dimensionless]. */
    double SCALING = 1;
    
    /** Number of channels per sample. */
    int BANDS;
    /** Centre of projection, latitude [degrees]. */
    double CENTER_LATITUDE;
    /** Centre of projection, longitude [degrees]. */
    double CENTER_LONGITUDE;
    /** Offset from samples origin, LINES direction [pixels]. */
    double LINE_PROJECTION_OFFSET;
    /** Offset from samples origin, SAMPLES direction [pixels]. */
    double SAMPLE_PROJECTION_OFFSET;
    /** 
     * Rotation of projection [degrees] - note we assume this is zero. The 
     * conversion will give dodgy results if not.
     */
    double MAP_PROJECTION_ROTATION;
    
    /** 
     * Equatorial radius of reference ellipsoid model, usually set equal to the
     * local radius at the centre of projection to avoid confusion [km].
     */
    double A_AXIS_RADIUS;
    
    double B_AXIS_RADIUS;
    /** 
     * Polar radius of reference ellipsoid model, usually set equal to the
     * local radius at the centre of projection to avoid confusion [km].
     */
    double C_AXIS_RADIUS;
    /** Physical scale of map [metres/pixel] */
    double MAP_SCALE;
    /** Angular scale of pixels [pixel/degree] */
    double MAP_RESOLUTION;

    
    static enum PDS_SAMPLE_TYPE
    {
        PC_REAL, UNSIGNED_INTEGER, MSB_INTEGER, LSB_INTEGER;
    }
    
    static enum PDS_PROJECTION_TYPE
    {
        EQUIRECTANGULAR, SINUSOIDAL, SIMPLE_CYLINDRICAL, POLAR_STEREOGRAPHIC;
    }
    
    public static void main(String[] args) throws IOException
    {
        processHiRISE();
//        processHRSC();
//        processMOLA();
//        processLOLA();
//        processLROC();
    }
    
    
    static void processHiRISE()
    {
        
        // Latitude and longitide of MSL landing site. Extracted topography will
        // be centred on this point, but all vertices are transformed to Mars or
        // MSL_TOPO frame so origin should lie at landing site and not at the
        // point specified here.
        double lat = Math.toRadians(-4.5895);
        double lon = Math.toRadians(137.4417);
        
        File parent = new File("/home/nickrowell/Projects/PANGU_4/Documents/TN01/"
                + "Planetary_Science/Available_Data/Shape_Models/Mars/HiRISE/Gale_Crater/"
                + "mrohr_0001_0001/dtm/psp/orb_010500_010599/psp_010573_1755_psp_010639_1755");
        
        String name = "dteec_010573_1755_010639_1755_u01";
        
        IMGtoPPM imgToPPM = new IMGtoPPM(new File(parent, name.concat(".img")),
                                         new File(parent, name.concat(".pxn")));
        
        // Set parameters fo HiRISE data
        imgToPPM.RECORD_BYTES             = 28744;
        imgToPPM.IMAGE_POINTER            = 2;
        imgToPPM.LINES                    = 16392;
        imgToPPM.LINE_SAMPLES             = 7186;
        imgToPPM.SAMPLE_BITS              = 32;                 // Multiple of 8
        imgToPPM.SAMPLE_TYPE              = PDS_SAMPLE_TYPE.PC_REAL;
        imgToPPM.ENDIANNESS               = ByteOrder.LITTLE_ENDIAN;
        imgToPPM.MISSING_DATA             = 0xFF7FFFFB;
        imgToPPM.BANDS                    = 1;
        
        imgToPPM.MAP_PROJECTION_TYPE      = PDS_PROJECTION_TYPE.EQUIRECTANGULAR;
        imgToPPM.CENTER_LATITUDE          = 0.0;                // [degrees]
        imgToPPM.CENTER_LONGITUDE         = 137.41;             // [degrees]
        imgToPPM.LINE_PROJECTION_OFFSET   = -255116.5;
        imgToPPM.SAMPLE_PROJECTION_OFFSET = 3827.5;
        imgToPPM.MAP_PROJECTION_ROTATION  = 0.0;                // Only zero supported
        imgToPPM.A_AXIS_RADIUS            = 3396.19;
        imgToPPM.B_AXIS_RADIUS            = 3396.19;
        imgToPPM.C_AXIS_RADIUS            = 3396.19;
        imgToPPM.MAP_SCALE                = 1.0117648673518;    // [metres/pixel]
        imgToPPM.MAP_RESOLUTION           = 58585.44750467;     // [pixels/degree]

        // Height values are relative to the Mars equipotential surface. We require
        // that the landing site is at an altitude of zero (the MSL_LANDING_SITE
        // is at 0,0,0 in the MSL_TOPO frame) so we apply the following offset to
        // achieve this. This was determined by starting with the mean Mars radius
        // of 3396190 and subtracting the residual altitude at the origin of the
        // resulting PANGU model.
        imgToPPM.OFFSET = 3395636.66; // metres
        
        // Transformation from Mars body frame to MSL_TOPO frame
        Vector3d origin = new Vector3d(-2489864.49476611, 2286205.60053228, -271345.826044048);
        
        Matrix R = new Matrix(new double[][]{{-0.0596355988672989,  0.676340057309038,    -0.734171452882054},
                                             { 0.0547576948095413,  0.736589524008594,     0.67411977272584},
                                             { 0.99671720673764,   -1.22062653688634e-16, -0.0809617798292272}});
        Quaternion q = new Quaternion(R);
        
        Pose IAU_MARS_2_MSL_TOPO = new Pose(origin, q);
        
        // Parse IMG model from file.
        try
        {
//            imgToPPM.extractPPM(lat, lon, true);
            imgToPPM.extractPXN(lat, lon, IAU_MARS_2_MSL_TOPO, true);
        }
        catch(IOException ioe)
        {
            System.err.println("Error extracting IMG topography: "+ioe.getMessage());
        }
    
    }
    
    
    static void processHRSC()
    {
        
        // Latitude and longitide of MSL landing site (extracted DEM will be
        // centred on this point).
        double lat = Math.toRadians(-4.5895);
        double lon = Math.toRadians(137.4417);
        
        File parent = new File("/home/nickrowell/Projects/PANGU_4/Documents/TN01/"
                + "Planetary_Science/Available_Data/Shape_Models/Mars/HRSC/mexhrs_2001/data/1927");
        
        String name = "h1927_0000_nd4";
        
        IMGtoPPM imgToPPM = new IMGtoPPM(new File(parent, name.concat(".img")),
                                         new File(parent, name.concat(".ppm")));
        
        // Set parameters fo HRSC data
        imgToPPM.RECORD_BYTES             = 7160;
        imgToPPM.IMAGE_POINTER            = 6;
        imgToPPM.LINES                    = 28824;
        imgToPPM.LINE_SAMPLES             = 7160;
        imgToPPM.SAMPLE_BITS              = 8;                 // Multiple of 8
        imgToPPM.SAMPLE_TYPE              = PDS_SAMPLE_TYPE.UNSIGNED_INTEGER;
        imgToPPM.ENDIANNESS               = ByteOrder.LITTLE_ENDIAN;
        imgToPPM.MISSING_DATA             = 0xFF7FFFFB;
        imgToPPM.BANDS                    = 1;
        
        imgToPPM.MAP_PROJECTION_TYPE      = PDS_PROJECTION_TYPE.SINUSOIDAL;
        imgToPPM.CENTER_LATITUDE          = 0.0;                // [degrees]
        imgToPPM.CENTER_LONGITUDE         = 138.0;              // [degrees]
        imgToPPM.LINE_PROJECTION_OFFSET   = -9297.5;
        imgToPPM.SAMPLE_PROJECTION_OFFSET = 4689.1;
        imgToPPM.MAP_PROJECTION_ROTATION  = 0.0;                // Only zero supported
        imgToPPM.A_AXIS_RADIUS            = 3396.0;
        imgToPPM.B_AXIS_RADIUS            = 3396.0;
        imgToPPM.C_AXIS_RADIUS            = 3396.0;
        imgToPPM.MAP_SCALE                = 12.5;               // [metres/pixel]
        imgToPPM.MAP_RESOLUTION           = 4741.71043093333;   // [pixels/degree]

        // Parse IMG model from file.
        try
        {
            imgToPPM.extractPPM(lat, lon, true);
        }
        catch(IOException ioe)
        {
            System.err.println("Error extracting IMG topography: "+ioe.getMessage());
        }
    
    }
    
    
    
    static void processMOLA()
    {
        
        // Latitude and longitide of MSL landing site (extracted DEM will be
        // centred on this point).
        double lat = Math.toRadians(-4.5895);
        double lon = Math.toRadians(137.4417);
        
        File parent = new File("/home/nickrowell/Projects/PANGU_4/Documents/TN01/"
                + "Planetary_Science/Available_Data/Shape_Models/Mars/MOLA/MEGR00N090HB");
        
        String name = "megr00n090hb";
        
        IMGtoPPM imgToPPM = new IMGtoPPM(new File(parent, name.concat(".img")),
                                         new File(parent, name.concat(".ppm")));
        
        // Set parameters fo MOLA data
        imgToPPM.RECORD_BYTES             = 11520;
        imgToPPM.IMAGE_POINTER            = 1;
        imgToPPM.LINES                    = 5632;
        imgToPPM.LINE_SAMPLES             = 11520;
        imgToPPM.SAMPLE_BITS              = 16;                 // Multiple of 8
        imgToPPM.SAMPLE_TYPE              = PDS_SAMPLE_TYPE.MSB_INTEGER;
        imgToPPM.ENDIANNESS               = ByteOrder.LITTLE_ENDIAN;
        imgToPPM.MISSING_DATA             = Integer.MAX_VALUE;
        imgToPPM.BANDS                    = 1;
        imgToPPM.OFFSET                   = 3396000;
        imgToPPM.SCALING                  = 1;
        
        imgToPPM.MAP_PROJECTION_TYPE      = PDS_PROJECTION_TYPE.SIMPLE_CYLINDRICAL;
        imgToPPM.CENTER_LATITUDE          = 0.0;                // [degrees]
        imgToPPM.CENTER_LONGITUDE         = 180.0;              // [degrees]
        imgToPPM.LINE_PROJECTION_OFFSET   = 0.5;
        imgToPPM.SAMPLE_PROJECTION_OFFSET = 11520.5;
        imgToPPM.MAP_PROJECTION_ROTATION  = 0.0;                // Only zero supported
        imgToPPM.A_AXIS_RADIUS            = 3396.0;
        imgToPPM.B_AXIS_RADIUS            = 3396.0;
        imgToPPM.C_AXIS_RADIUS            = 3396.0;
        imgToPPM.MAP_SCALE                = 463.0;              // [metres/pixel]
        imgToPPM.MAP_RESOLUTION           = 128.0;              // [pixels/degree]

        // Parse IMG model from file.
        try
        {
            imgToPPM.extractPPM(lat, lon, true);
        }
        catch(IOException ioe)
        {
            System.err.println("Error extracting IMG topography: "+ioe.getMessage());
        }
    
    }
    
    static void processLOLA()
    {
        
        // Latitude and longitide (extracted DEM will be
        // centred on this point).
        double lat = Math.toRadians(-80.76);
        double lon = Math.toRadians(326.26);
                
        File parent = new File("/home/nickrowell/Projects/PANGU_4/Documents/TN01/"
                + "Planetary_Science/Available_Data/Shape_Models/Moon/LRO/75S");
        
        String name = "LDEM_75S_60M";
        
        IMGtoPPM imgToPPM = new IMGtoPPM(new File(parent, name.concat(".IMG")),
                                         new File(parent, name.concat(".pxn")));
        
        // Set parameters fo MOLA data
        imgToPPM.RECORD_BYTES             = 30496;
        imgToPPM.IMAGE_POINTER            = 1;
        imgToPPM.LINES                    = 15248;
        imgToPPM.LINE_SAMPLES             = 15248;
        imgToPPM.SAMPLE_BITS              = 16;                 // Multiple of 8
        imgToPPM.SAMPLE_TYPE              = PDS_SAMPLE_TYPE.LSB_INTEGER;
        imgToPPM.ENDIANNESS               = ByteOrder.LITTLE_ENDIAN;
        imgToPPM.MISSING_DATA             = Integer.MAX_VALUE;
        imgToPPM.BANDS                    = 1;
        imgToPPM.OFFSET                   = 1737400;
        imgToPPM.SCALING                  = 0.5;
        
        imgToPPM.MAP_PROJECTION_TYPE      = PDS_PROJECTION_TYPE.POLAR_STEREOGRAPHIC;
        imgToPPM.CENTER_LATITUDE          = -90.0;              // [degrees]
        imgToPPM.CENTER_LONGITUDE         =   0.0;              // [degrees]
        imgToPPM.LINE_PROJECTION_OFFSET   = 7623.5;
        imgToPPM.SAMPLE_PROJECTION_OFFSET = 7623.5;
        imgToPPM.MAP_PROJECTION_ROTATION  = 0.0;                // Only zero supported
        imgToPPM.A_AXIS_RADIUS            = 1737.400;
        imgToPPM.B_AXIS_RADIUS            = 1737.400;
        imgToPPM.C_AXIS_RADIUS            = 1737.400;
        imgToPPM.MAP_SCALE                = 60.0;              // [metres/pixel]
        imgToPPM.MAP_RESOLUTION           = 505.389;           // [pixels/degree]

        // Parse IMG model from file.
        try
        {
//            imgToPPM.extractPPM(lat, lon, true);
            imgToPPM.extractPXN(lat, lon, new Pose(), true);
        }
        catch(IOException ioe)
        {
            System.err.println("Error extracting IMG topography: "+ioe.getMessage());
        }
    
    }
    
    static void processLROC()
    {
        
        File parent = new File("/home/nickrowell/Projects/PANGU_4/Documents/TN01/"
                + "Planetary_Science/Available_Data/Shape_Models/Moon/LRO/LROC/Einmarta");
        
        String name = "NAC_DTM_EIMMARTA_E240N0660";
        
        IMGtoPPM imgToPPM = new IMGtoPPM(new File(parent, name.concat(".IMG")),
                                         new File(parent, name.concat(".pxn")));
        
        // Set parameters fo MOLA data
        imgToPPM.RECORD_BYTES             = 12336;
        imgToPPM.IMAGE_POINTER            = 2;
        imgToPPM.LINES                    = 6610;
        imgToPPM.LINE_SAMPLES             = 3084;
        imgToPPM.SAMPLE_BITS              = 32;                 // Multiple of 8
        imgToPPM.SAMPLE_TYPE              = PDS_SAMPLE_TYPE.PC_REAL;
        imgToPPM.ENDIANNESS               = ByteOrder.LITTLE_ENDIAN;
        imgToPPM.MISSING_DATA             = 0xFF7FFFFB;
        imgToPPM.BANDS                    = 1;
        imgToPPM.OFFSET                   = 1737400;
        imgToPPM.SCALING                  = 0.5;
        
        imgToPPM.MAP_PROJECTION_TYPE      = PDS_PROJECTION_TYPE.EQUIRECTANGULAR;
        imgToPPM.CENTER_LATITUDE          =  24.0;              // [degrees]
        imgToPPM.CENTER_LONGITUDE         = 180.0;              // [degrees]
        imgToPPM.LINE_PROJECTION_OFFSET   = 149569.5;
        imgToPPM.SAMPLE_PROJECTION_OFFSET = 635219.5;
        imgToPPM.MAP_PROJECTION_ROTATION  = 0.0;                // Only zero supported
        imgToPPM.A_AXIS_RADIUS            = 1737.400;
        imgToPPM.B_AXIS_RADIUS            = 1737.400;
        imgToPPM.C_AXIS_RADIUS            = 1737.400;
        imgToPPM.MAP_SCALE                = 5.0;              // [metres/pixel]
        imgToPPM.MAP_RESOLUTION           = 6064.670;         // [pixels/degree]

        // Parse IMG model from file.
        try
        {
//            imgToPPM.extractPPM(lat, lon, true);
            imgToPPM.extractPXN(1542, 3305, new Pose(), true);
        }
        catch(IOException ioe)
        {
            System.err.println("Error extracting IMG topography: "+ioe.getMessage());
        }
    
    }
    
    public IMGtoPPM(File pimg_file, File pppm_file)
    {
        img_file  = pimg_file;
        ppm_file  = pppm_file;
    }
    
    /**
     * Read IMG file and extract a region around the selected site, writing the
     * result to a PPM image that can be converted to a PANGU DEM by dem2pan.
     * @param lat
     * @param lon
     * @param verbose
     * @throws IOException 
     */
    private void extractPPM(double lat, double lon, boolean verbose) throws IOException
    {
        
        // Get a Projection object used to interpret map.
        Projection projection = null;
        
        switch(MAP_PROJECTION_TYPE)
        {
            case EQUIRECTANGULAR:
            {
                // Check that all necessary parameters have been set
                projection = new Equirectangular(LINE_PROJECTION_OFFSET,
                                                 SAMPLE_PROJECTION_OFFSET,
                                                 MAP_SCALE,
                                                 CENTER_LONGITUDE,
                                                 CENTER_LATITUDE,
                                                 A_AXIS_RADIUS,
                                                 C_AXIS_RADIUS);
                break;
            }
            case SINUSOIDAL:
            {
                // Check that all necessary parameters have been set
                projection = new Sinusoidal(LINE_PROJECTION_OFFSET,
                                            SAMPLE_PROJECTION_OFFSET,
                                            MAP_RESOLUTION,
                                            CENTER_LONGITUDE,
                                            CENTER_LATITUDE);
                break;
            }
            case SIMPLE_CYLINDRICAL:
            {
                // Check that all necessary parameters have been set
                projection = new SimpleCylindrical(LINE_PROJECTION_OFFSET,
                                                   SAMPLE_PROJECTION_OFFSET,
                                                   MAP_RESOLUTION,
                                                   CENTER_LONGITUDE,
                                                   CENTER_LATITUDE);
            }
            case POLAR_STEREOGRAPHIC:
            {
                // Check that all necessary parameters have been set
                projection = new PolarStereographic(MAP_SCALE,
                                                    A_AXIS_RADIUS,
                                                    LINES,
                                                    PolarStereographic.HEMI.SOUTH);
            }
                
            default:
            {
                
            }
        }
        
        // Open reader on input file
        FileInputStream freader = new FileInputStream(img_file);
        
        // Skip past header to start of data section
        int bytes_so_far = 0;
        while(bytes_so_far < RECORD_BYTES*(IMAGE_POINTER-1)) { freader.read(); bytes_so_far++;}
        
        // Get indices in data array of sample closest to specfied point.
        int[] ls = projection.getLineSampleIndices(lon, lat);
        int s0   = ls[0];
        int l0   = ls[1];
        
        
        
        
        if(verbose)
            System.out.println("Centering extracted DTM on "+l0+", "+s0);
                
        // Set radius of region to extract (if it's odd we get a central pixel)
        int WIDTH  = 4097;
        int HEIGHT = 4097;
        
        // Stats used to digitize height values:
    
        // Maximum valid height.
        float hmax = 0.0f;
        // Minimum valid height.
        float hmin = 0.0f;
        // Height difference between consecutive levels in tga16.
        float hres = 0.0f;
    
        // Temporary storage for each sample
        int nbytes = SAMPLE_BITS/8;
        
        byte[] bsample = new byte[nbytes];
        
        if(verbose)
            System.out.println("First pass over the data: reading IMG and "
                    + "computing height stats...");
        
        // Record max & min heights (to allow later digitizing to 16-bit TGA)
        boolean initialised_heights = false;
        
        // Read all samples into a 2D array
        float[][] samples = new float[LINES][LINE_SAMPLES];
        
        // Count number of valid & missing samples in full grid (sanity check...)
        long n_valid   = 0;
        long n_invalid = 0;
        
        // Read in the data in samples of 32 bits (4 bytes)
        for(int l=0; l<LINES; l++)
        {
            for(int s=0; s<LINE_SAMPLES; s++)
            {
                
                int nread = freader.read(bsample);        // Read one sample
                assert(nread==bsample.length);
                ByteBuffer bb = ByteBuffer.wrap(bsample); // Create a byte buffer and wrap the array
                bb.order(ENDIANNESS);                     // Set endian-ness
                
                boolean isMissing = false;                // Stores valid/invalid status of sample
                float   sample;                           // Stores converted floating point value
                
                switch(SAMPLE_TYPE)
                {
                        case PC_REAL:
                        {
                            // 32-bit samples: read 4 bytes into an integer
                            int data = bb.getInt();
                            isMissing = (data==MISSING_DATA);
                            sample = Float.intBitsToFloat(data);
                            
                            break;
                        }
                        case UNSIGNED_INTEGER:
                        {
                            // 8-bit samples: read one byte
                            byte data = bb.get();
                            isMissing = (data==MISSING_DATA);
                            // Careful! Java bytes are signed.
                            sample = (float)(data & 0xFF);
                            break;
                        }
                        case MSB_INTEGER:
                        {
                            // 16-bit samples: read two bytes into an integer.
                            // Need to reverse bit order for MSB_INTEGER.
                            short data = Short.reverseBytes(bb.getShort());
                            isMissing = (data==MISSING_DATA);
                            sample = (float)(data);
                            break;
                        }
                        case LSB_INTEGER:
                        {
                            // 16-bit samples: read two bytes into an integer.
                            short data = bb.getShort();
                            isMissing = (data==MISSING_DATA);
                            sample = (float)(data);
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
                    n_valid++;
                    
                }
                else
                {
                    samples[l][s] = Float.NaN;
                    n_invalid++;
                }
                
            }
        }
        
        // Verify that end of file has been reached
        int n_tail = 0;
        while(freader.read() != -1) n_tail++;
        if(n_tail!=0)
            System.err.println("Encountered "+n_tail+" additional bytes after "
                                + "end of data section!");
        
        // Scale height to restrict to range +/-32768 when digitised.
        float scale_factor = 32768/Math.max(Math.abs(hmax),Math.abs(hmin));
        
        System.out.println("hmax = "+hmax);
        System.out.println("hmin = "+hmin);
        System.out.println("Encountered "+n_valid+" valid data (over full DTM)");
        System.out.println("Encountered "+n_invalid+" invalid data (over full DTM)");
        
        // Metres-per-digital-level
        hres = 1.0f/scale_factor;
        
        freader.close();
        
        if(verbose)
            System.out.println("Second pass over data: digitizing samples and"
                    + " writing PPM file...");
        
        // Open writer for final PPM file
        FileOutputStream ppmwriter = new FileOutputStream(ppm_file);
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
                if(l<0 || l>=LINES || s<0 || s>=LINE_SAMPLES)
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
        
        
        // Form .pan filename from input
        String img = img_file.getName();
        String ppm = ppm_file.getName();
        String pan = img.substring(0, img.lastIndexOf('.')).concat(".pan");
        
        // Print dem2pan command...
        System.out.println("/opt/pangu/source/pangu_testing/bin/linux_x86/dem2pan \\");
        System.out.println("-o "+pan+" \\");
        System.out.println("-z "+hres+" -s "+MAP_SCALE+" \\");
        System.out.println(ppm);

    }
    
    
    
    private Projection getProjection()
    {
        switch(MAP_PROJECTION_TYPE)
        {
            case EQUIRECTANGULAR:
            {
                // Check that all necessary parameters have been set
                return new Equirectangular(LINE_PROJECTION_OFFSET,
                                           SAMPLE_PROJECTION_OFFSET,
                                           MAP_SCALE,
                                           CENTER_LONGITUDE,
                                           CENTER_LATITUDE,
                                           A_AXIS_RADIUS,
                                           C_AXIS_RADIUS);
            }
            case SINUSOIDAL:
            {
                // Check that all necessary parameters have been set
                return new Sinusoidal(LINE_PROJECTION_OFFSET,
                                      SAMPLE_PROJECTION_OFFSET,
                                      MAP_RESOLUTION,
                                      CENTER_LONGITUDE,
                                      CENTER_LATITUDE);
            }
            case SIMPLE_CYLINDRICAL:
            {
                // Check that all necessary parameters have been set
                return new SimpleCylindrical(LINE_PROJECTION_OFFSET,
                                             SAMPLE_PROJECTION_OFFSET,
                                             MAP_RESOLUTION,
                                             CENTER_LONGITUDE,
                                             CENTER_LATITUDE);
            }
            case POLAR_STEREOGRAPHIC:
            {
                // Check that all necessary parameters have been set
                return new PolarStereographic(MAP_SCALE,
                                              A_AXIS_RADIUS,
                                              LINES,
                                              PolarStereographic.HEMI.SOUTH);
            }
                
            default:
            {
                return null;
            }
        }
    }
    
    
    
    /**
     * Read the IMG file and extract a region surrounding the selected site,
     * writing the results to a PXN file in absolute position units, ie in
     * the planet coordinate frame.
     * @param lat
     * @param lon
     * @param verbose
     * @throws IOException 
     */
    
    private void extractPXN(double lat, double lon, Pose pose, boolean verbose) throws IOException
    {
        
        
        // Get a Projection object used to interpret map.
        Projection projection = getProjection();
        // Get indices in data array of sample closest to specfied point.
        int[] ls = projection.getLineSampleIndices(lon, lat);
        int I_0 = ls[0];
        int J_0 = ls[1];
        
        if(verbose)
        {
            System.out.println("Extracting region surrounding long/lat: "+
                    Math.toDegrees(lon)+", "+Math.toDegrees(lat));
            System.out.println("Centering extracted DTM on pixel "+J_0+", "+I_0);
            double[] longlat = projection.getLongLat(I_0, J_0);
            System.out.println("Converts back to long/lat: "+
                    Math.toDegrees(longlat[0])+", "+Math.toDegrees(longlat[1]));
            
        }
        
        
        extractPXN(I_0, J_0, pose, verbose);
    }
    
    private void extractPXN(int I_0, int J_0, Pose pose, boolean verbose) throws IOException
    {
        // Set radius of region to extract (if it's odd we get a central pixel)
        int WIDTH  = 2049;
        int HEIGHT = 2049;
        
        // Get a Projection object used to interpret map.
        Projection projection = getProjection();
        
        // Open reader on input file
        FileInputStream freader = new FileInputStream(img_file);
        
        // Skip past header to start of data section
        int bytes_so_far = 0;
        while(bytes_so_far < RECORD_BYTES*(IMAGE_POINTER-1)) { freader.read(); bytes_so_far++;}
        
        // Read all samples into a 2D array
        Vector3d[][] samples = new Vector3d[WIDTH][HEIGHT];
        
        // Initialise all samples to the height offset. When we read over the DTM,
        // any samples found inside the DTM area will have their values reset.
        for(int j=0; j<WIDTH; j++)
        {
            for(int i=0; i<HEIGHT; i++)
            {
                // Corresponding line & sample indices in DTM grid
                int J = (int)Math.rint(J_0 - (HEIGHT-1)/2 + j);
                int I = (int)Math.rint(I_0 - (WIDTH-1)/2 + i);
                
                double[] latlon = projection.getLongLat(I, J);
                double slon = latlon[0];
                double slat = latlon[1];
                    
                // Convert to a 3D point in the frame of the model
                double x = OFFSET*Math.cos(slat)*Math.cos(slon);
                double y = OFFSET*Math.cos(slat)*Math.sin(slon);
                double z = OFFSET*Math.sin(slat);
                
                Vector3d model = new Vector3d(x,y,z);
                
                // Now apply any transformation
                samples[j][i] = pose.toBodyFrame(model);
            }
        }
        
        
        // Temporary storage for each sample
        int nbytes = SAMPLE_BITS/8;
        
        byte[] bsample = new byte[nbytes];
        
        if(verbose)
            System.out.println("First pass over the data: reading IMG...");
        
        // Read in the data in samples
        for(int J=0; J<LINES; J++)
        {
            for(int I=0; I<LINE_SAMPLES; I++)
            {
                
                int nread = freader.read(bsample);        // Read one sample
                
                // If this sample lies in the region of interest, we convert it
                // to a Vector3d
                if(Math.abs(I-I_0) <= (WIDTH-1)/2 && Math.abs(J-J_0) <= (HEIGHT-1)/2)
                {
                    
                    // parse the sample bits to a number
                    assert(nread==bsample.length);
                    ByteBuffer bb = ByteBuffer.wrap(bsample); // Create a byte buffer and wrap the array
                    bb.order(ENDIANNESS);                     // Set endian-ness
                    
                    boolean isMissing = false;                // Stores valid/invalid status of sample
                    float   height;                           // Stores converted floating point value
                    
                    switch(SAMPLE_TYPE)
                    {
                            case PC_REAL:
                            {
                                // 32-bit samples: read 4 bytes into an integer
                                int data = bb.getInt();
                                isMissing = (data==MISSING_DATA);
                                height = Float.intBitsToFloat(data);
                                break;
                            }
                            case UNSIGNED_INTEGER:
                            {
                                // 8-bit samples: read one byte
                                byte data = bb.get();
                                isMissing = (data==MISSING_DATA);
                                // Careful! Java bytes are signed.
                                height = (float)(data & 0xFF);
                                break;
                            }
                            case MSB_INTEGER:
                            {
                                // 16-bit samples: read two bytes into an integer.
                                // Need to reverse bit order for MSB_INTEGER.
                                short data = Short.reverseBytes(bb.getShort());
                                isMissing = (data==MISSING_DATA);
                                height = (float)(data);
                                break;
                            }
                            case LSB_INTEGER:
                            {
                                // 16-bit samples: read two bytes into an integer.
                                short data = bb.getShort();
                                isMissing = (data==MISSING_DATA);
                                height = (float)(data);
                                break;
                            }

                            default:
                                height = 0;
                    }
                    
                    double r = OFFSET;
                    
                    if(!isMissing)
                    {
                        // Get the radius of the point
                        r = (height*SCALING)+ OFFSET;
                    }
                    
                    // Convert the sample number to a latitude/longitude
                    double[] latlon = projection.getLongLat(I, J);
                    double slon = latlon[0];
                    double slat = latlon[1];
                    
                    // Convert to a 3D point in the body frame
                    double x = r*Math.cos(slat)*Math.cos(slon);
                    double y = r*Math.cos(slat)*Math.sin(slon);
                    double z = r*Math.sin(slat);
                    
                    Vector3d model = new Vector3d(x,y,z);
                
                    // Now apply any transformation
                    samples[J-J_0+(HEIGHT-1)/2][I-I_0+(WIDTH-1)/2] = pose.toBodyFrame(model);
                    
                }
                
            }
        }
        
        // Verify that end of file has been reached
        int n_tail = 0;
        while(freader.read() != -1) n_tail++;
        if(n_tail!=0)
            System.err.println("Encountered "+n_tail+" additional bytes after "
                                + "end of data section!");
        
        System.out.println("Parsed TAB; writing PXN...");
        
        // Open writer on output file
        BufferedWriter out = new BufferedWriter(new FileWriter(ppm_file));
        
        // Write PXN header
        out.write("<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
        out.write("<pangu_model ver=\"0\">\n"); 
        
        // Cease recursion and write out mesh for this node
        out.write("<mesh>\n");
        
        // Write out vertex positions
        out.write("<positions>\n");
        for(int J=0; J<HEIGHT; J++)
        {
            for(int I=0; I<WIDTH; I++)
            {
                Vector3d pixel = samples[J][I];
                out.write("<vec3> "+pixel.getX()+" "+pixel.getY()+" "+pixel.getZ()+" </vec3>\n");
            }
        }
        out.write("</positions>\n");
        // Write out surface normals
        out.write("<normals>\n");
        for(int J=0; J<HEIGHT; J++)
        {
            for(int I=0; I<WIDTH; I++)
            {
                //
                // Use the 4 nearest neighbours. These vertices are laid out:
                //
                //           l-->
                //
                //   s        E
                //   |    D   A   B
                //   V        C
                //
                Vector3d norm = new Vector3d();
                
                // Four corner points - only one facet can be constructed
                if(J==0 && I==0)
                {
                    Vector3d A = samples[J][I];
                    Vector3d B = samples[J+1][I];
                    Vector3d C = samples[J][I+1];
                    norm = Vector3d.getClockwiseSurfaceNormal(A, B, C);
                    
                    if(Double.isNaN(norm.getX()) || Double.isNaN(norm.getY()) || Double.isNaN(norm.getZ()))
                    {
                        System.out.println("NaN on I,J = "+I+","+J);
                        System.out.println("A = "+A.toString());
                        System.out.println("B = "+B.toString());
                        System.out.println("C = "+C.toString());
                    }
                    
                }
                else if(J==0 && I==(samples[J].length-1))
                {
                    Vector3d A = samples[J][I];
                    Vector3d B = samples[J+1][I];
                    Vector3d E = samples[J][I-1];
                    norm = Vector3d.getClockwiseSurfaceNormal(A, E, B);
                    
                    if(Double.isNaN(norm.getX()) || Double.isNaN(norm.getY()) || Double.isNaN(norm.getZ()))
                    {
                        System.out.println("NaN on I,J = "+I+","+J);
                        System.out.println("A = "+A.toString());
                        System.out.println("B = "+B.toString());
                        System.out.println("E = "+E.toString());
                    }
                }
                else if(J==(samples.length-1) && I==0)
                {
                    Vector3d A = samples[J][I];
                    Vector3d C = samples[J][I+1];
                    Vector3d D = samples[J-1][I];
                    norm = Vector3d.getClockwiseSurfaceNormal(A, C, D);
                    
                    if(Double.isNaN(norm.getX()) || Double.isNaN(norm.getY()) || Double.isNaN(norm.getZ()))
                    {
                        System.out.println("NaN on I,J = "+I+","+J);
                        System.out.println("A = "+A.toString());
                        System.out.println("C = "+C.toString());
                        System.out.println("D = "+D.toString());
                    }
                }
                else if(J==(samples.length-1) && I==(samples[J].length-1))
                {
                    Vector3d A = samples[J][I];
                    Vector3d D = samples[J-1][I];
                    Vector3d E = samples[J][I-1];
                    norm = Vector3d.getClockwiseSurfaceNormal(A, D, E);
                    
                    if(Double.isNaN(norm.getX()) || Double.isNaN(norm.getY()) || Double.isNaN(norm.getZ()))
                    {
                        System.out.println("NaN on I,J = "+I+","+J);
                        System.out.println("A = "+A.toString());
                        System.out.println("D = "+D.toString());
                        System.out.println("E = "+E.toString());
                    }
                }
                // Edge sections: two facets can be constructed
                else if(J==0)
                {
                    
                    Vector3d A = samples[J][I];
                    Vector3d B = samples[J+1][I];
                    Vector3d C = samples[J][I+1];
                    Vector3d E = samples[J][I-1];

                    Vector3d ABC = Vector3d.getClockwiseSurfaceNormal(A, B, C);
                    Vector3d AEB = Vector3d.getClockwiseSurfaceNormal(A, E, B);

                    norm = (ABC.add(AEB)).mult(0.5);
                    
                    if(Double.isNaN(norm.getX()) || Double.isNaN(norm.getY()) || Double.isNaN(norm.getZ()))
                    {
                        System.out.println("NaN on I,J = "+I+","+J);
                        System.out.println("A = "+A.toString());
                        System.out.println("B = "+B.toString());
                        System.out.println("C = "+C.toString());
                        System.out.println("E = "+E.toString());
                        System.out.println("ABC = "+ABC.toString());
                        System.out.println("AEB = "+AEB.toString());
                        
                    }
                }
                else if(J==(samples.length-1))
                {
                    Vector3d A = samples[J][I];
                    Vector3d C = samples[J][I+1];
                    Vector3d D = samples[J-1][I];
                    Vector3d E = samples[J][I-1];
                    Vector3d ACD = Vector3d.getClockwiseSurfaceNormal(A, C, D);
                    Vector3d ADE = Vector3d.getClockwiseSurfaceNormal(A, D, E);

                    norm = (ACD.add(ADE)).mult(0.5);
                    
                    if(Double.isNaN(norm.getX()) || Double.isNaN(norm.getY()) || Double.isNaN(norm.getZ()))
                    {
                        System.out.println("NaN on I,J = "+I+","+J);
                        System.out.println("A = "+A.toString());
                        System.out.println("C = "+C.toString());
                        System.out.println("D = "+D.toString());
                        System.out.println("E = "+E.toString());
                        System.out.println("ACD = "+ACD.toString());
                        System.out.println("ADE = "+ADE.toString());
                    }
                }
                else if(I==0)
                {
                    Vector3d A = samples[J][I];
                    Vector3d B = samples[J+1][I];
                    Vector3d C = samples[J][I+1];
                    Vector3d D = samples[J-1][I];

                    Vector3d ABC = Vector3d.getClockwiseSurfaceNormal(A, B, C);
                    Vector3d ACD = Vector3d.getClockwiseSurfaceNormal(A, C, D);

                    norm = (ABC.add(ACD)).mult(0.5);
                    
                    if(Double.isNaN(norm.getX()) || Double.isNaN(norm.getY()) || Double.isNaN(norm.getZ()))
                    {
                        System.out.println("NaN on I,J = "+I+","+J);
                        System.out.println("A = "+A.toString());
                        System.out.println("B = "+B.toString());
                        System.out.println("C = "+C.toString());
                        System.out.println("D = "+D.toString());
                        System.out.println("ABC = "+ABC.toString());
                        System.out.println("ACD = "+ACD.toString());
                        
                    }
                }
                else if(I==(samples[J].length-1))
                {
                    Vector3d A = samples[J][I];
                    Vector3d B = samples[J+1][I];
                    Vector3d D = samples[J-1][I];
                    Vector3d E = samples[J][I-1];

                    Vector3d ADE = Vector3d.getClockwiseSurfaceNormal(A, D, E);
                    Vector3d AEB = Vector3d.getClockwiseSurfaceNormal(A, E, B);

                    norm = (ADE.add(AEB)).mult(0.5);
                    
                    if(Double.isNaN(norm.getX()) || Double.isNaN(norm.getY()) || Double.isNaN(norm.getZ()))
                    {
                        System.out.println("NaN on I,J = "+I+","+J);
                        System.out.println("A = "+A.toString());
                        System.out.println("B = "+B.toString());
                        System.out.println("D = "+D.toString());
                        System.out.println("E = "+E.toString());
                        System.out.println("ADE = "+ADE.toString());
                        System.out.println("AEB = "+AEB.toString());
                        
                    }
                }
                // Everywhere else: can construct 4 facets
                else
                {
                    
                
                    Vector3d A = samples[J][I];
                    Vector3d B = samples[J+1][I];
                    Vector3d C = samples[J][I+1];
                    Vector3d D = samples[J-1][I];
                    Vector3d E = samples[J][I-1];

                    // Calculate normals of each facet
                    Vector3d ABC = Vector3d.getClockwiseSurfaceNormal(A, B, C);
                    Vector3d ACD = Vector3d.getClockwiseSurfaceNormal(A, C, D);
                    Vector3d ADE = Vector3d.getClockwiseSurfaceNormal(A, D, E);
                    Vector3d AEB = Vector3d.getClockwiseSurfaceNormal(A, E, B);

                    // Take average
                    norm = (ABC.add(ACD).add(ADE).add(AEB)).mult(0.25);
                    
                    
                if(Double.isNaN(norm.getX()) || Double.isNaN(norm.getY()) || Double.isNaN(norm.getZ()))
                    {
                        System.out.println("NaN on I,J = "+I+","+J);
                        System.out.println("A = "+A.toString());
                        System.out.println("B = "+B.toString());
                        System.out.println("C = "+C.toString());
                        System.out.println("D = "+D.toString());
                        System.out.println("E = "+E.toString());
                        
                        System.out.println("ABC = "+ABC.toString());
                        System.out.println("ACD = "+ACD.toString());
                        System.out.println("ADE = "+ADE.toString());
                        System.out.println("AEB = "+AEB.toString());
                        
                    }
                }
                
                // In southern hemisphere, need to reverse normals
                norm.multEquals(1);
                
                
                
                out.write("<vec3> "+norm.getX()+" "+norm.getY()+" "+norm.getZ()+" </vec3>\n");
            }
        }
        out.write("</normals>\n");

        // Write out triangular strips
        out.write("<tristrips>\n");
        
        
        for(int J=0; J<HEIGHT-1; J++)
        {
            out.write("<tristrip>\t");
            
            for(int I=0; I<WIDTH; I++)
            {
                out.write((J*WIDTH + I) + " ");
                out.write(((J+1)*WIDTH + I) + " ");
            }
            out.write("\t</tristrip>\n");
        }
        out.write("</tristrips>\n");
        
        out.write("</mesh>\n");
        
        out.write("</pangu_model>\n");  // Closes pangu model
        out.flush();
        out.close();
        

    }
    
}