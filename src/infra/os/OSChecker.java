/**
 * 
 * Name:
 *  OSChecker.java
 * 
 * Purpose:
 *  This class detects what OS we are currently running in.
 * 
 *  On windows/linux machines, main consideration is correct use of 
 *  newline characters in script files and path separation characters.
 *  Also consider typical fonts available to Gnuplot. Note that \n will
 *  always work as a newline when used in JTextArea. The newline is
 *  only important when we write stuff to disk.
 * 
 * Language:
 * Java
 *
 * Author:
 * Nicholas Rowell
 * 
 */
package infra.os;

import java.io.File;

public class OSChecker 
{
    /** Get OS name. */
    private static String os = System.getProperty("os.name").toLowerCase();
    
    /** Get file path separator symbol for this OS. */
    public static String pathSep = File.separator;
    
    /** Get newline symbol for this OS. */
    public static String newline = System.getProperty("line.separator");
    
    
    /** Enums to represent the different OS types. */
    public enum PLATFORM{WINDOWS,MAC,UNIX,SOLARIS,UNKNOWN}
    
    /** Get the OS type. */
    public static PLATFORM getOS()
    {
    
        if(os.indexOf("win") >= 0) return PLATFORM.WINDOWS;
        if(os.indexOf("mac") >= 0) return PLATFORM.MAC;
        if(os.indexOf("nix") >= 0 || os.indexOf("nux") >= 0 || os.indexOf("aix") > 0)
            return PLATFORM.UNIX;
        if(os.indexOf("sunos") >= 0) return PLATFORM.SOLARIS;
        
        return PLATFORM.UNKNOWN;
    }
    
    /** Check that user's platform is supported. */
    public static void checkOS()
    {
    
        // Check that platform is supported...
        switch(OSChecker.getOS())
        {
            case UNIX: break;
            case WINDOWS: break;
            case MAC: break;
            case SOLARIS:
            case UNKNOWN: 
                System.out.println("This software has not been "
                   + "tested on "+OSChecker.getOS().toString()+", exiting...");
                System.exit(1);
        }   
    
    
    }
    
    /**
     * Get the font to be used in plot scripts. Linux and Windows tend to have
     * a different set of basic fonts.
     * 
     * @return 
     */
    public static String getFont(){
        
        switch(getOS())
        {
            case UNIX:
            case MAC: return "helvetica";    
            default:   return "arial";
        }
    }
    
}
