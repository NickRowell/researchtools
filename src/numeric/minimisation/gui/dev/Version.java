/**
 * 
 * Name:
 * Version.java
 * 
 * Purpose:
 * Used to provide a means of tracking release versions. The values specified
 * in this class define the version number of the software: this should match
 * the numbers given in EDITS.txt log file so we can track changes. This class 
 * provides methods to extract & format the version number string, which is 
 * useful for setting appropriate titles etc and means the version number for a 
 * given build of the software can be extracted.
 * 
 * Language:
 * Java
 *
 * Author:
 * Nicholas Rowell
 * 
 */
package numeric.minimisation.gui.dev;

/**
 * Used to record details of the current version of this software. Clients
 * may extract the individual version number components themselves or they
 * can use the supplied {@link #string() string()} method.
 * <p>
 * There are three types of release and they control the way that the version
 * number must be interpreted.
 * <p>
 * Development versions [<code>printf("%d.%02d(%d)", maj, min, edit)</code>]:
 * <ul>
 *    <li>new functionality being added and tested
 *    <li>patch level is zero
 *    <li>edit level is non-zero
 * </ul>
 * <p>
 * Release versions [<code>printf("%d.%02d", maj, min)</code>]:
 * <ul>
 *    <li>all changes frozen
 *    <li>patch level is zero
 *    <li>edit level is zero
 * </ul>
 * <p>
 * Patched versions [<code>printf("%d.%02dp%d", maj, min, patch)</code>]:
 * <ul>
 *    <li>changes made to released versions (eg bug fixes)
 *    <li>patch level is non-zero
 *    <li>edit level is zero
 * </ul>
 * <p>
 * Example version numbers:
 * <p>
 * <ul>
 *    <li>1.00    : release version (major = 1, minor =  0, edit = 0, patch = 0)
 *    <li>1.00p2  : second patch    (major = 1, minor =  0, edit = 0, patch = 2)
 *    <li>1.12(3) : development     (major = 1, minor = 12, edit = 3, patch = 0)
 * </ul>
 * <p>
 * Example of use:                    (mj, mn, ed, pa) : label
 * <ul>
 *   <li>development begins           ( 0,  0,  0,  0) : 0.00(0)
 *   <li>first change made            ( 0,  0,  1,  0) : 0.00(1)
 *   <li>more changes made            ( 0,  0, 23,  0) : 0.00(23)
 *   <li>stable so make minor release ( 0,  1,  0,  0) : 0.01
 *   <li>some more changes made       ( 0,  1,  4,  0) : 0.01(4)
 *   <li>first major release          ( 1,  0,  0,  0) : 1.00
 *   <li>more changes                 ( 1,  0,  7,  0) : 1.00(7)
 *   <li>fixed bug in release 1.00    ( 1,  0,  0,  1) : 1.00p1
 *   <li>add bug fix to dev version   ( 1,  0,  8,  0) : 1.00(8)
 *   <li>another bug in releas 1.00   ( 1,  0,  0,  2) : 1.00p2
 *   <li>add bug fix to dev version   ( 1,  0,  9,  0) : 1.00(9)
 *   <li>minor release with bug fixes ( 1,  1,  0,  0) : 1.01
 * </ul>
 */
public class Version {


	/**
	 * Used to mark major releases (significant milestones)
	 */
	public static final int MajorVersion = 0;

	/**
	 * Used to mark minor releases (less-significant milestones)
	 */
	public static final int MinorVersion = 0;

	/**
	 * Used to track each change made to development versions
	 */
	public static final int EditLevel    = 0;

	/**
	 * Used to denote patched released versions (bug fixes)
	 */
	public static final int PatchLevel   = 0;

	/**
	 * Used to format the version number in a standard way.
	 *
	 * @return the formatted version number.
	 */
	public static String string()
	{
		return string(MajorVersion,MinorVersion,EditLevel,PatchLevel);
	}

	/**
	 * Used to format the version number in a standard way.
	 *
	 * @return the formatted version number.
	 */
	public static String string(int major, int minor, int edit, int patch)
	{
		// All version numbers start with major.minor format.
		String result = "";
		result += major + ".";
		if (minor < 10) result += "0";
		result += minor;

		// Pick the ending of according to the precise version.
		if (edit != 0)
			result += "(" + edit + ")";
		else if (patch != 0)
			result += "p" + patch;
		return result;
	}


}
