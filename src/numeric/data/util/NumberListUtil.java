package numeric.data.util;

import java.awt.image.BufferedImage;
import java.io.IOException;

import infra.io.Gnuplot;
import infra.os.OSChecker;
import numeric.data.FloatList;

public class NumberListUtil {
	
	/**
	 * Produces a histogram plot of the {@link FloatList} contents in the output location specified.
	 * @param data
	 * 	The {@link FloatList} to plot.
	 * @param output
	 * 	The output file
	 * @param xAxisLabel
	 * 	The label to use for the X axis
	 * @throws IOException
	 * 	If there's a problem writing the plot file
	 */
	public static BufferedImage plotHistogram(FloatList data, String xAxisLabel) throws IOException {

		// Get the median and MAD values
		float[] stats = data.getMAD();
		float median = stats[0];
		float mad = stats[1];
		
		// Initial axis range from spread of data
		float xmin = median - 4.0f * mad;
		float xmax = median + 4.0f * mad;
		
		// Bin width for histogram -> fit 50 bins in the core of the distribution
		float binWidth = (xmax - xmin)/50.0f;
		
		// Figure out which side to put the label on...
		boolean placeLabelOnLeft = true;
		if(Math.abs(xmin) > Math.abs(xmax)) {
			// Symmetric range is obtained by adding blank space on the right.
			// Put the label on this side.
			placeLabelOnLeft = false;
		}
		
		// ...now enforce symmetric axis range
		float max = Math.max(Math.abs(xmax), Math.abs(xmin));
		xmin = -1.0f * max;
		xmax = max;
		
		StringBuilder script = new StringBuilder();
		
		script.append("set terminal pngcairo font 'Helvetica' enhanced").append(OSChecker.newline);
		script.append("set key off").append(OSChecker.newline);
		
		// Set x axis
		script.append("set xrange ["+xmin+":"+xmax+"]").append(OSChecker.newline);
		script.append("set xtics font ',12' out nomirror offset 0.0,0.0").append(OSChecker.newline);
		script.append("set mxtics 5").append(OSChecker.newline);
		script.append("set xlabel font ',16' '"+xAxisLabel+"' offset 0.0,0.0").append(OSChecker.newline);
		
		// Set y axis
		script.append("set yrange [0:*]").append(OSChecker.newline);
		script.append("set ytics font ',12' out").append(OSChecker.newline);
		script.append("set mytics 5").append(OSChecker.newline);
		script.append("set ylabel font ',16' 'N' offset 0.0,0.0").append(OSChecker.newline);
		
		// Histogram parameters
		script.append("binwidth="+binWidth).append(OSChecker.newline);
		script.append("bin(x,width)=width*floor(x/width)").append(OSChecker.newline);
		
		// Grid lines
		script.append("set grid").append(OSChecker.newline);
		
		// Label stating the median & MAD
		script.append("median="+median).append(OSChecker.newline);
		script.append("mad="+mad).append(OSChecker.newline);
		script.append("LABEL = \""+String.format("Median = %8.6f\\n    MAD = %8.6f", median, mad)+"\"").append(OSChecker.newline);
		
		if(placeLabelOnLeft) {
			script.append("set obj 10 rect at screen 0.33,0.86 size char 22, char 3").append(OSChecker.newline);
			script.append("set obj 10 fillstyle solid border -1 front").append(OSChecker.newline);
			script.append("set label 10 at screen 0.33,0.9 LABEL front center font ',16'").append(OSChecker.newline);
		}
		else {
			script.append("set obj 10 rect at screen 0.75,0.86 size char 22, char 3").append(OSChecker.newline);
			script.append("set obj 10 fillstyle solid border -1 front").append(OSChecker.newline);
			script.append("set label 10 at screen 0.75,0.9 LABEL front center font ',16'").append(OSChecker.newline);
		}
		
		// Lines marking median & MAD
		script.append("set style arrow 4 nohead ls 1 lc rgb 'blue'").append(OSChecker.newline);
		script.append("set style arrow 5 nohead ls 1 lc rgb 'red'").append(OSChecker.newline);
		script.append("set arrow from median, graph 0.0 to median, graph 1.0 as 4").append(OSChecker.newline);
		script.append("set arrow from (median-mad), graph 0.0 to (median-mad), graph 1.0 as 5").append(OSChecker.newline);
		script.append("set arrow from (median+mad), graph 0.0 to (median+mad), graph 1.0 as 5").append(OSChecker.newline);
		
		script.append("plot '-' using (bin($1,binwidth)):(1.0) smooth freq with l lw 2 lc rgb 'dark-green' notitle").append(OSChecker.newline);

		for(float dAl : data) {
			script.append(dAl).append(OSChecker.newline);
		}
		script.append("e").append(OSChecker.newline);
		
		BufferedImage plot = Gnuplot.executeScript(script.toString());
		
		return plot;
	}
}
