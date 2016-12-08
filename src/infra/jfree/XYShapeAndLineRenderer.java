package infra.jfree;

import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;

/**
 * This class is a simple extension of the JFreeChart XYLineAndShapeRenderer in order to
 * overcome one of the features/limitations of that class: the behaviour of drawing all lines
 * followed by all shapes, regardless of the order of the series themselves. The effect that this
 * has is to always draw shapes on top of lines on plots, which is hopeless if we're plotting a
 * regression line fit to dense data points (for example), as the line will be hidden by the
 * shapes used to mark the data.
 *
 *
 * @author nrowell
 * @version $Id$
 */
public class XYShapeAndLineRenderer extends XYLineAndShapeRenderer {
	
	/**
	 * The serial version UID
	 */
	private static final long serialVersionUID = -4217561880731969064L;

	/**
	 * {@inheritDoc}}
	 */
	protected boolean isLinePass(int pass) {
	    return pass == 1;
	}
	
	/**
	 * {@inheritDoc}}
	 */
	protected boolean isItemPass(int pass) {
	    return pass == 0;
	}
	
}
