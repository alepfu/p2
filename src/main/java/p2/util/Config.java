package p2.util;

import java.awt.Color;
import java.awt.Paint;

import org.jfree.chart.plot.DefaultDrawingSupplier;

public class Config {

	public static int numClusters = 2;
	
	public static int numPointsCluster = 1000;
	
	public static int numNoisePoints = 0;
	
	public static int numPoints = (numClusters * numPointsCluster) + numNoisePoints;
	
	public static int numDimPerType = 2;
	
	public static long seed = 7;
	
	public static String workDir = "/home/alepfu/Desktop/P2";
	
	public static boolean verfiyDataset = false;
	
	public static boolean plotClusters = true;
	
	public static boolean plotClustersHighlighted = false;
	
	public static boolean displayPlots = true;
	
	
	
	
	public static DefaultDrawingSupplier customDrawingSupplier = new DefaultDrawingSupplier(
			new Paint[] { 
					new Color(31,120,180),
					new Color(227,26,28),
					new Color(51,160,44),
					new Color(106,61,154),
					new Color(177,89,40),
					new Color(255,127,0),
					new Color(178,223,138),
					new Color(251,154,153),
					new Color(166,206,227),
					new Color(253,191,111),
					new Color(202,178,214),
					new Color(255,255,153),
			},
			DefaultDrawingSupplier.DEFAULT_OUTLINE_PAINT_SEQUENCE,
			DefaultDrawingSupplier.DEFAULT_STROKE_SEQUENCE,
			DefaultDrawingSupplier.DEFAULT_OUTLINE_STROKE_SEQUENCE,
			DefaultDrawingSupplier.DEFAULT_SHAPE_SEQUENCE);
}
