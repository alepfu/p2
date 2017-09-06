package p2.util;

/**
 * Copyright (c) 2017 Alexander Pfundner
 * 
 * Integration of Density-based and Partitioning-based Clustering Methods
 * 
 */

import java.awt.Color;
import java.awt.Paint;

import org.jfree.chart.plot.DefaultDrawingSupplier;

/**
 * Holds parameters used by the classes DatasetGenerator and IntegratedRunner.
 *
 */
public class Config {

	/**
	 * Number of clusters to be generated.
	 */
	public static int numClusters = 2;
	
	/**
	 * Number of generated points per cluster. 
	 */
	public static int numPointsCluster = 1000;
	
	/**
	 * Number of noise points.
	 */
	public static int numNoisePoints = 100;
	
	/**
	 * Total number of points.
	 */
	public static int numPoints = (numClusters * numPointsCluster) + numNoisePoints;
	
	/**
	 * Number of features each cluster type has.
	 */
	public static int numDimPerType = 2;
	
	/**
	 * Seed for the random number generator.
	 */
	public static long seed = 7;
	
	/**
	 * Working directory, without trailing slash.
	 */
	public static String workDir = "/home/alepfu/Desktop/P2";
	
	/**
	 * Flag that enables the generation of verfication plots.
	 */
	public static boolean verfiyDataset = false;

	/**
	 * Flag that enables the generation of cluster plots.
	 */
	public static boolean plotClusters = true;
	
	/**
	 * Flag that enables the generation of plots that highlight the points that are easy or hard to cluster.
	 */
	public static boolean plotClustersHighlighted = false;
	
	/**
	 * Flag that enables the showing of plots directly.
	 */
	public static boolean displayPlots = false;
	
}
