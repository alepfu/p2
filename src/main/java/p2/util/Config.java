package p2.util;

public class Config {

	public static int numClusters = 2;
	
	public static int numPointsCluster = 1000;
	
	public static int numNoisePoints = 0;
	
	public static int numPoints = (numClusters * numPointsCluster) + numNoisePoints;
	
	public static int numDimPerType = 2;
	
	public static long seed = 7;
	
	public static String workDir = "/home/alepfu/Desktop/P2";
	
}
