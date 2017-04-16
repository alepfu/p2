package p2.datagenerator;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class DataGenerator {
	
	static long SEED = System.currentTimeMillis();
	static Random RANDOM = new Random(SEED);
	
	public static List<double[]> getGaussianCluster(int numPointsPerCluster, int numDimensions) {
		
		List<double[]> cluster = new ArrayList<double[]>();
		
		double deviation = 0.1 ;				     		
		double[] mean = new double[numDimensions];
		
		for (int j = 0; j < numDimensions; j++) 
			mean[j] = RANDOM.nextDouble();
		
		for (int x = 0; x < numPointsPerCluster; x++) {

			double[] dataPoint = new double[numDimensions];
			for (int k = 0; k < numDimensions; k++) 
				dataPoint[k] = RANDOM.nextGaussian() * deviation + mean[k];
			
			cluster.add(dataPoint);
		}

		return cluster;
	}
	
	public static List<double[]> getDensityCluster(int numPoints, double epsilon, int numDimensions, double[] firstPoint) {
	
		List<double[]> cluster = new ArrayList<double[]>();
	
		cluster.add(firstPoint.clone());
		
		for (int i = 0; i < (numPoints - 1); i++) {
			
			double[] nextPoint = new double[numDimensions]; 
			
			for (int j = 0; j < numDimensions; j++)
				nextPoint[j] = firstPoint[j] + ((epsilon * (RANDOM.nextBoolean() ? 1 : -1)) * RANDOM.nextDouble());
			cluster.add(nextPoint.clone());
			
			firstPoint = nextPoint;
		}
		
		return cluster;
	}
	
	public static List<double[]> mergeDataPoints(int numPoints, int numDimensions, List<double[]> gaussianDataPoints, List<double[]> densityDataPoints) {
		
		List<double[]> dataPoints = new ArrayList<double[]>();
		
		for (int i = 0; i < numPoints; i++) {
			
			double[] gaussianPoint = gaussianDataPoints.get(i);
			double[] densityPoint = densityDataPoints.get(i);
		
			double[] values = new double[numDimensions * 2];
			
			for (int j = 0; j < numDimensions; j++) { 
				values[j] = gaussianPoint[j];
				values[j + numDimensions] = densityPoint[j];
			}
			
			dataPoints.add(values);
		}
		
		return dataPoints;
	}
		
	
}
