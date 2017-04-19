package p2.datagenerator;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import p2.plot.ScatterPlot;

public class DataGenerator {
	
	static long SEED = 1234;
	static Random RANDOM = new Random(SEED);
	
	public static void main(String[] args) {
		
		int numDimensions = 2;
		int numCluster = 2;
		int numPointsPerCluster = 500;
		int numNoisePoints = 20;
	
		// Generate gaussian clusters
		
		List<double[]> gaussianDataPoints = new ArrayList<double[]>();
		
		for (int i = 0; i < numCluster; i++)
			gaussianDataPoints.addAll(DataGenerator.getGaussianCluster(numPointsPerCluster, numDimensions));
		
		for (int i = 0; i < numNoisePoints; i++) {
			
			double[] noisePoint = new double[numDimensions];
			for (int k = 0; k < numDimensions; k++) 
				noisePoint[k] = RANDOM.nextDouble() * 1.5 * (RANDOM.nextBoolean() ? 1 : -1);
			
			gaussianDataPoints.add(noisePoint);
		}
		
		ScatterPlot.plot(gaussianDataPoints, "plots/gauss.jpeg");  
		
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter("data/data_gauss_noisy.csv"));
			for (double[] row : gaussianDataPoints)
				writer.write(Arrays.toString(row).replace(" ", "").replace("[", "").replace("]", "\n"));
			writer.flush();
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	
		// Generate density-based clusters
		
		double epsilon = 1.0;
		double spacing = numPointsPerCluster * 0.05;  //Spacing between clusters
		
		List<double[]> densityDataPoints = new ArrayList<double[]>();
		
		double[] firstPoint = new double[numDimensions];
		Arrays.fill(firstPoint, 0);
		
		for (int i = 0; i < numCluster; i++) {
			
			densityDataPoints.addAll(DataGenerator.getDensityCluster(numPointsPerCluster, epsilon, numDimensions, firstPoint));
	
			//Ensure that clusters don't overlap 
			for (int j = 0; j < numDimensions; j++) 
				firstPoint[j] = firstPoint[j] + (RANDOM.nextDouble() * epsilon * spacing * (RANDOM.nextBoolean() ? 1 : -1));
		}
		
		for (int i = 0; i < numNoisePoints; i++) {
			
			double[] noisePoint = new double[numDimensions];
			for (int k = 0; k < numDimensions; k++) 
				noisePoint[k] = RANDOM.nextDouble() * spacing * (RANDOM.nextBoolean() ? 1 : -1);
			
			densityDataPoints.add(noisePoint);
		}
		
		ScatterPlot.plot(densityDataPoints, "plots/density.jpeg");
		
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter("data/data_density_noisy.csv"));
			for (double[] row : densityDataPoints)
				writer.write(Arrays.toString(row).replace(" ", "").replace("[", "").replace("]", "\n"));
			writer.flush();
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	
		//Merge data points together
		
		List<double[]> mergedDataPoints = DataGenerator.mergeDataPoints(numCluster * numPointsPerCluster, numDimensions, gaussianDataPoints, densityDataPoints);
	
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter("data/data_noisy.csv"));
			for (double[] row : mergedDataPoints)
				writer.write(Arrays.toString(row).replace(" ", "").replace("[", "").replace("]", "\n"));
			writer.flush();
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

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
