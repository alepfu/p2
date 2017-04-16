package p2.datagenerator;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import p2.plot.ScatterPlot;

public class RunDataGenerator { 

	static long SEED = System.currentTimeMillis();
	static Random RANDOM = new Random(SEED);
	
	public static void main(String[] args) {
		
		int numDimensions = 2;
		int numCluster = 2;
		int numPointsPerCluster = 500;
		int numNoisePoints = 0;
		
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
		
		ScatterPlot.show(gaussianDataPoints);
		
		
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
		
		ScatterPlot.show(densityDataPoints);
		
		
		//Merge data points together and save to file
		
		List<double[]> mergedDataPoints = DataGenerator.mergeDataPoints(numCluster * numPointsPerCluster, numDimensions, gaussianDataPoints, densityDataPoints);

		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter("/home/alepfu/Desktop/data.txt"));
			for (double[] row : mergedDataPoints)
				writer.write(Arrays.toString(row).replace(" ", "").replace("[", "").replace("]", "\n"));
			writer.flush();
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
