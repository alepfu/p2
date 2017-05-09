package p2.datagenerator;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class DataGenerator {
	
	static long seed = 1234;
	static Random random = new Random(seed);
	
	/**
	 * Parameters with default values
	 */
	static int numDimensions = 10;
	static int numClusters = 2;
	static int numPointsPerCluster = 1000;
	
	
	public static void main(String[] args) {
		
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();
		options.addOption("nd", "num-dimensions", true, "Number of dimensions.");
		options.addOption("nc", "num-clusters", true, "Number of clusters.");
		options.addOption("np", "num-points-per-cluster", true, "Number of points per cluster.");
		
		try {
			CommandLine line = parser.parse(options, args);
			
			if (line.hasOption("nd"))	
				numDimensions = Integer.parseInt(line.getOptionValue("nd"));
			
			if (line.hasOption("nc"))	
				numClusters = Integer.parseInt(line.getOptionValue("nc"));
			
			if (line.hasOption("np"))	
				numPointsPerCluster = Integer.parseInt(line.getOptionValue("np"));
			
		} catch (ParseException e) {
			System.out.println( "Error on parsing command line arguments.");
			e.printStackTrace();
		}
		
		
		/**
		 * Generate gaussian clusters
		 */
		List<double[]> gaussianDataPoints = generateGaussianClusters();
		//saveDataPointsToFile(gaussianDataPoints, "data/gauss.csv");
		//plot(gaussianDataPoints, "plots/gauss.jpeg", true);  
		
		
		/** 
		 * Generate density-based clusters
		 */
		List<double[]> densityDataPoints = generateDensityClusters();
		//saveDataPointsToFile(densityDataPoints, "data/density.csv");
		//plot(densityDataPoints, "plots/density.jpeg", true);
		
		
		/**
		 * Merge data points together
		 */
		List<double[]> mergedDataPoints = mergeDataPoints(gaussianDataPoints, densityDataPoints);
		saveDataPointsToFile(mergedDataPoints, "data/merged.csv");
		
	}

	private static List<double[]> generateGaussianClusters() {
		
		List<double[]> gaussianDataPoints = new ArrayList<double[]>();
		
		for (int i = 0; i < numClusters; i++)
			gaussianDataPoints.addAll(getGaussianCluster());
		
		return gaussianDataPoints;
	}

	private static List<double[]> getGaussianCluster() {
		
		List<double[]> cluster = new ArrayList<double[]>();
		
		double deviation = 0.1 ;				     		
		double[] mean = new double[numDimensions];
		
		for (int j = 0; j < numDimensions; j++) 
			mean[j] = random.nextDouble();
		
		for (int x = 0; x < numPointsPerCluster; x++) {

			double[] dataPoint = new double[numDimensions];
			for (int k = 0; k < numDimensions; k++) 
				dataPoint[k] = random.nextGaussian() * deviation + mean[k];
			
			cluster.add(dataPoint);
		}

		return cluster;
	}
	
	private static List<double[]> generateDensityClusters() {
		double epsilon = 1.0;
		double spacing = numPointsPerCluster * 0.05;  //Spacing between clusters
		
		List<double[]> densityDataPoints = new ArrayList<double[]>();
		
		double[] firstPoint = new double[numDimensions];
		Arrays.fill(firstPoint, 0);
		
		for (int i = 0; i < numClusters; i++) {
			
			densityDataPoints.addAll(getDensityCluster(epsilon, firstPoint));
	
			//Ensure that clusters don't overlap 
			for (int j = 0; j < numDimensions; j++) 
				firstPoint[j] = firstPoint[j] + (random.nextDouble() * epsilon * spacing * (random.nextBoolean() ? 1 : -1));
		}
		
		return densityDataPoints;
	}

	private static List<double[]> getDensityCluster(double epsilon, double[] firstPoint) {
	
		List<double[]> cluster = new ArrayList<double[]>();
	
		cluster.add(firstPoint.clone());
		
		for (int i = 0; i < (numPointsPerCluster - 1); i++) {
			
			double[] nextPoint = new double[numDimensions]; 
			
			for (int j = 0; j < numDimensions; j++)
				nextPoint[j] = firstPoint[j] + ((epsilon * (random.nextBoolean() ? 1 : -1)) * random.nextDouble());
			cluster.add(nextPoint.clone());
			
			firstPoint = nextPoint;
		}
		
		return cluster;
	}
	
	/**
	 * Merges 2 sets of data points together, appending a point and cluster id to every row.
	 * 
	 * @param gaussianDataPoints
	 * @param densityDataPoints
	 * @return List of double[] of the form "gaussian dimensions, density dimensions, point id, cluster id"
	 */
	private static List<double[]> mergeDataPoints(List<double[]> gaussianDataPoints, List<double[]> densityDataPoints) {
		
		int numPoints = numClusters * numPointsPerCluster;
		int clusterId = 0;
		
		List<double[]> dataPoints = new ArrayList<double[]>();
		
		for (int i = 0; i < numPoints; i++) {
			
			if ((i % numPointsPerCluster) == 0)
				clusterId++;
			
			double[] gaussianPoint = gaussianDataPoints.get(i);
			double[] densityPoint = densityDataPoints.get(i);
		
			double[] values = new double[(numDimensions * 2) + 2];
			
			for (int j = 0; j < numDimensions; j++) { 
				values[j] = gaussianPoint[j];
				values[j + numDimensions] = densityPoint[j];
			}

			//Append point and cluster ids
			values[values.length - 2] = i;
			values[values.length - 1] = clusterId;
			
			dataPoints.add(values);
		}
		
		return dataPoints;
	}
	
	private static void saveDataPointsToFile(List<double[]> dataPoints, String file) {
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(file));
			
			//Add header holding the set parameter values
			writer.write("#numClusters=" + numClusters + "\n" +
						 "#numDimensions=" + numDimensions + "\n" +
						 "#numPointsPerCluster=" + numPointsPerCluster + "\n");
			
			for (double[] row : dataPoints)
				writer.write(Arrays.toString(row).replace(" ", "")
						.replace("[", "").replace("]", "\n"));
			writer.flush();
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private static void plot(List<double[]> dataPoints, String filename, boolean show) {

		XYSeriesCollection dataset = new XYSeriesCollection();
		XYSeries series = new XYSeries("Cluster");

		for (double[] dataPoint : dataPoints)
			series.add(dataPoint[0], dataPoint[1]);

		dataset.addSeries(series);
		
		JFreeChart chart = ChartFactory.createScatterPlot("", "", "", dataset, 
				PlotOrientation.VERTICAL, false, false, false);
		
		if (show) {
			ChartFrame frame = new ChartFrame(filename, chart);
			frame.pack();
			frame.setVisible(true);
		}
		
		int width = 660;
		int height = 420;

		try {
			ChartUtilities.saveChartAsJPEG(new File(filename), chart, width, height);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
