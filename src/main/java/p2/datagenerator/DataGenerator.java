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
	
	/**
	 * Seed for random number generation
	 */
	static long seed = 12345;
	
	/**
	 * Random number generator.
	 */
	static Random random = new Random(seed);
	
	/**
	 * Number of dimensions generated for each cluster type (gaussian and density).
	 */
	static int numDimensions = 2;
	
	/**
	 * Number of clusters to be generated.
	 */
	static int numClusters = 5;
	
	/**
	 * Number of generated points per cluster, meaningful values range from 100 to 1,000,000. 
	 */
	static int numPointsPerCluster = 10000;
	
	/**
	 * Flag indicating if cluster are overlapping or not.
	 */
	static boolean overlapping = false;
	
	/**
	 * Main method making use of command line arguments.
	 * Generates data points and exports them to a file.
	 * Generates 2D scatter plots for each cluster type.
	 * @param args
	 */
	public static void main(String[] args) {
		
		//Initialze argument parser
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();
		options.addOption("d", "num-dimensions", true, "Number of dimensions.");
		options.addOption("c", "num-clusters", true, "Number of clusters.");
		options.addOption("p", "num-points-per-cluster", true, "Number of points per cluster.");
		options.addOption("o", "overlap", false, "Overlapping clusters");
		options.addOption("s", "seed", true, "Seed for random number generation.");
		
		//Parse arguments
		try {
			CommandLine line = parser.parse(options, args);
			
			if (line.hasOption("d"))	
				numDimensions = Integer.parseInt(line.getOptionValue("d"));
			
			if (line.hasOption("c"))	
				numClusters = Integer.parseInt(line.getOptionValue("c"));
			
			if (line.hasOption("p"))	
				numPointsPerCluster = Integer.parseInt(line.getOptionValue("p"));
			
			if (line.hasOption("o"))				
				overlapping = true;
			
			if (line.hasOption("s"))	
				seed = Integer.parseInt(line.getOptionValue("s"));
			
		} catch (ParseException e) {
			System.out.println( "Error on parsing command line arguments.");
			e.printStackTrace();
		}
		
		//Generate gaussian data points
		List<double[]> gaussianDataPoints = getGaussianClusters(); 
		
		//Generate densitiy data points
		List<double[]> densityDataPoints = getDensityClusters();
		
		//Merge data points together
		List<MergedDataPoint> mergedDataPoints = mergeDataPoints(gaussianDataPoints, densityDataPoints);
		
		//Export data points to file
		saveDataPointsToFile(mergedDataPoints, "data/merged.csv");
		
		//Plotting
		try {
			plot(mergedDataPoints, true);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		System.out.println("Finished.");
	}
	
	/**
	 * Generates the data points for the gaussian clusters.
	 * Standard deviations of clusters are choosen randomly.
	 * The 1st cluster is positioned at the origin.
	 * Separation of clusters is handeld by a multiple of the standard deviation (68–95–99.7 rule). 
	 * @return A list of data points, where a data point is a double[] with size numDimensions.
	 */
	private static List<double[]> getGaussianClusters() {
		
		System.out.println("Generate gaussian clusters ...");
		
		List<double[]> points = new ArrayList<double[]>(numClusters*numPointsPerCluster);
		
		double mean = 0;  //The 1st cluster is positioned at the origin
		double deviation = random.nextDouble();
		
		int numPoints = numClusters * numPointsPerCluster;
		
		//Get the needed factor for the separation of clusters
		double factor = 3;
		if (numPoints >= 1000000)
			factor = 5;
		if ((1000000 > numPoints) && (numPoints >= 100000))
			factor = 4.5;
		if ((100000 > numPoints) && (numPoints >= 10000))
			factor = 4;
		if ((10000 > numPoints) && (numPoints >= 1000))
			factor = 3.5;
		
		//For overlapping clusters we just cut factor in half
		if (overlapping)
			factor /= 2;

		//Generate data points
		for (int i = 0; i < numClusters; i++) {
			
			for (int j = 0; j < numPointsPerCluster; j++) {
			
				double[] point = new double[numDimensions];
				
				for (int k = 0; k < numDimensions; k++) 
					point[k] = random.nextGaussian() * deviation + mean;
				
				points.add(point);
			}
			
			//Set deviation and mean for the next cluster
			double prevDeviation = deviation;
			deviation = random.nextDouble();
			mean += (prevDeviation + deviation) * factor;
		}
		
		return points;
	}
	
	/**
	 * NEW METHOD
	 * 
	 * 
	 * 
	 * 
	 * 
	 * 
	 * 
	 * @return
	 */
	private static List<double[]> getDensityClusters() {
		
		System.out.println("Generate density clusters ...");
		
		double epsilon = 1.0;
		double space = 3 * epsilon;
		
		
		List<double[]> points = new ArrayList<double[]>(numClusters*numPointsPerCluster);
		
		double[] firstPoint = new double[numDimensions];
		Arrays.fill(firstPoint, 0);
		
		for (int i = 0; i < numClusters; i++) {
			
			points.add(firstPoint);
			
			for (int j = 1; j < numPointsPerCluster; j++) {

				
				//TODO wir müssen uns von Punkt zu Punkt handeln und nicht immer nur vom ersten weg!
				
				double[] point = new double[numDimensions];
				double val = epsilon * random.nextDouble();
				Arrays.fill(point, val);
				point[0] = (space * i) + val + random.nextDouble();
				
				points.add(point);
			}
		
			firstPoint = new double[numDimensions];
			Arrays.fill(firstPoint, 0);
			firstPoint[0] = space * (i + 1);
		}
		
		
		
		return points;
	}
	
	
	
	
	
	
	private static List<double[]> generateDensityClusters() {
		System.out.println("Generate density-based clusters ...");
		
		double epsilon = 1.0;		
		
		double spacing = numPointsPerCluster * epsilon;
		
		if (overlapping)
			spacing = numPointsPerCluster * epsilon * 0.035;  //TODO taylored to seed 1234567 with 2 dimensions and 3 clusters with each 1000 points
		
		List<double[]> densityDataPoints = new ArrayList<double[]>();
		
		double[] firstPoint = new double[numDimensions];
		Arrays.fill(firstPoint, 0);
		
		for (int i = 0; i < numClusters; i++) {
			
			densityDataPoints.addAll(getDensityCluster(epsilon, firstPoint));
	
			//Apply spacing
			for (int j = 0; j < numDimensions; j++) 
				firstPoint[j] = firstPoint[j] + (random.nextDouble() * epsilon * spacing);
		}
		
		return densityDataPoints;
	}

	private static List<double[]> getDensityCluster(double epsilon, double[] firstPoint) {
	
		List<double[]> cluster = new ArrayList<double[]>();
	
		cluster.add(firstPoint.clone());
		
		for (int i = 0; i < (numPointsPerCluster - 1); i++) {
			
			double[] nextPoint = new double[numDimensions]; 
			
			int direction = random.nextBoolean() ? 1 : -1;
			
			for (int j = 0; j < numDimensions; j++)
				nextPoint[j] = firstPoint[j] + ((epsilon * direction) * random.nextDouble());
				
			cluster.add(nextPoint.clone());
			
			firstPoint = nextPoint;
		}
		
		return cluster;
	}
	
	/**
	 * Merges sets of gaussian and density data points together, adding point and cluster ids to every row.
	 * @param gaussianDataPoints
	 * @param densityDataPoints
	 * @return List of MergedDataPoint
	 */
	private static List<MergedDataPoint> mergeDataPoints(List<double[]> gaussianDataPoints, List<double[]> densityDataPoints) {
		
		System.out.println("Merge data points together ...");
		
		List<MergedDataPoint> dataPoints = new ArrayList<MergedDataPoint>();
		int numPoints = numClusters * numPointsPerCluster;
		int clusterId = -1;
		
		for (int i = 0; i < numPoints; i++) {
			if ((i % numPointsPerCluster) == 0)
				clusterId++;
			
			MergedDataPoint p = new MergedDataPoint(gaussianDataPoints.get(i), densityDataPoints.get(i), i, clusterId);
			dataPoints.add(p);
		}
		
		return dataPoints;
	}
	
	/**
	 * Exports data points to file.
	 * Adds header information.
	 * @param dataPoints A list of data points.
	 * @param file A filename.
	 */
	private static void saveDataPointsToFile(List<MergedDataPoint> dataPoints, String file) {
		
		System.out.println("Export data points to file ...");
		
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(file));
			
			//Add header holding the set parameter values
			writer.write("#numClusters=" + numClusters + "\n" +
						 "#numDimensions=" + numDimensions + "\n" +
						 "#numPointsPerCluster=" + numPointsPerCluster + "\n" +
						 "#spacingType=" + overlapping + "\n");
			
			for (MergedDataPoint p : dataPoints) {
				writer.write(p.toString());
				writer.newLine();
			}
			writer.flush();
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Generates and displays scatter plots for the gaussian and density clusters.
	 * @param mergedDataPoints The full set of generated data points with point and cluster ids.
	 * @param showPlot Wether the generated plot image should be shown or not.
	 */
	private static void plot(List<MergedDataPoint> mergedDataPoints, boolean showPlot) throws Exception {
		
		System.out.println("Plotting gaussian clusters ...");
		String filename = "plots/gauss.jpeg";
		
		XYSeriesCollection datasetGauss = new XYSeriesCollection();
		for (int i = 0; i < numClusters; i++) {
			XYSeries series = new XYSeries("Cluster #" + i);
				for (MergedDataPoint dataPoint : mergedDataPoints)
					if (dataPoint.getClusterId() == i)
						series.add(dataPoint.getGaussFeatures()[0], dataPoint.getGaussFeatures()[1]);
				datasetGauss.addSeries(series);
		}
		
		JFreeChart chartGauss = ChartFactory.createScatterPlot("Gauss clusters", "x1", "x2", datasetGauss, 
				PlotOrientation.VERTICAL, false, false, false);
		
		if (showPlot) {
			ChartFrame frame = new ChartFrame(filename, chartGauss);
			frame.pack();
			frame.setVisible(true);
		}
		
		ChartUtilities.saveChartAsJPEG(new File(filename), chartGauss, 660, 420);
		
		System.out.println("Plotting density clusters ...");
		filename = "plots/density.jpeg";
		
		XYSeriesCollection datasetDensity = new XYSeriesCollection();
		for (int i = 0; i < numClusters; i++) {
			XYSeries series = new XYSeries("Cluster #" + i);
			for (MergedDataPoint dataPoint : mergedDataPoints)
				if (dataPoint.getClusterId() == i)
					series.add(dataPoint.getDensityFeatures()[0], dataPoint.getDensityFeatures()[1]);
				datasetDensity.addSeries(series);
		}
		
		JFreeChart chartDensity = ChartFactory.createScatterPlot("Density clusters", "x1", "x2", datasetDensity, 
				PlotOrientation.VERTICAL, false, false, false);
		
		if (showPlot) {
			ChartFrame frame = new ChartFrame(filename, chartDensity);
			frame.pack();
			frame.setVisible(true);
		}
		
		ChartUtilities.saveChartAsJPEG(new File(filename), chartGauss, 660, 420);
	}
}
