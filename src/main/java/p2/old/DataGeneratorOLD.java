package p2.old;
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

import p2.datagenerator.DataPoint;

public class DataGeneratorOLD {
	
	/**
	 * Seed for random number generation
	 */
	private static long seed = 7;
	
	/**
	 * Random number generator.
	 */
	private static Random rand = new Random(seed);
	
	/**
	 * Number of dimensions generated for each cluster type (gaussian and density).
	 */
	private static int numDimensions = 5;
	
	/**
	 * Number of clusters to be generated.
	 */
	private static int numClusters = 3;
	
	/**
	 * Number of generated points per cluster, meaningful values range from 100 to 1,000,000. 
	 */
	private static int numPointsCluster = 1000;
	
	/**
	 * Number of points in total.
	 */
	private static int numPoints = numClusters * numPointsCluster;
	
	/**
	 * Flag indicating if cluster are overlapping or not.
	 */
	private static boolean overlapping = true;
	
	/**
	 * Directory for exported files.
	 */
	private static String exportDir = "/home/alepfu/Desktop/P2_export";
	//private static String exportDir = "";
	
	/**
	 * Display plots.
	 */
	private static boolean showPlots = false;
	
	/**
	 * Main method using command line arguments.
	 * Generates data points and exports them to a file.
	 * Generates 2D scatter plots for each cluster type.
	 * @param args
	 */
	public static void main(String[] args) {
		
		//Long timeStart = System.currentTimeMillis();
		
		//Initialze argument parser
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();
		options.addOption("nd", "num-dimensions", true, "Number of dimensions.");
		options.addOption("nc", "num-clusters", true, "Number of clusters.");
		options.addOption("np", "num-points-cluster", true, "Number of points per cluster.");
		options.addOption("o", "overlap", false, "Overlapping clusters");
		options.addOption("s", "seed", true, "Seed for random number generation.");
		options.addOption("ed", "export-directory", true, "Directory for saving generated files.");
		options.addOption("p", "show-plots", false, "Display plots.");
		
		//Parse arguments
		try {
			CommandLine line = parser.parse(options, args);
			
			if (line.hasOption("nd"))	
				numDimensions = Integer.parseInt(line.getOptionValue("nd"));
			
			if (line.hasOption("nc"))	
				numClusters = Integer.parseInt(line.getOptionValue("nc"));
			
			if (line.hasOption("np"))	
				numPointsCluster = Integer.parseInt(line.getOptionValue("np"));
			
			numPoints = numClusters * numPointsCluster;
			
			if (line.hasOption("o"))				
				overlapping = true;
			
			if (line.hasOption("s")) {
				seed = Integer.parseInt(line.getOptionValue("s"));
				rand = new Random(seed);
			}
			
			if (line.hasOption("ed")) {
				exportDir = line.getOptionValue("ed");
			}
			
			if (line.hasOption("p"))				
				showPlots = true;
			
		} catch (ParseException e) {
			System.out.println( "Error on parsing command line arguments.");
			e.printStackTrace();
		}
		
		//Generate gaussian data points
		List<double[]> gaussianDataPoints = getGaussianClusters(); 
		
		//Generate densitiy data points
		List<double[]> densityDataPoints = getDensityClusters();
		
		//Merge data points together
		List<DataPoint> dataPoints = getDataPoints(gaussianDataPoints, densityDataPoints);
		
		//Export data points to file
		saveDataPointsToFile(dataPoints, exportDir + "/data_" + getFileName() + ".csv");
		
		//Plotting
		try {

			plotGaussianClusters(dataPoints, exportDir + "/gauss_" + getFileName() + ".jpeg");
			plotDensityClusters(dataPoints, exportDir + "/density_" + getFileName() + ".jpeg");
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		System.out.println("Finished.");
	}
	
	private static String getFileName() {
		
		return numDimensions + "d_" + numClusters + "c_" + numPointsCluster + "p_" + (overlapping ? "overlap_" : "separated_") + seed + "seed"; 
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
		
		List<double[]> points = new ArrayList<double[]>();
		
		double mean = 0;  //The 1st cluster is positioned at the origin
		double scaling = 100;
		double deviation = rand.nextDouble() * scaling;
		
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
		
		//For overlapping clusters we just cut factor
		if (overlapping)
			factor /= 3;

		//Generate data points
		for (int i = 0; i < numClusters; i++) {
			
			for (int j = 0; j < numPointsCluster; j++) {
			
				double[] point = new double[numDimensions];
				
				for (int k = 0; k < numDimensions; k++) 
					point[k] = rand.nextGaussian() * deviation + mean;
				
				points.add(point);
			}
			
			//Set deviation and mean for the next cluster
			double prevDeviation = deviation;
			deviation = rand.nextDouble() * scaling;
			mean += (prevDeviation + deviation) * factor;
		}
		
		return points;
	}
	
	/**
	 * Generates data points for density-connected clusters.
	 * Generates some core points with a random walk, then sets remaining points arround the core points.
	 * The first random walk starts at the origin.
	 * Separation of clusters is handeld by a fixed spacing. 
	 * @return
	 */
	private static List<double[]> getDensityClusters() {
		
		System.out.println("Generate density-based clusters ...");
		
		double epsilon = 1.0;		
		double spacing = epsilon * numPointsCluster * 0.03;
		
		if (overlapping) {
			spacing = epsilon * numPointsCluster * 0.015;
		
			//Empirical stuff
			if (numClusters == 2)
				spacing = epsilon * numPointsCluster * 0.005;
			
		}
		
		int numCorePoints = numPoints / 10;
		int numPointsCorePoint = (numPoints - numCorePoints) / numCorePoints;
		int numCorePointsPerCluster = numCorePoints / numClusters;
		
		List<double[]> points = new ArrayList<double[]>();
		
		for (int i = 0; i < numCorePoints; i++) {
			
			//Generate core point
			double[] corePoint = new double[numDimensions]; 
			
			if (i == 0) {
				Arrays.fill(corePoint, 0);
			} else {

				int direction = rand.nextBoolean() ? 1 : -1;
				
				for (int k = 0; k < numDimensions; k++)
					if ((i % numCorePointsPerCluster) != 0) 
						corePoint[k] = points.get(points.size() - 1)[k] + (epsilon * direction * rand.nextDouble());
					else
						corePoint[k] = points.get(points.size() - 1)[k] + (epsilon * direction * rand.nextDouble()) + spacing;
			}
			
			points.add(corePoint);
			
			//Add points arround core point
			for (int j = 0; j < numPointsCorePoint; j++) {
				
				double[] point = new double[numDimensions];
				
				for (int k = 0; k < numDimensions; k++) 					
					point[k] = corePoint[k] + (epsilon * rand.nextDouble() * (rand.nextBoolean() ? 1 : -1));
					
				points.add(point);
			}
			
		}
		
		return points;
	}
	
	/**
	 * Merges sets of gaussian and density data points together, adding cluster labels.
	 * @param gaussianDataPoints
	 * @param densityDataPoints
	 * @return List of MergedDataPoint
	 */
	private static List<DataPoint> getDataPoints(List<double[]> gaussianDataPoints, List<double[]> densityDataPoints) {
		
		System.out.println("Merge data points together ...");
		
		List<DataPoint> dataPoints = new ArrayList<DataPoint>();
		int numPoints = numClusters * numPointsCluster;
		int clusterId = -1;
		
		for (int i = 0; i < numPoints; i++) {
			if ((i % numPointsCluster) == 0)
				clusterId++;
			
			DataPoint p = new DataPoint(gaussianDataPoints.get(i), densityDataPoints.get(i), "c" + clusterId);
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
	private static void saveDataPointsToFile(List<DataPoint> dataPoints, String filename) {
		
		System.out.println("Export data points to file ...");
		
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(filename));
			
			//Add header holding the set parameter values
			writer.write("#numClusters=" + numClusters + "\n" +
						 "#numDimensions=" + numDimensions + "\n" +
						 "#numPointsCluster=" + numPointsCluster + "\n" +
						 "#overlapping=" + overlapping + "\n" +
						 "#seed=" + seed + "\n");
			
			for (DataPoint p : dataPoints) {
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
	 * Generates and displays scatter plots for the gaussian clusters.
	 * @param mergedDataPoints The full set of generated data points with point and cluster ids.
	 * @param showPlot Wether the generated plot image should be shown or not.
	 */
	private static void plotGaussianClusters(List<DataPoint> mergedDataPoints, String filename) throws Exception {
		
		System.out.println("Plotting gaussian clusters ...");
		
		XYSeriesCollection datasetGauss = new XYSeriesCollection();
		for (int i = 0; i < numClusters; i++) {
			XYSeries series = new XYSeries("Cluster #" + i);
				for (DataPoint dataPoint : mergedDataPoints)
					if (dataPoint.getClusterLabel().equals("c" + i))
						series.add(dataPoint.getGaussFeatures()[0], dataPoint.getGaussFeatures()[1]);
				datasetGauss.addSeries(series);
		}
		
		JFreeChart chartGauss = ChartFactory.createScatterPlot("Gauss clusters", "x1", "x2", datasetGauss, 
				PlotOrientation.VERTICAL, false, false, false);
		
		//Display plot
		if (showPlots) {
			ChartFrame frame = new ChartFrame(filename, chartGauss);
			frame.pack();
			frame.setVisible(true);
		}
		
		//Export plot to file
		ChartUtilities.saveChartAsJPEG(new File(filename), chartGauss, 660, 420);
	}
	
	/**
	 * Generates and displays scatter plots for the density clusters.
	 * @param mergedDataPoints The full set of generated data points with point and cluster ids.
	 * @param showPlot Wether the generated plot image should be shown or not.
	 */
	private static void plotDensityClusters(List<DataPoint> mergedDataPoints, String filename) throws Exception {
		
		System.out.println("Plotting density clusters ...");
		
		XYSeriesCollection datasetDensity = new XYSeriesCollection();
		for (int i = 0; i < numClusters; i++) {
			XYSeries series = new XYSeries("Cluster #" + i);
			for (DataPoint dataPoint : mergedDataPoints)
				if (dataPoint.getClusterLabel().equals("c" + i))
					series.add(dataPoint.getDensityFeatures()[0], dataPoint.getDensityFeatures()[1]);
				datasetDensity.addSeries(series);
		}
		
		JFreeChart chartDensity = ChartFactory.createScatterPlot("Density clusters", "x1", "x2", datasetDensity, 
				PlotOrientation.VERTICAL, false, false, false);
		
		//Display plot
		if (showPlots) {
			ChartFrame frame = new ChartFrame(filename, chartDensity);
			frame.pack();
			frame.setVisible(true);
		}
		
		//Export plot to file
		ChartUtilities.saveChartAsJPEG(new File(filename), chartDensity, 660, 420);
	}
}
