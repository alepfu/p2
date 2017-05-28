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
	//static Random random = new Random(seed);
	static Random random = new Random();  //seedless
	
	/**
	 * Number of dimensions generated for each cluster type (gaussian and density).
	 */
	static int numDimensions = 2;
	
	/**
	 * Number of clusters to be generated.
	 */
	static int numClusters = 3;
	
	/**
	 * Number of generated points per cluster, meaningful values range from 100 to 1,000,000. 
	 */
	static int numPointsPerCluster = 1000;
	
	/**
	 * Flag indicating if cluster are overlapping or not.
	 */
	static boolean overlapping = false;
	
	/**
	 * Filename for exporting data points.
	 */
	static String filename = "data/datapoints.csv";
	
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
		
		//TODO add filename as parameter
		
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
		List<DataPoint> dataPoints = getDataPoints(gaussianDataPoints, densityDataPoints);
		
		//Export data points to file
		saveDataPointsToFile(dataPoints, filename);
		
		//Plotting
		try {
			plot(dataPoints, true);   //TODO in final version remove showing of plots
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
		double scaling = 100;
		double deviation = random.nextDouble() * scaling;
		
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
			deviation = random.nextDouble() * scaling;
			mean += (prevDeviation + deviation) * factor;
		}
		
		return points;
	}
	
	/**
	 * DRAUF GESCHISSEN
	 * 
	 * ICH MACH DAS NICHT MEHR GSCHEIT !!!!!
	 * 
	 * SOLL HEISSEN EIN RICHTIGES OVERLAPPING GIBTS NUR BEIM GAUS,
	 * BEIM DBSCAN LÖSE ICH DAS SPACING BRUTE FORCE UND PASSE DIE PARAMETER DEMENTSPRECHEND AN !!!!!!!
	 * 
	 * 
	 * @return
	 */
	private static List<double[]> getDensityClusters() {
		
		System.out.println("Generate density-based clusters ...");
		
		double epsilon = 1.0;		
		
		double spacing = epsilon * numPointsPerCluster * 0.05;

		if (overlapping)
			spacing = epsilon * numPointsPerCluster * 0.02;
		
		int numPoints = numClusters * numPointsPerCluster;
		
		List<double[]> points = new ArrayList<double[]>();
		
		for (int i = 0; i < numPoints; i++) {
			
			double[] point = new double[numDimensions]; 
			
			if (i == 0) {
				Arrays.fill(point, 0);
			} else {
				
				/*
				double x = random.nextDouble();
				int direction = -1;
				if (x > 0.25)
					direction = 1;
				*/

				int direction = random.nextBoolean() ? 1 : -1;
				
				for (int k = 0; k < numDimensions; k++)
					if ((i % numPointsPerCluster) != 0) 
						point[k] = points.get(points.size() - 1)[k] + (epsilon * direction * random.nextDouble());
					else
						point[k] = points.get(points.size() - 1)[k] + (epsilon * direction * random.nextDouble()) + spacing;
			}
			
			points.add(point);
		}
		
		return points;
	}
	
	/**
	 * Merges sets of gaussian and density data points together, adding point and cluster ids to every row.
	 * @param gaussianDataPoints
	 * @param densityDataPoints
	 * @return List of MergedDataPoint
	 */
	private static List<DataPoint> getDataPoints(List<double[]> gaussianDataPoints, List<double[]> densityDataPoints) {
		
		System.out.println("Merge data points together ...");
		
		List<DataPoint> dataPoints = new ArrayList<DataPoint>();
		int numPoints = numClusters * numPointsPerCluster;
		int clusterId = -1;
		
		for (int i = 0; i < numPoints; i++) {
			if ((i % numPointsPerCluster) == 0)
				clusterId++;
			
			DataPoint p = new DataPoint(gaussianDataPoints.get(i), densityDataPoints.get(i), i, clusterId);
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
	private static void saveDataPointsToFile(List<DataPoint> dataPoints, String file) {
		
		System.out.println("Export data points to file ...");
		
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(file));
			
			//Add header holding the set parameter values
			writer.write("#numClusters=" + numClusters + "\n" +
						 "#numDimensions=" + numDimensions + "\n" +
						 "#numPointsPerCluster=" + numPointsPerCluster + "\n" +
						 "#spacingType=" + overlapping + "\n");
			
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
	 * Generates and displays scatter plots for the gaussian and density clusters.
	 * @param mergedDataPoints The full set of generated data points with point and cluster ids.
	 * @param showPlot Wether the generated plot image should be shown or not.
	 */
	private static void plot(List<DataPoint> mergedDataPoints, boolean showPlot) throws Exception {
		
		System.out.println("Plotting gaussian clusters ...");
		String filename = "plots/gauss.jpeg";
		
		XYSeriesCollection datasetGauss = new XYSeriesCollection();
		for (int i = 0; i < numClusters; i++) {
			XYSeries series = new XYSeries("Cluster #" + i);
				for (DataPoint dataPoint : mergedDataPoints)
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
			for (DataPoint dataPoint : mergedDataPoints)
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
