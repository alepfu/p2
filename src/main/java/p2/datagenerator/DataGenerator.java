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
	
	static long seed = 73;
	static Random random = new Random(seed);
	
	/**
	 * Parameters
	 */
	static int numDimensions = 2;
	static int numClusters = 2;
	static int numPointsPerCluster = 1000;
	static String spacingType = "o";
	static boolean isInconsistent = false;
	
	public static void main(String[] args) {
		
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();
		options.addOption("d", "num-dimensions", true, "Number of dimensions.");
		options.addOption("c", "num-clusters", true, "Number of clusters.");
		options.addOption("p", "num-points-per-cluster", true, "Number of points per cluster.");
		options.addOption("s", "spacing", true, "[s] separated, [o] overlapping");
		options.addOption("i", "inconsistent", false, "Generates an inconsitent data set");
		
		try {
			CommandLine line = parser.parse(options, args);
			
			if (line.hasOption("d"))	
				numDimensions = Integer.parseInt(line.getOptionValue("d"));
			
			if (line.hasOption("c"))	
				numClusters = Integer.parseInt(line.getOptionValue("c"));
			
			if (line.hasOption("p"))	
				numPointsPerCluster = Integer.parseInt(line.getOptionValue("p"));
			
			if (line.hasOption("s")) {				
				spacingType = line.getOptionValue("s");
				
				if (!(spacingType.equals("s") || spacingType.equals("o") || spacingType.equals("i")))
					throw new ParseException("Spacing type must be one of s|o|i.");
			}
			
		} catch (ParseException e) {
			System.out.println( "Error on parsing command line arguments.");
			e.printStackTrace();
		}
		
		System.out.println("Generate gaussian clusters ...");
		List<double[]> gaussianDataPoints = generateGaussianClusters();
		
		System.out.println("Generate density-based clusters ...");
		List<double[]> densityDataPoints = generateDensityClusters();
		
		System.out.println("Merge data points together ...");
		List<double[]> mergedDataPoints = mergeDataPoints(gaussianDataPoints, densityDataPoints);
		saveDataPointsToFile(mergedDataPoints, "data/merged.csv");
		
		
		try {
			plot(mergedDataPoints, true);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		System.out.println("Finished.");
	}

	private static List<double[]> generateGaussianClusters() {
		
		List<double[]> gaussianDataPoints = new ArrayList<double[]>();
		
		for (int i = 0; i < numClusters; i++)
			gaussianDataPoints.addAll(getGaussianCluster());
		
		return gaussianDataPoints;
	}
	
	private static List<double[]> getGaussianCluster() {
		
		List<double[]> cluster = new ArrayList<double[]>();
		
		double deviation = 0;
		
		if (spacingType.equals("s"))
			deviation = 0.01;
		else if (spacingType.equals("o"))
			deviation = 0.06;
		
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
		
		double spacing = 0;
		
		if (spacingType.equals("s"))
			spacing = numPointsPerCluster * epsilon;
		else if (spacingType.equals("o"))
			spacing = epsilon;
		
		List<double[]> densityDataPoints = new ArrayList<double[]>();
		
		double[] firstPoint = new double[numDimensions];
		Arrays.fill(firstPoint, 0);
		
		for (int i = 0; i < numClusters; i++) {
			
			densityDataPoints.addAll(getDensityCluster(epsilon, firstPoint));
	
			//Apply spacing
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
			
			int direction = random.nextBoolean() ? 1 : -1;
			
			for (int j = 0; j < numDimensions; j++)
				nextPoint[j] = firstPoint[j] + ((epsilon * direction) * random.nextDouble());
				
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
		int clusterId = -1;
		
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
						 "#numPointsPerCluster=" + numPointsPerCluster + "\n" +
						 "#spacingType=" + spacingType + "\n" +
						 "#isInconsistent=" + isInconsistent + "\n");
			
			for (double[] row : dataPoints)
				writer.write(Arrays.toString(row).replace(" ", "")
						.replace("[", "").replace("]", "\n"));
			writer.flush();
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Generates and displays scatter plots for the gaussian and density clusters.
	 * Only useful for 2D data points.
	 * 
	 * @param mergedDataPoints The full set of generated data points with point and cluster ids.
	 * @param showPlot Wether the generated plot image should be shown or not.
	 */
	private static void plot(List<double[]> mergedDataPoints, boolean showPlot) throws Exception {
		
		System.out.println("Plotting gaussian clusters ...");
		String filename = "plots/gauss.jpeg";
		
		XYSeriesCollection datasetGauss = new XYSeriesCollection();
		for (int i = 0; i < numClusters; i++) {
			XYSeries series = new XYSeries("Cluster #" + i);
				for (double[] dataPoint : mergedDataPoints)
					if (dataPoint[dataPoint.length - 1] == i)
						series.add(dataPoint[0], dataPoint[1]);
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
				for (double[] dataPoint : mergedDataPoints)
					if (dataPoint[dataPoint.length - 1] == i)
						series.add(dataPoint[0 + numDimensions], dataPoint[1 + numDimensions]);
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
