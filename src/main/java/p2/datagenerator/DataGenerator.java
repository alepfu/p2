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
	static int numDimensions = 5;		//Number of dimensions generated for each type of clusters (gaussian and density)
	static int numClusters = 2;
	static int numPointsPerCluster = 1000;
	static String spacingType = "s";
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
		List<MergedDataPoint> mergedDataPoints = mergeDataPoints(gaussianDataPoints, densityDataPoints);
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
	 * Merges sets of gaussian and density data points together, adding point and cluster ids to every row.
	 * 
	 * @param gaussianDataPoints
	 * @param densityDataPoints
	 * @return List of MergedDataPoint
	 */
	private static List<MergedDataPoint> mergeDataPoints(List<double[]> gaussianDataPoints, List<double[]> densityDataPoints) {
		
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
	
	private static void saveDataPointsToFile(List<MergedDataPoint> dataPoints, String file) {
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(file));
			
			//Add header holding the set parameter values
			writer.write("#numClusters=" + numClusters + "\n" +
						 "#numDimensions=" + numDimensions + "\n" +
						 "#numPointsPerCluster=" + numPointsPerCluster + "\n" +
						 "#spacingType=" + spacingType + "\n" +
						 "#isInconsistent=" + isInconsistent + "\n");
			
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
	 * Only useful for 2D data points.
	 * 
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
