package p2.datagenerator;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class DataGenerator2D {

	/**
	 * Seed for random number generation
	 */
	private static long seed = 7;
	
	/**
	 * Random number generator.
	 */
	private static Random rand = new Random(seed);
	
	/**
	 * Number of clusters to be generated.
	 */
	private static int numClusters = 3;
	
	/**
	 * Number of dimensions generated for each cluster type.
	 */
	private static int numDimensions = 2;
	
	/**
	 * Number of generated points per cluster. 
	 */
	private static int numPointsCluster = 3000;
	
	/**
	 * Working directory
	 */
	private static String workDir = "/home/alepfu/Desktop/P2";
	
	/**
	 * Shuffle rows in output file
	 */
	private static boolean shuffleOutput = false;
	
	public static void main(String[] args) {
		
		//Generate gauss clusters
		System.out.println("Generate gauss clusters ...");
		List<double[]> gaussData = new ArrayList<double[]>();
		
		double[] mean1 = {0, 0};
		double dev1 = 100;
		Gaussian2DCluster g1 = new Gaussian2DCluster(numPointsCluster, mean1, dev1, rand);
		gaussData.addAll(g1.getPoints());
		
		double[] mean2 = {200, 200};
		double dev2 = 50;
		Gaussian2DCluster g2 = new Gaussian2DCluster(numPointsCluster, mean2, dev2, rand);
		gaussData.addAll(g2.getPoints());
		
		double[] mean3 = {-200, 200};
		double dev3 = 75;
		Gaussian2DCluster g3 = new Gaussian2DCluster(numPointsCluster, mean3, dev3, rand);
		gaussData.addAll(g3.getPoints());
		
		
		//Generate density clusters
		System.out.println("Generate density clusters ...");
		List<double[]> densityData = new ArrayList<double[]>();
		
		double[] dpos1 = {0, 0};		
		Density2DCluster d1 = new Density2DCluster(numPointsCluster, Density2DCluster.TYPE_ARC_UP, dpos1, rand);
		densityData.addAll(d1.getPoints());
		
		double[] dpos2 = {100, 95.08};
		Density2DCluster d2 = new Density2DCluster(numPointsCluster, Density2DCluster.TYPE_ARC_DOWN, dpos2, rand);
		densityData.addAll(d2.getPoints());
		
		double[] dpos3 = {-35, 25.2};
		Density2DCluster d3 = new Density2DCluster(numPointsCluster, Density2DCluster.TYPE_CIRCLE, dpos3, rand);
		densityData.addAll(d3.getPoints());
		
		//Merge data points together
		List<DataPoint> dataPoints = getDataPoints(gaussData, densityData);
		
		//Export data points and true clustering
		saveDataPointsToFile(dataPoints);
		saveTrueClusteringToFile(dataPoints);
		
		//Plotting
		plotGaussianClusters(dataPoints, workDir + "/true_gauss.jpeg");
		plotDensityClusters(dataPoints, workDir + "/true_density.jpeg");
		
		
		System.out.println("Finished.");
	}
	
	/**
	 * Merges sets of gaussian and density data points together, adding cluster labels.
	 * @param gaussianDataPoints
	 * @param densityDataPoints
	 * @return List of MergedDataPoint
	 */
	private static List<DataPoint> getDataPoints(List<double[]> gaussianDataPoints, List<double[]> densityDataPoints) {
		
		System.out.println("Merge and label data points ...");
		
		List<DataPoint> dataPoints = new ArrayList<DataPoint>();
		int numPoints = numClusters * numPointsCluster;
		int clusterId = -1;
		
		for (int i = 0; i < numPoints; i++) {
			if ((i % numPointsCluster) == 0)
				clusterId++;
			
			DataPoint p = new DataPoint(gaussianDataPoints.get(i), densityDataPoints.get(i), "c" + (clusterId+1));
			dataPoints.add(p);
		}
		
		if (shuffleOutput)
			Collections.shuffle(dataPoints);
		
		return dataPoints;
	}
	
	/**
	 * Exports data points to file.
	 * Adds header information.
	 * @param dataPoints A list of data points.
	 */
	private static void saveDataPointsToFile(List<DataPoint> dataPoints) {
		
		System.out.println("Export data points ...");
		
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(workDir + "/data.csv"));
			
			//Add header holding the set parameter values
			writer.write("#numClusters=" + numClusters + "\n" +
						 "#numDimensions=" + numDimensions + "\n" +
						 "#numPointsCluster=" + numPointsCluster + "\n" +
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
	 * Generates a file which holds the true clustering
	 * File format: e.g. "1 1 1 2 2" would say objects 1-3 are in cluster 1 and objects 4-5 are in cluster 2.
	 * @param filename A filename
	 * @param clustering A clustering 
	 */
	private static void saveTrueClusteringToFile(List<DataPoint> dataPoints) {
	
		System.out.println("Export true clustering ...");
		
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(workDir + "/true.csv"));
						
			for (DataPoint p : dataPoints) 
				writer.write(p.getClusterLabel().replace("c", "") + " ");
			writer.flush();
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private static void plotGaussianClusters(List<DataPoint> dataPoints, String filename) {
		
		System.out.println("Plotting gaussian clusters ...");
		
		XYSeriesCollection datasetGauss = new XYSeriesCollection();
		for (int i = 1; i <= numClusters; i++) {
			XYSeries series = new XYSeries("Cluster #" + i);
				for (DataPoint dataPoint : dataPoints)
					if (dataPoint.getClusterLabel().equals("c" + i))
						series.add(dataPoint.getGaussFeatures()[0], dataPoint.getGaussFeatures()[1]);
				datasetGauss.addSeries(series);
		}
		
		JFreeChart chartGauss = ChartFactory.createScatterPlot("", "", "", datasetGauss, 
				PlotOrientation.VERTICAL, false, false, false);
		
		try {
			ChartUtilities.saveChartAsJPEG(new File(filename), chartGauss, 660, 420);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		ChartFrame frame = new ChartFrame(filename, chartGauss);
		frame.pack();
		frame.setVisible(true);
	}
	
	private static void plotDensityClusters(List<DataPoint> mergedDataPoints, String filename) {
		
		System.out.println("Plotting density clusters ...");
		
		XYSeriesCollection datasetDensity = new XYSeriesCollection();
		for (int i = 1; i <= numClusters; i++) {
			XYSeries series = new XYSeries("Cluster #" + i);
			for (DataPoint dataPoint : mergedDataPoints)
				if (dataPoint.getClusterLabel().equals("c" + i))
					series.add(dataPoint.getDensityFeatures()[0], dataPoint.getDensityFeatures()[1]);
				datasetDensity.addSeries(series);
		}
		
		JFreeChart chartDensity = ChartFactory.createScatterPlot("", "", "", datasetDensity, 
				PlotOrientation.VERTICAL, false, false, false);
		
		try {
			ChartUtilities.saveChartAsJPEG(new File(filename), chartDensity, 660, 420);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		ChartFrame frame = new ChartFrame(filename, chartDensity);
		frame.pack();
		frame.setVisible(true);
	}

}
