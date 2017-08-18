package p2.datagenerator;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import p2.util.Config;

public class DatasetGenerator {

	/**
	 * Number of clusters to be generated
	 */
	private static int numClusters = Config.numClusters;
	
	/**
	 * Number of generated points per cluster. 
	 */
	private static int numPointsCluster = Config.numPointsCluster;
	
	/**
	 * Working directory
	 */
	private static String workDir = Config.workDir;
	
	/**
	 * Random number generator
	 */
	private static Random rand = new Random(Config.seed);
	
	public static void main(String[] args) {
		
		generateDataset1();
		
		//generateDataset2();
		
	}
	
	private static void generateDataset2() {
		
		//Generate gauss clusters
		double[] mean1 = {0, 0};
		double dev1 = 30;
		GaussCluster g1 = new GaussCluster(numPointsCluster, mean1, dev1, rand);
		double[] mean2 = {60, 60};
		double dev2 = 70;
		GaussCluster g2 = new GaussCluster(numPointsCluster, mean2, dev2, rand);
		
		
		//Generate density clusters with low density regions
		double[] dpos1 = {0, 0};		
		DensityCluster d1 = new DensityCluster(numPointsCluster, DensityCluster.TYPE_ARC_UP_1, dpos1, rand);
		double[] dpos2 = {-10, 58.2};
		DensityCluster d2 = new DensityCluster(numPointsCluster, DensityCluster.TYPE_BOX, dpos2, rand);
		
		
		//Get low/high density ids
		List<Integer> lowIds1 = d1.getLowIds();
		List<Integer> lowIds2 = d2.getLowIds();
		List<Integer> highIds1 = d1.getHighIds();
		List<Integer> highIds2 = d2.getHighIds();
			
		
		//Get good/bad points from Gauss clusters
		List<Integer> goodIds1 = getGoodIdsGauss(mean1, g1, 25.0);
		List<Integer> goodIds2 = getGoodIdsGauss(mean2, g2, 25.0);
		List<List<Integer>> twoLists = getBadIdsGauss(mean1, mean2, dev1, dev2, g1, g2, 1.1, 0.7);
		List<Integer> badIds1 = twoLists.get(0);
		List<Integer> badIds2 = twoLists.get(1);		
		

		//Plotting
//		plotGaussClusters(g1.getPoints(), g2.getPoints(), workDir + "/true_gauss.jpeg");
//		plotGaussClustersGoodBad(g1.getPoints(), g2.getPoints(), workDir + "/badgood_gauss.jpeg", goodIds1, goodIds2, badIds1, badIds2);
		plotDensityClusters(d1.getPoints(), d2.getPoints(), workDir + "/true_density.jpeg");
//		plotDensityClustersLow(d1.getPoints(), d2.getPoints(), workDir + "/low_density.jpeg", lowIds1, lowIds2); 				
			
		
		//Merge data points together
		List<DataPoint> dataPoints1 = mergePoints(g1, d1, goodIds1, badIds1, lowIds1, highIds1, "1");
		List<DataPoint> dataPoints2 = mergePoints(g2, d2, goodIds2, badIds2, lowIds2, highIds2, "2");
		
		
		//Export data points and true clustering
		List<DataPoint> allDataPoints = dataPoints1;
		allDataPoints.addAll(dataPoints2);
		saveDataPoints(allDataPoints);
		saveTrueClustering(allDataPoints);
		
		
		//Verify the merged dataset
		int sizeGoodLow1 = goodIds1.size() < lowIds1.size() ? goodIds1.size() : lowIds1.size();
		int sizeBadHigh1 = badIds1.size() < highIds1.size() ? badIds1.size() : highIds1.size();
		int sizeGoodLow2 = goodIds2.size() < lowIds2.size() ? goodIds2.size() : lowIds2.size();
		int sizeBadHigh2 = badIds2.size() < highIds2.size() ? badIds2.size() : highIds2.size();
//		plotVerifyGauss(allDataPoints, workDir + "/verify_gauss.jpeg", sizeGoodLow1, sizeBadHigh1, sizeGoodLow2, sizeBadHigh2);
//		plotVerifyDensity(allDataPoints, workDir + "/verify_density.jpeg", sizeGoodLow1, sizeBadHigh1, sizeGoodLow2, sizeBadHigh2);
	}
	
	private static void generateDataset1() {
		
		//Generate gauss clusters
		double[] mean1 = {0, 0};
		double dev1 = 30;
		GaussCluster g1 = new GaussCluster(numPointsCluster, mean1, dev1, rand);
		double[] mean2 = {60, 60};
		double dev2 = 70;
		GaussCluster g2 = new GaussCluster(numPointsCluster, mean2, dev2, rand);
		
		
		//Generate density clusters with low density regions
		double[] dpos1 = {0, 0};		
		DensityCluster d1 = new DensityCluster(numPointsCluster, DensityCluster.TYPE_ARC_UP_1, dpos1, rand);
		double[] dpos2 = {0, 200.9};
		DensityCluster d2 = new DensityCluster(numPointsCluster, DensityCluster.TYPE_ARC_DOWN_1, dpos2, rand);
		
		
		//Get low/high density ids
		List<Integer> lowIds1 = d1.getLowIds();
		List<Integer> lowIds2 = d2.getLowIds();
		List<Integer> highIds1 = d1.getHighIds();
		List<Integer> highIds2 = d2.getHighIds();
			
		
		//Get good/bad points from Gauss clusters
		List<Integer> goodIds1 = getGoodIdsGauss(mean1, g1, 25.0);
		List<Integer> goodIds2 = getGoodIdsGauss(mean2, g2, 25.0);
		List<List<Integer>> twoLists = getBadIdsGauss(mean1, mean2, dev1, dev2, g1, g2, 1.1, 0.7);
		List<Integer> badIds1 = twoLists.get(0);
		List<Integer> badIds2 = twoLists.get(1);		
		

		//Plotting
		plotGaussClusters(g1.getPoints(), g2.getPoints(), workDir + "/true_gauss.jpeg");
//		plotGaussClustersGoodBad(g1.getPoints(), g2.getPoints(), workDir + "/badgood_gauss.jpeg", goodIds1, goodIds2, badIds1, badIds2);
		plotDensityClusters(d1.getPoints(), d2.getPoints(), workDir + "/true_density.jpeg");
//		plotDensityClustersLow(d1.getPoints(), d2.getPoints(), workDir + "/low_density.jpeg", lowIds1, lowIds2); 				
			
		
		//Merge data points together
		List<DataPoint> dataPoints1 = mergePoints(g1, d1, goodIds1, badIds1, lowIds1, highIds1, "1");
		List<DataPoint> dataPoints2 = mergePoints(g2, d2, goodIds2, badIds2, lowIds2, highIds2, "2");
		List<DataPoint> allDataPoints = new ArrayList<DataPoint>();
		allDataPoints.addAll(dataPoints1);
		allDataPoints.addAll(dataPoints2);
		
		
		//Add noise
		// 
		// !!!!!!!!!!!!! Don't forget to set Config.numNoise Points !!!!!!!!!!!!
		//
//		List<DataPoint> noisePoints = getNoisePoints(Config.numNoisePoints, -400, 400, -800, 500, -600, 300, -100, 300);
//		plotGaussClustersNoise(dataPoints1, dataPoints2, noisePoints, workDir + "/noise_gauss.jpeg");
//		plotDensityClustersNoise(dataPoints1, dataPoints2, noisePoints, workDir + "/noise_gauss.jpeg");
//		allDataPoints.addAll(noisePoints);
		
		
		
		
		//Export data points and true clustering
		saveDataPoints(allDataPoints);
		saveTrueClustering(allDataPoints);
		
		
		//Verify the merged dataset
		int sizeGoodLow1 = goodIds1.size() < lowIds1.size() ? goodIds1.size() : lowIds1.size();
		int sizeBadHigh1 = badIds1.size() < highIds1.size() ? badIds1.size() : highIds1.size();
		int sizeGoodLow2 = goodIds2.size() < lowIds2.size() ? goodIds2.size() : lowIds2.size();
		int sizeBadHigh2 = badIds2.size() < highIds2.size() ? badIds2.size() : highIds2.size();
//		plotVerifyGauss(allDataPoints, workDir + "/verify_gauss.jpeg", sizeGoodLow1, sizeBadHigh1, sizeGoodLow2, sizeBadHigh2);
//		plotVerifyDensity(allDataPoints, workDir + "/verify_density.jpeg", sizeGoodLow1, sizeBadHigh1, sizeGoodLow2, sizeBadHigh2);
	}
	
	private static List<DataPoint> mergePoints(GaussCluster g, DensityCluster d, List<Integer> goodIdsGauss, List<Integer> badIdsGauss, List<Integer> lowIdsDensity, List<Integer> highIdsDensity, String clusterLabel) {
		
		List<DataPoint> dataPoints = new ArrayList<DataPoint>();
		
		List<double[]> pointsGauss = g.getPoints();
		List<double[]> pointsDensity = d.getPoints();
		
		List<double[]> usedGauss = new ArrayList<double[]>();
		List<double[]> usedDensity = new ArrayList<double[]>();
		
		//Merge good Gauss with low Density
		int nGoodLow = goodIdsGauss.size() < lowIdsDensity.size() ? goodIdsGauss.size() : lowIdsDensity.size();
		for (int i = 0; i < nGoodLow; i++) {
			double[] pointGauss = pointsGauss.get(goodIdsGauss.get(i));
			double[] pointDensity = pointsDensity.get(lowIdsDensity.get(i));
			dataPoints.add(new DataPoint(pointGauss, pointDensity, clusterLabel));
			usedGauss.add(pointGauss);
			usedDensity.add(pointDensity);
		}
		
		//Merge bad Gauss with high Density
		int nBadHigh = badIdsGauss.size() < highIdsDensity.size() ? badIdsGauss.size() : highIdsDensity.size();
		for (int i = 0; i < nBadHigh; i++) {
			double[] pointGauss = pointsGauss.get(badIdsGauss.get(i));
			double[] pointDensity = pointsDensity.get(highIdsDensity.get(i));
			dataPoints.add(new DataPoint(pointGauss, pointDensity, clusterLabel));
			usedGauss.add(pointGauss);
			usedDensity.add(pointDensity);
		}
		
		//Remove used Ids
		pointsGauss.removeAll(usedGauss);
		pointsDensity.removeAll(usedDensity);
		
		//Merge residual points
		for (int i = 0; i < pointsGauss.size(); i++)
			dataPoints.add(new DataPoint(pointsGauss.get(i), pointsDensity.get(i), clusterLabel));
		
		return dataPoints;
	}
	
	/**
	 * Exports data points to file.
	 * @param dataPoints A list of data points.
	 */
	private static void saveDataPoints(List<DataPoint> dataPoints) {
		
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(workDir + "/data.csv"));
			
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
	private static void saveTrueClustering(List<DataPoint> dataPoints) {
	
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(workDir + "/true.csv"));
						
			for (DataPoint p : dataPoints) 
				writer.write(p.getClusterLabel() + " ");
			writer.flush();
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private static void plotGaussClustersGoodBad(List<double[]> points1, List<double[]> points2, String filename, List<Integer> goodIds1, List<Integer> goodIds2, List<Integer> badIds1, List<Integer> badIds2) {
		
		XYSeriesCollection dataset = new XYSeriesCollection();
		XYSeries series1 = new XYSeries("1");
		XYSeries series2 = new XYSeries("2");
		XYSeries bad = new XYSeries("bad");
		XYSeries good = new XYSeries("good");
		dataset.addSeries(good);
		dataset.addSeries(bad);
		dataset.addSeries(series1);
		dataset.addSeries(series2);
		
		for (int id = 0; id < points1.size(); id++) {
			double[] p = points1.get(id);
			if (goodIds1.contains(id))
				good.add(p[0], p[1]);
			else if (badIds1.contains(id))
				bad.add(p[0], p[1]);
			else
				series1.add(p[0], p[1]);
		}
		
		for (int id = 0; id < points2.size(); id++) {
			double[] p = points2.get(id);
			if (goodIds2.contains(id))
				good.add(p[0], p[1]);
			else if (badIds2.contains(id))
				bad.add(p[0], p[1]);
			else
				series2.add(p[0], p[1]);
		}
		
		JFreeChart chart = ChartFactory.createScatterPlot("", "", "", dataset, 
				PlotOrientation.VERTICAL, true, false, false);
		
		try {
			ChartUtilities.saveChartAsJPEG(new File(filename), chart, 660, 420);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		ChartFrame frame = new ChartFrame(filename, chart);
		frame.pack();
		frame.setVisible(true);
	}
	
	private static void plotGaussClusters(List<double[]> dataPoints1, List<double[]> dataPoints2, String filename) {
		
		XYSeriesCollection dataset = new XYSeriesCollection();
		XYSeries series1 = new XYSeries("1");
		XYSeries series2 = new XYSeries("2");
		dataset.addSeries(series1);
		dataset.addSeries(series2);
		
		for (double[] p : dataPoints1)
			series1.add(p[0], p[1]);
		for (double[] p : dataPoints2)
			series2.add(p[0], p[1]);
		
		JFreeChart chart = ChartFactory.createScatterPlot("", "", "", dataset, 
				PlotOrientation.VERTICAL, true, false, false);
		
		try {
			ChartUtilities.saveChartAsJPEG(new File(filename), chart, 660, 420);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		ChartFrame frame = new ChartFrame(filename, chart);
		frame.pack();
		frame.setVisible(true);
	}
	
	private static void plotDensityClustersLow(List<double[]> points1, List<double[]> points2, String filename, List<Integer> lowIds1, List<Integer> lowIds2) {
		
		XYSeriesCollection datasetDensity = new XYSeriesCollection();
		XYSeries low = new XYSeries("low");
		XYSeries series1 = new XYSeries("1");
		XYSeries series2 = new XYSeries("2");
		datasetDensity.addSeries(low);
		datasetDensity.addSeries(series1);
		datasetDensity.addSeries(series2);
		
		for (int id = 0; id < points1.size(); id++) {
			double[] p = points1.get(id);
			if (lowIds1.contains(id))
				low.add(p[0], p[1]);
			else
				series1.add(p[0], p[1]);
		}
		
		for (int id = 0; id < points2.size(); id++) {
			double[] p = points2.get(id);
			if (lowIds2.contains(id))
				low.add(p[0], p[1]);
			else
				series2.add(p[0], p[1]);
		}
		
		JFreeChart chartDensity = ChartFactory.createScatterPlot("", "", "", datasetDensity, 
				PlotOrientation.VERTICAL, true, false, false);
		
		try {
			ChartUtilities.saveChartAsJPEG(new File(filename), chartDensity, 660, 420);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		ChartFrame frame = new ChartFrame(filename, chartDensity);
		frame.pack();
		frame.setVisible(true);
	}
	
	private static void plotDensityClusters(List<double[]> dataPoints1, List<double[]> dataPoints2, String filename) {
		
		XYSeriesCollection dataset = new XYSeriesCollection();
		XYSeries series1 = new XYSeries("1");
		XYSeries series2 = new XYSeries("2");
		dataset.addSeries(series1);
		dataset.addSeries(series2);
		
		for (double[] p : dataPoints1)
			series1.add(p[0], p[1]);
		for (double[] p : dataPoints2)
			series2.add(p[0], p[1]);
		
		JFreeChart chartDensity = ChartFactory.createScatterPlot("", "", "", dataset, 
				PlotOrientation.VERTICAL, true, false, false);
		
		try {
			ChartUtilities.saveChartAsJPEG(new File(filename), chartDensity, 660, 420);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		ChartFrame frame = new ChartFrame(filename, chartDensity);
		frame.pack();
		frame.setVisible(true);
	}
	
	private static List<List<Integer>> getBadIdsGauss(double[] mean1, double[] mean2, double dev1, double dev2, GaussCluster g1, GaussCluster g2, double distFactor1, double distFactor2) {
		
		//Cluster 1
		double distMeans = Math.sqrt(Math.pow(mean1[0] - mean2[0], 2) + Math.pow(mean1[1] - mean2[1], 2));
		List<Integer> badIds1 = new ArrayList<Integer>();
		for(int i = 0; i < g1.getPoints().size(); i++) {
			double[] p = g1.getPoints().get(i);
			double distToMean1 = Math.sqrt(Math.pow(mean1[0] - p[0], 2) + Math.pow(mean1[1] - p[1], 2));
			if (distToMean1 > (dev1 * distFactor1)) {
				double distToMean2 = Math.sqrt(Math.pow(mean2[0] - p[0], 2) + Math.pow(mean2[1] - p[1], 2));
				if(distToMean2 < (distMeans * distFactor2))			
					badIds1.add(i);
			}
		}
		
		//Cluster 2
		List<Integer> badIds2 = new ArrayList<Integer>();
		for(int i = 0; i < g2.getPoints().size(); i++) {
			double[] p = g2.getPoints().get(i);
			double distToMean2 = Math.sqrt(Math.pow(mean2[0] - p[0], 2) + Math.pow(mean2[1] - p[1], 2));
			if (distToMean2 > (dev2 * distFactor1)) {
				double distToMean1 = Math.sqrt(Math.pow(mean1[0] - p[0], 2) + Math.pow(mean1[1] - p[1], 2));
				if(distToMean1 < (distMeans * distFactor2))			
					badIds2.add(i);
			}
		}
		
		List<List<Integer>> twoLists = new ArrayList<List<Integer>>();
		twoLists.add(badIds1);
		twoLists.add(badIds2);
		
		return twoLists;
	}

	private static List<Integer> getGoodIdsGauss(double[] mean, GaussCluster g, double maxDist) {
		
		List<Integer> goodIds = new ArrayList<Integer>();
		for(int i = 0; i < g.getPoints().size(); i++) {
			double[] p = g.getPoints().get(i);
			double distToMean = Math.sqrt(Math.pow(mean[0] - p[0], 2) + Math.pow(mean[1] - p[1], 2));
			if (distToMean < maxDist)
				goodIds.add(i); 
		}
		
		return goodIds;
	}
	
	private static void plotVerifyGauss(List<DataPoint> dataPoints, String filename, int sizeGoodLow1, int sizeBadHigh1, int sizeGoodLow2, int sizeBadHigh2) {
		
		XYSeriesCollection dataset = new XYSeriesCollection();
		
		XYSeries goodlow1 = new XYSeries("goodlow1");
		for (int i = 0; i < sizeGoodLow1; i++)
			goodlow1.add(dataPoints.get(i).getGaussFeatures()[0], dataPoints.get(i).getGaussFeatures()[1]);
		dataset.addSeries(goodlow1);
		XYSeries badhigh1 = new XYSeries("badhigh1");
		for (int i = 0; i < sizeBadHigh1; i++)
			badhigh1.add(dataPoints.get(i + sizeGoodLow1).getGaussFeatures()[0], dataPoints.get(i + sizeGoodLow1).getGaussFeatures()[1]);	
		dataset.addSeries(badhigh1);
		
		XYSeries goodlow2 = new XYSeries("goodlow2");
		for (int i = 0; i < sizeGoodLow2; i++)
			goodlow2.add(dataPoints.get(numPointsCluster + i).getGaussFeatures()[0], dataPoints.get(numPointsCluster + i).getGaussFeatures()[1]);
		dataset.addSeries(goodlow2);
		XYSeries badhigh2 = new XYSeries("badhigh2");
		for (int i = 0; i < sizeBadHigh2; i++)
			badhigh2.add(dataPoints.get(numPointsCluster + i + sizeGoodLow2).getGaussFeatures()[0], dataPoints.get(numPointsCluster + i + sizeGoodLow2).getGaussFeatures()[1]);	
		dataset.addSeries(badhigh2);
		
		JFreeChart chart = ChartFactory.createScatterPlot("", "", "", dataset, 
				PlotOrientation.VERTICAL, true, false, false);
		
		try {
			ChartUtilities.saveChartAsJPEG(new File(filename), chart, 660, 420);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		ChartFrame frame = new ChartFrame(filename, chart);
		frame.pack();
		frame.setVisible(true);
		
	}
	
	private static void plotVerifyDensity(List<DataPoint> dataPoints, String filename, int sizeGoodLow1, int sizeBadHigh1, int sizeGoodLow2, int sizeBadHigh2) {
		
		XYSeriesCollection dataset = new XYSeriesCollection();
		
		XYSeries goodlow1 = new XYSeries("goodlow1");
		for (int i = 0; i < sizeGoodLow1; i++)
			goodlow1.add(dataPoints.get(i).getDensityFeatures()[0], dataPoints.get(i).getDensityFeatures()[1]);
		dataset.addSeries(goodlow1);
		XYSeries badhigh1 = new XYSeries("badhigh1");
		for (int i = 0; i < sizeBadHigh1; i++)
			badhigh1.add(dataPoints.get(i + sizeGoodLow1).getDensityFeatures()[0], dataPoints.get(i + sizeGoodLow1).getDensityFeatures()[1]);	
		dataset.addSeries(badhigh1);
		
		XYSeries goodlow2 = new XYSeries("goodlow2");
		for (int i = 0; i < sizeGoodLow2; i++)
			goodlow2.add(dataPoints.get(numPointsCluster + i).getDensityFeatures()[0], dataPoints.get(numPointsCluster + i).getDensityFeatures()[1]);
		dataset.addSeries(goodlow2);
		XYSeries badhigh2 = new XYSeries("badhigh2");
		for (int i = 0; i < sizeBadHigh2; i++)
			badhigh2.add(dataPoints.get(numPointsCluster + i + sizeGoodLow2).getDensityFeatures()[0], dataPoints.get(numPointsCluster + i + sizeGoodLow2).getDensityFeatures()[1]);	
		dataset.addSeries(badhigh2);
		
		//Keep here because of dumb auto-scaling of JFreeChart
		for (int i = 1; i <= numClusters; i++) {
			XYSeries series = new XYSeries(""+ i);
			for (DataPoint dataPoint : dataPoints)
				if (dataPoint.getClusterLabel().equals(i))
					series.add(dataPoint.getDensityFeatures()[0], dataPoint.getDensityFeatures()[1]);
			dataset.addSeries(series);
		}
		
		JFreeChart chart = ChartFactory.createScatterPlot("", "", "", dataset, 
				PlotOrientation.VERTICAL, true, false, false);
		
		try {
			ChartUtilities.saveChartAsJPEG(new File(filename), chart, 660, 420);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		ChartFrame frame = new ChartFrame(filename, chart);
		frame.pack();
		frame.setVisible(true);
		
	}
	
	private static List<DataPoint> getNoisePoints(int numNoisePoints, double f1, double t1, double f2, double t2, double f3, double t3, double f4, double t4) {
		
		List<DataPoint> noisePoints = new ArrayList<DataPoint>(numNoisePoints);
		
		for (int i = 0; i < numNoisePoints; i++) {
			
			double[] gauss = new double[2];
			gauss[0] = f1 + (t1 - f1) * rand.nextDouble();
			gauss[1] = f2 + (t2 - f2) * rand.nextDouble();
			
			double[] density = new double[2];
			density[0] = f3 + (t3 - f3) * rand.nextDouble();
			density[1] = f4 + (t4 - f4) * rand.nextDouble();
			
			noisePoints.add(new DataPoint(gauss, density, "3"));
		}
		
		return noisePoints;
	}
	
	private static void plotGaussClustersNoise(List<DataPoint> dataPoints1, List<DataPoint> dataPoints2, List<DataPoint> noisePoints, String filename) {
		
		XYSeriesCollection dataset = new XYSeriesCollection();
		XYSeries series1 = new XYSeries("1");
		XYSeries series2 = new XYSeries("2");
		XYSeries series0 = new XYSeries("0");
		dataset.addSeries(series1);
		dataset.addSeries(series2);
		dataset.addSeries(series0);
		
		for (DataPoint p : dataPoints1)
			series1.add(p.getGaussFeatures()[0], p.getGaussFeatures()[1]);
		for (DataPoint p : dataPoints2)
			series2.add(p.getGaussFeatures()[0], p.getGaussFeatures()[1]);
		for (DataPoint p : noisePoints)
			series0.add(p.getGaussFeatures()[0], p.getGaussFeatures()[1]);
		
		JFreeChart chart = ChartFactory.createScatterPlot("", "", "", dataset, 
				PlotOrientation.VERTICAL, true, false, false);
		
		try {
			ChartUtilities.saveChartAsJPEG(new File(filename), chart, 660, 420);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		ChartFrame frame = new ChartFrame(filename, chart);
		frame.pack();
		frame.setVisible(true);
	}
	
	private static void plotDensityClustersNoise(List<DataPoint> dataPoints1, List<DataPoint> dataPoints2, List<DataPoint> noisePoints, String filename) {
		
		XYSeriesCollection dataset = new XYSeriesCollection();
		XYSeries series1 = new XYSeries("1");
		XYSeries series2 = new XYSeries("2");
		XYSeries series0 = new XYSeries("0");
		dataset.addSeries(series1);
		dataset.addSeries(series2);
		dataset.addSeries(series0);
		
		for (DataPoint p : dataPoints1)
			series1.add(p.getDensityFeatures()[0], p.getDensityFeatures()[1]);
		for (DataPoint p : dataPoints2)
			series2.add(p.getDensityFeatures()[0], p.getDensityFeatures()[1]);
		for (DataPoint p : noisePoints)
			series0.add(p.getDensityFeatures()[0], p.getDensityFeatures()[1]);
		
		JFreeChart chartDensity = ChartFactory.createScatterPlot("", "", "", dataset, 
				PlotOrientation.VERTICAL, true, false, false);
		
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