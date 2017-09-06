package p2.datagenerator;

/**
 * Copyright (c) 2017 Alexander Pfundner
 * 
 * Integration of Density-based and Partitioning-based Clustering Methods
 * 
 */

import java.awt.Color;
import java.awt.Paint;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.geom.Ellipse2D;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.plot.DefaultDrawingSupplier;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.graphics2d.svg.SVGGraphics2D;
import org.jfree.graphics2d.svg.SVGUtils;

import p2.util.Config;

/**
 * Generation of datasets for use with integrative workflow.
 * 
 */
public class DatasetGenerator {

	/**
	 * Random number generator.
	 */
	private static Random rand = new Random(Config.seed);
	
	/**
	 * Generates 2 differnet datasets.
	 * Does not take arguments, parameters are set via class Config.
	 * 
	 */
	public static void main(String[] args) {
		
		System.out.println("Generating dataset 1 ...");
		generateDataset1("dataset_1");
		
		System.out.println("Generating dataset 2 ...");
		generateDataset2("dataset_2");
		
		System.out.println("\nFinished.");
	}
	
	/**
	 * Generation of dataset 2.
	 * @param dir The subfolder for the generated files.
	 */
	private static void generateDataset2(String dir) {
		
		//Generate gauss clusters
		double[] mean1 = {0, 0};
		double dev1 = 30;
		GaussCluster g1 = new GaussCluster(Config.numPointsCluster, mean1, dev1, rand);
		double[] mean2 = {60, 60};
		double dev2 = 70;
		GaussCluster g2 = new GaussCluster(Config.numPointsCluster, mean2, dev2, rand);
		
		
		//Generate density clusters with low density regions
		double[] dpos1 = {0, 0};		
		DensityCluster d1 = new DensityCluster(Config.numPointsCluster, DensityCluster.TYPE_ARC_UP_1, dpos1, rand);
		double[] dpos2 = {-10, 58.2};
		DensityCluster d2 = new DensityCluster(Config.numPointsCluster, DensityCluster.TYPE_BOX, dpos2, rand);
		
		
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
		

		//Plot clusters
		if (Config.plotClusters) {
			plotGaussClusters(g1.getPoints(), g2.getPoints(), Config.workDir + "/" + dir + "/true_gauss");
			plotDensityClusters(d1.getPoints(), d2.getPoints(), Config.workDir + "/" + dir + "/true_density");
		}
		
		
		//Plot clusters with good/bad and high/low regions highlighted
		if (Config.plotClustersHighlighted) {
			plotGaussClustersGoodBad(g1.getPoints(), g2.getPoints(), Config.workDir + "/" + dir + "/badgood_gauss", goodIds1, goodIds2, badIds1, badIds2);
			plotDensityClustersLow(d1.getPoints(), d2.getPoints(), Config.workDir + "/" + dir + "/low_density", lowIds1, lowIds2); 				
		}		
			
		
		//Merge data points together
		List<DataPoint> dataPoints1 = mergePoints(g1, d1, goodIds1, badIds1, lowIds1, highIds1, "1");
		List<DataPoint> dataPoints2 = mergePoints(g2, d2, goodIds2, badIds2, lowIds2, highIds2, "2");
		
		
		//Export data points and true clustering
		List<DataPoint> allDataPoints = dataPoints1;
		allDataPoints.addAll(dataPoints2);
		saveDataPoints(allDataPoints, dir);
		saveTrueClustering(allDataPoints, dir);
		
		
		//Verify the merged dataset
		if (Config.verfiyDataset)
			verifyDataset(allDataPoints, dir, goodIds1, badIds1, goodIds2, badIds2, lowIds1, highIds1, lowIds2, highIds2);
	}
	
	/**
	 * Generation of dataset 1.
	 * @param dir The subfolder for the generated files.
	 */	
	private static void generateDataset1(String dir) {
		
		//Generate gauss clusters
		double[] mean1 = {0, 0};
		double dev1 = 30;
		GaussCluster g1 = new GaussCluster(Config.numPointsCluster, mean1, dev1, rand);
		double[] mean2 = {60, 60};
		double dev2 = 70;
		GaussCluster g2 = new GaussCluster(Config.numPointsCluster, mean2, dev2, rand);
		
		
		//Generate density clusters with low density regions
		double[] dpos1 = {0, 0};		
		DensityCluster d1 = new DensityCluster(Config.numPointsCluster, DensityCluster.TYPE_ARC_UP_1, dpos1, rand);
		double[] dpos2 = {0, 200.9};
		DensityCluster d2 = new DensityCluster(Config.numPointsCluster, DensityCluster.TYPE_ARC_DOWN_1, dpos2, rand);
		
		
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
		

		//Plot clusters
		if (Config.plotClusters) {
			plotGaussClusters(g1.getPoints(), g2.getPoints(), Config.workDir + "/" + dir + "/true_gauss");
			plotDensityClusters(d1.getPoints(), d2.getPoints(), Config.workDir + "/" + dir + "/true_density");
		}
		
		
		//Plot clusters with good/bad and high/low regions highlighted
		if (Config.plotClustersHighlighted) {
			plotGaussClustersGoodBad(g1.getPoints(), g2.getPoints(), Config.workDir + "/" + dir + "/badgood_gauss", goodIds1, goodIds2, badIds1, badIds2);
			plotDensityClustersLow(d1.getPoints(), d2.getPoints(), Config.workDir + "/" + dir + "/low_density", lowIds1, lowIds2); 				
		}
		
		
		//Merge data points together
		List<DataPoint> dataPoints1 = mergePoints(g1, d1, goodIds1, badIds1, lowIds1, highIds1, "1");
		List<DataPoint> dataPoints2 = mergePoints(g2, d2, goodIds2, badIds2, lowIds2, highIds2, "2");
		List<DataPoint> allDataPoints = new ArrayList<DataPoint>();
		allDataPoints.addAll(dataPoints1);
		allDataPoints.addAll(dataPoints2);
		
		
		//Add noise
		if (Config.numNoisePoints > 0) {
			List<DataPoint> noisePoints = getNoisePoints(Config.numNoisePoints, -400, 400, -800, 500, -600, 300, -100, 300);
			
			if (Config.plotClusters) {
				plotGaussClustersNoise(dataPoints1, dataPoints2, noisePoints, Config.workDir + "/" + dir + "/noise_gauss");
				plotDensityClustersNoise(dataPoints1, dataPoints2, noisePoints, Config.workDir + "/" + dir + "/noise_density");
			}
			
			allDataPoints.addAll(noisePoints);
		}
		
		
		//Export data points and true clustering
		saveDataPoints(allDataPoints, dir);
		saveTrueClustering(allDataPoints, dir);
		
		
		//Verify the merged dataset
		if (Config.verfiyDataset)
			verifyDataset(allDataPoints, dir, goodIds1, badIds1, goodIds2, badIds2, lowIds1, highIds1, lowIds2, highIds2);
	}
	
	/**
	 * Merges the generated clusters together.
	 * @param g A gauss cluster.
	 * @param d A density cluster. 
	 * @param goodIdsGauss List of all easy to cluster gauss Ids.
	 * @param badIdsGauss List of all hard to cluster gauss Ids.
	 * @param lowIdsDensity List of all hard to cluster density Ids.
	 * @param highIdsDensity List of all easy to cluster density Ids.
	 * @param clusterLabel The cluster label to use.
	 * @return A list of datapoints.
	 */
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
	 * @param dir The output directory.
	 */
	private static void saveDataPoints(List<DataPoint> dataPoints, String dir) {
		
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(Config.workDir + "/"+ dir + "/data.csv"));
			
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
	 * @param dir The output directory
	 */
	private static void saveTrueClustering(List<DataPoint> dataPoints, String dir) {
	
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(Config.workDir + "/"+ dir + "/true.csv"));
						
			for (DataPoint p : dataPoints) 
				writer.write(p.getClusterLabel() + " ");
			writer.flush();
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Generation of noise points.
	 * @param numNoisePoints Amount of noise points to be generated.
	 * @param f1 From range for noise points in feature 1.
	 * @param t1 To range for noise points in feature 1.
	 * @param f2 From range for noise points in feature 2.
	 * @param t2 To range for noise points in feature 2.
	 * @param f3 From range for noise points in feature 3.
	 * @param t3 To range for noise points in feature 3.
	 * @param f4 From range for noise points in feature 4.
	 * @param t4 To range for noise points in feature 4.
	 * @return A list of datapoints.
	 */
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
	
	/**
	 * Determines gauss datapoints that are hard to cluster.
	 * @param mean1
	 * @param mean2
	 * @param dev1
	 * @param dev2
	 * @param g1
	 * @param g2
	 * @param distFactor1
	 * @param distFactor2
	 * @return List of ids.
	 */
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

	/**
	 * Determines gauss datapoints that are easy to cluster.
	 * @param mean
	 * @param g
	 * @param maxDist
	 * @return A list of ids.
	 */
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
	
	/**
	 * Plots the gauss clusters.
	 * @param dataPoints1
	 * @param dataPoints2
	 * @param filename
	 */
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
		
		JFreeChart chart = ChartFactory.createScatterPlot("Gauss Dimensions", "f1", "f2", dataset, PlotOrientation.VERTICAL, true, false, false);
		XYPlot plot = (XYPlot)chart.getPlot();
		plot.setDrawingSupplier(getCustomDrawingSupplier());
		plot.setBackgroundPaint(new Color(210, 210, 210));
	    plot.setOutlinePaint(Color.white);
	    NumberAxis range = (NumberAxis)plot.getRangeAxis();
        range.setTickUnit(new NumberTickUnit(100));
        NumberAxis domain = (NumberAxis)plot.getDomainAxis();
        domain.setTickUnit(new NumberTickUnit(100));
        plot.setDomainGridlinesVisible(false);
        plot.setRangeGridlinesVisible(false);
        
		try {
			SVGGraphics2D svg = new SVGGraphics2D(600, 600);
	        Rectangle area = new Rectangle(0, 0, 600, 600);
	        chart.draw(svg, area);
	        SVGUtils.writeToSVG(new File(filename + ".svg"), svg.getSVGElement());
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		ChartFrame frame = new ChartFrame(filename, chart);
		frame.pack();
		frame.setBounds(0, 0, 600, 600);
		frame.setVisible(Config.displayPlots);
	}
	
	/**
	 * Plots the density clusters
	 * @param dataPoints1
	 * @param dataPoints2
	 * @param filename
	 */
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
		
		JFreeChart chart = ChartFactory.createScatterPlot("Density Dimensions", "\u03C1" + "1", "\u03C1" + "2", dataset, PlotOrientation.VERTICAL, true, false, false);
		XYPlot plot = (XYPlot)chart.getPlot();
		plot.setDrawingSupplier(getCustomDrawingSupplier());
		plot.setBackgroundPaint(new Color(210, 210, 210));
	    plot.setOutlinePaint(Color.white);
	    NumberAxis range = (NumberAxis)plot.getRangeAxis();
        range.setTickUnit(new NumberTickUnit(50));
        NumberAxis domain = (NumberAxis)plot.getDomainAxis();
        domain.setTickUnit(new NumberTickUnit(50));
        plot.setDomainGridlinesVisible(false);
        plot.setRangeGridlinesVisible(false);
        
		try {
			SVGGraphics2D svg = new SVGGraphics2D(600, 600);
	        Rectangle area = new Rectangle(0, 0, 600, 600);
	        chart.draw(svg, area);
	        SVGUtils.writeToSVG(new File(filename + ".svg"), svg.getSVGElement());
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		ChartFrame frame = new ChartFrame(filename, chart);
		frame.pack();
		frame.setBounds(0, 0, 600, 600);
		frame.setVisible(Config.displayPlots);
	}
	
	/**
	 * Plots the gaus clusters with highlighted hard and easy to cluster datapoints.
	 * @param points1
	 * @param points2
	 * @param filename
	 * @param goodIds1
	 * @param goodIds2
	 * @param badIds1
	 * @param badIds2
	 */
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
		
		JFreeChart chart = ChartFactory.createScatterPlot("Gauss Highlighted", "f1", "f2", dataset, PlotOrientation.VERTICAL, true, false, false);
		XYPlot plot = (XYPlot)chart.getPlot();
		plot.setDrawingSupplier(getCustomDrawingSupplier());
		plot.setBackgroundPaint(new Color(210, 210, 210));
	    plot.setOutlinePaint(Color.white);
	    NumberAxis range = (NumberAxis)plot.getRangeAxis();
        range.setTickUnit(new NumberTickUnit(100));
        NumberAxis domain = (NumberAxis)plot.getDomainAxis();
        domain.setTickUnit(new NumberTickUnit(100));
        plot.setDomainGridlinesVisible(false);
        plot.setRangeGridlinesVisible(false);
        
		try {
			SVGGraphics2D svg = new SVGGraphics2D(600, 600);
	        Rectangle area = new Rectangle(0, 0, 600, 600);
	        chart.draw(svg, area);
	        SVGUtils.writeToSVG(new File(filename + ".svg"), svg.getSVGElement());
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		ChartFrame frame = new ChartFrame(filename, chart);
		frame.pack();
		frame.setBounds(0, 0, 600, 600);
		frame.setVisible(Config.displayPlots);
	}
	
	/**
	 * Plots the density clusters with highlighted low density areas.
	 * @param points1
	 * @param points2
	 * @param filename
	 * @param lowIds1
	 * @param lowIds2
	 */
	private static void plotDensityClustersLow(List<double[]> points1, List<double[]> points2, String filename, List<Integer> lowIds1, List<Integer> lowIds2) {
		
		XYSeriesCollection dataset = new XYSeriesCollection();
		XYSeries low = new XYSeries("low");
		XYSeries series1 = new XYSeries("1");
		XYSeries series2 = new XYSeries("2");
		dataset.addSeries(low);
		dataset.addSeries(series1);
		dataset.addSeries(series2);
		
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
		
		JFreeChart chart = ChartFactory.createScatterPlot("Density Highlighted", "\u03C1" + "1", "\u03C1" + "2", dataset, PlotOrientation.VERTICAL, true, false, false);
		XYPlot plot = (XYPlot)chart.getPlot();
		plot.setDrawingSupplier(getCustomDrawingSupplier());
		plot.setBackgroundPaint(new Color(210, 210, 210));
	    plot.setOutlinePaint(Color.white);
	    NumberAxis range = (NumberAxis)plot.getRangeAxis();
        range.setTickUnit(new NumberTickUnit(50));
        NumberAxis domain = (NumberAxis)plot.getDomainAxis();
        domain.setTickUnit(new NumberTickUnit(50));
        plot.setDomainGridlinesVisible(false);
        plot.setRangeGridlinesVisible(false);
        
		try {
			SVGGraphics2D svg = new SVGGraphics2D(600, 600);
	        Rectangle area = new Rectangle(0, 0, 600, 600);
	        chart.draw(svg, area);
	        SVGUtils.writeToSVG(new File(filename + ".svg"), svg.getSVGElement());
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		ChartFrame frame = new ChartFrame(filename, chart);
		frame.pack();
		frame.setBounds(0, 0, 600, 600);
		frame.setVisible(Config.displayPlots);
	}
	
	/**
	 * Plots a verification of the gauss clusters.
	 * @param dataPoints
	 * @param filename
	 * @param sizeGoodLow1
	 * @param sizeBadHigh1
	 * @param sizeGoodLow2
	 * @param sizeBadHigh2
	 */
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
			goodlow2.add(dataPoints.get(Config.numPointsCluster + i).getGaussFeatures()[0], dataPoints.get(Config.numPointsCluster + i).getGaussFeatures()[1]);
		dataset.addSeries(goodlow2);
		XYSeries badhigh2 = new XYSeries("badhigh2");
		for (int i = 0; i < sizeBadHigh2; i++)
			badhigh2.add(dataPoints.get(Config.numPointsCluster + i + sizeGoodLow2).getGaussFeatures()[0], dataPoints.get(Config.numPointsCluster + i + sizeGoodLow2).getGaussFeatures()[1]);	
		dataset.addSeries(badhigh2);
		
		JFreeChart chart = ChartFactory.createScatterPlot("Gauss Verification", "f1", "f2", dataset, PlotOrientation.VERTICAL, true, false, false);
		XYPlot plot = (XYPlot)chart.getPlot();
		plot.setDrawingSupplier(getCustomDrawingSupplier());
		plot.setBackgroundPaint(new Color(210, 210, 210));
	    plot.setOutlinePaint(Color.white);
	    NumberAxis range = (NumberAxis)plot.getRangeAxis();
        range.setTickUnit(new NumberTickUnit(50));
        NumberAxis domain = (NumberAxis)plot.getDomainAxis();
        domain.setTickUnit(new NumberTickUnit(50));
        plot.setDomainGridlinesVisible(false);
        plot.setRangeGridlinesVisible(false);
        
		try {
			SVGGraphics2D svg = new SVGGraphics2D(600, 600);
	        Rectangle area = new Rectangle(0, 0, 600, 600);
	        chart.draw(svg, area);
	        SVGUtils.writeToSVG(new File(filename + ".svg"), svg.getSVGElement());
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		ChartFrame frame = new ChartFrame(filename, chart);
		frame.pack();
		frame.setBounds(0, 0, 600, 600);
		frame.setVisible(Config.displayPlots);
		
	}
	
	/**
	 * Plots a verification of the density clusters.
	 * @param dataPoints
	 * @param filename
	 * @param sizeGoodLow1
	 * @param sizeBadHigh1
	 * @param sizeGoodLow2
	 * @param sizeBadHigh2
	 */
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
			goodlow2.add(dataPoints.get(Config.numPointsCluster + i).getDensityFeatures()[0], dataPoints.get(Config.numPointsCluster + i).getDensityFeatures()[1]);
		dataset.addSeries(goodlow2);
		XYSeries badhigh2 = new XYSeries("badhigh2");
		for (int i = 0; i < sizeBadHigh2; i++)
			badhigh2.add(dataPoints.get(Config.numPointsCluster + i + sizeGoodLow2).getDensityFeatures()[0], dataPoints.get(Config.numPointsCluster + i + sizeGoodLow2).getDensityFeatures()[1]);	
		dataset.addSeries(badhigh2);
		
		for (int i = 1; i <= Config.numClusters; i++) {
			XYSeries series = new XYSeries(""+ i);
			for (DataPoint dataPoint : dataPoints)
				if (dataPoint.getClusterLabel().equals(i))
					series.add(dataPoint.getDensityFeatures()[0], dataPoint.getDensityFeatures()[1]);
			dataset.addSeries(series);
		}
		
		JFreeChart chart = ChartFactory.createScatterPlot("Density Verification", "\u03C1" + "1", "\u03C1" + "2", dataset, PlotOrientation.VERTICAL, true, false, false);
		XYPlot plot = (XYPlot)chart.getPlot();
		plot.setDrawingSupplier(getCustomDrawingSupplier());
		plot.setBackgroundPaint(new Color(210, 210, 210));
	    plot.setOutlinePaint(Color.white);
	    NumberAxis range = (NumberAxis)plot.getRangeAxis();
        range.setTickUnit(new NumberTickUnit(50));
        NumberAxis domain = (NumberAxis)plot.getDomainAxis();
        domain.setTickUnit(new NumberTickUnit(50));
        plot.setDomainGridlinesVisible(false);
        plot.setRangeGridlinesVisible(false);
        
		try {
			SVGGraphics2D svg = new SVGGraphics2D(600, 600);
	        Rectangle area = new Rectangle(0, 0, 600, 600);
	        chart.draw(svg, area);
	        SVGUtils.writeToSVG(new File(filename + ".svg"), svg.getSVGElement());
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		ChartFrame frame = new ChartFrame(filename, chart);
		frame.pack();
		frame.setBounds(0, 0, 600, 600);
		frame.setVisible(Config.displayPlots);
		
	}
	
	/**
	 * Plots the gauss cluster with noise.
	 * @param dataPoints1
	 * @param dataPoints2
	 * @param noisePoints
	 * @param filename
	 */
	private static void plotGaussClustersNoise(List<DataPoint> dataPoints1, List<DataPoint> dataPoints2, List<DataPoint> noisePoints, String filename) {
		
		XYSeriesCollection dataset = new XYSeriesCollection();
		XYSeries series1 = new XYSeries("1");
		XYSeries series2 = new XYSeries("2");
		XYSeries series0 = new XYSeries("Noise");
		dataset.addSeries(series1);
		dataset.addSeries(series2);
		dataset.addSeries(series0);
		
		for (DataPoint p : dataPoints1)
			series1.add(p.getGaussFeatures()[0], p.getGaussFeatures()[1]);
		for (DataPoint p : dataPoints2)
			series2.add(p.getGaussFeatures()[0], p.getGaussFeatures()[1]);
		for (DataPoint p : noisePoints)
			series0.add(p.getGaussFeatures()[0], p.getGaussFeatures()[1]);
		
		JFreeChart chart = ChartFactory.createScatterPlot("Gauss Dimensions with Noise", "f1", "f2", dataset, PlotOrientation.VERTICAL, true, false, false);
		XYPlot plot = (XYPlot)chart.getPlot();
		plot.setDrawingSupplier(getCustomDrawingSupplier());
		plot.setBackgroundPaint(new Color(210, 210, 210));
	    plot.setOutlinePaint(Color.white);
	    NumberAxis range = (NumberAxis)plot.getRangeAxis();
        range.setTickUnit(new NumberTickUnit(200));
        NumberAxis domain = (NumberAxis)plot.getDomainAxis();
        domain.setTickUnit(new NumberTickUnit(200));
        plot.setDomainGridlinesVisible(false);
        plot.setRangeGridlinesVisible(false);
        
		try {
			SVGGraphics2D svg = new SVGGraphics2D(600, 600);
	        Rectangle area = new Rectangle(0, 0, 600, 600);
	        chart.draw(svg, area);
	        SVGUtils.writeToSVG(new File(filename + ".svg"), svg.getSVGElement());
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		ChartFrame frame = new ChartFrame(filename, chart);
		frame.pack();
		frame.setBounds(0, 0, 600, 600);
		frame.setVisible(Config.displayPlots);
	}
	
	/**
	 * Plots density clusters with noise.
	 * @param dataPoints1
	 * @param dataPoints2
	 * @param noisePoints
	 * @param filename
	 */
	private static void plotDensityClustersNoise(List<DataPoint> dataPoints1, List<DataPoint> dataPoints2, List<DataPoint> noisePoints, String filename) {
		
		XYSeriesCollection dataset = new XYSeriesCollection();
		XYSeries series1 = new XYSeries("1");
		XYSeries series2 = new XYSeries("2");
		XYSeries series0 = new XYSeries("Noise");
		dataset.addSeries(series1);
		dataset.addSeries(series2);
		dataset.addSeries(series0);
		
		for (DataPoint p : dataPoints1)
			series1.add(p.getDensityFeatures()[0], p.getDensityFeatures()[1]);
		for (DataPoint p : dataPoints2)
			series2.add(p.getDensityFeatures()[0], p.getDensityFeatures()[1]);
		for (DataPoint p : noisePoints)
			series0.add(p.getDensityFeatures()[0], p.getDensityFeatures()[1]);
		
		JFreeChart chart = ChartFactory.createScatterPlot("Density Dimensions with Noise", "\u03C1" + "1", "\u03C1" + "2", dataset, PlotOrientation.VERTICAL, true, false, false);
		XYPlot plot = (XYPlot)chart.getPlot();
		plot.setDrawingSupplier(getCustomDrawingSupplier());
		plot.setBackgroundPaint(new Color(210, 210, 210));
	    plot.setOutlinePaint(Color.white);
	    NumberAxis range = (NumberAxis)plot.getRangeAxis();
        range.setTickUnit(new NumberTickUnit(200));
        NumberAxis domain = (NumberAxis)plot.getDomainAxis();
        domain.setTickUnit(new NumberTickUnit(200));
        plot.setDomainGridlinesVisible(false);
        plot.setRangeGridlinesVisible(false);
        
		try {
			SVGGraphics2D svg = new SVGGraphics2D(600, 600);
	        Rectangle area = new Rectangle(0, 0, 600, 600);
	        chart.draw(svg, area);
	        SVGUtils.writeToSVG(new File(filename + ".svg"), svg.getSVGElement());
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		ChartFrame frame = new ChartFrame(filename, chart);
		frame.pack();
		frame.setBounds(0, 0, 600, 600);
		frame.setVisible(Config.displayPlots);
	}
	
	/**
	 * Generates verificatoin plots for easy and hard to cluster datapoints.
	 * 
	 * @param allDataPoints
	 * @param dir The output directory
	 * @param goodIds1 List of all easy to cluster gauss Ids.
	 * @param badIds1 List of all hard to cluster gauss Ids.
	 * @param goodIds2 List of all easy to cluster gauss Ids.
	 * @param badIds2 List of all hard to cluster gauss Ids.
	 * @param lowIds1 List of all hard to cluster density Ids.
	 * @param highIds1 List of all easy to cluster gauss Ids.
	 * @param lowIds2 List of all hard to cluster density Ids.
	 * @param highIds2 List of all easy to cluster gauss Ids.
	 */
	private static void verifyDataset(List<DataPoint> allDataPoints, String dir, List<Integer> goodIds1, List<Integer> badIds1, List<Integer> goodIds2, List<Integer> badIds2, List<Integer> lowIds1, List<Integer> highIds1, List<Integer> lowIds2, List<Integer> highIds2) {
		
		int sizeGoodLow1 = goodIds1.size() < lowIds1.size() ? goodIds1.size() : lowIds1.size();
		int sizeBadHigh1 = badIds1.size() < highIds1.size() ? badIds1.size() : highIds1.size();
		int sizeGoodLow2 = goodIds2.size() < lowIds2.size() ? goodIds2.size() : lowIds2.size();
		int sizeBadHigh2 = badIds2.size() < highIds2.size() ? badIds2.size() : highIds2.size();
		
		plotVerifyGauss(allDataPoints, Config.workDir + "/" + dir + "/verify_gauss", sizeGoodLow1, sizeBadHigh1, sizeGoodLow2, sizeBadHigh2);
		plotVerifyDensity(allDataPoints, Config.workDir + "/" + dir + "/verify_density", sizeGoodLow1, sizeBadHigh1, sizeGoodLow2, sizeBadHigh2);
	}
	
	/** 
	 * Color palette for plotting.
	 *
	 */
	private static DefaultDrawingSupplier getCustomDrawingSupplier() {
		return new DefaultDrawingSupplier(
				new Paint[] { 
						new Color(31,120,180),
						new Color(227,26,28),
						new Color(51,160,44),
						new Color(106,61,154),
						new Color(177,89,40),
						new Color(255,127,0),
						new Color(178,223,138),
						new Color(251,154,153),
						new Color(166,206,227),
						new Color(253,191,111),
						new Color(202,178,214),
						new Color(255,255,153),
				},
				DefaultDrawingSupplier.DEFAULT_OUTLINE_PAINT_SEQUENCE,
				DefaultDrawingSupplier.DEFAULT_STROKE_SEQUENCE,
				DefaultDrawingSupplier.DEFAULT_OUTLINE_STROKE_SEQUENCE,
				DefaultDrawingSupplier.DEFAULT_SHAPE_SEQUENCE);
	}
}
