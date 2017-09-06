package p2.clustering;

/**
 * Copyright (c) 2017 Alexander Pfundner
 * 
 * Integration of Density-based and Partitioning-based Clustering Methods
 * 
 */

import java.awt.Color;
import java.awt.Paint;
import java.awt.Rectangle;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.plot.DefaultDrawingSupplier;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.graphics2d.svg.SVGGraphics2D;
import org.jfree.graphics2d.svg.SVGUtils;

import de.lmu.ifi.dbs.elki.algorithm.KNNDistancesSampler;
import de.lmu.ifi.dbs.elki.algorithm.KNNDistancesSampler.KNNDistanceOrderResult;
import de.lmu.ifi.dbs.elki.algorithm.clustering.DBSCAN;
import de.lmu.ifi.dbs.elki.algorithm.clustering.kmeans.KMeansMacQueen;
import de.lmu.ifi.dbs.elki.data.Cluster;
import de.lmu.ifi.dbs.elki.data.Clustering;
import de.lmu.ifi.dbs.elki.data.DoubleVector;
import de.lmu.ifi.dbs.elki.data.NumberVector;
import de.lmu.ifi.dbs.elki.data.model.KMeansModel;
import de.lmu.ifi.dbs.elki.data.model.Model;
import de.lmu.ifi.dbs.elki.data.type.TypeUtil;
import de.lmu.ifi.dbs.elki.database.AbstractDatabase;
import de.lmu.ifi.dbs.elki.database.Database;
import de.lmu.ifi.dbs.elki.database.StaticArrayDatabase;
import de.lmu.ifi.dbs.elki.database.ids.DBIDIter;
import de.lmu.ifi.dbs.elki.database.ids.DBIDRange;
import de.lmu.ifi.dbs.elki.database.relation.Relation;
import de.lmu.ifi.dbs.elki.datasource.ArrayAdapterDatabaseConnection;
import de.lmu.ifi.dbs.elki.distance.distancefunction.minkowski.EuclideanDistanceFunction;
import de.lmu.ifi.dbs.elki.index.tree.spatial.rstarvariants.rstar.RStarTreeFactory;
import de.lmu.ifi.dbs.elki.math.geometry.XYCurve;
import de.lmu.ifi.dbs.elki.utilities.ClassGenericsUtil;
import de.lmu.ifi.dbs.elki.utilities.optionhandling.parameterization.ListParameterization;
import p2.util.Config;
import p2.util.DataLoaderUtil;
import p2.util.StatisticsUtil;

/**
 * Processes generated datasets with an integrative workflow of KMeans and DBSCAN.
 * 
 */
public class IntegratedRunner {
	
	/**
	 * Number of clusters
	 */
	private int numClusters = Config.numClusters;
	
	/**
	 * Data points
	 */
	private double[][] data;
	
	/**
	 * Gauss part of the data
	 */
	private double[][] dataGauss;
	
	/**
	 * Density part of the data
	 */
	private double[][] dataDensity;
	
	/**
	 * Main method, processing dataset 1 and 2.
	 * Does not take arguments, parameters are set via class Config.
	 * 
	 */
	public static void main(String[] args) {
	
		System.out.println("Processing dataset 1 ...");
		runWorkflow("dataset_1", 25, 1, 1);
		
		System.out.println("Processing dataset 2 ...");
		runWorkflow("dataset_2", 22, 1, 0.85);
		
		System.out.println("\nFinished.");
	}
	
	/**
	 * Runs the integrative workflow for a dataset.
	 * 
	 * @param dir The directory wher the data is.
	 * @param epsilonFull DBSCAN epsilon parameter for the full feature space.
	 * @param epsilonDensity DBSCAN epsilon parameter for the density feature space.
	 * @param epsilonDummy DBSCAN epsilon parameter for the extended feature space.
	 */
	public static void runWorkflow(String dir, double epsilonFull, double epsilonDensity, double epsilonDummy) {
		
		IntegratedRunner runner = new IntegratedRunner(Config.workDir + "/" + dir +  "/data.csv");
		
		int minPts = 4;
		
		//runner.estimateDBSCANEpsilon(runner.data, minPts, "full");
		//runner.estimateDBSCANEpsilon(runner.dataDensity, minPts, "density");
		
		// 
		// Initial performance analysis for full and partial data
		// ------------------------------------------------------------
		//
		
		//Full data
		ExtKMeansClustering kmeansFull = runner.runKMeans(runner.data, "kmeans_full", dir); 		
		ExtDBSCANClustering dbscanFull = runner.runDBSCAN(runner.data, minPts, epsilonFull, "dbscan_full", dir);
		
		//Gauss data
		ExtKMeansClustering kmeansGauss = runner.runKMeans(runner.dataGauss, "kmeans_gauss", dir); 
		
		//Density data	
		ExtDBSCANClustering dbscanDensity = runner.runDBSCAN(runner.dataDensity, minPts, epsilonDensity, "dbscan_density", dir);
		
		//Plot results
		if (Config.plotClusters) {
			plotKMeansClustering(kmeansGauss.getClustering(), runner.dataGauss, kmeansGauss.getIds(), Config.workDir + "/" + dir + "/kmeans_gauss");
			plotDBSCANClustering(dbscanDensity.getClustering(), runner.dataDensity, dbscanDensity.getIds(), Config.workDir + "/" + dir + "/dbscan_density");
		}

		
		//
		// Running DBSCAN and KMeans alternately
		// --------------------------------------------
		//
		 
		ExtKMeansClustering kmeans1 = runner.runKMeans(runner.dataGauss, "kmeans1", dir);
		//runner.estimateDBSCANEpsilon(runner.getExtData(runner.dataDensity, kmeans1.getDummy()), minPts, "dummy1");
		ExtDBSCANClustering dbscan1 = runner.runDBSCAN(runner.getExtData(runner.dataDensity, kmeans1.getDummy()), minPts, epsilonDummy, "dbscan1", dir);
		
		ExtKMeansClustering kmeans2 = runner.runKMeans(runner.getExtData(runner.dataGauss, dbscan1.getDummy()), "kmeans2", dir);
		ExtDBSCANClustering dbscan2 = runner.runDBSCAN(runner.getExtData(runner.dataDensity, kmeans2.getDummy()), minPts, epsilonDummy, "dbscan2", dir);
		
		ExtKMeansClustering kmeans3 = runner.runKMeans(runner.getExtData(runner.dataGauss, dbscan2.getDummy()), "kmeans3", dir);
		ExtDBSCANClustering dbscan3 = runner.runDBSCAN(runner.getExtData(runner.dataDensity, kmeans3.getDummy()), minPts, epsilonDummy, "dbscan3", dir);
		
		ExtKMeansClustering kmeans4 = runner.runKMeans(runner.getExtData(runner.dataGauss, dbscan3.getDummy()), "kmeans4", dir);
		ExtDBSCANClustering dbscan4 = runner.runDBSCAN(runner.getExtData(runner.dataDensity, kmeans4.getDummy()), minPts, epsilonDummy, "dbscan4", dir);
		
		ExtKMeansClustering kmeans5 = runner.runKMeans(runner.getExtData(runner.dataGauss, dbscan4.getDummy()), "kmeans5", dir);
		ExtDBSCANClustering dbscan5 = runner.runDBSCAN(runner.getExtData(runner.dataDensity, kmeans5.getDummy()), minPts, epsilonDummy, "dbscan5", dir);
		
		ExtKMeansClustering kmeans6 = runner.runKMeans(runner.getExtData(runner.dataGauss, dbscan5.getDummy()), "kmeans6", dir);
		ExtDBSCANClustering dbscan6 = runner.runDBSCAN(runner.getExtData(runner.dataDensity, kmeans6.getDummy()), minPts, epsilonDummy, "dbscan6", dir);
		
		ExtKMeansClustering kmeans7 = runner.runKMeans(runner.getExtData(runner.dataGauss, dbscan6.getDummy()), "kmeans7", dir);
		ExtDBSCANClustering dbscan7 = runner.runDBSCAN(runner.getExtData(runner.dataDensity, kmeans7.getDummy()), minPts, epsilonDummy, "dbscan7", dir);
		
		ExtKMeansClustering kmeans8 = runner.runKMeans(runner.getExtData(runner.dataGauss, dbscan7.getDummy()), "kmeans8", dir);
		ExtDBSCANClustering dbscan8 = runner.runDBSCAN(runner.getExtData(runner.dataDensity, kmeans8.getDummy()), minPts, epsilonDummy, "dbscan8", dir);
		
		ExtKMeansClustering kmeans9 = runner.runKMeans(runner.getExtData(runner.dataGauss, dbscan8.getDummy()), "kmeans9", dir);
		ExtDBSCANClustering dbscan9 = runner.runDBSCAN(runner.getExtData(runner.dataDensity, kmeans9.getDummy()), minPts, epsilonDummy, "dbscan9", dir);
		
		ExtKMeansClustering kmeans10 = runner.runKMeans(runner.getExtData(runner.dataGauss, dbscan9.getDummy()), "kmeans10", dir);
		ExtDBSCANClustering dbscan10 = runner.runDBSCAN(runner.getExtData(runner.dataDensity, kmeans10.getDummy()), minPts, epsilonDummy, "dbscan10", dir);
	}
	
	/**
	 * Constructor, loads data and header information for a given file.
	 * @param file The full filename.
	 */
	public IntegratedRunner(String file) {

		DataLoaderUtil dataUtil = new DataLoaderUtil(file);	
		data = dataUtil.loadData();
		dataGauss = dataUtil.getGaussData();
		dataDensity = dataUtil.getDensityData();
	}
	
	/**
	 * Runs KMeans clustering using data extended by DBSCAN cluster labels in dummy encoding.
	 * @param iterLabel The name of the iteration
	 * @param dir The output directory.
	 * @return A clustering result with additional dummy encoding.
	 */
	private ExtKMeansClustering runKMeans(double[][] data, String iterLabel, String dir) {
		
		System.out.println("Running KMeans ...");
		
		ArrayAdapterDatabaseConnection conn = new ArrayAdapterDatabaseConnection(data);
		ListParameterization dbParams = new ListParameterization();
		dbParams.addParameter(AbstractDatabase.Parameterizer.DATABASE_CONNECTION_ID, conn);
		Database db = ClassGenericsUtil.parameterizeOrAbort(StaticArrayDatabase.class, dbParams);
		db.initialize();
		
		Relation<NumberVector> rel = db.getRelation(TypeUtil.NUMBER_VECTOR_FIELD);
	    DBIDRange ids = (DBIDRange) rel.getDBIDs();
		
	    ListParameterization kmeansParams = new ListParameterization();
		kmeansParams.addParameter(KMeansMacQueen.K_ID, numClusters);
		kmeansParams.addParameter(KMeansMacQueen.DISTANCE_FUNCTION_ID, EuclideanDistanceFunction.class);
		
		//Run KMeans
		KMeansMacQueen<DoubleVector> kmeans = ClassGenericsUtil.parameterizeOrAbort(KMeansMacQueen.class, kmeansParams);
		Clustering<KMeansModel> clu = kmeans.run(db);
		
		//Generate dummy encoding
		double[][] dummy = new double[Config.numPoints][Config.numClusters];		
		int clusterID = 0;
		for (Cluster<KMeansModel> c : clu.getAllClusters()) {
			for (DBIDIter it = c.getIDs().iter(); it.valid(); it.advance())
				dummy[ids.getOffset(it)][clusterID] = 1.0;
			++clusterID;
		}
		
		//Save clustering to file
		saveClusteringToFile(Config.workDir + "/" + dir + "/" + iterLabel + ".csv", clu, ids);
		
		return new ExtKMeansClustering(clu, dummy, ids);
	}
	
	/**
	 * Runs DBSCAN a single time for given minPts and epsilon.
	 * @param data The data to process.
	 * @param minPts The minPts parameter of the DBSCAN algorithm.
	 * @param epsilon The epsilon parameter of the DBSCAN algorithm.
	 * @param iterLabel The name of the iteration
	 * @param dir The output directory.
	 * @return A clustering result with additional dummy encoding.
	 */
	private ExtDBSCANClustering runDBSCAN(double[][] data, int minPts, double epsilon, String iterLabel, String dir) {
		
		System.out.println("Running DBSCAN ...");
		
		ArrayAdapterDatabaseConnection conn = new ArrayAdapterDatabaseConnection(data);
		ListParameterization dbParams = new ListParameterization();
		dbParams.addParameter(AbstractDatabase.Parameterizer.DATABASE_CONNECTION_ID, conn);
		dbParams.addParameter(AbstractDatabase.Parameterizer.INDEX_ID, RStarTreeFactory.class);
		Database db = ClassGenericsUtil.parameterizeOrAbort(StaticArrayDatabase.class, dbParams);
		db.initialize();
		
		Relation<DoubleVector> rel = db.getRelation(TypeUtil.NUMBER_VECTOR_FIELD);
	    DBIDRange ids = (DBIDRange) rel.getDBIDs();
		
		ListParameterization dbscanParams = new ListParameterization();
		dbscanParams.addParameter(DBSCAN.Parameterizer.EPSILON_ID, epsilon);
		dbscanParams.addParameter(DBSCAN.Parameterizer.MINPTS_ID, minPts);
		dbscanParams.addParameter(DBSCAN.DISTANCE_FUNCTION_ID, EuclideanDistanceFunction.class);
		
		//Run DBSCAN
		DBSCAN<DoubleVector> dbscan = ClassGenericsUtil.parameterizeOrAbort(DBSCAN.class, dbscanParams);
		Clustering<Model> clu = dbscan.run(db);
		
		//Strip 0-element clusters		
		Clustering<Model> stripped = new Clustering<Model>("DBSCAN Clustering", "dbscan-clustering");
		for (Cluster<Model> c : clu.getAllClusters())
			if (c.size() > 0)
				stripped.addToplevelCluster(c);
		clu = stripped;
		
		//Generate dummy encoding
		double[][] dummy = new double[Config.numPoints][clu.getAllClusters().size()];
		int clusterID = 0;
		for (Cluster<Model> c : clu.getAllClusters()) {
			for (DBIDIter it = c.getIDs().iter(); it.valid(); it.advance()) 
				dummy[ids.getOffset(it)][clusterID] = 1.0;
			++clusterID;
		}
		
		//Save clustering to file
		saveClusteringToFile(Config.workDir + "/" + dir + "/" + iterLabel + ".csv", clu, ids);
		
		return new ExtDBSCANClustering(clu, dummy, ids);
	}
	
	/**
	 * Exports a clustering to a file.
	 * Use for evaluation with Matlab implementation of Adjusted Mutual Information.
	 * E.g. "1 1 1 2 2" would say objects 1-3 are in cluster 1 and objects 4-5 are in cluster 2.
	 * @param filename The file to write to.
	 * @param clu The clustering to export to the file.
	 * @param ids The ids of the loaded dataset.
	 */
	private void saveClusteringToFile(String filename, Clustering<? extends Model> clu, DBIDRange ids) {
		
		int[] labels = new int[Config.numPoints];

		int cluId = 1;
		for (Cluster<? extends Model> c : clu.getAllClusters()) { 
			for (DBIDIter it = c.getIDs().iter(); it.valid(); it.advance())
				labels[ids.getOffset(it)] = cluId;
			++cluId;
		}
		StringBuilder sb = new StringBuilder();			
		for (int l : labels)
			sb.append(l + " ");
		
		try {
			FileWriter writer = new FileWriter(new File(filename));
			writer.write(sb.toString());
			writer.flush();
			writer.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Concats data with the normalized dummy encoding.
	 * 
	 * @param data The data as a double array.
	 * @param dummy The dummy encoded cluster labels.
	 * @return Concat between data and dummy encoding.
	 */
	public double[][] getExtData(double[][] data, double[][] dummy) {
		
		if (data.length != dummy.length)
			throw new IllegalStateException("data.length != dummy.length");
		
		int numRows = data.length;
		int numColsData = data[0].length;
		int numColsDummy = dummy[0].length;
		
		//Normalization of dummy
		StatisticsUtil statUtil = new StatisticsUtil();
		double sumVarData = statUtil.calcSumVar(data, numColsData, numRows);
		double sumVarDummy = statUtil.calcSumVar(dummy, numColsDummy, numRows);
		double factor = Math.sqrt(sumVarData / sumVarDummy);
		for (int row = 0; row < numRows; row++)
			for (int col = 0; col < numColsDummy; col++)
				dummy[row][col] *= factor;

		//Concat data and dummy
		int numColsConcat = data[0].length + dummy[0].length;
		double[][] extData = new double[numRows][numColsConcat];
		for (int row = 0; row < numRows; row++)
			for (int col = 0; col < numColsData; col++)
				extData[row][col] = data[row][col];
		for (int row = 0; row < numRows; row++)
			for (int col = numColsData; col < numColsConcat; col++)
				extData[row][col] = dummy[row][col - numColsData];		
		
		return extData;
	}
	
	/**
	 * Plots clustering obtained by KMeans.
	 * @param clustering
	 * @param plotData
	 * @param ids
	 * @param filename
	 */
	private static void plotKMeansClustering(Clustering<KMeansModel> clustering, double[][] plotData, DBIDRange ids, String filename) {
		
		XYSeriesCollection collection = new XYSeriesCollection();
		for (int i = 0; i < clustering.getAllClusters().size(); i++) {
			Cluster<KMeansModel> cluster = clustering.getAllClusters().get(i);
			XYSeries series = new XYSeries(i+1);
			for (DBIDIter it = cluster.getIDs().iter(); it.valid(); it.advance()) {
				int objectID = ids.getOffset(it);
				double[] row = plotData[objectID];
				series.add(row[0], row[1]);
			}
			collection.addSeries(series);
		}
		
		JFreeChart chart = ChartFactory.createScatterPlot("k-Means on Gauss Dimensions", "d1", "d2", collection, PlotOrientation.VERTICAL, true, false, false);
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
	 * Plots clustering obtained by DBSCAN.
	 * @param clustering
	 * @param plotData
	 * @param ids
	 * @param filename
	 */
	private static void plotDBSCANClustering(Clustering<Model> clustering, double[][] plotData, DBIDRange ids, String filename) {
		
		XYSeriesCollection collection = new XYSeriesCollection();
		for (int i = 0; i < clustering.getAllClusters().size(); i++) {
			Cluster<Model> cluster = clustering.getAllClusters().get(i);
			XYSeries series = new XYSeries(i+1);
			for (DBIDIter it = cluster.getIDs().iter(); it.valid(); it.advance()) {
				int objectID = ids.getOffset(it);
				double[] row = plotData[objectID];
				series.add(row[0], row[1]);
			}
			collection.addSeries(series);
		}
		
		JFreeChart chart = ChartFactory.createScatterPlot("DBSCAN on Density Dimensions", "d1", "d2", collection, PlotOrientation.VERTICAL, true, false, false);
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
	 * Plots the KNN-distances for estimating DBSCAN epsilon parameter.
	 * @param data
	 * @param minPts
	 * @param label
	 */
	public void estimateDBSCANEpsilon(double[][] data, int minPts, String label) {
		
		ArrayAdapterDatabaseConnection conn = new ArrayAdapterDatabaseConnection(data);
		ListParameterization dbParams = new ListParameterization();
		dbParams.addParameter(AbstractDatabase.Parameterizer.DATABASE_CONNECTION_ID, conn);
		dbParams.addParameter(AbstractDatabase.Parameterizer.INDEX_ID, RStarTreeFactory.class);
		Database db = ClassGenericsUtil.parameterizeOrAbort(StaticArrayDatabase.class, dbParams);
		db.initialize();
		
		Relation<DoubleVector> rel = db.getRelation(TypeUtil.NUMBER_VECTOR_FIELD);
		
		ListParameterization paramsKNN = new ListParameterization();
		paramsKNN.addParameter(KNNDistancesSampler.DISTANCE_FUNCTION_ID, EuclideanDistanceFunction.class);
		paramsKNN.addParameter(KNNDistancesSampler.Parameterizer.K_ID, minPts - 1);
		paramsKNN.addParameter(KNNDistancesSampler.Parameterizer.SAMPLING_ID, 1.0);
		paramsKNN.addParameter(KNNDistancesSampler.Parameterizer.SEED_ID, System.currentTimeMillis());
		
		KNNDistancesSampler<DoubleVector> knn = ClassGenericsUtil.parameterizeOrAbort(KNNDistancesSampler.class, paramsKNN);
		KNNDistanceOrderResult result = knn.run(db, rel);
		
		XYSeriesCollection coll = new XYSeriesCollection();
		XYSeries series = new XYSeries(label);
		coll.addSeries(series);
		for (XYCurve.Itr it = result.iterator(); it.valid(); it.advance()) 
			series.add(it.getX(), it.getY());
		JFreeChart chart = ChartFactory.createXYLineChart("", "", "", coll, PlotOrientation.VERTICAL, true, false, false);
		ChartFrame frame = new ChartFrame("KNN Distances", chart);
		frame.pack();
		frame.setVisible(true);

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
