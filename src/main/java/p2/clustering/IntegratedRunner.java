package p2.clustering;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.spi.NumberFormatProvider;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

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
import p2.util.DataLoaderUtil;
import p2.util.StatisticsUtil;

public class IntegratedRunner {

	/**
	 * Number of dimensions for each cluster type.
	 */
	private int numDimPerType;
	
	/**
	 * Number of clusters.
	 */
	private int numClusters;
	
	/**
	 * Number of points per cluster. 
	 */
	private int numPointsCluster;
	
	/**
	 * Total number of points. 
	 */
	private int numPoints;
	
	/**
	 * Data points
	 */
	private double[][] data;
	
	/**
	 * Gaussian part of the data.
	 */
	private double[][] dataGauss;
	
	/**
	 * Density part of the data.
	 */
	private double[][] dataDensity;
	
	/**
	 * Working directory
	 */
	private final static String workDir = "/home/alepfu/Desktop/P2";
	
	
	public static void main(String[] args) {
		
		IntegratedRunner runner = new IntegratedRunner(workDir + "/data.csv");
		
		/**
		 * 
		 * 
		 * Initial performance analysis for full and partial data
		 * =======================================================
		 * 
		 */
		//KMeans runs
		runner.runKMeans(runner.data, "kmeansfull");
		runner.runDBSCAN(runner.data, 10, 60, true, "dbscanfull");
		
//		saveClusteringToFile(workDir + "/kmeans_full.csv", kmeans1.getClustering());
//		ExtKMeansClustering kmeans2 = runner.runKMeans(runner.dataGauss);
//		saveClusteringToFile(workDir + "/kmeans_gauss.csv", kmeans2.getClustering());
//		plotKMeansClustering(kmeans2.getClustering(), runner.dataGauss, kmeans2.getIds(), workDir + "/kmeans_gauss.jpeg");
//		ExtKMeansClustering kmeans3 = runner.runKMeans(runner.dataDensity);
//		saveClusteringToFile(workDir + "/kmeans_density.csv", kmeans3.getClustering());
//		plotKMeansClustering(kmeans3.getClustering(), runner.dataDensity, kmeans3.getIds(), workDir + "/kmeans_density.jpeg");
		
		//DBSCAN runs
//		ExtDBSCANClustering dbscan1 = runner.runMultipleDBSCAN(runner.data, 40.0, 1.0, 5);
//		saveClusteringToFile(workDir + "/dbscan_full.csv", dbscan1.getClustering());
//		ExtDBSCANClustering dbscan2 = runner.runSingleDBSCAN(runner.dataGauss, 5, 2.0);
//		saveClusteringToFile(workDir + "/dbscan_gauss.csv", dbscan2.getClustering());
//		plotDBSCANClustering(dbscan2.getClustering(), runner.dataGauss, dbscan2.getIds(), workDir + "/dbscan_gauss.jpeg");
//		ExtDBSCANClustering dbscan3 = runner.runSingleDBSCAN(runner.dataDensity, 10, 1.75);
//		saveClusteringToFile(workDir + "/dbscan_density.csv", dbscan3.getClustering());
//		plotDBSCANClustering(dbscan3.getClustering(), runner.dataDensity, dbscan3.getIds(), workDir + "/dbscan_density.jpeg");
		
		
		/**
		 * 
		 * 
		 * Running DBSCAN and KMeans alternately
		 * =====================================
		 * 
		 */
		int minPts = 4;
		double epsilon = 7.5;
		boolean doEpsilonEstimation = false;
		
		ExtKMeansClustering kmeans1 = runner.runKMeans(runner.dataGauss, "kmeans1");
		
		ExtDBSCANClustering dbscan1 = runner.runDBSCAN(runner.getExtData(runner.dataDensity, kmeans1.getDummy()), minPts, epsilon, doEpsilonEstimation, "dbscan1");
		
		ExtKMeansClustering kmeans2 = runner.runKMeans(runner.getExtData(runner.dataGauss, dbscan1.getDummy()), "kmeans2");
		
		ExtDBSCANClustering dbscan2 = runner.runDBSCAN(runner.getExtData(runner.dataDensity, kmeans2.getDummy()), minPts, epsilon, doEpsilonEstimation, "dbscan2");
		
		ExtKMeansClustering kmeans3 = runner.runKMeans(runner.getExtData(runner.dataGauss, dbscan2.getDummy()), "kmeans3");
		
		ExtDBSCANClustering dbscan3 = runner.runDBSCAN(runner.getExtData(runner.dataDensity, kmeans3.getDummy()), minPts, epsilon, doEpsilonEstimation, "dbscan3");
		
		
		plotKMeansClustering(kmeans3.getClustering(), runner.dataGauss, kmeans3.getIds(), workDir + "/kmeans_final.jpeg");
		plotDBSCANClustering(dbscan3.getClustering(), runner.dataDensity, dbscan3.getIds(), workDir + "/dbscan_final.jpeg");
		
		
		/*int numRuns = 8;

		ExtKMeansClustering kmeans = null;
		ExtDBSCANClustering dbscan = null;
		
		int minPts = 10;
		boolean doEpsilonEstimation = false;
		
		for (int r = 1; r <= numRuns; r++) {
	
			//KMeans on gauss features with dummy encoding (for r > 1) from DBSCAN
			kmeans = runner.runKMeans(r == 1 ? runner.dataGauss : runner.getExtData(runner.dataGauss, dbscan.getDummy()));
			saveClusteringToFile(workDir + "/run_" + (r++) + ".csv", kmeans.getClustering());
			
			//DBSCAN on density features with dummy encoding from KMeans
			dbscan = runner.runSingleDBSCAN(runner.getExtData(runner.dataDensity, kmeans.getDummy()), minPts, 16, doEpsilonEstimation);
			saveClusteringToFile(workDir + "/run_" + r + ".csv", dbscan.getClustering());
			
		}*/
		


		
		System.out.println("\nFinished.");
	}
	
	/**
	 * Constructor
	 * Loads data and header information for a given file, generates ground truth clustering.
	 * @param file The full filename.
	 */
	public IntegratedRunner(String file) {

		System.out.println("Loading data ...");
		DataLoaderUtil dataUtil = new DataLoaderUtil(file);	
		numDimPerType = dataUtil.getNumDimPerType();
		numClusters = dataUtil.getNumClusters();
		numPointsCluster = dataUtil.getNumPointsCluster();
		numPoints = numClusters * numPointsCluster;
		data = dataUtil.loadData();
		dataGauss = dataUtil.getGaussianData();
		dataDensity = dataUtil.getDensityData();
	}
	
	/**
	 * Runs KMeans clustering using data extended by DBSCAN cluster labels in dummy encoding.
	 * @return A clustering result with additional dummy encoding.
	 */
	private ExtKMeansClustering runKMeans(double[][] data, String iterLabel) {
		
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
		double[][] dummy = new double[numPoints][numClusters];		
		int clusterID = 0;
		for (Cluster<KMeansModel> c : clu.getAllClusters()) {
			for (DBIDIter it = c.getIDs().iter(); it.valid(); it.advance())
				dummy[ids.getOffset(it)][clusterID] = 1.0;
			++clusterID;
		}
		
		//Save clustering to file
		saveClusteringToFile(workDir + "/" + iterLabel + ".csv", clu, ids);
		
		return new ExtKMeansClustering(clu, dummy, ids);
	}
	
	/**
	 * Runs DBSCAN a single time for given minPts and epsilon.
	 * @param data The data to process.
	 * @param minPts The minPts parameter of the DBSCAN algorithm.
	 * @param epsilon The epsilon parameter of the DBSCAN algorithm.
	 * @return A clustering result with additional dummy encoding.
	 */
	private ExtDBSCANClustering runDBSCAN(double[][] data, int minPts, double epsilon, boolean doEpsilonEstimation, String iterLabel) {
		
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
		
		//Sample some KNN-distances for estimating a good epsilon
		if (doEpsilonEstimation)
			estimateDBSCANEpsilonParameter(minPts, db, rel);
		
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
		double[][] dummy = new double[numPoints][clu.getAllClusters().size()];
		int clusterID = 0;
		for (Cluster<Model> c : clu.getAllClusters()) {
			for (DBIDIter it = c.getIDs().iter(); it.valid(); it.advance()) 
				dummy[ids.getOffset(it)][clusterID] = 1.0;
			++clusterID;
		}
		
		//Save clustering to file
		saveClusteringToFile(workDir + "/" + iterLabel + ".csv", clu, ids);
		
		return new ExtDBSCANClustering(clu, dummy, ids);
	}
	
	/**
	 * Exports a clustering to a file.
	 * Use for evaluation with Matlab implementation of Adjusted Mutual Information.
	 * E.g. "1 1 1 2 2" would say objects 1-3 are in cluster 1 and objects 4-5 are in cluster 2.
	 * @param filename The file to write to.
	 * @param clustering The clustering to export to the file.
	 */
	private void saveClusteringToFile(String filename, Clustering<? extends Model> clu, DBIDRange ids) {
		
		int[] labels = new int[numPoints];

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
		
		//TODO make logging configurable
		
		//Log the result
		/*NumberFormat nf = NumberFormat.getInstance();
		nf.setMinimumFractionDigits(2);
		nf.setMaximumFractionDigits(2);
		StringBuilder log = new StringBuilder("\nExtended data:\n");
		for (int row = 0; row < 3; row++) {
			for (int col = 0; col < numColsConcat; col++)
				log.append(nf.format(extData[row][col]) + " ");
			log.append("\n");
		}
		log.append("...\n");
		System.out.println(log);*/
		
		return extData;
	}
	
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
		
		JFreeChart chartDensity = ChartFactory.createScatterPlot("", "", "", collection, 
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
		
		JFreeChart chartDensity = ChartFactory.createScatterPlot("", "", "", collection, 
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
	
	public double estimateDBSCANEpsilonParameter(int minPts, Database db, Relation<DoubleVector> rel) {
		
		ListParameterization params = new ListParameterization();
		params.addParameter(KNNDistancesSampler.DISTANCE_FUNCTION_ID, EuclideanDistanceFunction.class);
		params.addParameter(KNNDistancesSampler.Parameterizer.K_ID, minPts - 1);
		params.addParameter(KNNDistancesSampler.Parameterizer.SAMPLING_ID, 0.2);
		params.addParameter(KNNDistancesSampler.Parameterizer.SEED_ID, 1234);
		
		KNNDistancesSampler<DoubleVector> knn = ClassGenericsUtil.parameterizeOrAbort(KNNDistancesSampler.class, params);
		KNNDistanceOrderResult result = knn.run(db, rel);
		
		
		XYSeriesCollection coll = new XYSeriesCollection();
		XYSeries series = new XYSeries("knn distances");
		coll.addSeries(series);
		for (XYCurve.Itr it = result.iterator(); it.valid(); it.advance()) 
			series.add(it.getX(), it.getY());
		JFreeChart chart = ChartFactory.createXYLineChart("", "", "", coll, PlotOrientation.VERTICAL, true, false, false);
		ChartFrame frame = new ChartFrame("knn distances", chart);
		frame.pack();
		frame.setVisible(true);
			 
		
		return 0;   //TODO work in progress
	}
}
