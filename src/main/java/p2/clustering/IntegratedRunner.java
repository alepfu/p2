package p2.clustering;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

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
import p2.util.Config;
import p2.util.DataLoaderUtil;
import p2.util.StatisticsUtil;

public class IntegratedRunner {
	
	/**
	 * Number of clusters
	 */
	private int numClusters = Config.numClusters;
	
	/**
	 * Total number of points
	 */
	private int numPoints = Config.numPointsCluster * Config.numClusters;
	
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
	 * Working directory
	 */
	private final static String workDir = Config.workDir;
	
	/**
	 * Diplay plots
	 */
	private static boolean showPlots = false;
	
	
	public static void main(String[] args) {
		
		IntegratedRunner runner = new IntegratedRunner(workDir + "/data.csv");
		
		int minPts = 4;
		runner.estimateDBSCANEpsilon(runner.data, minPts, "full");
		runner.estimateDBSCANEpsilon(runner.dataGauss, minPts, "gauss");
		runner.estimateDBSCANEpsilon(runner.dataDensity, minPts, "density");
		
		double epsilonFull = 45;
		double epsilonGauss= 40;
		double epsilonDensity = 4.5;
		
		/**
		 * 
		 * 
		 * Initial performance analysis for full and partial data
		 * =======================================================
		 * 
		 */
		
		//Full data
		ExtKMeansClustering kmeansFull = runner.runKMeans(runner.data, "kmeans_full"); 
		plotKMeansClustering(kmeansFull.getClustering(), runner.data, kmeansFull.getIds(), workDir + "/kmeans_full.jpeg");		
		ExtDBSCANClustering dbscanFull = runner.runDBSCAN(runner.data, minPts, epsilonFull, "dbscan_full");
		plotDBSCANClustering(dbscanFull.getClustering(), runner.data, dbscanFull.getIds(), workDir + "/dbscan_full.jpeg");
		
		//Gauss data
		ExtKMeansClustering kmeansGauss = runner.runKMeans(runner.dataGauss, "kmeans_gauss"); 
		plotKMeansClustering(kmeansGauss.getClustering(), runner.dataGauss, kmeansGauss.getIds(), workDir + "/kmeans_gauss.jpeg");		
		ExtDBSCANClustering dbscanGauss = runner.runDBSCAN(runner.dataGauss, minPts, epsilonGauss, "dbscan_gauss");
		plotDBSCANClustering(dbscanGauss.getClustering(), runner.dataGauss, dbscanGauss.getIds(), workDir + "/dbscan_gauss.jpeg");
		
		//Density data
		ExtKMeansClustering kmeansDensity = runner.runKMeans(runner.dataDensity, "kmeans_density"); 
		plotKMeansClustering(kmeansDensity.getClustering(), runner.dataDensity, kmeansDensity.getIds(), workDir + "/kmeans_density.jpeg");		
		ExtDBSCANClustering dbscanDensity = runner.runDBSCAN(runner.dataDensity, minPts, epsilonDensity, "dbscan_density");
		plotDBSCANClustering(dbscanDensity.getClustering(), runner.dataDensity, dbscanDensity.getIds(), workDir + "/dbscan_density.jpeg");
	
		
		/**
		 * 
		 * 
		 * Running DBSCAN and KMeans alternately
		 * =====================================
		 * 
		 */
		
		double epsilonDummy = 6;
		
		ExtKMeansClustering kmeans1 = runner.runKMeans(runner.dataGauss, "kmeans1");
		runner.estimateDBSCANEpsilon(runner.getExtData(runner.dataDensity, kmeans1.getDummy()), minPts, "dummy1");
		ExtDBSCANClustering dbscan1 = runner.runDBSCAN(runner.getExtData(runner.dataDensity, kmeans1.getDummy()), minPts, epsilonDummy, "dbscan1");
		ExtKMeansClustering kmeans2 = runner.runKMeans(runner.getExtData(runner.dataGauss, dbscan1.getDummy()), "kmeans2");
		ExtDBSCANClustering dbscan2 = runner.runDBSCAN(runner.getExtData(runner.dataDensity, kmeans2.getDummy()), minPts, epsilonDummy, "dbscan2");
		ExtKMeansClustering kmeans3 = runner.runKMeans(runner.getExtData(runner.dataGauss, dbscan2.getDummy()), "kmeans3");
		ExtDBSCANClustering dbscan3 = runner.runDBSCAN(runner.getExtData(runner.dataDensity, kmeans3.getDummy()), minPts, epsilonDummy, "dbscan3");
		ExtKMeansClustering kmeans4 = runner.runKMeans(runner.getExtData(runner.dataGauss, dbscan3.getDummy()), "kmeans4");
		ExtDBSCANClustering dbscan4 = runner.runDBSCAN(runner.getExtData(runner.dataDensity, kmeans4.getDummy()), minPts, epsilonDummy, "dbscan4");
		ExtKMeansClustering kmeans5 = runner.runKMeans(runner.getExtData(runner.dataGauss, dbscan4.getDummy()), "kmeans5");
		ExtDBSCANClustering dbscan5 = runner.runDBSCAN(runner.getExtData(runner.dataDensity, kmeans5.getDummy()), minPts, epsilonDummy, "dbscan5");
		
		
//		int numRuns = 10;
//		ExtKMeansClustering kmeans = null;
//		ExtDBSCANClustering dbscan = null;
//		
//		for (int r = 1; r <= numRuns;) {
//
//			kmeans = runner.runKMeans(r == 1 ? runner.dataGauss : runner.getExtData(runner.dataGauss, dbscan.getDummy()), "run_" + (r++));
//			dbscan = runner.runDBSCAN(runner.getExtData(runner.dataDensity, kmeans.getDummy()), minPts, epsilonDensity, "run_" + (r++));
//		}
		
		System.out.println("\nFinished.");
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
	private ExtDBSCANClustering runDBSCAN(double[][] data, int minPts, double epsilon, String iterLabel) {
		
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
		
		if (showPlots) {
			ChartFrame frame = new ChartFrame(filename, chartDensity);
			frame.pack();
			frame.setVisible(true);
		}
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
		
		if (showPlots) {
			ChartFrame frame = new ChartFrame(filename, chartDensity);
			frame.pack();
			frame.setVisible(true);
		}
	}
	
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
		paramsKNN.addParameter(KNNDistancesSampler.Parameterizer.SAMPLING_ID, 0.2);
		paramsKNN.addParameter(KNNDistancesSampler.Parameterizer.SEED_ID, 1234);
		
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
}
