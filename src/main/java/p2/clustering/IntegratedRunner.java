package p2.clustering;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import de.lmu.ifi.dbs.elki.algorithm.clustering.DBSCAN;
import de.lmu.ifi.dbs.elki.algorithm.clustering.kmeans.KMeansMacQueen;
import de.lmu.ifi.dbs.elki.algorithm.clustering.trivial.ByLabelClustering;
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
import de.lmu.ifi.dbs.elki.datasource.FileBasedDatabaseConnection;
import de.lmu.ifi.dbs.elki.datasource.filter.FixedDBIDsFilter;
import de.lmu.ifi.dbs.elki.datasource.filter.ObjectFilter;
import de.lmu.ifi.dbs.elki.datasource.parser.NumberVectorLabelParser;
import de.lmu.ifi.dbs.elki.distance.distancefunction.minkowski.EuclideanDistanceFunction;
import de.lmu.ifi.dbs.elki.evaluation.clustering.ClusterContingencyTable;
import de.lmu.ifi.dbs.elki.index.tree.spatial.rstarvariants.rstar.RStarTreeFactory;
import de.lmu.ifi.dbs.elki.utilities.ClassGenericsUtil;
import de.lmu.ifi.dbs.elki.utilities.optionhandling.parameterization.ListParameterization;
import p2.old.StatisticsUtil;
import p2.util.DataLoaderUtil;

public class IntegratedRunner {

	/**
	 * Number of dimensions for each cluster type.
	 */
	private int numDimensions;
	
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
	 * Working direcktory
	 */
	private final static String workDir = "/home/alepfu/Desktop/P2";
	
	
	public static void main(String[] args) {
		
		IntegratedRunner runner = new IntegratedRunner(workDir + "/data.csv");
		
/*		//KMeans runs on the whole data set, the gaussian features and the density features
		ExtKMeansClustering kmeans1 = runner.runKMeans(runner.data);
		saveClusteringToFile(workDir + "/kmeans_whole.csv", kmeans1.getClustering());
		ExtKMeansClustering kmeans2 = runner.runKMeans(runner.dataGauss);
		saveClusteringToFile(workDir + "/kmeans_gauss.csv", kmeans2.getClustering());
		ExtKMeansClustering kmeans3 = runner.runKMeans(runner.dataDensity);
		saveClusteringToFile(workDir + "/kmeans_density.csv", kmeans3.getClustering());
		
		//DBSCAN runs on the whole data set, the gaussian features and the density features
		double stepsize = 0.1;
		double initEpsilon = 1.5;
		int minPts = 5;
		ExtDBSCANClustering dbscan1 = runner.runMultipleDBSCAN(runner.data, 3.0, 1.0, minPts);
		saveClusteringToFile(workDir + "/dbscan_whole.csv", dbscan1.getClustering());
		ExtDBSCANClustering dbscan2 = runner.runMultipleDBSCAN(runner.dataGauss, 33.0, 1.0, 3);
		saveClusteringToFile(workDir + "/dbscan_gauss.csv", dbscan2.getClustering());
		*/ExtDBSCANClustering dbscan3 = runner.runMultipleDBSCAN(runner.dataDensity, 0.001, 0.1, 5);
		saveClusteringToFile(workDir + "/dbscan_density.csv", dbscan3.getClustering());

		//TODO parameter für guten run dbscan finden
		
		
	/*	
		
		//MinPTS = 2 * numDimensions - 1
		
		
		//Integrated run
		int numRuns = 48;
		double[][] extData;
		ExtKMeansClustering kmeans = null;
		ExtDBSCANClustering dbscan = null;
		for (int r = 1; r <= numRuns; r++) {
			
			//DBSCAN
			dbscan = runner.runMultipleDBSCAN(r == 1 ? runner.dataDensity : runner.getExtData(runner.dataDensity, kmeans.getDummy()),
					1.5, 0.5, 5);  
			saveClusteringToFile(workDir + "/run_" + r + ".csv", dbscan.getClustering());
			
			++r;
	
			//KMeans
			kmeans = runner.runKMeans(runner.getExtData(runner.dataGauss, dbscan.getDummy()));
			saveClusteringToFile(workDir + "/run_" + r + ".csv", kmeans.getClustering());
		}
	*/		
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
		numDimensions = dataUtil.getNumDimensions();
		numClusters = dataUtil.getNumClusters();
		numPointsCluster = dataUtil.getNumPointsCluster();
		numPoints = numClusters * numPointsCluster;
		data = dataUtil.loadData();
		dataGauss = dataUtil.loadGaussianData();
		dataDensity = dataUtil.loadDensityData();
	}
	
	/**
	 * Runs KMeans clustering using data extended by DBSCAN cluster labels in dummy encoding.
	 * @return A clustering result with additional dummy encoding.
	 */
	private ExtKMeansClustering runKMeans(double[][] data) {
		
		System.out.println("Running KMeans ...");
		
		ArrayAdapterDatabaseConnection gaussDBConn = new ArrayAdapterDatabaseConnection(data);
		ListParameterization gaussDBParams = new ListParameterization();
		gaussDBParams.addParameter(AbstractDatabase.Parameterizer.DATABASE_CONNECTION_ID, gaussDBConn);
		Database gaussDB = ClassGenericsUtil.parameterizeOrAbort(StaticArrayDatabase.class, gaussDBParams);
		gaussDB.initialize();
		
		Relation<NumberVector> gausRel = gaussDB.getRelation(TypeUtil.NUMBER_VECTOR_FIELD);
	    DBIDRange gaussIDs = (DBIDRange) gausRel.getDBIDs();
		
	    ListParameterization kmeansParams = new ListParameterization();
		kmeansParams.addParameter(KMeansMacQueen.K_ID, numClusters);
		kmeansParams.addParameter(KMeansMacQueen.DISTANCE_FUNCTION_ID, EuclideanDistanceFunction.class);
		
		//Run KMeans
		KMeansMacQueen<DoubleVector> kmeans = ClassGenericsUtil.parameterizeOrAbort(KMeansMacQueen.class, kmeansParams);
		Clustering<KMeansModel> kmeansClustering = kmeans.run(gaussDB);
		
		//Generate dummy encoding
		double[][] kmeansDummy = new double[numPoints][numClusters];		
		int clusterID = 0;
		for (Cluster<KMeansModel> c : kmeansClustering.getAllClusters()) {
			for (DBIDIter it = c.getIDs().iter(); it.valid(); it.advance())
				kmeansDummy[gaussIDs.getOffset(it)][clusterID] = 1.0;
			++clusterID;
		}
		
		return new ExtKMeansClustering(kmeansClustering, kmeansDummy, gaussIDs);
	}
	
	/**
	 * Runs DBSCAN a single time for given minPts and epsilon.
	 * @param data The data to process.
	 * @param minPts The minPts parameter of the DBSCAN algorithm.
	 * @param epsilon The epsilon parameter of the DBSCAN algorithm.
	 * @return A clustering result with additional dummy encoding.
	 */
	private ExtDBSCANClustering runSingleDBSCAN(double[][] data, int minPts, double epsilon) {
		
		System.out.println("Running Single-DBSCAN ...");
		
		ArrayAdapterDatabaseConnection densityDBConn = new ArrayAdapterDatabaseConnection(data);
		ListParameterization densityDBParams = new ListParameterization();
		densityDBParams.addParameter(AbstractDatabase.Parameterizer.DATABASE_CONNECTION_ID, densityDBConn);
		densityDBParams.addParameter(AbstractDatabase.Parameterizer.INDEX_ID, RStarTreeFactory.class);
		Database densityDB = ClassGenericsUtil.parameterizeOrAbort(StaticArrayDatabase.class, densityDBParams);
		densityDB.initialize();
		
		Relation<NumberVector> densityRel = densityDB.getRelation(TypeUtil.NUMBER_VECTOR_FIELD);
	    DBIDRange densityIDs = (DBIDRange) densityRel.getDBIDs();
		
		ListParameterization dbscanParams = new ListParameterization();
		dbscanParams.addParameter(DBSCAN.Parameterizer.EPSILON_ID, epsilon);
		dbscanParams.addParameter(DBSCAN.Parameterizer.MINPTS_ID, minPts);
		dbscanParams.addParameter(DBSCAN.DISTANCE_FUNCTION_ID, EuclideanDistanceFunction.class);
		
		//Run DBSCAN
		DBSCAN<DoubleVector> dbscan = ClassGenericsUtil.parameterizeOrAbort(DBSCAN.class, dbscanParams);
		Clustering<Model> dbscanClustering = dbscan.run(densityDB);
		
		//Generate dummy encoding
		double[][] dbscanDummy = new double[numPoints][dbscanClustering.getAllClusters().size()];
		int clusterID = 0;
		for (Cluster<Model> c : dbscanClustering.getAllClusters()) {
			for (DBIDIter it = c.getIDs().iter(); it.valid(); it.advance()) 
				dbscanDummy[densityIDs.getOffset(it)][clusterID] = 1.0;
			++clusterID;
		}
		
		return new ExtDBSCANClustering(dbscanClustering, dbscanDummy, densityIDs);
	}
	
	
	/**
	 * Runs DBSCAN multiple times until k clusters are found by increasing epsilon successivley.
	 * Since we know k (the number of clusters) we do multiple runs of DBSCAN with increasing epsilon until k clusters are found.
	 * @param data The data to process.
	 * @param initEpsilon The initial value for the parameter epsilon of the DBCAN algorithm.
	 * @param stepsize The stepsize epsilon is increased every iteration.
	 * @return A clustering result with additional dummy encoding.
	 */
	private ExtDBSCANClustering runMultipleDBSCAN(double[][] data, double initEpsilon, double stepsize, int minPts) {
		
		System.out.println("Running Multiple-DBSCAN ...");
		
		ArrayAdapterDatabaseConnection densityDBConn = new ArrayAdapterDatabaseConnection(data);
		ListParameterization densityDBParams = new ListParameterization();
		densityDBParams.addParameter(AbstractDatabase.Parameterizer.DATABASE_CONNECTION_ID, densityDBConn);
		densityDBParams.addParameter(AbstractDatabase.Parameterizer.INDEX_ID, RStarTreeFactory.class);
		Database densityDB = ClassGenericsUtil.parameterizeOrAbort(StaticArrayDatabase.class, densityDBParams);
		densityDB.initialize();
		
		Relation<NumberVector> densityRel = densityDB.getRelation(TypeUtil.NUMBER_VECTOR_FIELD);
	    DBIDRange densityIDs = (DBIDRange) densityRel.getDBIDs();
		
	    //Increase epsilon until DBSCAN finds k clusters. 
		int numFoundClusters = 0;
		double epsilon = initEpsilon;					
		Clustering<Model> dbscanClustering;
		do {
			ListParameterization dbscanParams = new ListParameterization();
			dbscanParams.addParameter(DBSCAN.Parameterizer.EPSILON_ID, epsilon);
			dbscanParams.addParameter(DBSCAN.Parameterizer.MINPTS_ID, minPts);
			dbscanParams.addParameter(DBSCAN.DISTANCE_FUNCTION_ID, EuclideanDistanceFunction.class);
			
			DBSCAN<DoubleVector> dbscan = ClassGenericsUtil.parameterizeOrAbort(DBSCAN.class, dbscanParams);
			dbscanClustering = dbscan.run(densityDB);
			
			//System.out.println("eps = " + epsilon);
			
			epsilon += stepsize;
			
			//Don't account for 0-element clusters
			numFoundClusters = 0;
			for (Cluster<Model> c : dbscanClustering.getAllClusters())
				if (c.size() > 0)	
					++numFoundClusters;

		} while (numFoundClusters != numClusters);

		//Strip 0-element clusters
		Clustering<Model> stripped = new Clustering<Model>("DBSCAN Clustering", "dbscan-clustering");
		if (numFoundClusters != dbscanClustering.getAllClusters().size()) {
			for (Cluster<Model> c : dbscanClustering.getAllClusters())
				if (c.size() > 0)
					stripped.addToplevelCluster(c);
			dbscanClustering = stripped;
		}
		
		//Generate dummy encoding
		double[][] dbscanDummy = new double[numPoints][numFoundClusters];
		int clusterID = 0;
		for (Cluster<Model> c : dbscanClustering.getAllClusters()) {
			for (DBIDIter it = c.getIDs().iter(); it.valid(); it.advance()) 
				dbscanDummy[densityIDs.getOffset(it)][clusterID] = 1.0;
			++clusterID;
		}
		
		return new ExtDBSCANClustering(dbscanClustering, dbscanDummy, densityIDs);
	}
	
	/**
	 * Exports a clustering to a file.
	 * Use for evaluation with Matlab implementation of Adjusted Mutual Information.
	 * E.g. "1 1 1 2 2" would say objects 1-3 are in cluster 1 and objects 4-5 are in cluster 2.
	 * @param filename The file to write to.
	 * @param clustering The clustering to export to the file.
	 */
	private static void saveClusteringToFile(String filename, Clustering<? extends Model> clustering) {
		
		StringBuilder log = new StringBuilder();
		
		int clusterLabel = 1;
		for (Cluster<? extends Model> cluster : clustering.getAllClusters()) {
			for (int i = 0; i < cluster.getIDs().size(); i++)
				log.append(clusterLabel + " ");
			++clusterLabel;
		}
		
		try {
			FileWriter writer = new FileWriter(new File(filename));
			writer.write(log.toString());
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
		NumberFormat nf = NumberFormat.getInstance();
		nf.setMinimumFractionDigits(2);
		nf.setMaximumFractionDigits(2);
		StringBuilder log = new StringBuilder("\nExtended data:\n");
		for (int row = 0; row < 3; row++) {
			for (int col = 0; col < numColsConcat; col++)
				log.append(nf.format(extData[row][col]) + " ");
			log.append("\n");
		}
		log.append("...\n");
		System.out.println(log);
		
		return extData;
	}
	
//	/**
//	 * TODO
//	 */
//	public double getMinDistance(double[][] data) {
//		
//		int n = data.length;
//		DescriptiveStatistics stats = new DescriptiveStatistics();
//		for (int i = 0; i < (n - 1); i++) {
//			for (int j = i + 1; j < n; j++) {
//			
//				double[] p = data[i];
//				double[] q = data[j];
//				double sqsum = 0;
//				for (int d = 0; d < numDimensions; d++) 
//					sqsum += Math.pow(p[d] - q[d], 2);
//				stats.addValue(Math.sqrt(sqsum));
//			}
//		}
//		return stats.getMin();
//	}
}
