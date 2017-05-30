package p2.clustering;

import java.util.ArrayList;
import java.util.List;

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
import de.lmu.ifi.dbs.elki.database.ids.integer.SimpleDBIDFactory;
import de.lmu.ifi.dbs.elki.database.relation.Relation;
import de.lmu.ifi.dbs.elki.datasource.ArrayAdapterDatabaseConnection;
import de.lmu.ifi.dbs.elki.distance.distancefunction.minkowski.EuclideanDistanceFunction;
import de.lmu.ifi.dbs.elki.index.tree.spatial.rstarvariants.rstar.RStarTreeFactory;
import de.lmu.ifi.dbs.elki.utilities.ClassGenericsUtil;
import de.lmu.ifi.dbs.elki.utilities.optionhandling.parameterization.ListParameterization;
import p2.elki.ClusterContingencyTable;

public class IntegratedRunner {

	private int numDimensions;
	private int numClusters;
	private int numPointsPerCluster;
	private int numPoints;
	private double[][] dataGauss;
	private double[][] dataDensity;
	
	private Clustering<Model> gtClustering;	

	public static void main(String[] args) {
		
		IntegratedRunner runner = new IntegratedRunner("/home/alepfu/Desktop/P2_export/data_1496098805253.csv");
		
		ExtData extData;
		ExtKMeansClustering extKMeansClustering;
		ExtDBSCANClustering extDBSCANClustering;
		
		//KMeans
		extData = new ExtData(runner.dataGauss, null);  //1st run has no dummy information yet
		extKMeansClustering = runner.runKMeans(extData);  
		runner.evaluateClusteringToGroundTruth(extKMeansClustering.getClustering());
		
		//DBSCAN
		extData = new ExtData(runner.dataDensity, extKMeansClustering.getDummy());
		extDBSCANClustering = runner.runDBSCANSuccessively(extData);
		runner.evaluateClusteringToGroundTruth(extDBSCANClustering.getClustering());
				
		//KMeans
		extData = new ExtData(runner.dataGauss, extDBSCANClustering.getDummy());
		extKMeansClustering = runner.runKMeans(extData);
		runner.evaluateClusteringToGroundTruth(extKMeansClustering.getClustering());
		
		//DBSCAN
		/*
		extDensityData = new ExtDensityData(runner.dataDensity, extKMeansClustering.getDummy());
		extDBSCANClustering = runner.runDBSCAN(extDensityData);
		System.out.println(runner.evaluateClustering(extDBSCANClustering.getClustering()));
*/		
		
		System.out.println("\nFinished.");
	}
	
	/**
	 * Constructor
	 * Loads data and header information for a given file, generates ground truth clustering.
	 * @param file The full filename.
	 */
	public IntegratedRunner(String file) {

		DataLoaderUtil dataUtil = new DataLoaderUtil(file);	
		
		numDimensions = dataUtil.getNumDimensions();
		numClusters = dataUtil.getNumClusters();
		numPointsPerCluster = dataUtil.getNumPointsPerCluster();
		numPoints = numClusters * numPointsPerCluster;
		dataGauss = dataUtil.loadGaussianData();
		dataDensity = dataUtil.loadDensityData();
		
		gtClustering = getGroundTruthClustering(numClusters, numPointsPerCluster);
	}
	
	/**
	 * Runs KMeans clustering using data extended by DBSCAN cluster labels in dummy encoding.
	 * @return A KMeans clustering result including the found cluster labels in dummy encoding.
	 */
	private ExtKMeansClustering runKMeans(ExtData extDataGauss) {
		
		ArrayAdapterDatabaseConnection gaussDBConn = new ArrayAdapterDatabaseConnection(extDataGauss.getExtData());
		ListParameterization gaussDBParams = new ListParameterization();
		gaussDBParams.addParameter(AbstractDatabase.Parameterizer.DATABASE_CONNECTION_ID, gaussDBConn);
		Database gaussDB = ClassGenericsUtil.parameterizeOrAbort(StaticArrayDatabase.class, gaussDBParams);
		gaussDB.initialize();
		
		Relation<NumberVector> gausRel = gaussDB.getRelation(TypeUtil.NUMBER_VECTOR_FIELD);
	    DBIDRange gaussIDs = (DBIDRange) gausRel.getDBIDs();
		
	    ListParameterization kmeansParams = new ListParameterization();
		kmeansParams.addParameter(KMeansMacQueen.K_ID, numClusters);
		kmeansParams.addParameter(KMeansMacQueen.DISTANCE_FUNCTION_ID, EuclideanDistanceFunction.class);
		
		KMeansMacQueen<DoubleVector> kmeans = ClassGenericsUtil.parameterizeOrAbort(KMeansMacQueen.class, kmeansParams);
		Clustering<KMeansModel> kmeansClustering = kmeans.run(gaussDB);
		
		double[][] kmeansDummy = new double[numPoints][numClusters];
		
		int clusterID = 0;
		for (Cluster<KMeansModel> c : kmeansClustering.getAllClusters()) {
			for (DBIDIter it = c.getIDs().iter(); it.valid(); it.advance())
				kmeansDummy[gaussIDs.getOffset(it)][clusterID] = 1.0;
			++clusterID;
		}
		
		ExtKMeansClustering extClustering = new ExtKMeansClustering(kmeansClustering, kmeansDummy, gaussIDs);
		
		return extClustering;
	}
	
	/**
	 * Just a single DBSCAN run on the extended data.
	 * 
	 * @param extDensityData
	 * @return
	 */
	private ExtDBSCANClustering runDBSCANSingle(ExtData extData) {
		
		ArrayAdapterDatabaseConnection densityDBConn = new ArrayAdapterDatabaseConnection(extData.getExtData());
		ListParameterization densityDBParams = new ListParameterization();
		densityDBParams.addParameter(AbstractDatabase.Parameterizer.DATABASE_CONNECTION_ID, densityDBConn);
		densityDBParams.addParameter(AbstractDatabase.Parameterizer.INDEX_ID, RStarTreeFactory.class);
		Database densityDB = ClassGenericsUtil.parameterizeOrAbort(StaticArrayDatabase.class, densityDBParams);
		densityDB.initialize();
		
		Relation<NumberVector> densityRel = densityDB.getRelation(TypeUtil.NUMBER_VECTOR_FIELD);
	    DBIDRange densityIDs = (DBIDRange) densityRel.getDBIDs();
		
		int minPts = 10;							
		double epsilon = 1.0;					
		
		ListParameterization dbscanParams = new ListParameterization();
		dbscanParams.addParameter(DBSCAN.Parameterizer.EPSILON_ID, epsilon);
		dbscanParams.addParameter(DBSCAN.Parameterizer.MINPTS_ID, minPts);
		dbscanParams.addParameter(DBSCAN.DISTANCE_FUNCTION_ID, EuclideanDistanceFunction.class);

		DBSCAN<DoubleVector> dbscan = ClassGenericsUtil.parameterizeOrAbort(DBSCAN.class, dbscanParams);
		Clustering<Model> dbscanClustering = dbscan.run(densityDB);
		
		int numFoundClusters = dbscanClustering.getAllClusters().size();
		
		double[][] dbscanDummy = new double[numPoints][numFoundClusters];
		
		int clusterID = 0;
		for (Cluster<Model> c : dbscanClustering.getAllClusters()) {
			for (DBIDIter it = c.getIDs().iter(); it.valid(); it.advance()) 
				dbscanDummy[densityIDs.getOffset(it)][clusterID] = 1.0;
			++clusterID;
		}
		
		ExtDBSCANClustering extClustering = new ExtDBSCANClustering(dbscanClustering, dbscanDummy, densityIDs);
		
		return extClustering;
	}
	
	
	/**
	 * Runs DBSCAN clustering using data extended by KMeans cluster labels in dummy encoding.
	 * Since we know k (the number of clusters) we do multiple runs of DBSCAN with increasing epsilon until k clusters are found.
	 * 
	 * @return A DBSCAN clustering result including the found cluster labels in dummy encoding.
	 */
	private ExtDBSCANClustering runDBSCANSuccessively(ExtData extData) {
		
		ArrayAdapterDatabaseConnection densityDBConn = new ArrayAdapterDatabaseConnection(extData.getExtData());
		ListParameterization densityDBParams = new ListParameterization();
		densityDBParams.addParameter(AbstractDatabase.Parameterizer.DATABASE_CONNECTION_ID, densityDBConn);
		densityDBParams.addParameter(AbstractDatabase.Parameterizer.INDEX_ID, RStarTreeFactory.class);
		Database densityDB = ClassGenericsUtil.parameterizeOrAbort(StaticArrayDatabase.class, densityDBParams);
		densityDB.initialize();
		
		Relation<NumberVector> densityRel = densityDB.getRelation(TypeUtil.NUMBER_VECTOR_FIELD);
	    DBIDRange densityIDs = (DBIDRange) densityRel.getDBIDs();
		
	    //Successively increase epsilon until DBSCAN finds k clusters. 
		double stepsize = 0.1;
		int numFoundClusters = 0;
		int minPts = 10;						//Suggested by Elki = 2 * numDimensions - 1
		double epsilon = 1.0;					//TODO Is this a good idea? Start big/small?
		Clustering<Model> dbscanClustering;
		
		do {
			ListParameterization dbscanParams = new ListParameterization();
			dbscanParams.addParameter(DBSCAN.Parameterizer.EPSILON_ID, epsilon);
			dbscanParams.addParameter(DBSCAN.Parameterizer.MINPTS_ID, minPts);
			dbscanParams.addParameter(DBSCAN.DISTANCE_FUNCTION_ID, EuclideanDistanceFunction.class);
			
			DBSCAN<DoubleVector> dbscan = ClassGenericsUtil.parameterizeOrAbort(DBSCAN.class, dbscanParams);
			dbscanClustering = dbscan.run(densityDB);
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
		
		
		double[][] dbscanDummy = new double[numPoints][numFoundClusters];
		
		int clusterID = 0;
		for (Cluster<Model> c : dbscanClustering.getAllClusters()) {
			for (DBIDIter it = c.getIDs().iter(); it.valid(); it.advance()) 
				dbscanDummy[densityIDs.getOffset(it)][clusterID] = 1.0;
			++clusterID;
		}
		
		ExtDBSCANClustering extClustering = new ExtDBSCANClustering(dbscanClustering, dbscanDummy, densityIDs);
		
		return extClustering;
	}
	
	/**
	 * Generates the ground truth clustering.
	 * Use for evaluation of clusterings.
	 * With knowing the number of clusters and the number of points per cluster, this reduces to a
	 * simple numbering procedure.
	 * @param numClusters The number of clusters.
	 * @param numPointsPerCluster The number of points per cluster.
	 * @return The ground truth clustering.
	 */
	private Clustering<Model> getGroundTruthClustering(int numClusters, int numPointsPerCluster) {
		
		Clustering<Model> gtClustering = new Clustering<Model>("Ground truth", "gd");
		SimpleDBIDFactory idFactory = new SimpleDBIDFactory();
		DBIDRange ids = idFactory.generateStaticDBIDRange(numClusters * numPointsPerCluster);

		//DEBUG log ground truth clustering
		StringBuilder log = new StringBuilder();
		log.append("Ground Truth Clustering:\n");
		
		for (int i = 0; i < numClusters; i++) {
			
			Cluster<Model> gtCluster = new Cluster<Model>(Integer.toString(i), ids.slice(i * numPointsPerCluster, (i + 1) * numPointsPerCluster));
			gtClustering.addToplevelCluster(gtCluster);
			
			//DEBUG log ground truth clustering
			log.append("#" + i + " [" + gtCluster.size() + "]");
			for (DBIDIter it = gtCluster.getIDs().iter(); it.valid(); it.advance())
				log.append(" " + ids.getOffset(it));
			log.append("\n");
		}
		
		//DEBUG log ground truth clustering
		System.out.print(log);

		return gtClustering;
	}
	
	/**
	 * Uses extracted stuff from Elki!
	 * 
	 * Compares a clustering with the ground truth clutering using mutual information.
	 * 
	 * High mutual information indicates a large reduction in uncertainty;
	 * low mutual information indicates a small reduction;
	 * and zero mutual information between two random variables means the variables are independent. 
	 * 
	 * @param c Clustering to compare to ground truth
	 * @return The mutual information measure
	 */
	private void evaluateClusteringToGroundTruthElki(Clustering<?> clustering) {
		
		ClusterContingencyTable table = new ClusterContingencyTable(false, false);
		table.process(gtClustering, clustering);
		
		//DEBUG log
		System.out.println(table.toString());
		
		double mi = table.getEntropy().entropyMutualInformation();
		
		System.out.println("Mutual information: " +  mi);
	}
	
	/**
	 * My own implementation of MI between two clusterings.
	 * 
	 */
	private void evaluateClusteringToGroundTruth(Clustering<?> clustering) {
		
		//TODO implement :D
		
		//See paper "Comparing Clusterings - An Overview" for tips on implementation!
		
	}

}
