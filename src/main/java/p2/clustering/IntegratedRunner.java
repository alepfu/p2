package p2.clustering;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;

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

public class IntegratedRunner {

	private int numDimensions;
	private int numClusters;
	private int numPointsCluster;
	private int numPoints;
	private double[][] dataGauss;
	private double[][] dataDensity;
	private double[][] data;
	private Clustering<Model> trueClustering;	
	
	//private final static String nr = "1496272536499";
	
	private final static String directory = "/home/alepfu/Desktop/P2_export";
	
	private final static String filename = directory + "/" + "data_5d_3c_1000p_overlap_7seed" + ".csv";
	
	
	
	/**
	 * Directory for exported files.
	 */
	static String exportDir = "/home/alepfu/Desktop/AMI";
	
	public static void main(String[] args) {
		
		IntegratedRunner runner = new IntegratedRunner(filename);
		
		ExtData extData;
		ExtKMeansClustering extKMeansClustering = null;
		ExtDBSCANClustering extDBSCANClustering = null;
		
		for (int r = 1; r <= 12; r++) {
			
			//KMeans
			if (r == 1)
				extData = new ExtData(runner.dataGauss, null);  //First run has no dummy information yet
			else
				extData = new ExtData(runner.dataGauss, extDBSCANClustering.getDummy());
			extKMeansClustering = runner.runKMeans(extData);  
			saveClusteringToFile(exportDir + "/run_" + r + ".csv", extKMeansClustering.getClustering());
			
			++r;
			
			//DBSCAN
			extData = new ExtData(runner.dataDensity, extKMeansClustering.getDummy());
			extDBSCANClustering = runner.runDBSCAN(extData);
			saveClusteringToFile(exportDir + "/run_" + r + ".csv", extDBSCANClustering.getClustering());
		}
			
		/*
		
		//KMeans
		extData = new ExtData(runner.dataGauss, null);  //First run has no dummy information yet
		extKMeansClustering = runner.runKMeans(extData);  
		saveClusteringToFile(exportDir + "/kmeans_1.csv", extKMeansClustering.getClustering());
				
		//DBSCAN
		extData = new ExtData(runner.dataDensity, extKMeansClustering.getDummy());
		extDBSCANClustering = runner.runDBSCAN(extData);
		saveClusteringToFile(exportDir + "/dbscan_1.csv", extDBSCANClustering.getClustering());
		
		//KMeans
		extData = new ExtData(runner.dataGauss, extDBSCANClustering.getDummy());
		extKMeansClustering = runner.runKMeans(extData);
		saveClusteringToFile(exportDir + "/kmeans_2.csv", extKMeansClustering.getClustering());
		
		//DBSCAN
		extData = new ExtData(runner.dataDensity, extKMeansClustering.getDummy());
		extDBSCANClustering = runner.runDBSCAN(extData);
		saveClusteringToFile(exportDir + "/dbscan_2.csv", extDBSCANClustering.getClustering());
		
		//KMeans
		extData = new ExtData(runner.dataGauss, extDBSCANClustering.getDummy());
		extKMeansClustering = runner.runKMeans(extData);
		saveClusteringToFile(exportDir + "/kmeans_3.csv", extKMeansClustering.getClustering());
		
		//DBSCAN
		extData = new ExtData(runner.dataDensity, extKMeansClustering.getDummy());
		extDBSCANClustering = runner.runDBSCAN(extData);
		saveClusteringToFile(exportDir + "/dbscan_3.csv", extDBSCANClustering.getClustering());
		
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
		
		System.out.println("Loading data ...");
		
		numDimensions = dataUtil.getNumDimensions();
		numClusters = dataUtil.getNumClusters();
		numPointsCluster = dataUtil.getNumPointsCluster();
		numPoints = numClusters * numPointsCluster;
		data = dataUtil.loadData();
		dataGauss = dataUtil.loadGaussianData();
		dataDensity = dataUtil.loadDensityData();
		
		try {
			
			System.out.println("Generate true clustering ...");
			
			trueClustering = getTrueClustering();
			saveClusteringToFile(exportDir + "/true.csv", trueClustering);
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	
	/**
	 * Runs KMeans clustering using data extended by DBSCAN cluster labels in dummy encoding.
	 * @return A KMeans clustering result including the found cluster labels in dummy encoding.
	 */
	private ExtKMeansClustering runKMeans(ExtData extDataGauss) {
		
		System.out.println("Running KMeans ...");
		
		double[][] dataWithDummy = extDataGauss.getExtData();
		
		ArrayAdapterDatabaseConnection gaussDBConn = new ArrayAdapterDatabaseConnection(dataWithDummy);
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
		
		//Generate dummy encoding
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
	 * Runs DBSCAN clustering using data extended by KMeans cluster labels in dummy encoding.
	 * Since we know k (the number of clusters) we do multiple runs of DBSCAN with increasing epsilon until k clusters are found.
	 * 
	 * @return A DBSCAN clustering result including the found cluster labels in dummy encoding.
	 */
	private ExtDBSCANClustering runDBSCAN(ExtData extData) {
		
		System.out.println("Running DBSCAN ...");
		
		double[][] dataWithDummy = extData.getExtData();
		
		ArrayAdapterDatabaseConnection densityDBConn = new ArrayAdapterDatabaseConnection(dataWithDummy);
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
		int minPts = 2 * numDimensions - 1;		//Suggested by Elki = 2 * numDimensions - 1
		double epsilon = 0.5;					
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
		
		//Generate dummy encoding
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
	 * Loads the true clustering from the data file
	 * Use for evaluation of clusterings.
	 * @return The true clustering.
	 */
	private Clustering<Model> getTrueClustering() throws FileNotFoundException {

		List<ObjectFilter> filterlist = new ArrayList<>();
		filterlist.add(new FixedDBIDsFilter(1));
		
		NumberVectorLabelParser<DoubleVector> parser = new NumberVectorLabelParser<>(DoubleVector.FACTORY);
		FileBasedDatabaseConnection dbc = new FileBasedDatabaseConnection(filterlist, parser, filename);
		
		ListParameterization params = new ListParameterization();
		params.addParameter(AbstractDatabase.Parameterizer.DATABASE_CONNECTION_ID, dbc);
		
		Database db = ClassGenericsUtil.parameterizeOrAbort(StaticArrayDatabase.class, params);
		db.initialize();
		
		Relation<?> rel = db.getRelation(TypeUtil.ANY);
		DBIDRange ids = (DBIDRange) rel.getDBIDs();

		ByLabelClustering algo = new ByLabelClustering();
		Clustering<Model> clustering = algo.run(db);
		
		//DEBUG log clustering
		/*StringBuilder log = new StringBuilder();
		log.append("\nTrue Clustering:\n");
		int clusterId = 1;
		for (Cluster<Model> c : clustering.getAllClusters()) {
			log.append("#" + clusterId + " [" + c.size() + "]");
			
			List<Integer> idList = new ArrayList<Integer>();
			for (DBIDIter it = c.getIDs().iter(); it.valid(); it.advance())
				idList.add(ids.getOffset(it));
			
			Collections.sort(idList);
			for (Integer id : idList)
				log.append(" " + id);
			
			log.append("\n");
			++clusterId;
		}
		System.out.println(log);*/
		
		return clustering;
	}

	/**
	 * Compares a clustering with the true clutering using mutual information.
	 * 
	 * High mutual information indicates a large reduction in uncertainty;
	 * low mutual information indicates a small reduction;
	 * and zero mutual information between two random variables means the variables are independent. 
	 * 
	 * @param clustering Clustering to compare to true clustering
	 * @return The mutual information measure
	 */
	private double evaluateClustering(Clustering<?> clustering) {
		
		ClusterContingencyTable ct = new ClusterContingencyTable(false, false);
		ct.process(trueClustering, clustering);
		
		return ct.getEntropy().entropyMutualInformation();
	}
	
	/**
	 * Exports a clustering to a file.
	 * Use for evaluation with Matlab implementation of Adjusted Mutual Information.
	 * File format: e.g. "1 1 1 2 2" would say objects 1-3 are in cluster 1 and objects 4-5 are in cluster 2.
	 * @param filename A filename
	 * @param clustering A clustering 
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
	

}
