package p2.clustering;

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

public class IntegratedRunner {

	private int numDimensions;
	private int numClusters;
	private int numPointsPerCluster;
	private int numPoints;
	private double[][] dataGauss;
	private double[][] dataDensity;
	
	private Clustering<Model> gtClustering;	

	public static void main(String[] args) {
		
		IntegratedRunner runner = new IntegratedRunner("/home/alepfu/Desktop/P2_export/data_1496004163759.csv");
		
		ExtGaussData extDataGauss;
		ExtDensityData extDensityData;
		ExtKMeansClustering extKMeansClustering;
		ExtDBSCANClustering extDBSCANClustering = null;
		
		//KMeans
		extDataGauss = new ExtGaussData(runner.dataGauss, null);  //1st run has no dummy information yet
		extKMeansClustering = runner.runKMeans(extDataGauss);  
		System.out.println(runner.evaluateClustering(extKMeansClustering.getClustering()));
		
		//DBSCAN
		extDensityData = new ExtDensityData(runner.dataDensity, extKMeansClustering.getDummy());
		extDBSCANClustering = runner.runDBSCAN(extDensityData);
		System.out.println(runner.evaluateClustering(extDBSCANClustering.getClustering()));
		
		//KMeans
		/*		extDataGauss = new ExtGaussData(runner.dataGauss, extDBSCANClustering.getDummy());
		extKMeansClustering = runner.runKMeans(extDataGauss);
		System.out.println(runner.evaluateClustering(extKMeansClustering.getClustering()));
		
		//DBSCAN
		extDensityData = new ExtDensityData(runner.dataDensity, extKMeansClustering.getDummy());
		extDBSCANClustering = runner.runDBSCAN(extDensityData);
		System.out.println(runner.evaluateClustering(extDBSCANClustering.getClustering()));
*/		
		
		System.out.println("Finished.");
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
		
		//DEBUG: log ground truth clustering
		StringBuilder log = new StringBuilder();
		log.append("Ground Truth Clustering:\n");
		int clusterId = 0;
		for (Cluster<Model> c : gtClustering.getAllClusters()) {
			if (c.size() > 0) {
				log.append("#" + clusterId + " [" + c.size() + "]");
				log.append("\n");
				++clusterId;
			}
		}
		System.out.println(log);
	}
	
	/**
	 * Runs KMeans clustering using data extended by DBSCAN cluster labels in dummy encoding.
	 * @return A KMeans clustering result including the found cluster labels in dummy encoding.
	 */
	private ExtKMeansClustering runKMeans(ExtGaussData extDataGauss) {
		
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
	 * Runs DBSCAN clustering using data extended by KMeans cluster labels in dummy encoding.
	 * Since we know k (the number of clusters) we do multiple runs of DBSCAN with increasing epsilon until k clusters are found.
	 * @return A DBSCAN clustering result including the found cluster labels in dummy encoding.
	 */
	private ExtDBSCANClustering runDBSCAN(ExtDensityData extDensityData) {
		
		ArrayAdapterDatabaseConnection densityDBConn = new ArrayAdapterDatabaseConnection(extDensityData.getExtData());
		ListParameterization densityDBParams = new ListParameterization();
		densityDBParams.addParameter(AbstractDatabase.Parameterizer.DATABASE_CONNECTION_ID, densityDBConn);
		densityDBParams.addParameter(AbstractDatabase.Parameterizer.INDEX_ID, RStarTreeFactory.class);
		Database densityDB = ClassGenericsUtil.parameterizeOrAbort(StaticArrayDatabase.class, densityDBParams);
		densityDB.initialize();
		
		Relation<NumberVector> densityRel = densityDB.getRelation(TypeUtil.NUMBER_VECTOR_FIELD);
	    DBIDRange densityIDs = (DBIDRange) densityRel.getDBIDs();
		
	    //Successively increase epsilon until DBSCAN finds k clusters. 
		double stepsize = 0.1;
		int nFoundClusters = 0;
		int minPts = 2 * numDimensions - 1;		//As suggested by Elki
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

			nFoundClusters = 0;
			for (Cluster<Model> c : dbscanClustering.getAllClusters())
				if (c.size() > 0)
					++nFoundClusters;

		} while (nFoundClusters > numClusters);
		
		//System.out.println("Epsilon = " + epsilon);
		
		double[][] dbscanDummy = new double[numPoints][numClusters];
		
		int clusterID = 0;
		for (Cluster<Model> c : dbscanClustering.getAllClusters()) {
			if (c.size() > 0) {
				for (DBIDIter it = c.getIDs().iter(); it.valid(); it.advance()) 
					dbscanDummy[densityIDs.getOffset(it)][clusterID] = 1.0;
				++clusterID;
			}
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
		
		Clustering<Model> gdClustering = new Clustering<Model>("Ground truth", "gd");
		SimpleDBIDFactory idFactory = new SimpleDBIDFactory();
		
		for (int i = 0; i < numClusters; i++) {
			DBIDRange ids = idFactory.generateStaticDBIDRange(i * numPointsPerCluster, numPointsPerCluster);
			Cluster<Model> gdCluster = new Cluster<Model>(Integer.toString(i), ids);
			gdClustering.addToplevelCluster(gdCluster);
		}		

		return gdClustering;
	}
	
	/**
	 * Compares a clustering with the ground truth clutering.
	 * @param c Clustering to compare to ground truth
	 * @return The mutual information measure
	 */
	private double evaluateClustering(Clustering<?> c) {
		
		ClusterContingencyTable table = new ClusterContingencyTable(false, false);
		table.process(gtClustering,  c);
		
		return new Entropy(table).normalizedVariationOfInformation();
	}

}
