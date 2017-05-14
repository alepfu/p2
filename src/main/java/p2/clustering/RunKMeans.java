package p2.clustering;

import de.lmu.ifi.dbs.elki.algorithm.clustering.kmeans.KMeansMacQueen;
import de.lmu.ifi.dbs.elki.data.Cluster;
import de.lmu.ifi.dbs.elki.data.Clustering;
import de.lmu.ifi.dbs.elki.data.DoubleVector;
import de.lmu.ifi.dbs.elki.data.NumberVector;
import de.lmu.ifi.dbs.elki.data.model.KMeansModel;
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
import de.lmu.ifi.dbs.elki.distance.similarityfunction.cluster.ClusteringAdjustedRandIndexSimilarityFunction;
import de.lmu.ifi.dbs.elki.utilities.ClassGenericsUtil;
import de.lmu.ifi.dbs.elki.utilities.optionhandling.parameterization.ListParameterization;

public class RunKMeans {

	private static final String FILE = "data/merged.csv";
	
	public static void main(String[] args) {

		//Load data and access header information
		LoadDataUtil util = new LoadDataUtil(FILE);	
		int numClusters = util.getNumClusters();
		int numPointsPerCluster = util.getNumPointsPerCluster();
		double[][] data = util.loadMergedData();
		
		//Initialize DB and clustering algorithm
		ArrayAdapterDatabaseConnection dbc = new ArrayAdapterDatabaseConnection(data);
		ListParameterization dbParams = new ListParameterization();
		dbParams.addParameter(AbstractDatabase.Parameterizer.DATABASE_CONNECTION_ID, dbc);
		Database db = ClassGenericsUtil.parameterizeOrAbort(StaticArrayDatabase.class, dbParams);
		db.initialize();
		Relation<NumberVector> rel = db.getRelation(TypeUtil.NUMBER_VECTOR_FIELD);
	    DBIDRange ids = (DBIDRange) rel.getDBIDs();
		ListParameterization kmeansParams = new ListParameterization();
		kmeansParams.addParameter(KMeansMacQueen.K_ID, numClusters);
		kmeansParams.addParameter(KMeansMacQueen.DISTANCE_FUNCTION_ID, EuclideanDistanceFunction.class);
		
		//Run clustering
		System.out.println("Running KMeans clustering ...");
		KMeansMacQueen<DoubleVector> kmeans = ClassGenericsUtil.parameterizeOrAbort(KMeansMacQueen.class, kmeansParams);
		Clustering<KMeansModel> kMeansClustering = kmeans.run(db);
		
		//Output the found clusters and objects
		int clusterId = 0;
		for (Cluster<KMeansModel> c : kMeansClustering.getAllClusters()) {
			System.out.print("#" + clusterId + " [" + c.size() + "]");
			for (DBIDIter it = c.getIDs().iter(); it.valid(); it.advance())
				System.out.print(" " + ids.getOffset(it));
			System.out.println();
			++clusterId;
		}
		
		//Evaluate the clustering
		System.out.println("Evaluating found clusters ...");
		Clustering<KMeansModel> gdClustering = getGroundTruthClustering(data, numClusters, numPointsPerCluster);		
		ClusteringAdjustedRandIndexSimilarityFunction ari = new ClusteringAdjustedRandIndexSimilarityFunction();
		double similarity = ari.similarity(kMeansClustering, gdClustering);
		System.out.println("ARI similarity to ground truth = " + similarity);
		
		
		
		System.out.println("Finished.");
	}
	
	private static Clustering<KMeansModel> getGroundTruthClustering(double[][] data, int numClusters, int numPointsPerCluster) {
		
		Clustering<KMeansModel> gdClustering = new Clustering<KMeansModel>("Ground truth", "gd");
		SimpleDBIDFactory idFactory = new SimpleDBIDFactory();
		
		for (int i = 0; i < numClusters; i++) {
			DBIDRange ids = idFactory.generateStaticDBIDRange(i * numPointsPerCluster, numPointsPerCluster);
			Cluster<KMeansModel> gdCluster = new Cluster<KMeansModel>(Integer.toString(i), ids);
			gdClustering.addToplevelCluster(gdCluster);
		}		

		return gdClustering;
	}

}
