package p2.clustering;

import de.lmu.ifi.dbs.elki.algorithm.clustering.DBSCAN;
import de.lmu.ifi.dbs.elki.data.Cluster;
import de.lmu.ifi.dbs.elki.data.Clustering;
import de.lmu.ifi.dbs.elki.data.DoubleVector;
import de.lmu.ifi.dbs.elki.data.NumberVector;
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
import de.lmu.ifi.dbs.elki.distance.similarityfunction.cluster.ClusteringAdjustedRandIndexSimilarityFunction;
import de.lmu.ifi.dbs.elki.index.tree.spatial.rstarvariants.rstar.RStarTreeFactory;
import de.lmu.ifi.dbs.elki.utilities.ClassGenericsUtil;
import de.lmu.ifi.dbs.elki.utilities.optionhandling.parameterization.ListParameterization;

public class DBSCANRunner {

	private static final String FILE = "data/merged.csv";
	
	public static void main(String[] args) {

		//Load data and access header information
		DataLoaderUtil util = new DataLoaderUtil(FILE);	
		int numDimensions = util.getNumDimensions();
		int numClusters = util.getNumClusters();
		int numPointsPerCluster = util.getNumPointsPerCluster();
		double[][] data = util.loadMergedData();
		
		//Initialize DB
		ArrayAdapterDatabaseConnection dbc = new ArrayAdapterDatabaseConnection(data);
		ListParameterization dbParams = new ListParameterization();
		dbParams.addParameter(AbstractDatabase.Parameterizer.DATABASE_CONNECTION_ID, dbc);
		dbParams.addParameter(AbstractDatabase.Parameterizer.INDEX_ID, RStarTreeFactory.class);
		Database db = ClassGenericsUtil.parameterizeOrAbort(StaticArrayDatabase.class, dbParams);
		db.initialize();
		Relation<NumberVector> rel = db.getRelation(TypeUtil.NUMBER_VECTOR_FIELD);
	    DBIDRange ids = (DBIDRange) rel.getDBIDs();
				
		//Repeat running the algorithm until the right number of clusters is found
		double stepsize = 0.1;
		int nFoundClusters = 0;
		int minPts = 2 * numDimensions - 1;		//As suggested by ELKI
		double epsilon = 1.0 * (numDimensions / 2.0);	//Empiric
		Clustering<Model> dbscanClustering;
		do {
			//Set clustering algorithm parameters
			ListParameterization dbscanParams = new ListParameterization();
			dbscanParams.addParameter(DBSCAN.Parameterizer.EPSILON_ID, epsilon);
			dbscanParams.addParameter(DBSCAN.Parameterizer.MINPTS_ID, minPts);
			dbscanParams.addParameter(DBSCAN.DISTANCE_FUNCTION_ID, EuclideanDistanceFunction.class);
	
			//Run clustering
			DBSCAN<DoubleVector> dbscan = ClassGenericsUtil.parameterizeOrAbort(DBSCAN.class, dbscanParams);
			dbscanClustering = dbscan.run(db);
			System.out.println("Running DBSCAN clustering (eps = " + epsilon + ") ...");
			epsilon += stepsize; 
			
			//Since this DBSCAN implementation counts 0-element-clusters as clusters we have to count for ourself
			nFoundClusters = 0;
			for (Cluster<Model> c : dbscanClustering.getAllClusters())
				if (c.size() > 0)
					++nFoundClusters;
		
		} while (nFoundClusters > numClusters);
			
		//Output the found clusters and objects, omit 0-element-clusters
		int clusterId = 0;
		for (Cluster<Model> c : dbscanClustering.getAllClusters()) {
			if (c.size() > 0) {
				System.out.print("#" + clusterId + " [" + c.size() + "]");
				for (DBIDIter it = c.getIDs().iter(); it.valid(); it.advance())
					System.out.print(" " + ids.getOffset(it));
				System.out.println();
				++clusterId; 
			}
		}
		
		//Evaluate the clustering
		System.out.println("Evaluating found clusters ...");
		Clustering<Model> gtClustering = getGroundTruthClustering(data, numClusters, numPointsPerCluster);		
		ClusteringAdjustedRandIndexSimilarityFunction ari = new ClusteringAdjustedRandIndexSimilarityFunction();
		double similarity = ari.similarity(dbscanClustering, gtClustering);
		System.out.println("ARI similarity to ground truth = " + similarity);
		
		
		
		
		
		
		System.out.println("Finished.");
	}
	
	private static Clustering<Model> getGroundTruthClustering(double[][] data, int numClusters, int numPointsPerCluster) {
		
		Clustering<Model> gtClustering = new Clustering<Model>("Ground truth", "gd");
		SimpleDBIDFactory idFactory = new SimpleDBIDFactory();
		
		for (int i = 0; i < numClusters; i++) {
			DBIDRange ids = idFactory.generateStaticDBIDRange(i * numPointsPerCluster, numPointsPerCluster);
			Cluster<Model> gdCluster = new Cluster<Model>(Integer.toString(i), ids);
			gtClustering.addToplevelCluster(gdCluster);
		}		

		return gtClustering;
	}
}
