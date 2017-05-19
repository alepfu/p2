package p2.clustering;

import de.lmu.ifi.dbs.elki.algorithm.clustering.DBSCAN;
import de.lmu.ifi.dbs.elki.algorithm.clustering.kmeans.KMeansMacQueen;
import de.lmu.ifi.dbs.elki.data.Clustering;
import de.lmu.ifi.dbs.elki.data.DoubleVector;
import de.lmu.ifi.dbs.elki.data.NumberVector;
import de.lmu.ifi.dbs.elki.data.model.KMeansModel;
import de.lmu.ifi.dbs.elki.data.model.Model;
import de.lmu.ifi.dbs.elki.data.type.TypeUtil;
import de.lmu.ifi.dbs.elki.database.AbstractDatabase;
import de.lmu.ifi.dbs.elki.database.Database;
import de.lmu.ifi.dbs.elki.database.StaticArrayDatabase;
import de.lmu.ifi.dbs.elki.database.ids.DBIDRange;
import de.lmu.ifi.dbs.elki.database.relation.Relation;
import de.lmu.ifi.dbs.elki.datasource.ArrayAdapterDatabaseConnection;
import de.lmu.ifi.dbs.elki.distance.distancefunction.minkowski.EuclideanDistanceFunction;
import de.lmu.ifi.dbs.elki.index.tree.spatial.rstarvariants.rstar.RStarTreeFactory;
import de.lmu.ifi.dbs.elki.utilities.ClassGenericsUtil;
import de.lmu.ifi.dbs.elki.utilities.optionhandling.parameterization.ListParameterization;

public class IntegratedRunner {

	private String file;
//	private double[] kmeansDummy;
//	private double[] dbscanDummy;
	private DataLoaderUtil dataUtil;	
	int numDimensions;
	int numClusters;
	int numPointsPerCluster;
	double[][] dataGauss;
	double[][] dataDensity;
	
	private void setup() {

		dataUtil = new DataLoaderUtil("data/merged.csv");	
		
		numDimensions = dataUtil.getNumDimensions();
		numClusters = dataUtil.getNumClusters();
		numPointsPerCluster = dataUtil.getNumPointsPerCluster();
		
		dataGauss = dataUtil.loadGaussianData();
		dataDensity = dataUtil.loadDensityData();
		
//		kmeansDummy = new double[numClusters];
//		dbscanDummy = new double[numClusters];	
		
	}

	private void run() {
		
		setup();
		
		//TODO run KMeans and DBSCAN alternatly until convergence
		
		// wir starten mit Kmeans und gehen dann mit der normalisierten dummy codierung in DBSCAN rein.
		
		
	}
	
	private void runKMeans() {
		
		ArrayAdapterDatabaseConnection gaussDBConn = new ArrayAdapterDatabaseConnection(dataGauss);
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
		Clustering<KMeansModel> kMeansClustering = kmeans.run(gaussDB);
		
		
		double[][] kmeansDummy = new double[numClusters][numPointsPerCluster];

		//TODO loop over clustering results and set 1/0 in the rows of dummy
		
		
		
		//TODO return ExtKMeansClustering
	}
	
	private void runDBSCAN() {
		
		ArrayAdapterDatabaseConnection densityDBConn = new ArrayAdapterDatabaseConnection(dataDensity);
		ListParameterization densityDBParams = new ListParameterization();
		densityDBParams.addParameter(AbstractDatabase.Parameterizer.DATABASE_CONNECTION_ID, densityDBConn);
		densityDBParams.addParameter(AbstractDatabase.Parameterizer.INDEX_ID, RStarTreeFactory.class);
		Database densityDB = ClassGenericsUtil.parameterizeOrAbort(StaticArrayDatabase.class, densityDBParams);
		densityDB.initialize();
		
		Relation<NumberVector> densityRel = densityDB.getRelation(TypeUtil.NUMBER_VECTOR_FIELD);
	    DBIDRange densityIDs = (DBIDRange) densityRel.getDBIDs();
		
	    ListParameterization dbscanParams = new ListParameterization();
		dbscanParams.addParameter(DBSCAN.Parameterizer.EPSILON_ID, 2 * numDimensions - 1);
		dbscanParams.addParameter(DBSCAN.Parameterizer.MINPTS_ID, 1.0 * (numDimensions / 2.0));
		dbscanParams.addParameter(DBSCAN.DISTANCE_FUNCTION_ID, EuclideanDistanceFunction.class);
		
		DBSCAN<DoubleVector> dbscan = ClassGenericsUtil.parameterizeOrAbort(DBSCAN.class, dbscanParams);
		Clustering<Model> dbscanClustering = dbscan.run(densityDB);
		
		
		double[][] dbscanDummy = new double[numClusters][numPointsPerCluster];
		
		//TODO loop over clustering results and set 1/0 in the rows of dummy
		
		
		
		//TODO return ExtDBSCANClustering
	}
	
	public static void main(String[] args) {
		
		IntegratedRunner runner = new IntegratedRunner();
		runner.run();
		
		//TODO print and process results
		
		
	}

}
