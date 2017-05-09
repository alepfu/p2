package p2.clustering;

import de.lmu.ifi.dbs.elki.algorithm.clustering.DBSCAN;
import de.lmu.ifi.dbs.elki.data.Cluster;
import de.lmu.ifi.dbs.elki.data.Clustering;
import de.lmu.ifi.dbs.elki.data.DoubleVector;
import de.lmu.ifi.dbs.elki.data.model.Model;
import de.lmu.ifi.dbs.elki.database.AbstractDatabase;
import de.lmu.ifi.dbs.elki.database.Database;
import de.lmu.ifi.dbs.elki.database.StaticArrayDatabase;
import de.lmu.ifi.dbs.elki.datasource.ArrayAdapterDatabaseConnection;
import de.lmu.ifi.dbs.elki.distance.distancefunction.minkowski.EuclideanDistanceFunction;
import de.lmu.ifi.dbs.elki.index.tree.spatial.rstarvariants.rstar.RStarTreeFactory;
import de.lmu.ifi.dbs.elki.utilities.ClassGenericsUtil;
import de.lmu.ifi.dbs.elki.utilities.optionhandling.parameterization.ListParameterization;

public class RunDBSCANOnGaussianDimensions {

	public static void main(String[] args) {

		int minPts = 5;
		double epsilon = 2.0;
		
		Util util = new Util();
		
		double[][] data = util.loadGaussianData("data/merged.csv");
		
		ArrayAdapterDatabaseConnection dbc = new ArrayAdapterDatabaseConnection(data);

		ListParameterization dbParams = new ListParameterization();
		dbParams.addParameter(AbstractDatabase.Parameterizer.DATABASE_CONNECTION_ID, dbc);
		dbParams.addParameter(AbstractDatabase.Parameterizer.INDEX_ID, RStarTreeFactory.class);

		Database db = ClassGenericsUtil.parameterizeOrAbort(StaticArrayDatabase.class, dbParams);
		db.initialize();

		ListParameterization dbscanParams = new ListParameterization();
		dbscanParams.addParameter(DBSCAN.Parameterizer.EPSILON_ID, epsilon);
		dbscanParams.addParameter(DBSCAN.Parameterizer.MINPTS_ID, minPts);
		dbscanParams.addParameter(DBSCAN.DISTANCE_FUNCTION_ID, EuclideanDistanceFunction.class);

		DBSCAN<DoubleVector> dbscan = ClassGenericsUtil.parameterizeOrAbort(DBSCAN.class, dbscanParams);
		
		Clustering<Model> c = dbscan.run(db);
		
		int i = 0;
		for (Cluster<Model> clu : c.getAllClusters())
			System.out.println("Cluster #" + i++ + "\n  Size: " + clu.size() + "\n");

	}
}
