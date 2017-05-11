package p2.clustering;

import java.util.Map;

import de.lmu.ifi.dbs.elki.algorithm.clustering.kmeans.KMeansMacQueen;
import de.lmu.ifi.dbs.elki.data.Cluster;
import de.lmu.ifi.dbs.elki.data.Clustering;
import de.lmu.ifi.dbs.elki.data.DoubleVector;
import de.lmu.ifi.dbs.elki.data.model.KMeansModel;
import de.lmu.ifi.dbs.elki.database.AbstractDatabase;
import de.lmu.ifi.dbs.elki.database.Database;
import de.lmu.ifi.dbs.elki.database.StaticArrayDatabase;
import de.lmu.ifi.dbs.elki.datasource.ArrayAdapterDatabaseConnection;
import de.lmu.ifi.dbs.elki.distance.distancefunction.minkowski.EuclideanDistanceFunction;
import de.lmu.ifi.dbs.elki.utilities.ClassGenericsUtil;
import de.lmu.ifi.dbs.elki.utilities.optionhandling.parameterization.ListParameterization;

public class RunKMeans {

	public static void main(String[] args) {

		String file = "data/merged.csv";
		Util util = new Util();	
		
		Map<String, String> header = util.getHeader(file);
		int k = Integer.parseInt(header.get("numClusters"));
		
		double[][] data = util.loadData(file);
		
		ArrayAdapterDatabaseConnection dbc = new ArrayAdapterDatabaseConnection(data);

		ListParameterization dbParams = new ListParameterization();
		dbParams.addParameter(AbstractDatabase.Parameterizer.DATABASE_CONNECTION_ID, dbc);

		Database db = ClassGenericsUtil.parameterizeOrAbort(StaticArrayDatabase.class, dbParams);
		db.initialize();

		ListParameterization kmeansParams = new ListParameterization();
		kmeansParams.addParameter(KMeansMacQueen.K_ID, k);
		kmeansParams.addParameter(KMeansMacQueen.DISTANCE_FUNCTION_ID, EuclideanDistanceFunction.class);
		
		KMeansMacQueen<DoubleVector> kmeans = ClassGenericsUtil.parameterizeOrAbort(KMeansMacQueen.class, kmeansParams);
		
		Clustering<KMeansModel> c = kmeans.run(db);

		int i = 0;
		for (Cluster<KMeansModel> clu : c.getAllClusters()) 
			System.out.println("Cluster #" + i++ + "\n  Size: " + clu.size() + "\n");
		

	}
	

}
