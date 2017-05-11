package p2.clustering;

import java.util.Map;

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

public class RunDBSCAN {

	public static void main(String[] args) {

		String file = "data/merged.csv";
		Util util = new Util();	
		
		Map<String, String> header = util.getHeader(file);
		int dim = Integer.parseInt(header.get("numDimensions"));
		int k = Integer.parseInt(header.get("numClusters"));
		
		int minPts = 2 * dim - 1;	//As suggested by ELKI
		double epsilon = 1.0;
		
		double[][] data = util.loadData(file);
		
		ArrayAdapterDatabaseConnection dbc = new ArrayAdapterDatabaseConnection(data);

		ListParameterization dbParams = new ListParameterization();
		dbParams.addParameter(AbstractDatabase.Parameterizer.DATABASE_CONNECTION_ID, dbc);
		dbParams.addParameter(AbstractDatabase.Parameterizer.INDEX_ID, RStarTreeFactory.class);

		Database db = ClassGenericsUtil.parameterizeOrAbort(StaticArrayDatabase.class, dbParams);
		db.initialize();

				
		//Repeat running the algorithm until the right number of clusters is found.
		
		double inc = 0.01;
		int nFoundClusters = 0;
		Clustering<Model> c;
		
		do {
			ListParameterization dbscanParams = new ListParameterization();
			dbscanParams.addParameter(DBSCAN.Parameterizer.EPSILON_ID, epsilon);
			dbscanParams.addParameter(DBSCAN.Parameterizer.MINPTS_ID, minPts);
			dbscanParams.addParameter(DBSCAN.DISTANCE_FUNCTION_ID, EuclideanDistanceFunction.class);
	
			DBSCAN<DoubleVector> dbscan = ClassGenericsUtil.parameterizeOrAbort(DBSCAN.class, dbscanParams);
			c = dbscan.run(db);
			
			epsilon += inc; 
			
			//Since this DBSCAN implementation counts 0-element-clusters as clusters we have to count for ourself.
			nFoundClusters = 0;
			for (Cluster<Model> clu : c.getAllClusters())
				if (clu.size() > 0)
					++nFoundClusters;
		
		} while (nFoundClusters > k);
			
		
		int i = 0;
		for (Cluster<Model> clu : c.getAllClusters())
			if (clu.size() > 0)
				System.out.println("Cluster #" + i++ + "\n  Size: " + clu.size() + "\n");
	}
}
