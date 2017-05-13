package p2.clustering;

import java.util.List;
import java.util.Map;

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
import de.lmu.ifi.dbs.elki.database.relation.Relation;
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
		
		double[][] data = util.loadMergedData(file);
		
		ArrayAdapterDatabaseConnection dbc = new ArrayAdapterDatabaseConnection(data);

		ListParameterization dbParams = new ListParameterization();
		dbParams.addParameter(AbstractDatabase.Parameterizer.DATABASE_CONNECTION_ID, dbc);

		Database db = ClassGenericsUtil.parameterizeOrAbort(StaticArrayDatabase.class, dbParams);
		db.initialize();
		
		Relation<NumberVector> rel = db.getRelation(TypeUtil.NUMBER_VECTOR_FIELD);
	    DBIDRange ids = (DBIDRange) rel.getDBIDs();
		
		ListParameterization kmeansParams = new ListParameterization();
		kmeansParams.addParameter(KMeansMacQueen.K_ID, k);
		kmeansParams.addParameter(KMeansMacQueen.DISTANCE_FUNCTION_ID, EuclideanDistanceFunction.class);
		
		KMeansMacQueen<DoubleVector> kmeans = ClassGenericsUtil.parameterizeOrAbort(KMeansMacQueen.class, kmeansParams);
		
		Clustering<KMeansModel> c = kmeans.run(db);
		
		int clusterId = 0;
		for (Cluster<KMeansModel> clu : c.getAllClusters()) {
			System.out.print("#" + clusterId + " [" + clu.size() + "]");
			for (DBIDIter it = clu.getIDs().iter(); it.valid(); it.advance())
				System.out.print(" " + ids.getOffset(it));
			System.out.println();
			++clusterId;
		}

		
		//Validate clustering, compare cluster objects of ground truth (input file) and clustering (DBIDs)
		//https://en.wikipedia.org/wiki/Cluster_analysis#Evaluation_and_assessment
		//Jaccard index wont't work, since clusters don't have an order!
		
		//Use Adjusted Rand Index from ELKI
		
		
		
	}
	

}
