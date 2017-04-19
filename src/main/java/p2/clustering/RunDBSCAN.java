package p2.clustering;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;

import de.lmu.ifi.dbs.elki.algorithm.clustering.DBSCAN;
import de.lmu.ifi.dbs.elki.data.Cluster;
import de.lmu.ifi.dbs.elki.data.Clustering;
import de.lmu.ifi.dbs.elki.data.DoubleVector;
import de.lmu.ifi.dbs.elki.data.NumberVector;
import de.lmu.ifi.dbs.elki.data.model.Model;
import de.lmu.ifi.dbs.elki.data.type.TypeUtil;
import de.lmu.ifi.dbs.elki.database.Database;
import de.lmu.ifi.dbs.elki.database.StaticArrayDatabase;
import de.lmu.ifi.dbs.elki.database.ids.DBIDRange;
import de.lmu.ifi.dbs.elki.database.relation.Relation;
import de.lmu.ifi.dbs.elki.datasource.ArrayAdapterDatabaseConnection;
import de.lmu.ifi.dbs.elki.datasource.DatabaseConnection;
import de.lmu.ifi.dbs.elki.distance.distancefunction.minkowski.EuclideanDistanceFunction;

public class RunDBSCAN {

	public static void main(String[] args) {

		double[][] data = loadData("data/data_density.csv");
		
		DatabaseConnection dbc = new ArrayAdapterDatabaseConnection(data);
		Database db = new StaticArrayDatabase(dbc, null);
		db.initialize();
		Relation<NumberVector> rel = db.getRelation(TypeUtil.NUMBER_VECTOR_FIELD);
		DBIDRange ids = (DBIDRange) rel.getDBIDs();

		EuclideanDistanceFunction dist = EuclideanDistanceFunction.STATIC;
		
		int minPts = 5;
		double epsilon = 2.0;
		
		DBSCAN<DoubleVector> dbscan = new DBSCAN<>(dist, epsilon, minPts);
		
		Clustering<Model> c = dbscan.run(db);
		 
		int i = 0;
		for (Cluster<Model> clu : c.getAllClusters()) {
			System.out.println("Cluster #" + i);
			System.out.println("Size: " + clu.size());
			System.out.println();
			++i;
		} 

	}

	public static double[][] loadData(String file) {
		
		ArrayList<double[]> rows = new ArrayList<double[]>();
		String line = "";

		try {
			BufferedReader br = new BufferedReader(new FileReader(file));
			while ((line = br.readLine()) != null) {

				double[] doubleValues = Arrays.stream(line.split(",")).mapToDouble(Double::parseDouble).toArray();
				rows.add(doubleValues);
			}

			br.close();

		} catch (Exception e) {
			e.printStackTrace();
		}

		double[][] data = new double[rows.size()][rows.get(0).length];

		for (int i = 0; i < rows.size(); i++)
			data[i] = rows.get(i);

		return data;
	}
	
}
