package p2.clustering;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;

import de.lmu.ifi.dbs.elki.algorithm.clustering.kmeans.KMeansMacQueen;
import de.lmu.ifi.dbs.elki.algorithm.clustering.kmeans.initialization.RandomlyGeneratedInitialMeans;
import de.lmu.ifi.dbs.elki.data.Cluster;
import de.lmu.ifi.dbs.elki.data.Clustering;
import de.lmu.ifi.dbs.elki.data.DoubleVector;
import de.lmu.ifi.dbs.elki.data.NumberVector;
import de.lmu.ifi.dbs.elki.data.model.KMeansModel;
import de.lmu.ifi.dbs.elki.data.type.TypeUtil;
import de.lmu.ifi.dbs.elki.database.Database;
import de.lmu.ifi.dbs.elki.database.StaticArrayDatabase;
import de.lmu.ifi.dbs.elki.database.ids.DBIDRange;
import de.lmu.ifi.dbs.elki.database.relation.Relation;
import de.lmu.ifi.dbs.elki.datasource.ArrayAdapterDatabaseConnection;
import de.lmu.ifi.dbs.elki.datasource.DatabaseConnection;
import de.lmu.ifi.dbs.elki.distance.distancefunction.minkowski.EuclideanDistanceFunction;
import de.lmu.ifi.dbs.elki.math.random.RandomFactory;

public class RunKMeans {

	public static void main(String[] args) {

		double[][] data = loadData("data/data_gauss.csv");

		DatabaseConnection dbc = new ArrayAdapterDatabaseConnection(data);
		Database db = new StaticArrayDatabase(dbc, null);
		db.initialize();
		Relation<NumberVector> rel = db.getRelation(TypeUtil.NUMBER_VECTOR_FIELD);
		DBIDRange ids = (DBIDRange) rel.getDBIDs();

		EuclideanDistanceFunction dist = EuclideanDistanceFunction.STATIC;
		RandomlyGeneratedInitialMeans init = new RandomlyGeneratedInitialMeans(RandomFactory.DEFAULT);

		KMeansMacQueen<DoubleVector> km = new KMeansMacQueen<>(dist, 2, 0, init);

		Clustering<KMeansModel> c = km.run(db);

		int i = 0;
		for (Cluster<KMeansModel> clu : c.getAllClusters()) {
			System.out.println("Cluster #" + i);
			System.out.println("  Size: " + clu.size());
			System.out.println("  Center: " + clu.getModel().getPrototype().toString());
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
