package p2.clustering;

import de.lmu.ifi.dbs.elki.data.Clustering;
import de.lmu.ifi.dbs.elki.data.model.Model;

public class ExtDBSCANClustering {
	
	private Clustering<Model> dbscanClustering;
	private double[][] dbscanDummy;
	
	public ExtDBSCANClustering(Clustering<Model> dbscanClustering, double[][] dbscanDummy) {
		this.dbscanClustering = dbscanClustering;
		this.dbscanDummy = dbscanDummy;
	}

	public Clustering<Model> getDbscanClustering() {
		return dbscanClustering;
	}
	
	public double[][] getDBSCANDummy() {
		return dbscanDummy;
	}
	
	
}
