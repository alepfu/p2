package p2.clustering;

import de.lmu.ifi.dbs.elki.data.Clustering;
import de.lmu.ifi.dbs.elki.data.model.KMeansModel;

public class ExtKMeansClustering {
	
	private Clustering<KMeansModel> kmeansClustering;
	private double[][] kmeansDummy;
	
	public ExtKMeansClustering(Clustering<KMeansModel> kmeansClustering, double[][] kmeansDummy) {
		this.kmeansClustering = kmeansClustering;
		this.kmeansDummy = kmeansDummy;
	}

	public Clustering<KMeansModel> getKmeansClustering() {
		return kmeansClustering;
	}

	public double[][] getKMeansDummy() {
		return kmeansDummy;
	}
}
