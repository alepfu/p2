package p2.clustering;

import de.lmu.ifi.dbs.elki.data.Cluster;
import de.lmu.ifi.dbs.elki.data.Clustering;
import de.lmu.ifi.dbs.elki.data.model.KMeansModel;
import de.lmu.ifi.dbs.elki.database.ids.DBIDIter;
import de.lmu.ifi.dbs.elki.database.ids.DBIDRange;

public class ExtKMeansClustering {
	
	private Clustering<KMeansModel> clustering;
	private double[][] dummy;
	private DBIDRange ids;
	
	public ExtKMeansClustering(Clustering<KMeansModel> clustering, double[][] dummy, DBIDRange ids) {
		this.clustering = clustering;
		this.dummy = dummy;
		this.ids = ids;
	}

	public Clustering<KMeansModel> getClustering() {
		return clustering;
	}

	public double[][] getDummy() {
		return dummy;
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		
		sb.append("KMeans Clustering:\n");
		int clusterId = 0;
		for (Cluster<KMeansModel> c : clustering.getAllClusters()) {
			sb.append("#" + clusterId + " [" + c.size() + "]");
			for (DBIDIter it = c.getIDs().iter(); it.valid(); it.advance())
				sb.append(" " + ids.getOffset(it));
			sb.append("\n");
			++clusterId;
		}

		sb.append("\nKMeans Dummy:\n");
		for (int row = 0; row < dummy.length; row++) { 
			String rowString = row + ": ";
			for (int col = 0; col < dummy[row].length; col++)
				rowString += dummy[row][col] + " ";
			sb.append(rowString + "\n");
		}
		
		return sb.toString();
	}
}
