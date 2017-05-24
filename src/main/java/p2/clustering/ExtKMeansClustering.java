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
		
		//DEBUG: log the clustering
		StringBuilder log = new StringBuilder();
		log.append("KMeans Clustering:\n");
		int clusterId = 0;
		for (Cluster<KMeansModel> c : clustering.getAllClusters()) {
			log.append("#" + clusterId + " [" + c.size() + "]");
			for (DBIDIter it = c.getIDs().iter(); it.valid(); it.advance())
				log.append(" " + ids.getOffset(it));
			log.append("\n");
			++clusterId;
		}
		System.out.println(log);

		//DEBUG: log dummy encoding
		/*log = new StringBuilder();
		log.append("\nKMeans Dummy:\n");
		for (int row = 0; row < dummy.length; row++) { 
			String rowString = row + ": ";
			for (int col = 0; col < dummy[row].length; col++)
				rowString += dummy[row][col] + " ";
			log.append(rowString + "\n");
		}
		System.out.println(log);*/
	}

	public Clustering<KMeansModel> getClustering() {
		return clustering;
	}

	public double[][] getDummy() {
		return dummy;
	}
}
