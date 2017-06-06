package p2.clustering;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

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

		//Log the clustering
		StringBuilder log = new StringBuilder();
		log.append("\nKMeans clustering:\n");
		int clusterId = 1;
		for (Cluster<KMeansModel> c : clustering.getAllClusters()) {
			log.append("#" + clusterId + " [" + c.size() + ", centroid " + c.getModel().getPrototype().toString() + "]");
			
			
			List<Integer> idList = new ArrayList<Integer>();
			for (DBIDIter it = c.getIDs().iter(); it.valid(); it.advance())
				idList.add(ids.getOffset(it));
			
			Collections.sort(idList);
			for (Integer id : idList)
				log.append(" " + id);
			
			
			log.append("\n");
			++clusterId;
		}
		System.out.println(log);

		//Log dummy encoding
		log = new StringBuilder();
		log.append("\nKMeans clustering dummy encoded:\n");
		for (int row = 0; row < 3; row++) { 
			String rowString = "";
			for (int col = 0; col < dummy[row].length; col++)
				rowString += dummy[row][col] + " ";
			log.append(rowString + "\n");
		}
		log.append("...\n");
		System.out.println(log);
	}

	public Clustering<KMeansModel> getClustering() {
		return clustering;
	}

	public double[][] getDummy() {
		return dummy;
	}

	public DBIDRange getIds() {
		return ids;
	}
	
}
