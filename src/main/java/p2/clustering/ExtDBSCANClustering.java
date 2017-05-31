package p2.clustering;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import de.lmu.ifi.dbs.elki.data.Cluster;
import de.lmu.ifi.dbs.elki.data.Clustering;
import de.lmu.ifi.dbs.elki.data.model.Model;
import de.lmu.ifi.dbs.elki.database.ids.DBIDIter;
import de.lmu.ifi.dbs.elki.database.ids.DBIDRange;

public class ExtDBSCANClustering {
	
	private Clustering<Model> clustering;
	private double[][] dummy;
	private DBIDRange ids;
	
	public ExtDBSCANClustering(Clustering<Model> clustering, double[][] dummy, DBIDRange ids) {
		this.clustering = clustering;
		this.dummy = dummy;
		this.ids = ids;
		
		//DEBUG log clustering
		/*StringBuilder log = new StringBuilder();
		log.append("\nDBSCAN Clustering:\n");
		int clusterId = 1;
		for (Cluster<Model> c : clustering.getAllClusters()) {
			log.append("#" + clusterId + " [" + c.size() + "]");
			
			List<Integer> idList = new ArrayList<Integer>();
			for (DBIDIter it = c.getIDs().iter(); it.valid(); it.advance())
				idList.add(ids.getOffset(it));
			
			Collections.sort(idList);
			for (Integer id : idList)
				log.append(" " + id);
			
			log.append("\n");
			++clusterId;
		}
		System.out.println(log);*/
		
		//DEBUG log dummy encoding
		/*log = new StringBuilder();
		log.append("\nDBSCAN Dummy:\n");
		for (int row = 0; row < dummy.length; row++) { 
			String rowString = row + ": ";
			for (int col = 0; col < dummy[row].length; col++)
				rowString += dummy[row][col] + " ";
			log.append(rowString + "\n");
		}
		System.out.println(log);*/
	}

	public Clustering<Model> getClustering() {
		return clustering;
	}

	public double[][] getDummy() {
		return dummy;
	}


	
	
}
