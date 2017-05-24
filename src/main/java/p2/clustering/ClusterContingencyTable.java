package p2.clustering;

import java.util.Iterator;
import java.util.List;

import de.lmu.ifi.dbs.elki.data.Cluster;
import de.lmu.ifi.dbs.elki.data.Clustering;
import de.lmu.ifi.dbs.elki.database.ids.DBIDUtil;
import de.lmu.ifi.dbs.elki.database.ids.DBIDs;
import de.lmu.ifi.dbs.elki.utilities.BitsUtil;


public class ClusterContingencyTable {
  /**
   * Noise cluster handling
   */
  protected boolean breakNoiseClusters = false;

  /**
   * Self pairing
   */
  protected boolean selfPairing = true;

  /**
   * Number of clusters.
   */
  protected int size1 = -1, size2 = -1;

  /**
   * Contingency matrix
   */
  protected int[][] contingency = null;

  /**
   * Noise flags
   */
  protected long[] noise1 = null, noise2 = null;

 

  /**
   * Constructor.
   * 
   * @param selfPairing Build self-pairs
   * @param breakNoiseClusters Break noise clusters into individual objects
   */
  public ClusterContingencyTable(boolean selfPairing, boolean breakNoiseClusters) {
    super();
    this.selfPairing = selfPairing;
    this.breakNoiseClusters = breakNoiseClusters;
  }

  /**
   * Process two clustering results.
   * 
   * @param result1 First clustering
   * @param result2 Second clustering
   */
  public void process(Clustering<?> result1, Clustering<?> result2) {
    // Get the clusters
    final List<? extends Cluster<?>> cs1 = result1.getAllClusters();
    final List<? extends Cluster<?>> cs2 = result2.getAllClusters();

    // Initialize
    size1 = cs1.size();
    size2 = cs2.size();
    contingency = new int[size1 + 2][size2 + 2];
    noise1 = BitsUtil.zero(size1);
    noise2 = BitsUtil.zero(size2);

    // Fill main part of matrix
    {
      final Iterator<? extends Cluster<?>> it2 = cs2.iterator();
      for(int i2 = 0; it2.hasNext(); i2++) {
        final Cluster<?> c2 = it2.next();
        if(c2.isNoise()) {
          BitsUtil.setI(noise2, i2);
        }
        contingency[size1 + 1][i2] = c2.size();
        contingency[size1 + 1][size2] += c2.size();
      }
    }
    final Iterator<? extends Cluster<?>> it1 = cs1.iterator();
    for(int i1 = 0; it1.hasNext(); i1++) {
      final Cluster<?> c1 = it1.next();
      if(c1.isNoise()) {
        BitsUtil.setI(noise1, i1);
      }
      final DBIDs ids = DBIDUtil.ensureSet(c1.getIDs());
      contingency[i1][size2 + 1] = c1.size();
      contingency[size1][size2 + 1] += c1.size();

      final Iterator<? extends Cluster<?>> it2 = cs2.iterator();
      for(int i2 = 0; it2.hasNext(); i2++) {
        final Cluster<?> c2 = it2.next();
        int count = DBIDUtil.intersectionSize(ids, c2.getIDs());
        contingency[i1][i2] = count;
        contingency[i1][size2] += count;
        contingency[size1][i2] += count;
        contingency[size1][size2] += count;
      }
    }
  }

  @Override
  public String toString() {
    StringBuilder buf = new StringBuilder();
    if(contingency != null) {
      for(int i1 = 0; i1 < size1 + 2; i1++) {
        if(i1 >= size1) {
          buf.append("------\n");
        }
        for(int i2 = 0; i2 < size2 + 2; i2++) {
          if(i2 >= size2) {
            buf.append("| ");
          }
          buf.append(contingency[i1][i2]).append(' ');
        }
        buf.append('\n');
      }
    }
    return buf.toString();
  }

}
