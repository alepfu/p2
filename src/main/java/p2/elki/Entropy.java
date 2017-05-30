package p2.elki;


public class Entropy {

	public double entropyFirst = -1.0;
	public double entropySecond = -1.0;
	public double entropyJoint = -1.0;

	public Entropy(ClusterContingencyTable table) {
		
		double norm = 1.0 / table.contingency[table.size1][table.size2];
		{
			entropyFirst = 0.0;
			for (int i1 = 0; i1 < table.size1; i1++) {
				if (table.contingency[i1][table.size2] > 0) {
					double probability = norm * table.contingency[i1][table.size2];
					entropyFirst -= probability * Math.log(probability);
				}
			}
		}
		{
			entropySecond = 0.0;
			for (int i2 = 0; i2 < table.size2; i2++) {
				if (table.contingency[table.size1][i2] > 0) {
					double probability = norm * table.contingency[table.size1][i2];
					entropySecond -= probability * Math.log(probability);
				}
			}
		}
		{
			entropyJoint = 0.0;
			for (int i1 = 0; i1 < table.size1; i1++) {
				for (int i2 = 0; i2 < table.size2; i2++) {
					if (table.contingency[i1][i2] > 0) {
						double probability = norm * table.contingency[i1][i2];
						entropyJoint -= probability * Math.log(probability);
					}
				}
			}
		}
	}

	public double entropyFirst() {
		return entropyFirst;
	}

	public double entropySecond() {
		return entropySecond;
	}

	public double entropyJoint() {
		return entropyJoint;
	}

	public double entropyConditionalFirst() {
		return (entropyJoint() - entropySecond());
	}

	public double entropyConditionalSecond() {
		return (entropyJoint() - entropyFirst());
	}

	public double entropyPowers() {
		return (2 * entropyJoint() / (entropyFirst() + entropySecond()) - 1);
	}

	public double entropyMutualInformation() {
		return (entropyFirst() + entropySecond() - entropyJoint());
	}

	public double entropyNMIJoint() {
		if (entropyJoint() == 0) {
			return 0;
		}
		return (entropyMutualInformation() / entropyJoint());
	}

	public double entropyNMIMin() {
		return (entropyMutualInformation() / Math.min(entropyFirst(), entropySecond()));
	}

	public double entropyNMIMax() {
		return (entropyMutualInformation() / Math.max(entropyFirst(), entropySecond()));
	}

	public double entropyNMISum() {
		return (2 * entropyMutualInformation() / (entropyFirst() + entropySecond()));
	}

	public double entropyNMISqrt() {
		if (entropyFirst() * entropySecond() <= 0) {
			return entropyMutualInformation();
		}
		return (entropyMutualInformation() / Math.sqrt(entropyFirst() * entropySecond()));
	}

	public double variationOfInformation() {
		return (2 * entropyJoint() - (entropyFirst() + entropySecond()));
	}

	public double normalizedVariationOfInformation() {
		return (1.0 - (entropyMutualInformation() / entropyJoint()));
	}
}