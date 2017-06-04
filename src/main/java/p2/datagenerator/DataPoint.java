package p2.datagenerator;

import java.util.Arrays;

public class DataPoint {
	
	private double[] gaussFeatures;
	private double[] densityFeatures;
	private String clusterLabel;
	
	public DataPoint(double[] gaussFeatures, double[] densityFeatures, String cluserLabel) {
		
		this.gaussFeatures = gaussFeatures;
		this.densityFeatures = densityFeatures;
		this.clusterLabel = cluserLabel;
	}
	
	public double[] getGaussFeatures() {
		return gaussFeatures;
	}
	
	public void setGaussFeatures(double[] gaussFeatures) {
		this.gaussFeatures = gaussFeatures;
	}
	
	public double[] getDensityFeatures() {
		return densityFeatures;
	}
	
	public void setDensityFeatures(double[] densityFeatures) {
		this.densityFeatures = densityFeatures;
	}
	
	public String getClusterLabel() {
		return clusterLabel;
	}

	public void setClusterLabel(String clusterLabel) {
		this.clusterLabel = clusterLabel;
	}

	public String toString() {
		
		StringBuilder sb = new StringBuilder();
		sb.append(Arrays.toString(gaussFeatures).replace("[", "").replace("]", "").replace(",", ""));
		sb.append(" ");
		sb.append(Arrays.toString(densityFeatures).replace("[", "").replace("]", "").replace(",", ""));
		sb.append(" ");
		sb.append(clusterLabel);
		
		return sb.toString();
	}
}
