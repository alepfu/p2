package p2.datagenerator;

import java.util.Arrays;

public class MergedDataPoint {
	
	double[] gaussFeatures;
	double[] densityFeatures;
	int id;
	int clusterId;
	
	public MergedDataPoint(double[] gaussFeatures, double[] densityFeatures, int id, int clusterId) {
		
		this.gaussFeatures = gaussFeatures;
		this.densityFeatures = densityFeatures;
		this.id = id;
		this.clusterId = clusterId;
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
	
	public int getId() {
		return id;
	}
	
	public void setId(int id) {
		this.id = id;
	}
	
	public int getClusterId() {
		return clusterId;
	}
	
	public void setClusterId(int clusterId) {
		this.clusterId = clusterId;
	}
	
	public String toString() {
		
		StringBuilder sb = new StringBuilder();
		sb.append(id).append(",");
		sb.append(Arrays.toString(gaussFeatures).replace(" ", "").replace("[", "").replace("]", "")).append(",");
		sb.append(Arrays.toString(densityFeatures).replace(" ", "").replace("[", "").replace("]", "")).append(",");
		sb.append(clusterId);
		
		return sb.toString();
	}
}
