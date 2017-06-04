package p2.util;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashMap;
import java.util.Map;

public class DataLoaderUtil {
	
	private String file;
	
	private int numDimensions;
	private int numPointsCluster;
	private int numClusters;
	private int numPoints;
	private long seed; 
	
	private double[][] data;
	
	public DataLoaderUtil(String file) {
		
		this.file = file;
		
		Map<String, String> header = new HashMap<String, String>();

		try {
			String line = "";
			BufferedReader br = new BufferedReader(new FileReader(file));
			while ((line = br.readLine()) != null) {

				if (line.startsWith("#")) {
					String[] split = line.split("=");
					header.put(split[0].substring(1), split[1]);
				}
				else
					break;
			}

			br.close();
			
			numDimensions = Integer.parseInt(header.get("numDimensions"));
			numPointsCluster = Integer.parseInt(header.get("numPointsCluster"));
			numClusters = Integer.parseInt(header.get("numClusters"));
			numPoints = numClusters * numPointsCluster;
			seed = Integer.parseInt(header.get("seed"));

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public double[][] loadData() {
		
		if (data != null)
			return data;
		
		data = new double[numPoints][numDimensions * 2];
		
		try {
			String line = "";
			int id = 0;
			BufferedReader br = new BufferedReader(new FileReader(file));
			while ((line = br.readLine()) != null) {
			
				if (!line.startsWith("#")) {
					String[] split = line.split(" ");
					for (int j = 0; j < (numDimensions * 2); j++)
						data[id][j] = Double.parseDouble(split[j]);
					++id;
				}
			}

			br.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return data;
	}
	
	public double[][] loadGaussianData() {
		
		if (data == null)
			data = this.loadData();
		
		//Copy gaussian data points to double[][] array
		double[][] gaussianData = new double[data.length][numDimensions];
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < numDimensions; j++)
				gaussianData[i][j] = data[i][j];
		
		return gaussianData; 
	}
	
	public double[][] loadDensityData() {
		
		if (data == null)
			data = this.loadData();
		
		//Copy density data points to double[][] array
		double[][] densityData = new double[data.length][numDimensions];
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < numDimensions; j++)
				densityData[i][j] = data[i][j +  numDimensions];
		
		return densityData; 
	}

	public int getNumDimensions() {
		return numDimensions;
	}

	public int getNumPointsCluster() {
		return numPointsCluster;
	}

	public int getNumClusters() {
		return numClusters;
	}

	public int getNumPoints() {
		return numPoints;
	}
	
}
