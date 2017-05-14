package p2.clustering;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashMap;
import java.util.Map;

public class LoadDataUtil {
	
	private String file;
	
	private int numDimensions;
	private int numPointsPerCluster;
	private int numClusters;
	private int numPoints;
	
	public LoadDataUtil(String file) {
		
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
			numPointsPerCluster = Integer.parseInt(header.get("numPointsPerCluster"));
			numClusters = Integer.parseInt(header.get("numClusters"));
			numPoints = numClusters * numPointsPerCluster;

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public double[][] loadMergedData() {
		
		double[][] data = new double[numPoints][numDimensions * 2];
		
		try {
			String line = "";
			BufferedReader br = new BufferedReader(new FileReader(file));
			while ((line = br.readLine()) != null) {
			
				if (!line.startsWith("#")) {
				
					String[] split = line.split(",");
					int id = Integer.parseInt(split[0]);
					
					for (int j = 0; j < (numDimensions * 2); j++)
						data[id][j] = Double.parseDouble(split[j + 1]);
				}
			}

			br.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return data;
	}
	
	public double[][] loadGaussianData() {
		
		double[][] mergedData = this.loadMergedData();
		
		//Copy gaussian data points to double[][] array
		double[][] gaussianData = new double[mergedData.length][numDimensions];
		for (int i = 0; i < mergedData.length; i++)
			for (int j = 0; j < numDimensions; j++)
				gaussianData[i][j] = mergedData[i][j];
		
		return gaussianData; 
	}
	
	public double[][] loadDensityData() {
		
		double[][] mergedData = this.loadMergedData();
		
		//Copy density data points to double[][] array
		double[][] densityData = new double[mergedData.length][numDimensions];
		for (int i = 0; i < mergedData.length; i++)
			for (int j = 0; j < numDimensions; j++)
				densityData[i][j] = mergedData[i][j +  numDimensions];
		
		return densityData; 
	}

	public int getNumDimensions() {
		return numDimensions;
	}

	public int getNumPointsPerCluster() {
		return numPointsPerCluster;
	}

	public int getNumClusters() {
		return numClusters;
	}

	public int getNumPoints() {
		return numPoints;
	}
	
	
}
