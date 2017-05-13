package p2.clustering;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Util {
	
	public double[][] loadMergedData(String file) {
		
		Map<String, String> header = this.getHeader(file);
		int numDimensions = Integer.parseInt(header.get("numDimensions"));
		int numPointsPerCluster = Integer.parseInt(header.get("numPointsPerCluster"));
		int numClusters = Integer.parseInt(header.get("numClusters"));
		int numPoints = numClusters * numPointsPerCluster;

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
	
	public double[][] loadGaussianData(String file) {
		
		double[][] mergedData = this.loadMergedData(file);
		
		Map<String, String> header = this.getHeader(file);
		int numDimensions = Integer.parseInt(header.get("numDimensions"));
		
		//Copy gaussian data points to double[][] array
		double[][] gaussianData = new double[mergedData.length][numDimensions];
		for (int i = 0; i < mergedData.length; i++)
			for (int j = 0; j < numDimensions; j++)
				gaussianData[i][j] = mergedData[i][j];
		
		return gaussianData; 
	}
	
	public double[][] loadDensityData(String file) {
		
		double[][] mergedData = this.loadMergedData(file);
		
		Map<String, String> header = this.getHeader(file);
		int numDimensions = Integer.parseInt(header.get("numDimensions"));
		
		//Copy density data points to double[][] array
		double[][] densityData = new double[mergedData.length][numDimensions];
		for (int i = 0; i < mergedData.length; i++)
			for (int j = 0; j < numDimensions; j++)
				densityData[i][j] = mergedData[i][j +  numDimensions];
		
		return densityData; 
	}

	public Map<String, String> getHeader(String file) {
		
		Map<String, String> header = new HashMap<String, String>();

		try {
			String line = "";
			BufferedReader br = new BufferedReader(new FileReader(file));
			while ((line = br.readLine()) != null) {

				if (line.startsWith("#")) {
					String[] split = line.split("=");
					header.put(split[0].substring(1), split[1]);
				}
			}

			br.close();

		} catch (Exception e) {
			e.printStackTrace();
		}

		return header;
	}
	
	public List<int[]> getClusterAssignments(String file) {
		
//		Map<String, String> header = this.getHeader(file);
//		int numPointsPerCluster = Integer.parseInt(header.get("numPointsPerCluster"));
//		int numClusters = Integer.parseInt(header.get("numClusters"));

		List<int[]> assignments = new ArrayList<int[]>();
		
		try {
			String line = "";
			BufferedReader br = new BufferedReader(new FileReader(file));
			while ((line = br.readLine()) != null) {

				if (!line.startsWith("#")) {
					String[] split = line.split(",");
					int[] a = new int[2];
					a[0] = Integer.parseInt(split[0]);
					a[1] = Integer.parseInt(split[split.length - 1]);
					assignments.add(a);
				}
			}

			br.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return assignments;
	}
}
