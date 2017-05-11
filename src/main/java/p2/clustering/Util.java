package p2.clustering;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class Util {
	
	public double[][] loadData(String file) {
		
		ArrayList<double[]> rows = new ArrayList<double[]>();
		String line = "";

		//Read lines from CSV file, omit header
		try {
			BufferedReader br = new BufferedReader(new FileReader(file));
			while ((line = br.readLine()) != null) {

				if (!line.startsWith("#")) {
					double[] doubleValues = Arrays.stream(line.split(",")).mapToDouble(Double::parseDouble).toArray();
					rows.add(doubleValues);
				}
			}

			br.close();

		} catch (Exception e) {
			e.printStackTrace();
		}

		//Copy values to double[][] array, omit the last 2 columns holding the point and cluster ids
		double[][] data = new double[rows.size()][rows.get(0).length - 2];
		for (int i = 0; i < rows.size(); i++)
			for (int j = 0; j < rows.get(0).length - 2; j++)
				data[i][j] = rows.get(i)[j];
		
		return data;
	}
	
	public double[][] loadGaussianData(String file) {
		
		double[][] mergedData = this.loadData(file);
		
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
		
		double[][] mergedData = this.loadData(file);
		
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
		
		String line = "";
		try {
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
}
