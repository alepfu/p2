package p2.util;

/**
 * Copyright (c) 2017 Alexander Pfundner
 * 
 * Integration of Density-based and Partitioning-based Clustering Methods
 * 
 */

import java.io.BufferedReader;
import java.io.FileReader;

public class DataLoaderUtil {
	
	private String file;
	
	private double[][] data;
	
	public DataLoaderUtil(String file) {
		
		this.file = file;
	}

	public double[][] loadData() {
		
		if (data != null)
			return data;
		
		data = new double[Config.numPoints][4];
		
		try {
			String line = "";
			int id = 0;
			BufferedReader br = new BufferedReader(new FileReader(file));
			while ((line = br.readLine()) != null) {
				
				String[] split = line.split(" ");
				for (int j = 0; j < (Config.numDimPerType * 2); j++)
					data[id][j] = Double.parseDouble(split[j]);
				++id;
			}

			br.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return data;
	}
	
	public double[][] getGaussData() {
		
		if (data == null)
			data = this.loadData();
		
		//Copy gauss data points to double[][] array
		double[][] gaussData = new double[data.length][Config.numDimPerType];
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < Config.numDimPerType; j++)
				gaussData[i][j] = data[i][j];
		
		return gaussData; 
	}
	
	public double[][] getDensityData() {
		
		if (data == null)
			data = this.loadData();
		
		//Copy density data points to double[][] array
		double[][] densityData = new double[data.length][Config.numDimPerType];
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < Config.numDimPerType; j++)
				densityData[i][j] = data[i][j + Config.numDimPerType];
		
		return densityData; 
	}
	
}
