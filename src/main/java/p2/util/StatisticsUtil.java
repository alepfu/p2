package p2.util;

/**
 * Copyright (c) 2017 Alexander Pfundner
 * 
 * Integration of Density-based and Partitioning-based Clustering Methods
 * 
 */

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

public class StatisticsUtil {
	
	/**
	 * Calculates the sum of all column-variances for a 2D array.
	 * @param a A 2D array
	 * @param numCols The number of columns in a.
	 * @param numRows The number of rows in a.
	 * @return The column-variances summed up.
	 */
	public double calcSumVar(double[][] a, int numCols, int numRows) {
		DescriptiveStatistics stats = new DescriptiveStatistics();
		
		double sumVar = 0;
		for (int col = 0; col < numCols; col++) {
			for (int row = 0; row < numRows; row++) 
				stats.addValue(a[row][col]);
			sumVar += stats.getVariance();
			stats = new DescriptiveStatistics();
		}
		
		return sumVar;
	}
	
	
}
