package p2.old;

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
		
		double sumVarData = 0;
		for (int col = 0; col < numCols; col++) {
			for (int row = 0; row < numRows; row++) 
				stats.addValue(a[row][col]);
			sumVarData += stats.getVariance();
			stats = new DescriptiveStatistics();
		}
		
		return sumVarData;
	}
}
