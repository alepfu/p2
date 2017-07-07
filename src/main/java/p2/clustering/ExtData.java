package p2.clustering;

import java.text.NumberFormat;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

public class ExtData {
	
	private double[][] data;
	private double[][] dummy;
	
	public ExtData(double[][] data, double[][] dummy) {
		this.data = data;
		this.dummy = dummy;
	}

	/**
	 * Normalizes the dummy encoding to the variance of the data, then appends the dummy encoding to the data.
	 * Handles the 1st run, which has no dummy encoding yet.
	 * @return Concat between data and dummy encoding.
	 */
	public double[][] getExtData() {
		
		if (dummy == null)  //For handling the 1st run
			return data;
		
		if (data.length != dummy.length)
			throw new IllegalStateException("data.length != dummy.length");
		
		int numRows = data.length;
		int numColsData = data[0].length;
		int numColsDummy = dummy[0].length;
		
		//Normalization of dummy
		double sumVarData = calcSumVar(data, numColsData, numRows);
		double sumVarDummy = calcSumVar(dummy, numColsDummy, numRows);
		double factor = Math.sqrt(sumVarData / sumVarDummy);
		for (int row = 0; row < numRows; row++)
			for (int col = 0; col < numColsDummy; col++)
				dummy[row][col] *= factor;

		//Concat data and dummy
		int numColsConcat = data[0].length + dummy[0].length;
		
		double[][] extData = new double[numRows][numColsConcat];
		
		//Insert values from data
		for (int row = 0; row < numRows; row++)
			for (int col = 0; col < numColsData; col++)
				extData[row][col] = data[row][col];
		
		//Insert values from dummy
		for (int row = 0; row < numRows; row++)
			for (int col = numColsData; col < numColsConcat; col++)
				extData[row][col] = dummy[row][col - numColsData];		
		
		
		//DEBUG log
		/*NumberFormat nf = NumberFormat.getInstance();
		nf.setMinimumFractionDigits(2);
		nf.setMaximumFractionDigits(2);
		StringBuilder log = new StringBuilder("\nExtended data:\n");
		for (int row = 0; row < 3; row++) {
			for (int col = 0; col < numColsConcat; col++)
				log.append(nf.format(extData[row][col]) + " ");
			log.append("\n");
		}
		log.append("...");
		System.out.println(log);*/
		
		return extData;
	}
	
	/**
	 * Calculates the sum of all column-variances for a 2D array.
	 * @param a A 2D array
	 * @param numCols The number of columns in a.
	 * @param numRows The number of rows in a.
	 * @return The column-variances summed up.
	 */
	private double calcSumVar(double[][] a, int numCols, int numRows) {
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
