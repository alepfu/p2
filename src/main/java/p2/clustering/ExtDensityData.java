package p2.clustering;

import java.text.NumberFormat;

public class ExtDensityData {
	
	private double[][] data;
	private double[][] dummy;
	
	//TODO In theory we don't need 2 separate ExtData classes, reduce to one
	
	public ExtDensityData(double[][] data, double[][] dummy) {
		this.data = data;
		this.dummy = dummy;
	}

	/**
	 * Normalizes the dummy encoding to the variance of the data, then appends the dummy encoding to the data.
	 * @return Concat between data and dummy encoding.
	 */
	public double[][] getExtData() {
		
		if (data.length != dummy.length)
			throw new IllegalStateException("data.length != dummy.length");
		
		int numRows = data.length;
		int numColsData = data[0].length;
		int numColsDummy = dummy[0].length;
		
		//Normalization of dummy
		StatisticsUtil statUtil = new StatisticsUtil();
		double sumVarData = statUtil.calcSumVar(data, numColsData, numRows);
		double sumVarDummy = statUtil.calcSumVar(dummy, numColsDummy, numRows);
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
		
		//DEBUG: print data with normalized dummy encoding attached
		/*NumberFormat nf = NumberFormat.getInstance();
		nf.setMinimumFractionDigits(2);
		nf.setMaximumFractionDigits(2);
		StringBuilder log = new StringBuilder("Extended Density Data:\n");
		for (int row = 0; row < numRows; row++) {
			for (int col = 0; col < numColsConcat; col++)
				log.append(nf.format(extData[row][col]) + " ");
			log.append("\n");
		}
		System.out.println(log);*/
		
		return extData;
	}
}
